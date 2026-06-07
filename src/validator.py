import os
import sys
import numpy as np
from rdkit import Chem
import mdtraj as md
import yaml
import logging

logger = logging.getLogger('simdock')

# OpenMM imports
try:
    from openmm.app import PDBFile, Modeller, ForceField, DCDReporter, PDBReporter, HBonds
    from openmm import LangevinIntegrator, MonteCarloBarostat, Platform, unit
    import openmm as mm
    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    OPENMM_AVAILABLE = True
except ImportError as e:
    import sys
    logger.warning(f"OpenMM or openmmforcefields import failed: {e}")
    OPENMM_AVAILABLE = False

def run_md_simulation(complex_pdb, ligand_pdb, config, quick_test=True):
    """
    Runs an OpenMM molecular dynamics simulation following the document protocol:
      1. Solvation: TIP3P water box with configurable padding (default 10Å)
      2. Energy Minimization: 500 steps
      3. NVT Equilibration: Short equilibration with velocity generation
      4. NPT Equilibration: Short equilibration with barostat for density relaxation
      5. Production Run: Langevin dynamics (10ns default, or 500 steps for quick_test)
    Falls back to mock trajectory if OpenMM or dependencies are unavailable.
    """
    out_dir = config['paths'].get('output_dir', './results')
    os.makedirs(out_dir, exist_ok=True)
    trajectory_dcd = os.path.join(out_dir, "trajectory.dcd")
    
    if not OPENMM_AVAILABLE:
        logger.warning("OpenMM or openmmforcefields not available. Falling back to mock trajectory.")
        return generate_mock_trajectory(complex_pdb, trajectory_dcd), True
        
    try:
        # 1. Load ligand and register GAFF parameters
        ligand_mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
        if ligand_mol is None:
            ligand_mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=True)
            
        ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        try:
            from openff.toolkit.topology import Molecule
            off_mol = Molecule.from_rdkit(ligand_mol, allow_undefined_stereo=True)
            gaff = GAFFTemplateGenerator(off_mol, forcefield='gaff-2.11')
        except ImportError:
            raise RuntimeError("openff-toolkit is not installed; small-molecule parameterization is unavailable.")
        except Exception as e:
            raise RuntimeError(f"Failed to parameterize ligand with GAFF: {e}")
            
        ff.registerTemplateGenerator(gaff.generator)
        
        # 2. Load complex and solvate with configurable padding
        padding_angstrom = config['md'].get('solvation_padding_A', 10.0)
        padding_nm = padding_angstrom / 10.0  # Convert Å to nm
        
        pdb = PDBFile(complex_pdb)
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addSolvent(ff, model='tip3p', padding=padding_nm * unit.nanometers)
        
        # Verify parameterization succeeded
        try:
            test_system = ff.createSystem(modeller.topology)
            for force in test_system.getForces():
                if hasattr(force, 'getNumParticles'):
                    assert force.getNumParticles() == modeller.topology.getNumAtoms(), \
                        'GAFF2 failed: atom count mismatch — some atoms were not parameterized'
        except Exception as e:
            raise RuntimeError(f"Ligand parameterization verification failed: {e}")
        
        # 3. Create OpenMM System
        system = ff.createSystem(modeller.topology, 
                                 nonbondedMethod=mm.app.PME,
                                 nonbondedCutoff=0.9 * unit.nanometer, 
                                 constraints=HBonds)
                                 
        # 4. Select compute platform
        try:
            platform = Platform.getPlatformByName('CUDA')
        except Exception:
            try:
                platform = Platform.getPlatformByName('OpenCL')
            except Exception:
                platform = Platform.getPlatformByName('CPU')
                
        logger.info(f"Running simulation on platform: {platform.getName()}")
        
        # =========================================================
        # STEP A: Energy Minimization
        # =========================================================
        integrator_min = LangevinIntegrator(300 * unit.kelvin, 
                                            1.0 / unit.picoseconds, 
                                            2.0 * unit.femtoseconds)
        
        simulation = mm.app.Simulation(modeller.topology, system, integrator_min, platform)
        simulation.context.setPositions(modeller.positions)
        
        min_steps = config['md'].get('minimization_steps', 500)
        logger.info(f"Minimizing system energy (max {min_steps} steps)...")
        simulation.minimizeEnergy(maxIterations=min_steps)
        
        # Save minimized positions
        minimized_positions = simulation.context.getState(getPositions=True).getPositions()
        
        if quick_test:
            # =========================================================
            # Quick test mode: abbreviated protocol
            # =========================================================
            # Short NVT (50 steps)
            logger.info("Running abbreviated NVT equilibration (50 steps)...")
            simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
            simulation.step(50)
            
            # Short NPT (50 steps)
            logger.info("Running abbreviated NPT equilibration (50 steps)...")
            barostat = MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin)
            system.addForce(barostat)
            simulation.context.reinitialize(preserveState=True)
            simulation.step(50)
            
            # Short production (500 steps)
            logger.info("Running abbreviated production MD (500 steps)...")
            report_interval = 100
            simulation.reporters.append(DCDReporter(trajectory_dcd, report_interval))
            simulation.step(500)
            
        else:
            # =========================================================
            # Full protocol (per architecture document)
            # =========================================================
            
            # Gradual Heating Ramp (from 0K to 300K in steps of 10K)
            logger.info("Gradually heating system from 0K to 300K...")
            heating_steps_per_increment = 250   # 0.5ps per 10K step (dt = 2fs)
            for temp_k in range(10, 310, 10):
                simulation.integrator.setTemperature(temp_k * unit.kelvin)
                simulation.context.setVelocitiesToTemperature(temp_k * unit.kelvin)
                simulation.step(heating_steps_per_increment)
            logger.info("Heating complete. Starting NVT equilibration.")
            
            # Ensure integrator is set to final temperature
            simulation.integrator.setTemperature(300 * unit.kelvin)
            
            # NVT Equilibration
            nvt_ps = config['md'].get('equilibration_nvt_ps', 100)
            nvt_steps = int(nvt_ps / 0.002)  # dt = 2fs
            logger.info(f"Running NVT equilibration ({nvt_ps}ps, {nvt_steps} steps)...")
            simulation.step(nvt_steps)
            
            # NPT Equilibration
            npt_ps = config['md'].get('equilibration_npt_ps', 100)
            npt_steps = int(npt_ps / 0.002)
            logger.info(f"Running NPT equilibration ({npt_ps}ps, {npt_steps} steps)...")
            barostat = MonteCarloBarostat(1.0 * unit.atmosphere, 300 * unit.kelvin)
            system.addForce(barostat)
            simulation.context.reinitialize(preserveState=True) # State is already preserved
            simulation.step(npt_steps)
            
            # Production Run
            production_ns = config['md'].get('production_ns', 10)
            production_steps = int((production_ns * 1000) / 0.002)  # ns → ps → steps
            report_interval = int(10 / 0.002)  # Report every 10ps
            logger.info(f"Running production MD ({production_ns}ns, {production_steps} steps)...")
            simulation.reporters.append(DCDReporter(trajectory_dcd, report_interval))
            simulation.step(production_steps)
        
        logger.info("MD simulation completed successfully!")
        return trajectory_dcd, False
        
    except Exception as e:
        logger.warning(f"OpenMM simulation failed with error: {e}. Falling back to mock trajectory.")
        return generate_mock_trajectory(complex_pdb, trajectory_dcd), True

def generate_mock_trajectory(complex_pdb, output_dcd):
    """
    Generates a mock DCD trajectory containing slightly perturbed coordinates
    of the starting complex to allow offline testing and validation of analysis code.
    """
    logger.info("Generating mock DCD trajectory...")
    traj = md.load(complex_pdb)
    
    frames = []
    r = np.random.RandomState(42)
    
    # Select ligand indices using correct MDTraj comparison operators
    ligand_idx = traj.topology.select("resname UNL or resname LIG or resname POL")
    
    for f in range(10):
        frame = md.Trajectory(np.copy(traj.xyz), traj.topology)
        if len(ligand_idx) > 0:
            noise = r.normal(0.0, 0.02, size=(1, len(ligand_idx), 3))
            frame.xyz[0, ligand_idx, :] += noise[0]
        frames.append(frame)
        
    mock_traj = md.join(frames)
    mock_traj.save_dcd(output_dcd)
    logger.info(f"Mock trajectory saved to {output_dcd}")
    return output_dcd

def analyze_trajectory(trajectory_dcd, complex_pdb, config, ligand_resname='UNL', enzyme_data=None):
    """
    Calculates ligand RMSD and catalytic triad attack distance over the MD trajectory.
    """
    traj = md.load(trajectory_dcd, top=complex_pdb)
    
    # 1. Superpose protein to align frames (align vs frame 0)
    protein_idx = traj.topology.select("protein")
    if len(protein_idx) > 0:
        traj.superpose(traj, frame=0, atom_indices=protein_idx)
        
    # 2. Compute ligand RMSD vs frame 0 (in Angstroms)
    ligand_idx = traj.topology.select(f"resname {ligand_resname} or resname LIG or resname POL")
    if len(ligand_idx) == 0:
        ligand_idx = traj.topology.select("not protein and not water")
        
    if len(ligand_idx) > 0:
        rmsd = md.rmsd(traj, traj, frame=0, atom_indices=ligand_idx) * 10.0 # nm to A
    else:
        rmsd = np.zeros(traj.n_frames)
        
    # 3. Compute catalytic distance over trajectory
    # Get nucleophile from enzyme_data
    nuc_idx = []
    if enzyme_data:
        nuc_res_num = enzyme_data.get('nucleophile_res_num', 160)
        nuc_res_name = enzyme_data.get('nucleophile_res_name', 'SER')
        nuc_atom = enzyme_data.get('nucleophile_atom_name', 'OG')
        nuc_idx = traj.topology.select(f"resSeq {nuc_res_num} and name {nuc_atom}")
        
    if len(nuc_idx) == 0:
        nuc_idx = traj.topology.select("resname SER and residue 160 and name OG")
    if len(nuc_idx) == 0:
        nuc_idx = traj.topology.select("resname SER and name OG")
    if len(nuc_idx) == 0:
        nuc_idx = traj.topology.select("name OG")
        
    # Identify scissile carbons in the ligand geometrically on the first frame
    ligand_carbons = [idx for idx in ligand_idx if traj.topology.atom(idx).element.symbol == 'C']
    ligand_oxygens = [idx for idx in ligand_idx if traj.topology.atom(idx).element.symbol == 'O']
    
    scissile_bond_type = enzyme_data.get('scissile_bond_type', 'ester_carbonyl') if enzyme_data else 'ester_carbonyl'
    scissile_bond_position = enzyme_data.get('scissile_bond_position', 'terminal') if enzyme_data else 'terminal'
    
    scissile_carbons = []
    xyz0 = traj.xyz[0]
    
    if scissile_bond_type == 'glycosidic_carbon':
        # Find carbon atoms with >= 2 oxygen neighbors within 1.5 A
        for c_idx in ligand_carbons:
            o_neighbors = 0
            for o_idx in ligand_oxygens:
                dist = np.linalg.norm(xyz0[c_idx] - xyz0[o_idx]) * 10.0
                if dist <= 1.5:
                    o_neighbors += 1
            if o_neighbors >= 2:
                scissile_carbons.append(c_idx)
    else:
        # Default: ester_carbonyl
        for c_idx in ligand_carbons:
            is_carbonyl = False
            for o_idx in ligand_oxygens:
                dist = np.linalg.norm(xyz0[c_idx] - xyz0[o_idx]) * 10.0
                if 1.15 <= dist <= 1.30:
                    is_carbonyl = True
                    break
            if is_carbonyl:
                scissile_carbons.append(c_idx)
                
        if scissile_bond_position == 'terminal':
            # Filter to keep only terminal ones
            terminal_carbons = []
            for c_idx in scissile_carbons:
                o_neighbors = [o_idx for o_idx in ligand_oxygens
                               if np.linalg.norm(xyz0[c_idx] - xyz0[o_idx]) * 10.0 <= 1.45]
                is_terminal = True
                for o_idx in o_neighbors:
                    other_c = [other_idx for other_idx in ligand_carbons
                               if other_idx != c_idx and
                               np.linalg.norm(xyz0[o_idx] - xyz0[other_idx]) * 10.0 <= 1.6]
                    if len(other_c) > 0:
                        is_terminal = False
                        break
                if is_terminal:
                    terminal_carbons.append(c_idx)
            if terminal_carbons:
                scissile_carbons = terminal_carbons
                
    distances = []
    if len(nuc_idx) > 0 and len(scissile_carbons) > 0:
        nuc_atom = nuc_idx[0]
        for frame in range(traj.n_frames):
            frame_xyz = traj.xyz[frame]
            p_nuc = frame_xyz[nuc_atom]
            frame_dists = [np.linalg.norm(frame_xyz[c_idx] - p_nuc) * 10.0 for c_idx in scissile_carbons]
            distances.append(min(frame_dists))
        distances = np.array(distances)
    else:
        distances = np.ones(traj.n_frames) * 99.9
        
    # Compute stability metrics
    cutoff_static = config['filters'].get('catalytic_cutoff', 4.5)
    cutoff_md = config['filters'].get('catalytic_cutoff_md', 5.0)
    min_rmsd_threshold = config['filters'].get('min_stability_rmsd', 3.0)
    
    rmsd_stable_fraction = np.mean(rmsd < min_rmsd_threshold)
    catalytic_stable_fraction = np.mean(distances < cutoff_md)
    
    verdict = 'STABLE' if rmsd_stable_fraction > 0.5 and catalytic_stable_fraction > 0.5 else 'UNSTABLE'
    
    return {
        'rmsd': rmsd.tolist(),
        'distances': distances.tolist(),
        'avg_rmsd': float(np.mean(rmsd)),
        'avg_distance': float(np.mean(distances)),
        'rmsd_stable_fraction': float(rmsd_stable_fraction),
        'catalytic_stable_fraction': float(catalytic_stable_fraction),
        'verdict': verdict
    }
