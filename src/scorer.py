import mdtraj as md
import numpy as np
from Bio.PDB import PDBParser
import os

# Standard Amber vdW parameters (R_min/2 in Angstroms, epsilon in kcal/mol)
VDW_PARAMS = {
    'H': (1.0, 0.0157),
    'C': (1.9, 0.1094),
    'N': (1.8, 0.1700),
    'O': (1.6, 0.2100),
    'S': (2.0, 0.2500),
    'P': (2.1, 0.2000),
    'DEFAULT': (2.0, 0.1500)
}

# Standard amino acid atom partial charges (Amber99SB)
# We fall back to neutral if not matched, but we can match common backbone atoms:
BACKBONE_CHARGES = {
    'N': -0.4157,
    'H': 0.2719,
    'CA': 0.0337,
    'HA': 0.0823,
    'C': 0.5973,
    'O': -0.5679
}

def get_atom_params(atom_symbol):
    """
    Returns (R_min_half, epsilon) for an atom symbol.
    """
    sym = atom_symbol.upper()
    return VDW_PARAMS.get(sym, VDW_PARAMS['DEFAULT'])

def compute_interaction_energy(complex_pdb, ligand_resname='UNL'):
    """
    Computes Lennard-Jones and electrostatic interaction energy between the
    protein and the ligand.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', complex_pdb)
    model = structure[0]
    
    protein_atoms = []
    ligand_atoms = []
    
    for residue in model.get_residues():
        res_name = residue.get_resname()
        res_id = residue.get_id()
        is_ligand = res_id[0].startswith('H_') or res_name in [ligand_resname, 'LIG', 'UNK', 'POL']
        
        for atom in residue.get_atoms():
            # Estimate partial charge
            charge = 0.0
            if is_ligand:
                # Mock ligand charges (small dipole estimates, or neutral for vdW dominant)
                if atom.element == 'O': charge = -0.4
                elif atom.element == 'C' and any(n.element == 'O' for n in atom.get_parent().get_atoms()): charge = 0.4
                ligand_atoms.append((atom.get_vector().get_array(), atom.element, charge))
            else:
                # Protein backbone charges
                atom_name = atom.get_name()
                charge = BACKBONE_CHARGES.get(atom_name, 0.0)
                # Side chain charge estimates (acidic/basic)
                if res_name in ['ASP', 'GLU'] and atom_name in ['CG', 'CD', 'OE1', 'OE2', 'OD1', 'OD2']:
                    charge = -0.5
                elif res_name in ['LYS', 'ARG'] and atom_name in ['NZ', 'NH1', 'NH2', 'CZ']:
                    charge = 0.5
                protein_atoms.append((atom.get_vector().get_array(), atom.element, charge))
                
    if not ligand_atoms or not protein_atoms:
        return 0.0
        
    e_lj = 0.0
    e_el = 0.0
    raw_lj = 0.0
    raw_el = 0.0
    
    # Cutoff distance for speed (12 Angstroms)
    cutoff = 12.0
    dielectric = 4.0 # distance-dependent implicit dielectric constant
    
    for l_pos, l_elem, l_q in ligand_atoms:
        R_l, eps_l = get_atom_params(l_elem)
        for p_pos, p_elem, p_q in protein_atoms:
            diff = l_pos - p_pos
            r = np.linalg.norm(diff)
            
            if r < cutoff and r > 0.1:
                # 1. Lennard-Jones term
                R_p, eps_p = get_atom_params(p_elem)
                R_min = R_l + R_p
                eps_ij = np.sqrt(eps_l * eps_p)
                
                # U_LJ = eps * [ (R_min/r)^12 - 2 * (R_min/r)^6 ]
                ratio = R_min / r
                ratio6 = ratio ** 6
                ratio12 = ratio6 ** 2
                pairwise_lj = eps_ij * (ratio12 - 2 * ratio6)
                raw_lj += pairwise_lj
                
                # Cap the repulsion to handle unminimized steric clashes gracefully
                if pairwise_lj > 10.0:
                    pairwise_lj = 10.0
                e_lj += pairwise_lj
                
                # 2. Electrostatic term (Coulomb's Law with implicit solvent dielectric)
                # k_e = 332.0637 kcal/mol * A / e^2
                if abs(l_q) > 1e-4 and abs(p_q) > 1e-4:
                    raw_el += (332.0637 * l_q * p_q) / (dielectric * r)
                    
                    pairwise_el = (332.0637 * l_q * p_q) / (dielectric * max(r, 1.5))
                    # Cap electrostatic contribution per pair to prevent numerical explosion
                    pairwise_el = max(-10.0, min(10.0, pairwise_el))
                    e_el += pairwise_el
                    
    print(f"Scoring Complex {complex_pdb}:")
    print(f"  Uncapped Raw LJ: {raw_lj:.2f} kcal/mol | Capped LJ: {e_lj:.2f} kcal/mol")
    print(f"  Uncapped Raw Elec: {raw_el:.2f} kcal/mol | Capped Elec: {e_el:.2f} kcal/mol")
    return e_lj + e_el

def score_binding(complex_pdb, trajectory_dcd=None, ligand_resname='UNL'):
    """
    Computes a hybrid MM-GBSA score for the complex PDB or trajectory.
    Score = E_interaction + G_nonpolar (SASA buried term)
    """
    # 1. Load trajectory/pdb using MDTraj
    if trajectory_dcd and os.path.exists(trajectory_dcd):
        traj = md.load(trajectory_dcd, top=complex_pdb)
    else:
        traj = md.load(complex_pdb)
        
    # 2. Select ligand atoms
    ligand_indices = traj.topology.select(f"resname {ligand_resname} or resname LIG or resname POL")
    if len(ligand_indices) == 0:
        # Fallback to heteroatoms
        ligand_indices = traj.topology.select("not protein and not water")
        
    if len(ligand_indices) == 0:
        raise ValueError("Could not locate ligand in PDB structure for scoring.")
        
    # 3. Compute buried SASA (nonpolar solvation contribution)
    # Standard surface tension parameter gamma = -0.025 kcal/mol/A^2
    gamma = -0.025
    
    # Calculate SASA for all atoms in the complex (frame 0 for single PDB)
    sasa_complex = md.shrake_rupley(traj, probe_radius=0.14, mode='atom')
    
    # Slice ligand to calculate its unbound SASA
    lig_traj = traj.atom_slice(ligand_indices)
    sasa_ligand = md.shrake_rupley(lig_traj, probe_radius=0.14, mode='atom')
    
    # Average over all trajectory frames
    avg_buried_sasa = 0.0
    num_frames = traj.n_frames
    for f in range(num_frames):
        s_complex_f = sasa_complex[f]
        s_lig_f = sasa_ligand[f]
        buried_f = np.sum(s_lig_f) - np.sum(s_complex_f[ligand_indices])
        avg_buried_sasa += buried_f
    avg_buried_sasa /= num_frames
    
    g_nonpolar = gamma * avg_buried_sasa
    
    # 4. Compute Lennard-Jones + electrostatic interaction energy
    # If a trajectory is provided, we compute the average over frames
    avg_e_int = 0.0
    if num_frames > 1 and trajectory_dcd:
        # For efficiency, we sample up to 10 frames from the trajectory
        step = max(1, num_frames // 10)
        frames_sampled = 0
        for f in range(0, num_frames, step):
            # Save frame to temporary PDB file
            temp_frame_pdb = f"results/temp_frame_{f}.pdb"
            traj[f].save_pdb(temp_frame_pdb)
            try:
                avg_e_int += compute_interaction_energy(temp_frame_pdb, ligand_resname)
            finally:
                if os.path.exists(temp_frame_pdb):
                    os.remove(temp_frame_pdb)
            frames_sampled += 1
        avg_e_int /= frames_sampled
    else:
        avg_e_int = compute_interaction_energy(complex_pdb, ligand_resname)
        
    # Final score (lower is better binding)
    final_score = avg_e_int + g_nonpolar
    return {
        'interaction_energy': avg_e_int,
        'nonpolar_solvation': g_nonpolar,
        'buried_sasa': avg_buried_sasa,
        'final_score': final_score
    }
