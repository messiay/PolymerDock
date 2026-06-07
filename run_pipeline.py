import argparse
import sys
import os
import shutil
import urllib.request
import json
import yaml
import numpy as np

# Add src/ to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'src')))

from src.builder import build_polymer, get_active_site_center, validate_input_structure
from src.docking import dock_anchor
from src.grower import grow_polymer
from src.scanner import scan_catalytic_viability
from src.scorer import score_binding
from src.validator import run_md_simulation, analyze_trajectory, OPENMM_AVAILABLE, OPENMM_ERROR
from src.utils import save_complex, setup_logging

def load_config():
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Config not found at {config_path}")
        sys.exit(1)
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_enzymes():
    enzymes_path = "data/enzymes.json"
    if not os.path.exists(enzymes_path):
        print(f"Error: Enzymes database not found at {enzymes_path}")
        sys.exit(1)
    with open(enzymes_path, 'r') as f:
        return json.load(f)

def generate_mock_anchor_poses(anchor_mol, center, num_poses=9):
    from rdkit import Chem
    import numpy as np
    poses = []
    conf = anchor_mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(anchor_mol.GetNumAtoms())])
    centroid = np.mean(coords, axis=0)
    translation = center - centroid
    
    # Base translated molecule
    base_mol = Chem.Mol(anchor_mol)
    b_conf = base_mol.GetConformer()
    for i in range(base_mol.GetNumAtoms()):
        pos = np.array(b_conf.GetAtomPosition(i))
        b_conf.SetAtomPosition(i, pos + translation)
        
    poses.append(base_mol)
    
    # Slightly perturbed versions for remaining modes
    r = np.random.RandomState(42)
    for p in range(1, num_poses):
        perturbed_mol = Chem.Mol(base_mol)
        p_conf = perturbed_mol.GetConformer()
        offset = r.normal(0.0, 0.8, size=3)
        for i in range(perturbed_mol.GetNumAtoms()):
            pos = np.array(p_conf.GetAtomPosition(i))
            p_conf.SetAtomPosition(i, pos + offset)
        poses.append(perturbed_mol)
    return poses

def main():
    parser = argparse.ArgumentParser(description="SimDock Polymer CLI Pipeline Driver")
    parser.add_argument("--enzyme", type=str, default="Amylase", help="Name of target enzyme (e.g. Amylase, PETase, Nylonase)")
    parser.add_argument("--smiles", type=str, default="OC1C(O)C(O)C(O)C(CO)O1", help="Monomer SMILES representation")
    parser.add_argument("--length", type=int, default=5, help="Target polymer chain length")
    parser.add_argument("--quick-test", type=str, default="True", help="True (500 steps for validation) or False (full production run)")
    parser.add_argument("--exhaustiveness", type=int, default=None, help="Vina docking exhaustiveness override")
    
    args = parser.parse_args()
    quick_test = args.quick_test.lower() in ("true", "1", "yes")
    
    print("\n" + "="*60)
    print("SimDock Polymer CLI Simulation Engine v2.0")
    print("="*60)
    
    config = load_config()
    setup_logging()
    enzymes_db = load_enzymes()
    
    if args.enzyme not in enzymes_db:
        print(f"Error: Enzyme '{args.enzyme}' not found in database. Available: {list(enzymes_db.keys())}")
        sys.exit(1)
        
    enzyme_data = enzymes_db[args.enzyme]
    print(f"\nTarget Enzyme: {args.enzyme} (PDB ID: {enzyme_data['pdb_id']})")
    print(f"Monomer SMILES: {args.smiles}")
    print(f"Chain Length: {args.length}-mer")
    
    # ===============================================================
    # PHASE 1: Structure Generation & Active Site Detection
    # ===============================================================
    print("\n--- PHASE 1: Structure Generation & Active Site Detection ---")
    scissile_type = enzyme_data.get('scissile_bond_type', 'ester_carbonyl')
    if 'ester' in scissile_type:
        linkage_type = 'ester'
    elif 'amide' in scissile_type:
        linkage_type = 'amide'
    elif 'glycosidic' in scissile_type:
        linkage_type = 'glycosidic'
    else:
        linkage_type = 'ester'
        
    print(f"Building full polymer ({linkage_type} linkage)...")
    polymer = build_polymer(args.smiles, args.length, config, linkage_type=linkage_type)
    success, checks = validate_input_structure(polymer)
    print(f"Validation checks: {checks}")
    if not success:
        print("Error: Input polymer failed validation.")
        sys.exit(1)
        
    # Download PDB if missing
    os.makedirs("data", exist_ok=True)
    pdb_file = f"data/{enzyme_data['pdb_id'].lower()}.pdb"
    if not os.path.exists(pdb_file):
        print(f"Enzyme structure {pdb_file} not found. Downloading from RCSB...")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{enzyme_data['pdb_id']}.pdb", pdb_file)
        print("Downloaded PDB successfully.")
        
    center = get_active_site_center(pdb_file, enzyme_data['catalytic_residues'])
    print(f"Active site center: [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}]")
    
    # ===============================================================
    # PHASE 2: Anchor Docking & Growth
    # ===============================================================
    print("\n--- PHASE 2: Anchor Docking & Growth Loop ---")
    anchor_length = config['docking'].get('anchor_length', 3)
    print(f"Building {anchor_length}-mer anchor fragment...")
    anchor_base = build_polymer(args.smiles, anchor_length, config, linkage_type=linkage_type)
    
    # Docking engine check
    engine_path = config['paths'].get('gnina_binary', 'gnina')
    has_docking_engine = (shutil.which(engine_path) is not None or 
                          shutil.which(config['paths'].get('vina_binary', 'vina')) is not None)
                          
    num_modes = config['docking'].get('num_modes', 9)
    if args.exhaustiveness:
        config['docking']['vina_exhaustiveness'] = args.exhaustiveness
        config['docking']['vina_exhaustiveness_colab'] = args.exhaustiveness
        
    if has_docking_engine:
        print(f"Real docking engine detected. Docking {anchor_length}-mer anchor pose candidates (num_modes={num_modes})...")
        anchor_poses = dock_anchor(anchor_base, pdb_file, center, config)
        if not anchor_poses:
            print("Error: Docking engine failed to yield poses.")
            sys.exit(1)
    else:
        print(f"No docking engine found on PATH. Generating {num_modes} mock anchor poses around active site center...")
        anchor_poses = generate_mock_anchor_poses(anchor_base, center, num_poses=num_modes)
        
    print(f"Successfully generated {len(anchor_poses)} anchor poses.")
    
    # Grow loop
    remaining_growth = args.length - anchor_length
    print(f"Growing {len(anchor_poses)} poses to {args.length}-mer (+{remaining_growth} monomers)...")
    grown_poses = []
    os.makedirs("results", exist_ok=True)
    
    for i, anchor_pose in enumerate(anchor_poses):
        print(f"  Pose {i+1}/{len(anchor_poses)}: Growing...", end="")
        if remaining_growth > 0:
            grown = grow_polymer(anchor_pose, args.smiles, args.length, pdb_file, config, linkage_type=linkage_type)
        else:
            grown = anchor_pose
            
        if grown is not None:
            print(f" OK: Growth succeeded ({grown.GetNumAtoms()} atoms)")
            grown_poses.append((i+1, grown))
        else:
            print(" FAIL: Growth failed (unresolvable clashes)")
            
    if not grown_poses:
        print("Error: All poses failed during growth. No candidates to filter.")
        sys.exit(1)
        
    print(f"Growth complete: {len(grown_poses)}/{len(anchor_poses)} poses survived.")
    
    # ===============================================================
    # PHASE 3: Catalytic Geometry Filter & Scoring
    # ===============================================================
    print("\n--- PHASE 3: Catalytic Geometry Filter & MM-GBSA Scoring ---")
    passing_poses = []
    
    for pose_num, grown_mol in grown_poses:
        ligand_pdb = f"results/grown_pose_{pose_num}.pdb"
        complex_pdb = f"results/complex_pose_{pose_num}.pdb"
        from rdkit import Chem
        Chem.MolToPDBFile(grown_mol, ligand_pdb)
        save_complex(pdb_file, grown_mol, complex_pdb)
        
        verdict, distance = scan_catalytic_viability(complex_pdb, enzyme_data, config)
        
        cutoff = config['filters'].get('catalytic_cutoff', 4.5)
        if 'catalytic_cutoff_override' in enzyme_data:
            cutoff = enzyme_data['catalytic_cutoff_override']
            
        if verdict == 'PASS':
            print(f"  Pose {pose_num}: PASS (distance {distance:.1f}A < {cutoff:.1f}A)")
            score_data = score_binding(complex_pdb, ligand_resname='UNL', config=config)
            print(f"            Score: {score_data['final_score']:.2f} kcal/mol (SASA: {score_data['buried_sasa']:.1f} A**2)")
            passing_poses.append({
                'pose_num': pose_num,
                'grown_mol': grown_mol,
                'ligand_pdb': ligand_pdb,
                'complex_pdb': complex_pdb,
                'distance': distance,
                'score_data': score_data
            })
        else:
            print(f"  Pose {pose_num}: REJECTED (distance {distance:.1f}A > {cutoff:.1f}A)")
            
    if not passing_poses:
        print("Error: No poses passed the catalytic geometry filter.")
        sys.exit(1)
        
    # Rank by binding energy
    passing_poses.sort(key=lambda p: p['score_data']['final_score'])
    print(f"\nRanked {len(passing_poses)} poses by binding energy:")
    for rank, p in enumerate(passing_poses):
        print(f"  Rank {rank+1}: Pose {p['pose_num']} -- {p['score_data']['final_score']:.2f} kcal/mol")
        
    # ===============================================================
    # PHASE 4: Molecular Dynamics Validation
    # ===============================================================
    print("\n--- PHASE 4: Molecular Dynamics Validation ---")
    if not OPENMM_AVAILABLE:
        print("[WARN] Warning: OpenMM or dependencies are not available in this Python interpreter.")
        print(f"   Detailed Import Error: {OPENMM_ERROR}")
        print("   Falling back to generating Mock DCD trajectories for verification.")
        
    # Validate top pose(s)
    n_to_validate = min(1, len(passing_poses))
    for rank in range(n_to_validate):
        p = passing_poses[rank]
        print(f"\nRunning MD simulation on Rank {rank+1} (Pose {p['pose_num']})...")
        traj_dcd, is_mock = run_md_simulation(p['complex_pdb'], p['ligand_pdb'], config, quick_test=quick_test)
        
        if is_mock:
            print("[WARN] OpenMM failed or was unavailable; generated MOCK trajectory.")
        else:
            print("[OK] Real OpenMM simulation completed successfully!")
            
        md_analysis = analyze_trajectory(traj_dcd, p['complex_pdb'], config, ligand_resname='UNL', enzyme_data=enzyme_data)
        print(f"MD Verdict: {md_analysis['verdict']}")
        print(f"  Average Ligand RMSD: {md_analysis['avg_rmsd']:.2f}A")
        print(f"  Average Catalytic Attack Distance: {md_analysis['avg_distance']:.2f}A")
        
        # Trajectory-averaged MM-GBSA
        if os.path.exists(traj_dcd):
            try:
                traj_score = score_binding(p['complex_pdb'], trajectory_dcd=traj_dcd, ligand_resname='UNL', config=config)
                print(f"  Trajectory-averaged Binding Score: {traj_score['final_score']:.2f} kcal/mol")
            except Exception as e:
                print(f"  Trajectory-averaged scoring failed: {e}")
                
    print("\n" + "="*60)
    print("SUCCESS: SimDock Pipeline execution finished successfully!")
    print("="*60 + "\n")

if __name__ == '__main__':
    main()
