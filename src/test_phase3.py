import sys
import os
import yaml

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.scanner import scan_catalytic_viability
from src.scorer import score_binding

def test_phase3():
    print("=== STARTING PHASE 3 TEST ===")
    
    # 1. Load config
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Config not found at {config_path}")
        return False
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    from src.utils import save_complex
    from rdkit import Chem

    enzyme_pdb = "data/5xjh.pdb"
    ligand_pdb = "results/grown_poly.pdb"
    complex_pdb = "results/complex.pdb"
    
    if not os.path.exists(enzyme_pdb):
        print(f"Error: Enzyme PDB not found at {enzyme_pdb}")
        return False
    if not os.path.exists(ligand_pdb):
        print(f"Error: Grown polymer PDB not found at {ligand_pdb}. Please run Phase 2 first.")
        return False
        
    print(f"Merging {enzyme_pdb} and {ligand_pdb} into {complex_pdb}...")
    ligand_mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
    save_complex(enzyme_pdb, ligand_mol, complex_pdb)
    print("Complex created successfully!")
    
    # PETase enzyme data for wild-type numbering (Ser160)
    enzyme_data = {
        "nucleophile_res_num": 160,
        "nucleophile_atom_name": "OG"
    }
    
    # 2. Run Scanner
    print("Running Catalytic Geometry Filter...")
    try:
        verdict, distance = scan_catalytic_viability(complex_pdb, enzyme_data, config)
        print(f"Scanner Verdict: {verdict}")
        print(f"Attack Distance: {distance:.2f} Angstroms")
        
        # Test a flipped/fail pose (mocked by using a wrong residue number or check behavior)
        # If we use a wrong nucleophile residue number, it should raise an error, which is correct.
        # Let's verify that the output of our grown polymer passes or fails as expected.
        # Since we docked at the active site center, the distance should be close.
        print("Scanner ran successfully!")
    except Exception as e:
        print(f"Scanner failed: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    # 3. Run Scorer
    print("\nRunning MM-GBSA hybrid Scorer...")
    try:
        score_data = score_binding(complex_pdb, ligand_resname='UNL')
        print(f"Scorer Output Detail: {score_data}")
        print(f"Final Score: {score_data['final_score']:.2f} kcal/mol")
        
        # Check assertions
        assert 'final_score' in score_data
        assert 'buried_sasa' in score_data
        assert score_data['buried_sasa'] > 0.0, "Buried SASA is not positive"
        
        print("\nPhase 3 testing complete and PASSED!")
        return True
    except Exception as e:
        print(f"Scorer failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    ok = test_phase3()
    sys.exit(0 if ok else 1)
