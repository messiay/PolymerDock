import sys
import os
import yaml
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.builder import build_polymer, get_active_site_center
from src.grower import grow_polymer

def test_phase2():
    print("=== STARTING PHASE 2 TEST ===")
    
    # 1. Load config and data
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Config not found at {config_path}")
        return False
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    enzyme_pdb = "data/5xjh.pdb"
    if not os.path.exists(enzyme_pdb):
        print(f"Error: Enzyme PDB not found at {enzyme_pdb}")
        return False
        
    # Get active site center
    catalytic_res = {'SER': 160, 'HIS': 237, 'ASP': 206}
    center = get_active_site_center(enzyme_pdb, catalytic_res)
    print(f"Active site center for growth: {center}")
    
    # 2. Generate a mock anchor pose
    smiles = '[C@@H](OC(=O)c1ccc(C(=O)O)cc1)'
    print("Generating mock anchor (1-mer)...")
    anchor = build_polymer(smiles, 1)
    
    # Translate anchor to the active site center to mimic docking
    conf = anchor.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(anchor.GetNumAtoms())])
    centroid = np.mean(coords, axis=0)
    translation = center - centroid
    
    for i in range(anchor.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i))
        conf.SetAtomPosition(i, pos + translation)
        
    print(f"Mock anchor translated. Centroid now at: {np.mean([conf.GetAtomPosition(i) for i in range(anchor.GetNumAtoms())], axis=0)}")
    
    # Save mock anchor to results/mock_anchor.pdb
    os.makedirs("results", exist_ok=True)
    Chem.MolToPDBFile(anchor, "results/mock_anchor.pdb")
    print("Saved mock anchor to results/mock_anchor.pdb")
    
    # 3. Grow polymer from length 1 to 6
    print("\nRunning polymer growth loop to 6-mer...")
    try:
        grown_poly = grow_polymer(anchor, smiles, 6, enzyme_pdb, config)
        if grown_poly is None:
            print("Error: Growth loop returned None (all poses clashed).")
            return False
            
        print("Growth loop completed successfully!")
        print(f"Grown polymer atom count: {grown_poly.GetNumAtoms()}")
        
        # Save output to results/grown_poly.pdb
        Chem.MolToPDBFile(grown_poly, "results/grown_poly.pdb")
        print("Saved grown polymer to results/grown_poly.pdb")
        
        # Check assertions
        assert grown_poly.GetNumAtoms() > anchor.GetNumAtoms(), "Grown polymer has fewer or equal atoms than anchor"
        assert grown_poly.GetNumConformers() == 1, "Grown polymer has incorrect conformer count"
        
        # Chemical sanitization check
        Chem.SanitizeMol(grown_poly)
        print("Sanitization check passed!")
        
        print("Phase 2 testing complete and PASSED!")
        return True
    except Exception as e:
        print(f"Error during growth test: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    ok = test_phase2()
    sys.exit(0 if ok else 1)
