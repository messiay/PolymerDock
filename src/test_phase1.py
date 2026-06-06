import sys
import os
import json
from rdkit import Chem

# Adjust path to find builder
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.builder import build_polymer, validate_input_structure, get_active_site_center

def test_phase1():
    print("=== STARTING PHASE 1 TEST ===")
    
    # 1. Test polymer generation
    smiles = '[C@@H](OC(=O)c1ccc(C(=O)O)cc1)'
    length = 6
    print(f"Building polymer from SMILES: {smiles} with length: {length}")
    try:
        polymer = build_polymer(smiles, length)
        print("Polymer build success!")
        print(f"Number of atoms: {polymer.GetNumAtoms()}")
        
        # Check conformers
        conf_count = polymer.GetNumConformers()
        print(f"Conformers generated: {conf_count}")
        assert conf_count == 1, "Conformer count is not 1"
        
        # Validate structure
        success, checks = validate_input_structure(polymer)
        print(f"Validation status: {success}")
        print(f"Validation checks detail: {checks}")
        assert success, "Valid polymer failed validation checks"
        
    except Exception as e:
        print(f"Error building/validating polymer: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    # 2. Test validation gate with bad structure (no conformer)
    print("\nTesting validation gate on bad structure (no conformer)...")
    bad_mol = Chem.MolFromSmiles(smiles)  # just 2D, no conformers embedded
    success, checks = validate_input_structure(bad_mol)
    print(f"Validation status for 2D mol: {success}")
    print(f"Validation checks detail: {checks}")
    assert not success, "Invalid (2D-only) molecule passed validation checks"
    
    # 3. Test active site auto-detection
    print("\nTesting active site auto-detection from PDB coordinates...")
    pdb_path = "data/5xjh.pdb"
    enzymes_path = "data/enzymes.json"
    
    if not os.path.exists(pdb_path):
        print(f"Error: PDB file not found at {pdb_path}")
        return False
    if not os.path.exists(enzymes_path):
        print(f"Error: Enzymes config not found at {enzymes_path}")
        return False
        
    with open(enzymes_path, 'r') as f:
        enzymes_db = json.load(f)
        
    petase_data = enzymes_db['PETase']
    catalytic_res = petase_data['catalytic_residues']
    print(f"PETase Catalytic residues: {catalytic_res}")
    
    try:
        center = get_active_site_center(pdb_path, catalytic_res)
        print(f"Computed Active Site Center: {center}")
        print("Phase 1 testing complete and PASSED!")
        return True
    except Exception as e:
        print(f"Error computing active site center: {e}")
        return False

if __name__ == '__main__':
    ok = test_phase1()
    sys.exit(0 if ok else 1)
