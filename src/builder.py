import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser
import os

def build_polymer(smiles_string, chain_length, config=None):
    """
    Generates a chemically valid 3D polymer molecule from a monomer SMILES string
    by chaining them together through ester-like single bonds.
    """
    monomer = Chem.MolFromSmiles(smiles_string)
    if monomer is None:
        raise ValueError(f"Invalid monomer SMILES string: {smiles_string}")
        
    # Find linking atoms
    # 1. Acid oxygen tail (part of C(=O)[OH])
    acid_pat = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    matches = monomer.GetSubstructMatches(acid_pat)
    if not matches:
        tail_idx = monomer.GetNumAtoms() - 1
    else:
        tail_idx = matches[0][2]
        
    # 2. Radical carbon head (carbon with radical electrons)
    head_idx = None
    for atom in monomer.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetNumRadicalElectrons() > 0:
            head_idx = atom.GetIdx()
            break
    if head_idx is None:
        head_idx = 0
        
    # Build polymer using RWMol
    rw_mol = Chem.RWMol(monomer)
    num_atoms = monomer.GetNumAtoms()
    active_tail = tail_idx
    
    for i in range(1, chain_length):
        offset = rw_mol.GetNumAtoms()
        rw_mol.InsertMol(monomer)
        new_head = head_idx + offset
        rw_mol.AddBond(active_tail, new_head, Chem.BondType.SINGLE)
        active_tail = tail_idx + offset
        
    # Adjust radical electrons for linked units
    for i in range(1, chain_length):
        h_idx = head_idx + i * num_atoms
        rw_mol.GetAtomWithIdx(h_idx).SetNumRadicalElectrons(0)
        
    for i in range(chain_length - 1):
        t_idx = tail_idx + i * num_atoms
        rw_mol.GetAtomWithIdx(t_idx).SetNumRadicalElectrons(0)
        
    polymer = rw_mol.GetMol()
    Chem.SanitizeMol(polymer)
    polymer = Chem.AddHs(polymer)
    
    # Embed molecule in 3D
    res = AllChem.EmbedMolecule(polymer, randomSeed=42, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
    if res == -1:
        res = AllChem.EmbedMolecule(polymer, randomSeed=42)
    if res == -1:
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        res = AllChem.EmbedMolecule(polymer, params)
        
    if res == -1:
        raise RuntimeError("Failed to generate 3D conformer for polymer.")
        
    # Optimize structure using UFF force field
    AllChem.UFFOptimizeMolecule(polymer)
    
    return polymer

def validate_input_structure(mol):
    """
    Validates that the molecule has 3D coordinates, passes chemical sanitization,
    and has no major self-clashing non-bonded atoms.
    """
    checks = {
        'has_3d_coords': mol.GetNumConformers() > 0,
        'correct_atom_count': mol.GetNumAtoms() > 0,
        'passes_sanitization': True,
        'no_clashing_atoms': True
    }
    
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        checks['passes_sanitization'] = False
        
    if checks['has_3d_coords']:
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        clash_threshold = 0.8  # angstroms
        
        # Build bond connectivity map to avoid checking bonded pairs
        bonded = {}
        for bond in mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonded[(i, j)] = True
            bonded[(j, i)] = True
            
        for i in range(num_atoms):
            pos_i = conf.GetAtomPosition(i)
            for j in range(i + 1, num_atoms):
                if not bonded.get((i, j), False):
                    dist = (pos_i - conf.GetAtomPosition(j)).Length()
                    if dist < clash_threshold:
                        checks['no_clashing_atoms'] = False
                        break
            if not checks['no_clashing_atoms']:
                break
                
    success = all(checks.values())
    return success, checks

def get_active_site_center(pdb_file, catalytic_residues):
    """
    Detects active site center coordinates dynamically by averaging CA atoms 
    of specified catalytic residues in the enzyme PDB.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('enzyme', pdb_file)
    coords = []
    
    model = structure[0]
    # Default to chain A, or pick the first available chain
    if 'A' in model:
        chain = model['A']
    else:
        chain = list(model.get_chains())[0]
        
    for res_name, res_num in catalytic_residues.items():
        if res_num in chain:
            residue = chain[res_num]
            if 'CA' in residue:
                coords.append(residue['CA'].get_vector().get_array())
            else:
                res_coords = [atom.get_vector().get_array() for atom in residue]
                if res_coords:
                    coords.append(np.mean(res_coords, axis=0))
        else:
            found = False
            for c in model.get_chains():
                if res_num in c:
                    residue = c[res_num]
                    if 'CA' in residue:
                        coords.append(residue['CA'].get_vector().get_array())
                        found = True
                        break
            if not found:
                raise ValueError(f"Residue {res_name} {res_num} not found in the structure.")
                
    if not coords:
        raise ValueError("No coordinates found for catalytic residues.")
    return np.mean(coords, axis=0)
