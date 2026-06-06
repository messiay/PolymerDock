import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBParser
import os
import logging

logger = logging.getLogger('simdock')

def get_enzyme_coords(enzyme_pdb):
    """
    Extracts all atom coordinates from the enzyme PDB to use for clash checking.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('enzyme', enzyme_pdb)
    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.get_vector().get_array())
    return np.array(coords)

def find_heavy_linking_atoms(mol):
    """
    Finds the head carbon and tail oxygen in a heavy-atom-only molecule.
    """
    # Find carboxylic acid oxygen (single-bonded to carbonyl C, and has degree 1)
    pat = Chem.MolFromSmarts('[CX3](=O)[OX2]')
    matches = mol.GetSubstructMatches(pat)
    tail_idx = None
    for match in matches:
        o_idx = match[2]
        if mol.GetAtomWithIdx(o_idx).GetDegree() == 1:
            tail_idx = o_idx
            break
    if tail_idx is None:
        tail_idx = mol.GetNumAtoms() - 1
        
    # Find head: Carbon with radical electrons (or fallback to first Carbon)
    head_idx = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetNumRadicalElectrons() > 0:
            head_idx = atom.GetIdx()
            break
            
    return head_idx, tail_idx

def find_active_tail(mol):
    """
    Dynamically finds the active tail oxygen (degree 1 carboxylic acid oxygen)
    in the growing polymer chain.
    """
    pat = Chem.MolFromSmarts('[CX3](=O)[OX2]')
    matches = mol.GetSubstructMatches(pat)
    for match in matches:
        o_idx = match[2]
        if mol.GetAtomWithIdx(o_idx).GetDegree() == 1:
            return o_idx
    # Fallback to the last atom
    return mol.GetNumAtoms() - 1

def attach_at_angle(current_mol, monomer_mol, angle, tail_idx, head_idx):
    """
    Attaches monomer_mol to current_mol at a specified torsion angle around the new bond,
    and performs a local UFF minimization to relax connection bond lengths and angles.
    """
    current = Chem.RWMol(current_mol)
    monomer = Chem.Mol(monomer_mol)
    
    n1 = current.GetNumAtoms()
    n2 = monomer.GetNumAtoms()
    
    # 1. Calculate positions for alignment
    conf1 = current.GetConformer()
    conf2 = monomer.GetConformer()
    
    P_tail = np.array(conf1.GetAtomPosition(tail_idx))
    P_head = np.array(conf2.GetAtomPosition(head_idx))
    
    # Find the atom bonded to the tail in current (carbonyl carbon)
    tail_atom = current.GetAtomWithIdx(tail_idx)
    neighbors = tail_atom.GetNeighbors()
    if not neighbors:
        raise ValueError("Tail atom has no neighbors to determine bond direction.")
    P_carbonyl = np.array(conf1.GetAtomPosition(neighbors[0].GetIdx()))
    
    # Direction vector extending from carbonyl C through tail O
    v = P_tail - P_carbonyl
    norm_v = np.linalg.norm(v)
    if norm_v < 1e-5:
        u = np.array([1.0, 0.0, 0.0])
    else:
        u = v / norm_v
        
    # Standard ester O-C bond length is ~1.43 A
    d = 1.43
    P_head_new = P_tail + d * u
    
    # 2. Translate monomer to align head with P_head_new
    translation = P_head_new - P_head
    for i in range(n2):
        pos = np.array(conf2.GetAtomPosition(i))
        new_pos = pos + translation
        conf2.SetAtomPosition(i, new_pos)
        
    # 3. Rotate monomer around the new bond axis (u) passing through P_tail by the specified angle
    theta = np.radians(angle)
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    ux, uy, uz = u
    
    cross_matrix = np.array([
        [0.0, -uz, uy],
        [uz, 0.0, -ux],
        [-uy, ux, 0.0]
    ])
    outer_matrix = np.outer(u, u)
    R = cos_t * np.eye(3) + sin_t * cross_matrix + (1 - cos_t) * outer_matrix
    
    for i in range(n2):
        pos = np.array(conf2.GetAtomPosition(i))
        pos_rel = pos - P_tail
        pos_rot = R.dot(pos_rel)
        new_pos = pos_rot + P_tail
        conf2.SetAtomPosition(i, new_pos)
        
    # 4. Combine molecules and form bond
    current.InsertMol(monomer)
    new_head_idx = head_idx + n1
    
    current.AddBond(tail_idx, new_head_idx, Chem.BondType.SINGLE)
    
    # Clean up radical electrons and set no implicit Hs on the connected oxygen
    current.GetAtomWithIdx(tail_idx).SetNumRadicalElectrons(0)
    current.GetAtomWithIdx(tail_idx).SetNoImplicit(True)
    current.GetAtomWithIdx(new_head_idx).SetNumRadicalElectrons(0)
    
    candidate = current.GetMol()
    
    # 5. Local UFF minimization to relax the newly formed bond geometry
    try:
        mol_opt = Chem.Mol(candidate)
        mol_opt.UpdatePropertyCache(strict=False)
        ff = AllChem.UFFGetMoleculeForceField(mol_opt)
        if ff:
            # Freeze the entire parent polymer (indices 0 to n1-1)
            for i in range(n1):
                ff.AddFixedPoint(i)
            # Minimize the newly attached monomer coordinates
            ff.Minimize(maxIts=300)
        candidate = mol_opt
    except Exception:
        return None
        
    try:
        Chem.SanitizeMol(candidate, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS)
    except Exception:
        return None
        
    return candidate

def check_clash(candidate_mol, new_atoms_range, enzyme_coords, config):
    """
    Checks if any of the newly added monomer atoms clash with either the enzyme
    or the rest of the growing polymer chain.
    """
    clash_threshold_enzyme = config['growth'].get('clash_threshold', 0.5)
    clash_threshold_self = 0.8 # internal polymer clash limit
    
    conf = candidate_mol.GetConformer()
    new_coords = np.array([conf.GetAtomPosition(i) for i in new_atoms_range])
    
    # Check clashes with enzyme
    if len(enzyme_coords) > 0:
        dists_enzyme = np.linalg.norm(new_coords[:, np.newaxis, :] - enzyme_coords[np.newaxis, :, :], axis=2)
        if np.any(dists_enzyme < clash_threshold_enzyme):
            return True
            
    # Check clashes with the rest of the polymer
    old_atoms_count = new_atoms_range[0]
    if old_atoms_count > 0:
        old_coords = np.array([conf.GetAtomPosition(i) for i in range(old_atoms_count)])
        
        bonded = {}
        for bond in candidate_mol.GetBonds():
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bonded[(i, j)] = True
            bonded[(j, i)] = True
            
        for idx_new, i in enumerate(new_atoms_range):
            for j in range(old_atoms_count):
                if not bonded.get((i, j), False):
                    dist = np.linalg.norm(new_coords[idx_new] - old_coords[j])
                    if dist < clash_threshold_self:
                        return True
                        
    return False

def score_energy(candidate_mol):
    """
    Computes a quick energy score for the candidate conformer using UFF force field.
    """
    try:
        mol_copy = Chem.Mol(candidate_mol)
        mol_copy.UpdatePropertyCache(strict=False)
        ff = AllChem.UFFGetMoleculeForceField(mol_copy)
        if ff:
            return ff.CalcEnergy()
    except Exception:
        pass
    return float('inf')

def optimize_hydrogens(mol):
    """
    Optimizes only the hydrogen atom positions using UFF while freezing heavy atoms.
    """
    try:
        ff = AllChem.UFFGetMoleculeForceField(mol)
        if ff:
            for i in range(mol.GetNumAtoms()):
                if mol.GetAtomWithIdx(i).GetAtomicNum() > 1: # Heavy atom
                    ff.AddFixedPoint(i)
            ff.Minimize(maxIts=500)
    except Exception:
        pass

def grow_polymer(anchor_pose, monomer_smiles, chain_length, enzyme_pdb, config):
    """
    Sequentially grows the anchor pose to the target chain length using 
    torsion-sampling and energy scoring.
    """
    # 1. Strip hydrogens from the anchor pose to work with heavy atoms
    current = Chem.RemoveHs(anchor_pose)
    num_anchor_heavy = current.GetNumAtoms()
    
    # 2. Prepare 3D monomer template (heavy atoms only)
    monomer_unit = Chem.MolFromSmiles(monomer_smiles)
    monomer_3d = Chem.Mol(monomer_unit)
    AllChem.EmbedMolecule(monomer_3d, randomSeed=42)
    AllChem.UFFOptimizeMolecule(monomer_3d)
    monomer_3d = Chem.RemoveHs(monomer_3d)
    
    # Find linking indices in the heavy monomer template
    monomer_head_idx, monomer_tail_idx = find_heavy_linking_atoms(monomer_3d)
    num_monomer_atoms = monomer_3d.GetNumAtoms()
    
    # Load enzyme coordinates
    enzyme_coords = get_enzyme_coords(enzyme_pdb)
    
    # 3. Growth loop
    for step in range(1, chain_length):
        # Dynamically find the active tail oxygen
        active_tail = find_active_tail(current)
        
        candidates = []
        num_samples = config['growth'].get('rotation_samples', 36)
        angles = np.linspace(0, 360, num_samples, endpoint=False)
        
        n_current = current.GetNumAtoms()
        new_atoms_range = range(n_current, n_current + num_monomer_atoms)
        
        for angle in angles:
            candidate = attach_at_angle(current, monomer_3d, angle, active_tail, monomer_head_idx)
            if candidate is not None:
                if not check_clash(candidate, new_atoms_range, enzyme_coords, config):
                    energy = score_energy(candidate)
                    candidates.append((candidate, energy))
                    
        if not candidates:
            logger.warning(f"Growth failed at step {step}: No clash-free conformations found.")
            return None
            
        current = min(candidates, key=lambda x: x[1])[0]
        
    # 4. Final sanitization and hydrogen restoration
    Chem.SanitizeMol(current)
    current = Chem.AddHs(current, addCoords=True)
    # Relax only the hydrogens to preserve coordinates of the clash-free heavy skeleton
    optimize_hydrogens(current)
    
    # Final hydrogen clash check against enzyme (on grown monomers only)
    if not validate_final_hydrogens(current, enzyme_coords, num_anchor_heavy, threshold=1.2):
        logger.warning("Final hydrogen clash check failed: optimized polymer hydrogens clash with the enzyme.")
        return None
        
    # 5. Compute Gasteiger partial charges for downstream scoring/docking
    try:
        AllChem.ComputeGasteigerCharges(current)
    except Exception:
        pass  # Non-fatal — charges are a bonus for scoring accuracy
    
    return current

def validate_final_hydrogens(mol, enzyme_coords, num_anchor_heavy, threshold=1.2):
    conf = mol.GetConformer()
    h_indices = []
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 1:
            neighbors = atom.GetNeighbors()
            if neighbors:
                parent_idx = neighbors[0].GetIdx()
                if parent_idx >= num_anchor_heavy:
                    h_indices.append(i)
                    
    if len(enzyme_coords) > 0 and len(h_indices) > 0:
        h_coords = np.array([conf.GetAtomPosition(i) for i in h_indices])
        dists = np.linalg.norm(
            h_coords[:, np.newaxis, :] - enzyme_coords[np.newaxis, :, :], axis=2)
        return not np.any(dists < threshold)
    return True
