from Bio.PDB import PDBParser
import os
import yaml

def scan_catalytic_viability(complex_pdb, enzyme_data, config):
    """
    Checks if the scissile carbonyl carbon of the polymer ligand is within 
    catalytic attack distance of the enzyme nucleophile (e.g. Ser160 OG).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', complex_pdb)
    model = structure[0]
    
    # Get nucleophile parameters
    nuc_res_num = enzyme_data.get("nucleophile_res_num", 160)
    nuc_atom_name = enzyme_data.get("nucleophile_atom_name", "OG")
    
    # 1. Locate the nucleophile atom
    nucleophile_atom = None
    # Look in Chain A first, or fallback to any chain
    for chain in model:
        if nuc_res_num in chain:
            res = chain[nuc_res_num]
            if nuc_atom_name in res:
                nucleophile_atom = res[nuc_atom_name]
                break
                
    if nucleophile_atom is None:
        # Detailed search across all chains/residues
        for residue in model.get_residues():
            if residue.id[1] == nuc_res_num:
                if nuc_atom_name in residue:
                    nucleophile_atom = residue[nuc_atom_name]
                    break
                    
    if nucleophile_atom is None:
        raise ValueError(f"Nucleophile atom ({nuc_atom_name} of residue {nuc_res_num}) not found in structure.")
        
    P_nuc = nucleophile_atom.get_vector()
    
    # 2. Extract ligand atoms (residues with non-blank hetero flags or names like UNL/LIG/POL)
    ligand_atoms = []
    for residue in model.get_residues():
        res_name = residue.get_resname()
        res_id = residue.get_id()
        # RDKit writes PDB blocks with residue name UNL by default.
        if res_id[0].startswith('H_') or res_name in ['UNL', 'LIG', 'UNK', 'POL']:
            for atom in residue.get_atoms():
                ligand_atoms.append(atom)
                
    if not ligand_atoms:
        raise ValueError("No ligand atoms (residue UNL/LIG/UNK/POL or HETATM) found in complex.")
        
    # 3. Detect carbonyl carbons in the ligand geometrically (C=O distance 1.15 to 1.30 A)
    carbons = [atom for atom in ligand_atoms if atom.element == 'C']
    oxygens = [atom for atom in ligand_atoms if atom.element == 'O']
    
    carbonyl_carbons = []
    for c_atom in carbons:
        for o_atom in oxygens:
            dist = (c_atom.get_vector() - o_atom.get_vector()).norm()
            if 1.15 <= dist <= 1.30:
                carbonyl_carbons.append(c_atom)
                break
                
    if not carbonyl_carbons:
        raise ValueError("No carbonyl carbons detected in the polymer ligand.")
        
    # 4. Measure minimum attack distance
    min_dist = float('inf')
    best_c_atom = None
    for c_atom in carbonyl_carbons:
        dist = (c_atom.get_vector() - P_nuc).norm()
        if dist < min_dist:
            min_dist = dist
            best_c_atom = c_atom
            
    cutoff = config['filters'].get('catalytic_cutoff', 4.5)
    verdict = 'PASS' if min_dist < cutoff else 'FAIL'
    
    return verdict, min_dist
