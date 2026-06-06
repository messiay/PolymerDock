from rdkit import Chem
import os

def save_complex(enzyme_pdb, ligand_mol, output_pdb):
    """
    Combines the enzyme PDB file and the RDKit ligand molecule into a single 
    complex PDB file, stripping redundant END/CONECT records.
    """
    if not os.path.exists(enzyme_pdb):
        raise FileNotFoundError(f"Enzyme PDB file not found at {enzyme_pdb}")
        
    with open(enzyme_pdb, 'r') as f:
        enzyme_lines = f.readlines()
        
    # Remove lines with END or CONECT to prevent PDB structure errors
    clean_enzyme_lines = []
    for line in enzyme_lines:
        if line.startswith("END") or line.startswith("CONECT"):
            continue
        clean_enzyme_lines.append(line)
        
    ligand_pdb_block = Chem.MolToPDBBlock(ligand_mol)
    
    # Write merged files
    with open(output_pdb, 'w') as f:
        f.writelines(clean_enzyme_lines)
        f.write(ligand_pdb_block)
        f.write("END\n")
        
    return output_pdb
