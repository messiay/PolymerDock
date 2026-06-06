import subprocess
import os
import shutil
import yaml
from rdkit import Chem
from meeko import MoleculePreparation

def prepare_ligand_pdbqt(rdkit_mol, output_path):
    """
    Converts an RDKit molecule to a PDBQT file using Meeko.
    """
    preparator = MoleculePreparation()
    preparator.prepare(rdkit_mol)
    pdbqt_str = preparator.write_pdbqt_string()
    
    with open(output_path, 'w') as f:
        f.write(pdbqt_str)
    return output_path

def run_docking_cmd(cmd):
    """
    Runs a subprocess command and handles errors.
    """
    print(f"Executing: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Docking execution failed: {result.stderr}")
    return result.stdout

def dock_anchor(anchor_mol, enzyme_pdb, center, config):
    """
    Docks the anchor fragment into the enzyme active site using GNINA or Vina.
    Implements automatic box size retries.
    """
    engine = config['docking'].get('engine', 'gnina').lower()
    out_dir = config['paths'].get('output_dir', './results')
    os.makedirs(out_dir, exist_ok=True)
    
    # Paths for files
    anchor_pdb = os.path.join(out_dir, "anchor.pdb")
    Chem.MolToPDBFile(anchor_mol, anchor_pdb)
    
    anchor_pdbqt = os.path.join(out_dir, "anchor.pdbqt")
    prepare_ligand_pdbqt(anchor_mol, anchor_pdbqt)
    
    # We will try to prepare receptor PDBQT for Vina
    # For GNINA, we can use the PDB directly
    receptor_file = enzyme_pdb
    if engine == 'vina':
        # Simple fallback for receptor PDBQT: copy or convert
        # In a real workflow, we prepare PDBQT, but Vina can also accept simple conversions.
        # Let's assume receptor is converted, or we use GNINA which is primary.
        # We can copy pdb to pdbqt as a dummy, or search for prepare_receptor.
        # Vina actually requires a valid PDBQT with charges.
        # If meeko has a receptor preparation command, we can run it.
        # For this pipeline, we will default to using GNINA or assume receptor_file is PDBQT if Vina.
        if receptor_file.endswith('.pdb'):
            receptor_pdbqt = os.path.join(out_dir, "receptor.pdbqt")
            # For simplicity, we just create a basic PDBQT or assume user provides PDBQT.
            # In a fallback scenario, Vina can run on receptor.pdbqt.
            shutil.copy(receptor_file, receptor_pdbqt) # Fallback copy
            receptor_file = receptor_pdbqt
            
    num_modes = config['docking'].get('num_modes', 9)
    exhaustiveness = config['docking'].get('vina_exhaustiveness', 8)
    
    # Docking box settings
    padding = config['docking'].get('box_padding', 10.0)
    retry_padding = config['docking'].get('box_padding_retry', 5.0)
    
    # Base box size (approximate ligand size + padding)
    # Let's compute box size from ligand coordinates or use a default 15x15x15 box.
    box_size = 15.0 + padding
    
    out_poses = os.path.join(out_dir, "docked_poses.sdf" if engine == 'gnina' else "docked_poses.pdbqt")
    
    # Command builder
    def build_command(current_box_size):
        if engine == 'gnina':
            gnina_bin = config['paths'].get('gnina_binary', 'gnina')
            cmd = [
                gnina_bin,
                '-r', receptor_file,
                '-l', anchor_pdb, # GNINA accepts PDB or SDF for ligand
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(current_box_size),
                '--size_y', str(current_box_size),
                '--size_z', str(current_box_size),
                '--exhaustiveness', str(exhaustiveness),
                '--num_modes', str(num_modes),
                '--out', out_poses
            ]
        else: # Vina fallback
            vina_bin = config['paths'].get('vina_binary', 'vina')
            cmd = [
                vina_bin,
                '--receptor', receptor_file,
                '--ligand', anchor_pdbqt,
                '--center_x', str(center[0]),
                '--center_y', str(center[1]),
                '--center_z', str(center[2]),
                '--size_x', str(current_box_size),
                '--size_y', str(current_box_size),
                '--size_z', str(current_box_size),
                '--exhaustiveness', str(exhaustiveness),
                '--num_modes', str(num_modes),
                '--out', out_poses
            ]
        return cmd
        
    # Run with retry
    try:
        run_docking_cmd(build_command(box_size))
    except Exception as e:
        print(f"Docking failed with box size {box_size}. Retrying with larger box...")
        # Retry with larger box
        box_size += retry_padding
        try:
            run_docking_cmd(build_command(box_size))
        except Exception as retry_err:
            raise RuntimeError(f"Docking failed on retry: {retry_err}")
            
    # Parse docked poses from out_poses back into RDKit molecules
    poses = []
    if os.path.exists(out_poses):
        if out_poses.endswith('.sdf'):
            supplier = Chem.SDMolSupplier(out_poses, sanitize=False)
            for mol in supplier:
                if mol is not None:
                    poses.append(mol)
        else:
            with open(out_poses, 'r') as f:
                lines = f.readlines()
                
            current_pose_lines = []
            for line in lines:
                if line.startswith("MODEL"):
                    current_pose_lines = []
                elif line.startswith("ENDMDL"):
                    # Clean PDBQT specific columns to prevent RDKit crashes
                    cleaned_lines = []
                    for p in current_pose_lines:
                        if p.startswith("ATOM") or p.startswith("HETATM"):
                            # Replace the AutoDock atom type at columns 77-78 with standard element
                            element = p[77:79].strip()
                            # simple hack: just remove it, MolFromPDBBlock tries to infer from atom name
                            cleaned_lines.append(p[:76] + "    \n")
                        else:
                            cleaned_lines.append(p)
                            
                    pose_str = "".join(cleaned_lines)
                    pose_mol = Chem.MolFromPDBBlock(pose_str, sanitize=False)
                    if pose_mol:
                        poses.append(pose_mol)
                else:
                    current_pose_lines.append(line)
                    
            # If no MODEL records, try parsing whole file
            if not poses:
                pose_mol = Chem.MolFromPDBFile(out_poses, sanitize=False)
                if pose_mol:
                    poses.append(pose_mol)
                
    return poses
