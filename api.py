import os
import uuid
import threading
import shutil
import urllib.request
import json
import yaml
import logging
from fastapi import FastAPI, BackgroundTasks, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel

# Import SimDock modules
from src.builder import build_polymer, get_active_site_center, validate_input_structure
from src.docking import dock_anchor
from src.grower import grow_polymer
from src.scanner import scan_catalytic_viability
from src.scorer import score_binding
from src.validator import run_md_simulation, analyze_trajectory, OPENMM_AVAILABLE, OPENMM_ERROR
from src.utils import save_complex

app = FastAPI(title="SimDock Polymer v2.0 API")

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# In-memory task store
SIMULATION_TASKS = {}
SIMULATION_TASKS_LOCK = threading.Lock()

# Custom logger for simdock
logger = logging.getLogger('simdock')
logger.setLevel(logging.INFO)

# Load configuration and enzymes database helper functions
def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def load_enzymes():
    with open("data/enzymes.json", "r") as f:
        return json.load(f)

class SimulationRequest(BaseModel):
    enzyme: str
    smiles: str
    length: int
    quick_test: bool = True

class TaskLogHandler(logging.Handler):
    def __init__(self, task_id):
        super().__init__()
        self.task_id = task_id
        self.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))

    def emit(self, record):
        log_entry = self.format(record)
        with SIMULATION_TASKS_LOCK:
            if self.task_id in SIMULATION_TASKS:
                SIMULATION_TASKS[self.task_id]['logs'].append(log_entry)

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

def run_pipeline_thread(task_id, request: SimulationRequest):
    # Setup directories
    task_dir = os.path.join("results", task_id)
    os.makedirs(task_dir, exist_ok=True)
    
    # Add custom logging handler for this task
    task_log_handler = TaskLogHandler(task_id)
    logger.addHandler(task_log_handler)
    
    # Add a file logging handler for backup
    file_handler = logging.FileHandler(os.path.join(task_dir, 'run.log'))
    file_handler.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(file_handler)
    
    try:
        config = load_config()
        enzymes_db = load_enzymes()
        enzyme_data = enzymes_db[request.enzyme]
        
        # Override output_dir in config path so that file outputs save inside the task_dir
        config['paths']['output_dir'] = task_dir
        
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['progress'] = 0.05
            SIMULATION_TASKS[task_id]['phase'] = "Phase 1: Input structure generation & validation"
            
        logger.info(f"═══ Phase 1: Input structure generation & validation ═══")
        scissile_type = enzyme_data.get('scissile_bond_type', 'ester_carbonyl')
        if 'ester' in scissile_type:
            linkage_type = 'ester'
        elif 'amide' in scissile_type:
            linkage_type = 'amide'
        elif 'glycosidic' in scissile_type:
            linkage_type = 'glycosidic'
        else:
            linkage_type = 'ester'
            
        logger.info(f"Generating 3D polymer coordinates ({linkage_type} linkage)...")
        polymer = build_polymer(request.smiles, request.length, config, linkage_type=linkage_type)
        success, checks = validate_input_structure(polymer)
        logger.info(f"Polymer validation checks: {checks}")
        if not success:
            raise RuntimeError(f"Input polymer structure failed validation checks: {checks}")
            
        pdb_file = f"data/{enzyme_data['pdb_id'].lower()}.pdb"
        if not os.path.exists(pdb_file):
            logger.info(f"Enzyme PDB {pdb_file} not found locally. Downloading...")
            os.makedirs("data", exist_ok=True)
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{enzyme_data['pdb_id']}.pdb", pdb_file)
            logger.info(f"Downloaded {pdb_file}")
            
        center = get_active_site_center(pdb_file, enzyme_data['catalytic_residues'])
        logger.info(f"Active site center: [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}]")
        
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['progress'] = 0.20
            SIMULATION_TASKS[task_id]['phase'] = "Phase 2: Anchor docking & growth loop (Funnel)"
            
        logger.info(f"═══ Phase 2: Anchor docking & growth loop (Funnel) ═══")
        anchor_length = config['docking'].get('anchor_length', 3)
        logger.info(f"Building {anchor_length}-mer anchor fragment ({linkage_type} linkage)...")
        anchor_base = build_polymer(request.smiles, anchor_length, config, linkage_type=linkage_type)
        
        engine_path = config['paths'].get('gnina_binary', 'gnina')
        has_docking_engine = (shutil.which(engine_path) is not None or 
                              shutil.which(config['paths'].get('vina_binary', 'vina')) is not None)
        
        num_modes = config['docking'].get('num_modes', 9)
        if has_docking_engine:
            logger.info(f"Real docking engine detected. Running GNINA/Vina for {num_modes} poses...")
            anchor_poses = dock_anchor(anchor_base, pdb_file, center, config)
        else:
            logger.info(f"No docking engine found on PATH. Generating {num_modes} mock anchor poses around active site center...")
            anchor_poses = generate_mock_anchor_poses(anchor_base, center, num_poses=num_modes)
            
        if not anchor_poses:
            raise RuntimeError("Docking returned 0 poses.")
            
        logger.info(f"Docking returned {len(anchor_poses)} anchor poses.")
        
        remaining_growth = request.length - anchor_length
        logger.info(f"Growing each of {len(anchor_poses)} anchor poses by {remaining_growth} monomers to {request.length}-mer...")
        grown_poses = []
        for i, anchor_pose in enumerate(anchor_poses):
            logger.info(f"  Pose {i+1}/{len(anchor_poses)}: Growing...")
            if remaining_growth > 0:
                grown = grow_polymer(anchor_pose, request.smiles, request.length, pdb_file, config, linkage_type=linkage_type)
            else:
                grown = anchor_pose
                
            if grown is not None:
                logger.info(f"  Pose {i+1}/{len(anchor_poses)}: ✓ Growth succeeded ({grown.GetNumAtoms()} atoms)")
                grown_poses.append((i+1, grown))
            else:
                logger.info(f"  Pose {i+1}/{len(anchor_poses)}: ✗ Growth FAILED (unresolvable clashes)")
                
        if not grown_poses:
            raise RuntimeError("All anchor poses failed to grow clash-free.")
            
        logger.info(f"Growth complete: {len(grown_poses)}/{len(anchor_poses)} poses survived.")
        
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['progress'] = 0.50
            SIMULATION_TASKS[task_id]['phase'] = "Phase 3: Catalytic Geometry Filter & MM-GBSA Scoring"
            
        logger.info(f"═══ Phase 3: Catalytic Geometry Filter & MM-GBSA Scoring ═══")
        passing_poses = []
        for pose_num, grown_mol in grown_poses:
            ligand_pdb = os.path.join(task_dir, f"grown_pose_{pose_num}.pdb")
            complex_pdb = os.path.join(task_dir, f"complex_pose_{pose_num}.pdb")
            from rdkit import Chem
            Chem.MolToPDBFile(grown_mol, ligand_pdb)
            save_complex(pdb_file, grown_mol, complex_pdb)
            
            verdict, distance = scan_catalytic_viability(complex_pdb, enzyme_data, config)
            cutoff = config['filters'].get('catalytic_cutoff', 4.5)
            if 'catalytic_cutoff_override' in enzyme_data:
                cutoff = enzyme_data['catalytic_cutoff_override']
                
            if verdict == 'PASS':
                logger.info(f"  Pose {pose_num}: ✓ PASS (distance {distance:.1f}Å < {cutoff:.1f}Å)")
                score_data = score_binding(complex_pdb, ligand_resname='UNL', config=config)
                logger.info(f"  Pose {pose_num}: Score = {score_data['final_score']:.2f} kcal/mol (SASA: {score_data['buried_sasa']:.1f} Å²)")
                passing_poses.append({
                    'pose_num': pose_num,
                    'grown_mol': grown_mol,
                    'ligand_pdb': ligand_pdb,
                    'complex_pdb': complex_pdb,
                    'distance': distance,
                    'score_data': score_data
                })
            else:
                logger.info(f"  Pose {pose_num}: ✗ REJECTED (distance {distance:.1f}Å > {cutoff:.1f}Å)")
                
        n_passed = len(passing_poses)
        logger.info(f"Filter complete: {n_passed}/{len(grown_poses)} poses passed catalytic geometry check.")
        if n_passed == 0:
            raise RuntimeError("No poses passed the catalytic geometry check.")
            
        # Rank by score (lower score is better)
        passing_poses.sort(key=lambda p: p['score_data']['final_score'])
        logger.info(f"Ranked {n_passed} poses by binding energy.")
        for rank, p in enumerate(passing_poses):
            logger.info(f"  Rank {rank+1}: Pose {p['pose_num']} — {p['score_data']['final_score']:.2f} kcal/mol")
            
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['progress'] = 0.70
            SIMULATION_TASKS[task_id]['phase'] = "Phase 4: Molecular Dynamics Validation (OpenMM)"
            
        logger.info(f"═══ Phase 4: Molecular Dynamics Validation (OpenMM) ═══")
        n_to_validate = min(3, n_passed)
        md_results = []
        is_mock_run = False
        
        for rank in range(n_to_validate):
            p = passing_poses[rank]
            logger.info(f"  Running MD on Rank {rank+1} (Pose {p['pose_num']})...")
            
            traj_dcd, is_mock = run_md_simulation(p['complex_pdb'], p['ligand_pdb'], config, quick_test=request.quick_test)
            if is_mock:
                is_mock_run = True
                logger.warning("OpenMM unavailable. MD results are MOCK DATA.")
                
            md_analysis = analyze_trajectory(traj_dcd, p['complex_pdb'], config, ligand_resname='UNL', enzyme_data=enzyme_data)
            logger.info(f"  Rank {rank+1}: {md_analysis['verdict']} (RMSD: {md_analysis['avg_rmsd']:.2f}Å, Cat.Dist: {md_analysis['avg_distance']:.2f}Å)")
            
            if os.path.exists(traj_dcd):
                try:
                    traj_score = score_binding(p['complex_pdb'], trajectory_dcd=traj_dcd, ligand_resname='UNL', config=config)
                    logger.info(f"  Rank {rank+1}: Trajectory-averaged score = {traj_score['final_score']:.2f} kcal/mol")
                    p['score_data'] = traj_score
                except Exception as e:
                    logger.warning(f"  Trajectory-averaged scoring failed: {e}")
                    
            md_results.append({
                **p,
                'md_analysis': md_analysis,
                'traj_dcd': traj_dcd
            })
            
            with SIMULATION_TASKS_LOCK:
                SIMULATION_TASKS[task_id]['progress'] = 0.70 + (0.20 * (rank + 1) / n_to_validate)
                
        # Find the best pose
        stable_results = [r for r in md_results if r['md_analysis']['verdict'] == 'STABLE']
        if stable_results:
            best = min(stable_results, key=lambda r: r['score_data']['final_score'])
        else:
            best = min(md_results, key=lambda r: r['score_data']['final_score'])
            
        # Save canonical outputs for this task
        best_complex = os.path.join(task_dir, "complex.pdb")
        best_ligand = os.path.join(task_dir, "grown_poly.pdb")
        shutil.copy(best['complex_pdb'], best_complex)
        shutil.copy(best['ligand_pdb'], best_ligand)
        if os.path.exists(best['traj_dcd']):
            shutil.copy(best['traj_dcd'], os.path.join(task_dir, "trajectory.dcd"))
            
        logger.info("═══ Pipeline Complete ═══")
        logger.info(f"Best pose: #{best['pose_num']} | Score: {best['score_data']['final_score']:.2f} | "
                    f"Geometry: {'PASS' if best['distance'] < config['filters']['catalytic_cutoff'] else 'FAIL'} | "
                    f"MD: {best['md_analysis']['verdict']}")
        
        # Prepare task result dictionary
        final_results = {
            "best_pose_num": best['pose_num'],
            "score": float(best['score_data']['final_score']),
            "interaction_energy": float(best['score_data']['interaction_energy']),
            "solvation_energy": float(best['score_data']['solvation_energy']),
            "sasa": float(best['score_data']['buried_sasa']),
            "attack_distance": float(best['distance']),
            "md_verdict": best['md_analysis']['verdict'],
            "avg_rmsd": float(best['md_analysis']['avg_rmsd']),
            "avg_distance": float(best['md_analysis']['avg_distance']),
            "is_mock": is_mock_run,
            "openmm_available": OPENMM_AVAILABLE,
            "openmm_error": OPENMM_ERROR,
            "all_poses_count": len(grown_poses),
            "passing_poses_count": n_passed,
            "linkage_type": linkage_type
        }
        
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['progress'] = 1.0
            SIMULATION_TASKS[task_id]['status'] = "completed"
            SIMULATION_TASKS[task_id]['phase'] = "Completed successfully"
            SIMULATION_TASKS[task_id]['results'] = final_results
            
    except Exception as e:
        logger.exception(f"Simulation pipeline failed with error: {e}")
        with SIMULATION_TASKS_LOCK:
            SIMULATION_TASKS[task_id]['status'] = "failed"
            SIMULATION_TASKS[task_id]['phase'] = "Failed"
            SIMULATION_TASKS[task_id]['error'] = str(e)
            
    finally:
        # Clean up logger handlers to prevent leak
        logger.removeHandler(task_log_handler)
        logger.removeHandler(file_handler)
        file_handler.close()

@app.get("/api/enzymes")
def get_enzymes():
    try:
        return load_enzymes()
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load enzymes: {e}")

@app.post("/api/simulate")
def start_simulation(request: SimulationRequest, background_tasks: BackgroundTasks):
    # Validate enzyme
    try:
        enzymes_db = load_enzymes()
        if request.enzyme not in enzymes_db:
            raise HTTPException(status_code=400, detail=f"Enzyme '{request.enzyme}' not found in database.")
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Enzymes database error: {e}")
        
    task_id = str(uuid.uuid4())
    
    with SIMULATION_TASKS_LOCK:
        SIMULATION_TASKS[task_id] = {
            "status": "running",
            "progress": 0.0,
            "phase": "Initialized",
            "logs": [],
            "results": None,
            "error": None
        }
        
    background_tasks.add_task(run_pipeline_thread, task_id, request)
    return {"task_id": task_id}

@app.get("/api/status/{task_id}")
def get_task_status(task_id: str):
    with SIMULATION_TASKS_LOCK:
        if task_id not in SIMULATION_TASKS:
            raise HTTPException(status_code=404, detail="Task not found")
        return SIMULATION_TASKS[task_id]

@app.get("/api/files/{task_id}/{file_type}")
def get_task_file(task_id: str, file_type: str):
    task_dir = os.path.join("results", task_id)
    if not os.path.exists(task_dir):
        raise HTTPException(status_code=404, detail="Task directory not found")
        
    if file_type == "complex":
        file_path = os.path.join(task_dir, "complex.pdb")
    elif file_type == "ligand":
        file_path = os.path.join(task_dir, "grown_poly.pdb")
    elif file_type == "trajectory":
        file_path = os.path.join(task_dir, "trajectory.dcd")
    else:
        raise HTTPException(status_code=400, detail="Invalid file type. Available: complex, ligand, trajectory")
        
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail=f"File {file_type} not generated or not found.")
        
    return FileResponse(file_path, filename=os.path.basename(file_path))

# Mount Frontend static files (built Vite React bundle) if the directory exists
frontend_dist = os.path.join("frontend", "dist")
if os.path.exists(frontend_dist):
    app.mount("/", StaticFiles(directory=frontend_dist, html=True), name="frontend")
else:
    @app.get("/")
    def index_fallback():
        return {"message": "SimDock Polymer API running. Please compile the React frontend inside frontend/dist to serve the UI."}
