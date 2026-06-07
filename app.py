import streamlit as st
import os
import yaml
import json
import numpy as np
import time
import shutil
import io
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from src.builder import build_polymer, get_active_site_center, validate_input_structure
from src.docking import dock_anchor
from src.grower import grow_polymer
from src.scanner import scan_catalytic_viability
from src.scorer import score_binding
from src.validator import run_md_simulation, analyze_trajectory, OPENMM_AVAILABLE
from src.utils import save_complex, setup_logging

# Setup page config
st.set_page_config(page_title="SimDock Polymer v2.0", layout="wide", initial_sidebar_state="expanded")

# Initialize session state variables
if "pipeline_running" not in st.session_state:
    st.session_state.pipeline_running = False
if "results" not in st.session_state:
    st.session_state.results = None

# ─── Premium Dark Glassmorphism Styling ───────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    .main {
        background-color: #0e1117;
        color: #c9d1d9;
        font-family: 'Inter', sans-serif;
    }
    .report-card {
        background: rgba(30, 30, 38, 0.65);
        backdrop-filter: blur(10px);
        border-radius: 16px;
        padding: 24px;
        box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.37);
        border: 1px solid rgba(255, 255, 255, 0.05);
        margin-bottom: 20px;
    }
    .card-title {
        font-size: 1.25rem;
        font-weight: 600;
        margin-bottom: 16px;
        color: #58a6ff;
    }
    .traffic-row {
        display: flex;
        align-items: center;
        margin-bottom: 14px;
    }
    .traffic-light {
        width: 18px;
        height: 18px;
        border-radius: 50%;
        margin-right: 14px;
        display: inline-block;
        flex-shrink: 0;
    }
    .light-green {
        background-color: #39ff14;
        box-shadow: 0 0 12px #39ff14;
    }
    .light-orange {
        background-color: #ff8c00;
        box-shadow: 0 0 12px #ff8c00;
    }
    .light-yellow {
        background-color: #ffcc00;
        box-shadow: 0 0 12px #ffcc00;
    }
    .light-red {
        background-color: #ff3333;
        box-shadow: 0 0 12px #ff3333;
    }
    .metric-value {
        font-size: 1.8rem;
        font-weight: bold;
        color: #ffffff;
    }
    .metric-label {
        font-size: 0.85rem;
        color: #8b949e;
    }
    .funnel-stat {
        background: rgba(30, 30, 38, 0.65);
        backdrop-filter: blur(10px);
        border-radius: 12px;
        padding: 14px 18px;
        border: 1px solid rgba(255, 255, 255, 0.05);
        margin-bottom: 8px;
        display: flex;
        justify-content: space-between;
        align-items: center;
    }
    .funnel-label {
        font-size: 0.9rem;
        color: #8b949e;
    }
    .funnel-value {
        font-size: 1.1rem;
        font-weight: 600;
        color: #58a6ff;
    }
</style>
""", unsafe_allow_html=True)


# ─── Helper Functions ─────────────────────────────────────────────────────────

def load_config():
    with open("config.yaml", "r") as f:
        return yaml.safe_load(f)

def load_enzymes():
    with open("data/enzymes.json", "r") as f:
        return json.load(f)

def render_3d_complex(complex_pdb_path):
    """Render py3Dmol viewer directly as an HTML component."""
    if not os.path.exists(complex_pdb_path):
        st.error("Complex structure file not found.")
        return
        
    with open(complex_pdb_path, 'r') as f:
        pdb_content = f.read().replace('\\', '\\\\').replace('\n', '\\n').replace('`', '\\`').replace("'", "\\'")
        
    html_code = f"""
    <div id="container" style="height: 450px; width: 100%; position: relative; border-radius: 12px; overflow: hidden; border: 1px solid #333;"></div>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <script>
      let viewer = $3Dmol.createViewer(document.getElementById('container'), {{defaultcolors: $3Dmol.rasmolElementColors}});
      let pdbData = '{pdb_content}';
      viewer.addModel(pdbData, "pdb");
      
      // Style protein
      viewer.setStyle({{protein: true}}, {{cartoon: {{color: 'spectrum', opacity: 0.85}}}});
      
      // Style ligand (UNL, LIG, POL, or HETATM)
      viewer.setStyle({{resname: 'UNL'}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}}}});
      viewer.setStyle({{resname: 'LIG'}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}}}});
      viewer.setStyle({{resname: 'POL'}}, {{stick: {{colorscheme: 'cyanCarbon', radius: 0.25}}}});
      
      // Highlight Catalytic Residues (Ser160)
      viewer.setStyle({{resseq: 160}}, {{stick: {{colorscheme: 'magentaCarbon', radius: 0.25}}, cartoon: {{color: 'magenta'}}}});
      
      viewer.zoomTo();
      viewer.render();
    </script>
    """
    st.components.v1.html(html_code, height=460)

def generate_mock_anchor_poses(anchor_mol, center, num_poses=9):
    """
    Generates multiple displaced/rotated anchor poses around the active site center.
    In a real docking run, GNINA/Vina produces these. For local testing without
    docking binaries, we create 9 poses by rotating and slightly displacing the anchor.
    """
    poses = []
    conf = anchor_mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(anchor_mol.GetNumAtoms())])
    centroid = np.mean(coords, axis=0)
    
    rng = np.random.RandomState(42)
    
    for pose_idx in range(num_poses):
        pose_mol = Chem.RWMol(anchor_mol)
        pose_conf = pose_mol.GetConformer()
        
        # Generate a random rotation matrix
        angle = (pose_idx / num_poses) * 2 * np.pi
        # Random rotation around z-axis + slight tilt
        cos_a, sin_a = np.cos(angle), np.sin(angle)
        tilt = rng.normal(0, 0.15)
        R = np.array([
            [cos_a, -sin_a, tilt],
            [sin_a,  cos_a, tilt],
            [0,      0,     1.0]
        ])
        
        # Small random displacement (up to 2 Å)
        displacement = rng.uniform(-2.0, 2.0, size=3)
        
        for i in range(pose_mol.GetNumAtoms()):
            pos = np.array(pose_conf.GetAtomPosition(i))
            # Center, rotate, displace, translate to active site
            pos_centered = pos - centroid
            pos_rotated = R.dot(pos_centered)
            pos_final = pos_rotated + center + displacement
            pose_conf.SetAtomPosition(i, pos_final)
            
        poses.append(pose_mol.GetMol())
    
    return poses

def geometry_traffic_light(distance, cutoff=4.5):
    """Returns (css_class, label) per the document spec."""
    if distance < cutoff:
        return 'light-green', 'Excellent'
    elif distance < cutoff * 1.5:
        return 'light-orange', 'Borderline'
    else:
        return 'light-red', 'Bad'

def energy_traffic_light(score):
    """Returns (css_class, label) per the document spec."""
    if score < -10:
        return 'light-green', 'Strong binding'
    elif score < 0:
        return 'light-yellow', 'Moderate'
    else:
        return 'light-red', 'Bad'

def reactivity_traffic_light(geom_pass, md_stable):
    """Returns (css_class, label) per the document spec."""
    if geom_pass and md_stable:
        return 'light-green', 'Ready to degrade'
    else:
        return 'light-red', 'Non-productive'

def generate_pdf_report(results, enzyme_name, polymer_type, chain_length):
    """Generates a one-page PDF summary report using matplotlib."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    buf = io.BytesIO()
    
    with PdfPages(buf) as pdf:
        fig, axes = plt.subplots(3, 1, figsize=(8.5, 11), gridspec_kw={'height_ratios': [1, 2, 1.5]})
        fig.patch.set_facecolor('#1a1a2e')
        
        # ── Header ──
        ax_header = axes[0]
        ax_header.set_facecolor('#1a1a2e')
        ax_header.axis('off')
        ax_header.text(0.5, 0.8, 'SimDock Polymer v2.0', fontsize=22, fontweight='bold',
                       color='#58a6ff', ha='center', va='center', transform=ax_header.transAxes)
        ax_header.text(0.5, 0.4, 'Catalytic Simulation Report', fontsize=14,
                       color='#c9d1d9', ha='center', va='center', transform=ax_header.transAxes)
        ax_header.text(0.5, 0.1, f'Enzyme: {enzyme_name}  |  Polymer: {polymer_type}  |  Chain: {chain_length}-mer',
                       fontsize=10, color='#8b949e', ha='center', va='center', transform=ax_header.transAxes)
        
        # ── Traffic Light Verdicts ──
        ax_tl = axes[1]
        ax_tl.set_facecolor('#16213e')
        ax_tl.axis('off')
        
        geom_color = '#39ff14' if results['geometry_verdict'] == 'PASS' else ('#ff8c00' if results['distance'] < 6.75 else '#ff3333')
        energy_color = '#39ff14' if results['binding_score'] < -10 else ('#ffcc00' if results['binding_score'] < 0 else '#ff3333')
        react_color = '#39ff14' if (results['geometry_verdict'] == 'PASS' and results['md_verdict'] == 'STABLE') else '#ff3333'
        
        verdicts = [
            ('Geometry Filter', f"{results['geometry_verdict']} — Attack distance: {results['distance']:.2f} Å", geom_color),
            ('Binding Energy', f"{results['binding_score']:.2f} kcal/mol (vdW: {results['interaction_energy']:.2f})", energy_color),
            ('MD Stability', f"{results['md_verdict']} — Avg RMSD: {results['md_rmsd']:.2f} Å", react_color),
        ]
        
        for i, (label, value, color) in enumerate(verdicts):
            y = 0.8 - i * 0.25
            ax_tl.add_patch(plt.Circle((0.08, y), 0.02, color=color, transform=ax_tl.transAxes))
            ax_tl.text(0.14, y, f'{label}:', fontsize=12, fontweight='bold', color='#ffffff',
                       va='center', transform=ax_tl.transAxes)
            ax_tl.text(0.14, y - 0.08, value, fontsize=10, color='#c9d1d9',
                       va='center', transform=ax_tl.transAxes)
        
        # ── Funnel Summary ──
        ax_funnel = axes[2]
        ax_funnel.set_facecolor('#1a1a2e')
        ax_funnel.axis('off')
        ax_funnel.text(0.5, 0.9, 'Pipeline Funnel Summary', fontsize=14, fontweight='bold',
                       color='#58a6ff', ha='center', va='center', transform=ax_funnel.transAxes)
        
        funnel_data = [
            ('Anchor Poses Generated', str(results.get('n_anchor_poses', '?'))),
            ('Grown Successfully', str(results.get('n_grown', '?'))),
            ('Passed Geometry Filter', str(results.get('n_passed_filter', '?'))),
            ('MD Validated', str(results.get('n_md_validated', '?'))),
        ]
        for i, (label, value) in enumerate(funnel_data):
            y = 0.65 - i * 0.18
            ax_funnel.text(0.15, y, label, fontsize=10, color='#8b949e',
                           va='center', transform=ax_funnel.transAxes)
            ax_funnel.text(0.85, y, value, fontsize=12, fontweight='bold', color='#58a6ff',
                           va='center', ha='right', transform=ax_funnel.transAxes)
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    
    buf.seek(0)
    return buf.getvalue()


# ─── Main App Layout ──────────────────────────────────────────────────────────

st.title("🔬 SimDock Polymer v2.0")
st.subheader("Universal Polymer-Enzyme Catalytic Simulation Engine")
st.write("---")

if not OPENMM_AVAILABLE:
    st.warning("⚠️ **Warning:** OpenMM or its forcefields are not available in this environment. "
               "The pipeline will fall back to generating **MOCK DATA** for the Molecular Dynamics (Phase 4) validation. "
               "Do not use these simulated trajectories or verdicts for publication.")

config = load_config()
logger = setup_logging()
enzymes_db = load_enzymes()

# ─── Sidebar: Screen 1 — Lab Bench (Setup) ───────────────────────────────────

st.sidebar.header("🧪 Lab Bench Controls")
enzyme_choice = st.sidebar.selectbox("Select Target Enzyme", list(enzymes_db.keys()))
enzyme_data = enzymes_db[enzyme_choice]

st.sidebar.markdown("**Enzyme Metadata:**")
st.sidebar.write(f"- PDB Code: `{enzyme_data['pdb_id']}`")
if 'description' in enzyme_data:
    st.sidebar.write(f"- {enzyme_data['description']}")
st.sidebar.write("- Catalytic Triad:")
for res, num in enzyme_data['catalytic_residues'].items():
    st.sidebar.write(f"  - {res}: `{num}`")
    
st.sidebar.write("---")

polymer_mode = st.sidebar.radio("Polymer Input Mode", ["Presets", "Custom SMILES"])
if polymer_mode == "Presets":
    poly_type = st.sidebar.selectbox("Select Polymer Repeat Unit", ["PET (Polyethylene Terephthalate)", "PLA (Polylactic Acid)"])
    if "PET" in poly_type:
        monomer_smiles = "[C@@H](OC(=O)c1ccc(C(=O)O)cc1)"
    else:
        monomer_smiles = "C(C(=O)O)O"  # PLA repeat unit analogue
else:
    monomer_smiles = st.sidebar.text_input("Monomer SMILES (with radical C end and acid end)", "[C@@H](OC(=O)c1ccc(C(=O)O)cc1)")

# Document specifies slider range 5 to 20 (Gap #5)
chain_length = st.sidebar.slider("Chain Length (Monomers)", 5, 20, 6)

st.sidebar.write("---")
run_btn = st.sidebar.button("🚀 Run Simulation", use_container_width=True)

# Execution State Controller
if run_btn:
    st.session_state.pipeline_running = True
    st.session_state.results = None
    
# Layout containers
status_container = st.container()
results_container = st.container()

# ─── Screen 2 — Simulation Status (Funnel Pipeline) ──────────────────────────

if st.session_state.pipeline_running:
    with status_container:
        st.header("🧬 Simulation Pipeline Progress")
        progress_bar = st.progress(0.0)
        status_text = st.empty()
        
        log_box = st.empty()
        logs = []
        
        def add_log(msg):
            logs.append(f"[{time.strftime('%H:%M:%S')}] {msg}")
            log_box.code("\n".join(logs), language="text")
            logger.info(msg)
            
        try:
            # ═══════════════════════════════════════════════════════════════
            # PHASE 1: Input Processing & Structure Generation
            # ═══════════════════════════════════════════════════════════════
            add_log("═══ Phase 1: Input structure generation & validation ═══")
            progress_bar.progress(0.05)
            
            monomer = Chem.MolFromSmiles(monomer_smiles)
            if monomer is None:
                raise ValueError("Invalid SMILES string entered.")
            
            # Get linkage type based on enzyme's scissile_bond_type
            scissile_type = enzyme_data.get('scissile_bond_type', 'ester_carbonyl')
            if 'ester' in scissile_type:
                linkage_type = 'ester'
            elif 'amide' in scissile_type:
                linkage_type = 'amide'
            elif 'glycosidic' in scissile_type:
                linkage_type = 'glycosidic'
            else:
                linkage_type = 'ester'
                
            # Build the full-length polymer for validation
            add_log(f"Generating 3D polymer coordinates ({linkage_type} linkage)...")
            polymer = build_polymer(monomer_smiles, chain_length, config, linkage_type=linkage_type)
            success, checks = validate_input_structure(polymer)
            add_log(f"Polymer validation checks: {checks}")
            if not success:
                raise RuntimeError("Polymer failed chemical validation checks.")
            
            # Locate active site center
            pdb_file = f"data/{enzyme_data['pdb_id'].lower()}.pdb"
            if not os.path.exists(pdb_file):
                add_log(f"Enzyme PDB {pdb_file} not found locally. Downloading...")
                import urllib.request
                urllib.request.urlretrieve(f"https://files.rcsb.org/download/{enzyme_data['pdb_id']}.pdb", pdb_file)
                add_log(f"Downloaded {pdb_file}")
                
            center = get_active_site_center(pdb_file, enzyme_data['catalytic_residues'])
            add_log(f"Active site center: [{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}]")
            
            # ═══════════════════════════════════════════════════════════════
            # PHASE 2: Anchor Docking & Growth Loop (FUNNEL ARCHITECTURE)
            # ═══════════════════════════════════════════════════════════════
            progress_bar.progress(0.15)
            add_log("═══ Phase 2: Anchor docking & growth loop (Funnel) ═══")
            
            # Build anchor fragment (first 3 monomers per document spec, Gap #13)
            anchor_length = config['docking'].get('anchor_length', 3)
            add_log(f"Building {anchor_length}-mer anchor fragment ({linkage_type} linkage)...")
            anchor_base = build_polymer(monomer_smiles, anchor_length, config, linkage_type=linkage_type)
            
            # Check if real docking engines are available
            engine_path = config['paths'].get('gnina_binary', 'gnina')
            has_docking_engine = (shutil.which(engine_path) is not None or 
                                  shutil.which(config['paths'].get('vina_binary', 'vina')) is not None)
            
            num_modes = config['docking'].get('num_modes', 9)
            
            if has_docking_engine:
                add_log(f"Real docking engine detected. Running GNINA/Vina for {num_modes} poses...")
                anchor_poses = dock_anchor(anchor_base, pdb_file, center, config)
                if not anchor_poses:
                    raise RuntimeError("Docking returned no poses.")
                add_log(f"Docking returned {len(anchor_poses)} anchor poses.")
            else:
                add_log(f"No docking engine on PATH. Generating {num_modes} mock anchor poses...")
                anchor_poses = generate_mock_anchor_poses(anchor_base, center, num_poses=num_modes)
                add_log(f"Generated {len(anchor_poses)} mock anchor poses around active site.")
            
            n_anchor_poses = len(anchor_poses)
            
            # ── Growth Loop: Extend each anchor to full chain length ──
            progress_bar.progress(0.25)
            remaining_growth = chain_length - anchor_length
            add_log(f"Growing each of {n_anchor_poses} anchor poses by {remaining_growth} monomers to {chain_length}-mer...")
            
            grown_poses = []
            os.makedirs("results", exist_ok=True)
            
            for i, anchor_pose in enumerate(anchor_poses):
                add_log(f"  Pose {i+1}/{n_anchor_poses}: Growing...")
                if remaining_growth > 0:
                    grown = grow_polymer(anchor_pose, monomer_smiles, chain_length, pdb_file, config)
                else:
                    grown = anchor_pose
                    
                if grown is None:
                    add_log(f"  Pose {i+1}/{n_anchor_poses}: ✗ Growth FAILED (unresolvable clashes)")
                else:
                    add_log(f"  Pose {i+1}/{n_anchor_poses}: ✓ Growth succeeded ({grown.GetNumAtoms()} atoms)")
                    grown_poses.append((i+1, grown))
                    
                # Update progress within growth phase
                frac = 0.25 + (0.20 * (i + 1) / n_anchor_poses)
                progress_bar.progress(min(frac, 0.45))
            
            n_grown = len(grown_poses)
            add_log(f"Growth complete: {n_grown}/{n_anchor_poses} poses survived.")
            
            if n_grown == 0:
                raise RuntimeError("All poses failed during growth. No candidates to filter.")
            
            # ═══════════════════════════════════════════════════════════════
            # PHASE 3: Catalytic Geometry Filter + MM-GBSA Scoring
            # ═══════════════════════════════════════════════════════════════
            progress_bar.progress(0.50)
            add_log("═══ Phase 3: Catalytic Geometry Filter & MM-GBSA Scoring ═══")
            
            passing_poses = []
            
            for pose_num, grown_mol in grown_poses:
                # Save grown pose and create complex
                ligand_pdb = f"results/grown_pose_{pose_num}.pdb"
                complex_pdb = f"results/complex_pose_{pose_num}.pdb"
                Chem.MolToPDBFile(grown_mol, ligand_pdb)
                save_complex(pdb_file, grown_mol, complex_pdb)
                
                # Run catalytic geometry filter
                verdict, distance = scan_catalytic_viability(complex_pdb, enzyme_data, config)
                
                if verdict == 'PASS':
                    add_log(f"  Pose {pose_num}: ✓ PASS (distance {distance:.1f}Å < {config['filters']['catalytic_cutoff']}Å)")
                    
                    # Score passing poses with MM-GBSA
                    score_data = score_binding(complex_pdb, ligand_resname='UNL', config=config)
                    add_log(f"  Pose {pose_num}: Score = {score_data['final_score']:.2f} kcal/mol (SASA: {score_data['buried_sasa']:.1f} Å²)")
                    
                    passing_poses.append({
                        'pose_num': pose_num,
                        'grown_mol': grown_mol,
                        'ligand_pdb': ligand_pdb,
                        'complex_pdb': complex_pdb,
                        'distance': distance,
                        'score_data': score_data,
                    })
                else:
                    add_log(f"  Pose {pose_num}: ✗ REJECTED (distance {distance:.1f}Å > {config['filters']['catalytic_cutoff']}Å)")
            
            n_passed_filter = len(passing_poses)
            add_log(f"Filter complete: {n_passed_filter}/{n_grown} poses passed catalytic geometry check.")
            
            if n_passed_filter == 0:
                add_log("🚨 No poses passed the geometry filter.")
                st.error("No poses passed the geometry filter.")
                st.stop()
            
            # Rank passing poses by binding score (lower = better)
            passing_poses.sort(key=lambda p: p['score_data']['final_score'])
            add_log(f"Ranked {len(passing_poses)} poses by binding energy.")
            for rank, p in enumerate(passing_poses):
                add_log(f"  Rank {rank+1}: Pose {p['pose_num']} — {p['score_data']['final_score']:.2f} kcal/mol")
            
            # ═══════════════════════════════════════════════════════════════
            # PHASE 4: MD Validation (Top 1-3 poses)
            # ═══════════════════════════════════════════════════════════════
            progress_bar.progress(0.70)
            add_log("═══ Phase 4: Molecular Dynamics Validation (OpenMM) ═══")
            
            # Validate top 1-3 poses
            n_to_validate = min(3, len(passing_poses))
            md_results = []
            is_mock_run = False
            
            for rank in range(n_to_validate):
                p = passing_poses[rank]
                add_log(f"  Running MD on Rank {rank+1} (Pose {p['pose_num']})...")
                
                traj_dcd, is_mock = run_md_simulation(p['complex_pdb'], p['ligand_pdb'], config, quick_test=True)
                if is_mock:
                    is_mock_run = True
                    st.error("⚠️ OpenMM unavailable. MD results are MOCK DATA.")
                md_analysis = analyze_trajectory(traj_dcd, p['complex_pdb'], config, ligand_resname='UNL', enzyme_data=enzyme_data)
                
                add_log(f"  Rank {rank+1}: {md_analysis['verdict']} (RMSD: {md_analysis['avg_rmsd']:.2f}Å, Cat.Dist: {md_analysis['avg_distance']:.2f}Å)")
                
                # Gap #11: Re-score with trajectory for MM-GBSA averaging
                if os.path.exists(traj_dcd):
                    try:
                        traj_score = score_binding(p['complex_pdb'], trajectory_dcd=traj_dcd, ligand_resname='UNL', config=config)
                        add_log(f"  Rank {rank+1}: Trajectory-averaged score = {traj_score['final_score']:.2f} kcal/mol")
                        p['score_data'] = traj_score  # Update with trajectory-averaged score
                    except Exception:
                        pass  # Keep single-frame score if trajectory scoring fails
                
                md_results.append({
                    **p,
                    'md_analysis': md_analysis,
                    'traj_dcd': traj_dcd,
                })
                
                frac = 0.70 + (0.20 * (rank + 1) / n_to_validate)
                progress_bar.progress(min(frac, 0.90))
            
            n_md_validated = len(md_results)
            
            # Select best pose (STABLE preferred, then best score)
            stable_results = [r for r in md_results if r['md_analysis']['verdict'] == 'STABLE']
            if stable_results:
                best = min(stable_results, key=lambda r: r['score_data']['final_score'])
            else:
                best = min(md_results, key=lambda r: r['score_data']['final_score'])
            
            # Copy best complex to canonical path
            best_complex = "results/complex.pdb"
            best_ligand = "results/grown_poly.pdb"
            shutil.copy(best['complex_pdb'], best_complex)
            shutil.copy(best['ligand_pdb'], best_ligand)
            
            # ═══════════════════════════════════════════════════════════════
            # FINALIZE
            # ═══════════════════════════════════════════════════════════════
            progress_bar.progress(1.0)
            add_log("═══ Pipeline Complete ═══")
            add_log(f"Best pose: #{best['pose_num']} | Score: {best['score_data']['final_score']:.2f} | "
                    f"Geometry: {'PASS' if best['distance'] < config['filters']['catalytic_cutoff'] else 'FAIL'} | "
                    f"MD: {best['md_analysis']['verdict']}")
            
            st.session_state.pipeline_running = False
            
            geom_verdict = 'PASS' if best['distance'] < config['filters']['catalytic_cutoff'] else 'FAIL'
            
            st.session_state.results = {
                'distance': best['distance'],
                'geometry_verdict': geom_verdict,
                'binding_score': best['score_data']['final_score'],
                'interaction_energy': best['score_data']['interaction_energy'],
                'buried_sasa': best['score_data']['buried_sasa'],
                'md_verdict': best['md_analysis']['verdict'],
                'md_rmsd': best['md_analysis']['avg_rmsd'],
                'md_rmsd_fraction': best['md_analysis']['rmsd_stable_fraction'],
                'md_catalytic_fraction': best['md_analysis']['catalytic_stable_fraction'],
                'md_avg_distance': best['md_analysis']['avg_distance'],
                'complex_pdb_path': best_complex,
                'ligand_pdb_path': best_ligand,
                'best_pose_num': best['pose_num'],
                'is_mock': is_mock_run,
                # Funnel statistics
                'n_anchor_poses': n_anchor_poses,
                'n_grown': n_grown,
                'n_passed_filter': n_passed_filter,
                'n_md_validated': n_md_validated,
            }
            
            st.rerun()
            
        except Exception as e:
            add_log(f"🚨 Pipeline error: {e}")
            st.session_state.pipeline_running = False
            st.error(f"Simulation failed: {e}")


# ─── Screen 3 — Report Card ──────────────────────────────────────────────────

if st.session_state.results is not None:
    res = st.session_state.results
    
    with results_container:
        st.header("📊 Simulation Analysis & Report Card")
        
        if res.get('is_mock', False):
            st.error("⚠️ **Warning:** OpenMM was unavailable or crashed during validation. "
                     "The MD stability metrics above are based on **MOCK DATA**. "
                     "Do NOT use these results for publication.")
            
        # Grid layout for report card and 3D viewer
        col1, col2 = st.columns([1, 1])
        
        with col1:
            # ── Traffic Light Verdict Display ──
            geom_class, geom_label = geometry_traffic_light(res['distance'])
            energy_class, energy_label = energy_traffic_light(res['binding_score'])
            react_class, react_label = reactivity_traffic_light(
                res['geometry_verdict'] == 'PASS',
                res['md_verdict'] == 'STABLE'
            )
            
            st.markdown(f"""
            <div class="report-card">
                <div class="card-title">🔬 Catalytic Triad Report Card — Pose #{res['best_pose_num']}</div>
                <div class="traffic-row">
                    <div class="traffic-light {geom_class}"></div>
                    <div>
                        <b>Geometry:</b> {geom_label} — Attack distance: {res['distance']:.2f} Å<br/>
                        <span style="font-size: 0.8rem; color: #8b949e;">Scissile C=O must be within {config['filters']['catalytic_cutoff']}Å of Ser160 OG nucleophile.</span>
                    </div>
                </div>
                <div class="traffic-row">
                    <div class="traffic-light {energy_class}"></div>
                    <div>
                        <b>Energy:</b> {energy_label} — {res['binding_score']:.2f} kcal/mol (vdW: {res['interaction_energy']:.2f})<br/>
                        <span style="font-size: 0.8rem; color: #8b949e;">MM-GBSA implicit solvent scoring. Negative = stable binding.</span>
                    </div>
                </div>
                <div class="traffic-row">
                    <div class="traffic-light {react_class}"></div>
                    <div>
                        <b>Reactivity:</b> {react_label}<br/>
                        <span style="font-size: 0.8rem; color: #8b949e; display: block; margin-top: 4px;">
                            • <b>RMSD Stability:</b> {'✅ STABLE' if res.get('md_rmsd_fraction', 1.0) > 0.5 else '❌ UNSTABLE'} 
                            (avg {res['md_rmsd']:.2f} Å, {res.get('md_rmsd_fraction', 1.0)*100:.0f}% of frames &lt; 3.0 Å)<br/>
                            • <b>Catalytic Distance:</b> {'✅ MAINTAINED' if res.get('md_catalytic_fraction', 1.0) > 0.5 else '❌ LOST'} 
                            (avg {res.get('md_avg_distance', 0.0):.2f} Å, {res.get('md_catalytic_fraction', 1.0)*100:.0f}% of frames &lt; {config['filters'].get('catalytic_cutoff_md', 5.0)} Å)
                        </span>
                    </div>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # ── Funnel Statistics ──
            st.markdown(f"""
            <div class="report-card">
                <div class="card-title">🔽 Pipeline Funnel Summary</div>
                <div class="funnel-stat">
                    <span class="funnel-label">Anchor Poses Generated</span>
                    <span class="funnel-value">{res['n_anchor_poses']}</span>
                </div>
                <div class="funnel-stat">
                    <span class="funnel-label">Grown Successfully</span>
                    <span class="funnel-value">{res['n_grown']}</span>
                </div>
                <div class="funnel-stat">
                    <span class="funnel-label">Passed Geometry Filter</span>
                    <span class="funnel-value">{res['n_passed_filter']}</span>
                </div>
                <div class="funnel-stat">
                    <span class="funnel-label">MD Validated</span>
                    <span class="funnel-value">{res['n_md_validated']}</span>
                </div>
            </div>
            """, unsafe_allow_html=True)
            
            # ── Metric cards ──
            m1, m2 = st.columns(2)
            with m1:
                st.markdown(f"""
                <div class="report-card" style="text-align: center;">
                    <div class="metric-label">BURIED SASA</div>
                    <div class="metric-value">{res['buried_sasa']:.1f} Å²</div>
                </div>
                """, unsafe_allow_html=True)
            with m2:
                final_outcome = 'REACTIVE' if (res['geometry_verdict'] == 'PASS' and res['md_verdict'] == 'STABLE') else 'INERT'
                outcome_color = '#39ff14' if final_outcome == 'REACTIVE' else '#ff3333'
                st.markdown(f"""
                <div class="report-card" style="text-align: center;">
                    <div class="metric-label">FINAL OUTCOME</div>
                    <div class="metric-value" style="color: {outcome_color};">
                        {final_outcome}
                    </div>
                </div>
                """, unsafe_allow_html=True)
                
            # ── Download Section ──
            st.subheader("💾 Export Results")
            
            dl_col1, dl_col2 = st.columns(2)
            
            with dl_col1:
                # PDB download
                with open(res['complex_pdb_path'], 'r') as f:
                    complex_pdb_data = f.read()
                st.download_button("📦 Download Complex PDB", complex_pdb_data, "complex_docked.pdb", 
                                   "text/plain", use_container_width=True)
            
            with dl_col2:
                # PDF summary report (Gap #6)
                try:
                    pdf_bytes = generate_pdf_report(
                        res, enzyme_choice,
                        poly_type if polymer_mode == "Presets" else "Custom",
                        chain_length
                    )
                    st.download_button("📄 Download PDF Report", pdf_bytes, "simdock_report.pdf",
                                       "application/pdf", use_container_width=True)
                except Exception as e:
                    st.warning(f"PDF generation failed: {e}")
            
        with col2:
            st.subheader("🖥️ Interactive 3D Active Site Pocket")
            render_3d_complex(res['complex_pdb_path'])
            st.write("**Magenta** stick residue = Ser160 nucleophile. **Cyan** stick structure = grown polymer chain.")
