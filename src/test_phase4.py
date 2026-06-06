import sys
import os
import yaml

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.validator import run_md_simulation, analyze_trajectory

def test_phase4():
    print("=== STARTING PHASE 4 TEST ===")
    
    # 1. Load config
    config_path = "config.yaml"
    if not os.path.exists(config_path):
        print(f"Error: Config not found at {config_path}")
        return False
        
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    complex_pdb = "results/complex.pdb"
    ligand_pdb = "results/grown_poly.pdb"
    
    if not os.path.exists(complex_pdb) or not os.path.exists(ligand_pdb):
        print(f"Error: Required PDB files not found in results/. Please run previous phases.")
        return False
        
    # 2. Run MD Simulation (quick 500-step test)
    print("Launching OpenMM simulation...")
    try:
        traj_dcd, is_mock = run_md_simulation(complex_pdb, ligand_pdb, config, quick_test=True)
        print(f"Trajectory DCD written to: {traj_dcd} (Is Mock: {is_mock})")
        assert os.path.exists(traj_dcd), "Trajectory DCD was not written"
    except Exception as e:
        print(f"MD Simulation run failed: {e}")
        import traceback
        traceback.print_exc()
        return False
        
    # 3. Analyze Trajectory
    print("\nAnalyzing DCD trajectory...")
    try:
        mock_enzyme_data = {
            "nucleophile_res_num": 160,
            "nucleophile_res_name": "SER",
            "nucleophile_atom_name": "OG",
            "scissile_bond_type": "ester_carbonyl",
            "scissile_bond_position": "terminal"
        }
        analysis = analyze_trajectory(traj_dcd, complex_pdb, config, ligand_resname='UNL', enzyme_data=mock_enzyme_data)
        print(f"Analysis Output Detail: {analysis}")
        print(f"Trajectory Verdict: {analysis['verdict']}")
        print(f"Average Ligand RMSD: {analysis['avg_rmsd']:.2f} Angstroms")
        print(f"Average Catalytic Distance: {analysis['avg_distance']:.2f} Angstroms")
        
        # Check assertions
        assert 'verdict' in analysis
        assert 'avg_rmsd' in analysis
        assert 'avg_distance' in analysis
        
        print("\nPhase 4 testing complete and PASSED!")
        return True
    except Exception as e:
        print(f"Trajectory analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    ok = test_phase4()
    sys.exit(0 if ok else 1)
