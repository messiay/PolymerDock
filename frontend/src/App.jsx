import React, { useState, useEffect, useRef } from 'react';
import { 
  Play, 
  Settings, 
  Terminal, 
  Activity, 
  Download, 
  CheckCircle, 
  XCircle, 
  Cpu, 
  Shield, 
  ChevronRight, 
  RotateCcw,
  Sparkles,
  Info
} from 'lucide-react';
import MoleculeViewer from './components/MoleculeViewer';

export default function App() {
  // Config & Database states
  const [enzymes, setEnzymes] = useState({});
  const [selectedEnzyme, setSelectedEnzyme] = useState('Amylase');
  const [smiles, setSmiles] = useState('OC1C(O)C(O)C(O)C(CO)O1');
  const [length, setLength] = useState(5);
  const [quickTest, setQuickTest] = useState(true);
  const [loadingEnzymes, setLoadingEnzymes] = useState(true);

  // Active Task States
  const [taskId, setTaskId] = useState(null);
  const [status, setStatus] = useState('idle'); // idle, running, completed, failed
  const [progress, setProgress] = useState(0);
  const [phase, setPhase] = useState('');
  const [logs, setLogs] = useState([]);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);

  const terminalEndRef = useRef(null);
  const pollIntervalRef = useRef(null);

  // Fetch enzymes on mount
  useEffect(() => {
    async function fetchEnzymes() {
      try {
        setLoadingEnzymes(true);
        const response = await fetch('/api/enzymes');
        if (!response.ok) throw new Error("Failed to load enzyme database.");
        const data = await response.ok ? await response.json() : {};
        setEnzymes(data);
        setLoadingEnzymes(false);
      } catch (err) {
        console.error("Error loading enzymes:", err);
        // Fallback placeholder database if backend fails temporarily
        setEnzymes({
          "Amylase": {
            "pdb_id": "1HNY",
            "description": "Human pancreatic alpha-amylase - degrades starch (alpha-1,4-glucosidic bonds).",
            "catalytic_residues": [197, 233, 300],
            "nucleophile_res_num": 197,
            "nucleophile_res_name": "ASP",
            "scissile_bond_type": "glycosidic_carbon"
          }
        });
        setLoadingEnzymes(false);
      }
    }
    fetchEnzymes();
  }, []);

  // Update default SMILES based on selected enzyme
  useEffect(() => {
    if (!enzymes[selectedEnzyme]) return;
    const bondType = enzymes[selectedEnzyme].scissile_bond_type || '';
    if (bondType.includes('glycosidic')) {
      setSmiles('OC1C(O)C(O)C(O)C(CO)O1'); // Glucose / Starch monomer
    } else if (selectedEnzyme.toLowerCase().includes('nylon')) {
      setSmiles('C(CCCC(=O)O)CN'); // 6-aminohexanoic acid (nylon monomer)
    } else if (selectedEnzyme.toLowerCase().includes('pet')) {
      setSmiles('O=C(O)c1ccc(C(=O)OCCO)cc1'); // MHET (PET monomer)
    } else {
      setSmiles('O=C(O)CCO'); // General ester monomer
    }
  }, [selectedEnzyme, enzymes]);

  // Scroll to bottom of terminal whenever logs update
  useEffect(() => {
    if (terminalEndRef.current) {
      terminalEndRef.current.scrollIntoView({ behavior: 'smooth' });
    }
  }, [logs]);

  // Clear polling interval on unmount
  useEffect(() => {
    return () => {
      if (pollIntervalRef.current) clearInterval(pollIntervalRef.current);
    };
  }, []);

  // Start the simulation run
  const handleStartSimulation = async () => {
    try {
      setStatus('running');
      setProgress(0);
      setPhase('Initializing task...');
      setLogs([]);
      setResults(null);
      setError(null);
      setTaskId(null);

      const response = await fetch('/api/simulate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          enzyme: selectedEnzyme,
          smiles: smiles,
          length: parseInt(length),
          quick_test: quickTest
        })
      });

      if (!response.ok) {
        const errData = await response.json();
        throw new Error(errData.detail || "Failed to start simulation.");
      }

      const data = await response.json();
      setTaskId(data.task_id);
      
      // Start polling status
      startPolling(data.task_id);
    } catch (err) {
      setStatus('failed');
      setError(err.message);
    }
  };

  // Poll status endpoint
  const startPolling = (tid) => {
    if (pollIntervalRef.current) clearInterval(pollIntervalRef.current);
    
    pollIntervalRef.current = setInterval(async () => {
      try {
        const response = await fetch(`/api/status/${tid}`);
        if (!response.ok) throw new Error("Status polling error.");
        
        const data = await response.json();
        
        // Update states
        setProgress(data.progress);
        setPhase(data.phase);
        setLogs(data.logs || []);
        
        if (data.status === 'completed') {
          setStatus('completed');
          setResults(data.results);
          clearInterval(pollIntervalRef.current);
        } else if (data.status === 'failed') {
          setStatus('failed');
          setError(data.error || "Simulation run encountered a fatal error.");
          clearInterval(pollIntervalRef.current);
        }
      } catch (err) {
        console.error(err);
      }
    }, 1000);
  };

  const resetState = () => {
    setStatus('idle');
    setProgress(0);
    setPhase('');
    setLogs([]);
    setResults(null);
    setError(null);
    setTaskId(null);
  };

  const getPhaseIndex = () => {
    if (phase.includes("Phase 1")) return 1;
    if (phase.includes("Phase 2")) return 2;
    if (phase.includes("Phase 3")) return 3;
    if (phase.includes("Phase 4") || phase.includes(" equilibration") || phase.includes("MD")) return 4;
    if (phase.includes("Complete")) return 5;
    return 0;
  };

  const currentEnzymeData = enzymes[selectedEnzyme] || {};

  return (
    <div style={{ minHeight: '100vh', display: 'flex', flexDirection: 'column' }}>
      {/* Top Header */}
      <header className="glass-card" style={{ 
        margin: '12px 20px', 
        padding: '12px 24px', 
        display: 'flex', 
        alignItems: 'center', 
        justifyContent: 'space-between',
        borderRadius: '12px'
      }}>
        <div style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
          <div style={{
            background: 'linear-gradient(135deg, #8b5cf6 0%, #3b82f6 100%)',
            padding: '8px',
            borderRadius: '10px',
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            boxShadow: '0 0 15px rgba(139, 92, 246, 0.4)'
          }}>
            <Activity size={24} color="white" />
          </div>
          <div>
            <h1 className="title-gradient" style={{ margin: 0, fontSize: '1.4rem', lineHeight: '1' }}>
              SimDock Polymer v2.0
            </h1>
            <span style={{ fontSize: '0.75rem', color: '#8b949e', letterSpacing: '0.5px' }}>
              UNIVERSAL CATALYTIC SIMULATION ENGINE
            </span>
          </div>
        </div>
        
        <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
          <div className="status-badge" style={{ 
            background: 'rgba(59, 130, 246, 0.1)', 
            border: '1px solid rgba(59, 130, 246, 0.2)',
            color: '#3b82f6'
          }}>
            <Cpu size={14} />
            <span>GPU Accelerated</span>
          </div>
        </div>
      </header>

      {/* Main Dashboard Layout */}
      <div className="dashboard-grid">
        {/* Left Sidebar: Controls */}
        <aside className="glass-card" style={{ padding: '20px', display: 'flex', flexDirection: 'column', gap: '20px', height: 'fit-content' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '8px', borderBottom: '1px solid rgba(255,255,255,0.06)', paddingBottom: '12px' }}>
            <Settings size={18} color="#8b5cf6" />
            <h2 className="font-outfit" style={{ fontSize: '1.1rem', margin: 0, fontWeight: '600' }}>Simulation Controls</h2>
          </div>

          {/* Enzyme Selector */}
          <div>
            <label style={{ fontSize: '0.85rem', color: '#9ca3af', fontWeight: '500', display: 'block', marginBottom: '8px' }}>
              Target Enzyme
            </label>
            {loadingEnzymes ? (
              <div style={{ fontSize: '0.9rem', color: '#8b949e' }}>Loading database...</div>
            ) : (
              <select 
                className="form-input"
                value={selectedEnzyme}
                onChange={(e) => setSelectedEnzyme(e.target.value)}
                disabled={status === 'running'}
              >
                {Object.keys(enzymes).map((name) => (
                  <option key={name} value={name}>{name} (PDB: {enzymes[name].pdb_id})</option>
                ))}
              </select>
            )}
            {currentEnzymeData.description && (
              <div style={{ marginTop: '6px', fontSize: '0.75rem', color: '#8b949e', lineHeight: '1.3' }}>
                {currentEnzymeData.description}
              </div>
            )}
          </div>

          {/* SMILES Input */}
          <div>
            <label style={{ fontSize: '0.85rem', color: '#9ca3af', fontWeight: '500', display: 'block', marginBottom: '8px' }}>
              Monomer SMILES String
            </label>
            <input 
              type="text" 
              className="form-input" 
              value={smiles}
              onChange={(e) => setSmiles(e.target.value)}
              disabled={status === 'running'}
            />
            <div style={{ marginTop: '6px', fontSize: '0.72rem', color: '#8b949e' }}>
              Valid chemical SMILES representing the building block.
            </div>
          </div>

          {/* Length Slider */}
          <div>
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.85rem', color: '#9ca3af', fontWeight: '500', marginBottom: '8px' }}>
              <span>Polymer Chain Length</span>
              <span style={{ color: '#8b5cf6', fontWeight: 'bold' }}>{length}-mer</span>
            </div>
            <input 
              type="range" 
              min="3" 
              max="10" 
              value={length} 
              onChange={(e) => setLength(e.target.value)}
              style={{ width: '100%', accentColor: '#8b5cf6', cursor: 'pointer' }}
              disabled={status === 'running'}
            />
            <div style={{ display: 'flex', justifyContent: 'space-between', fontSize: '0.7rem', color: '#8b949e', marginTop: '4px' }}>
              <span>3 monomers</span>
              <span>10 monomers</span>
            </div>
          </div>

          {/* MD Run Toggle */}
          <div>
            <label style={{ fontSize: '0.85rem', color: '#9ca3af', fontWeight: '500', display: 'block', marginBottom: '8px' }}>
              OpenMM MD Mode
            </label>
            <div style={{ display: 'flex', gap: '8px' }}>
              <button 
                type="button"
                onClick={() => setQuickTest(true)}
                style={{
                  flex: 1,
                  padding: '8px',
                  borderRadius: '6px',
                  border: '1px solid',
                  borderColor: quickTest ? '#3b82f6' : 'rgba(255,255,255,0.06)',
                  background: quickTest ? 'rgba(59, 130, 246, 0.15)' : 'transparent',
                  color: quickTest ? '#3b82f6' : '#9ca3af',
                  fontSize: '0.8rem',
                  fontWeight: '600',
                  cursor: 'pointer'
                }}
                disabled={status === 'running'}
              >
                Quick Test (500s)
              </button>
              <button 
                type="button"
                onClick={() => setQuickTest(false)}
                style={{
                  flex: 1,
                  padding: '8px',
                  borderRadius: '6px',
                  border: '1px solid',
                  borderColor: !quickTest ? '#8b5cf6' : 'rgba(255,255,255,0.06)',
                  background: !quickTest ? 'rgba(139, 92, 246, 0.15)' : 'transparent',
                  color: !quickTest ? '#a78bfa' : '#9ca3af',
                  fontSize: '0.8rem',
                  fontWeight: '600',
                  cursor: 'pointer'
                }}
                disabled={status === 'running'}
              >
                Production (10ns)
              </button>
            </div>
          </div>

          {/* Action Buttons */}
          {status !== 'running' ? (
            <button 
              className="btn-primary"
              onClick={handleStartSimulation}
              disabled={loadingEnzymes}
              style={{ marginTop: '10px' }}
            >
              <Play size={16} fill="white" />
              <span>Launch Engine Run</span>
            </button>
          ) : (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '10px', marginTop: '10px' }}>
              <div className="btn-primary" style={{ cursor: 'wait' }}>
                <div style={{
                  width: '16px',
                  height: '16px',
                  border: '2px solid rgba(255,255,255,0.3)',
                  borderTop: '2px solid white',
                  borderRadius: '50%',
                  animation: 'spin 0.8s linear infinite'
                }} />
                <span>Running Pipeline...</span>
              </div>
              <style dangerouslySetInnerHTML={{__html: `@keyframes spin { 0% { transform: rotate(0deg); } 100% { transform: rotate(360deg); } }`}} />
            </div>
          )}

          {status !== 'idle' && status !== 'running' && (
            <button 
              onClick={resetState}
              style={{
                width: '100%',
                background: 'rgba(255, 255, 255, 0.05)',
                border: '1px solid rgba(255,255,255,0.08)',
                color: '#f3f4f6',
                borderRadius: '8px',
                padding: '10px',
                fontSize: '0.85rem',
                fontWeight: '600',
                cursor: 'pointer',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                gap: '8px'
              }}
            >
              <RotateCcw size={16} />
              <span>Reset Configuration</span>
            </button>
          )}
        </aside>

        {/* Right Area: Logs, 3D structure and output metrics card */}
        <main style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
          
          {/* Phase Progress Bar for Simulation Run */}
          {status === 'running' && (
            <div className="glass-card pulse-border-glow" style={{ padding: '20px' }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px', fontSize: '0.9rem' }}>
                <span className="font-outfit" style={{ fontWeight: '600', color: '#a78bfa' }}>
                  {phase}
                </span>
                <span style={{ fontWeight: 'bold' }}>{Math.round(progress * 100)}%</span>
              </div>
              
              {/* Actual Progress Bar */}
              <div style={{ width: '100%', height: '6px', background: 'rgba(255,255,255,0.05)', borderRadius: '3px', overflow: 'hidden' }}>
                <div style={{ 
                  width: `${progress * 100}%`, 
                  height: '100%', 
                  background: 'linear-gradient(90deg, #8b5cf6 0%, #3b82f6 100%)',
                  boxShadow: '0 0 8px #8b5cf6',
                  transition: 'width 0.4s ease'
                }} />
              </div>

              {/* Progress Steps Indicators */}
              <div style={{ display: 'flex', justifyContent: 'space-between', marginTop: '16px' }}>
                {[
                  { id: 1, name: "P1: Build" },
                  { id: 2, name: "P2: Dock" },
                  { id: 3, name: "P3: Filter" },
                  { id: 4, name: "P4: OpenMM" }
                ].map((step) => {
                  const pIdx = getPhaseIndex();
                  const isDone = pIdx > step.id;
                  const isActive = pIdx === step.id;
                  return (
                    <div key={step.id} style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                      <div style={{
                        width: '18px',
                        height: '18px',
                        borderRadius: '50%',
                        background: isDone ? '#10b981' : isActive ? '#8b5cf6' : 'rgba(255,255,255,0.05)',
                        border: '1px solid',
                        borderColor: isDone ? '#10b981' : isActive ? '#a78bfa' : 'rgba(255,255,255,0.1)',
                        display: 'flex',
                        alignItems: 'center',
                        justifyContent: 'center',
                        fontSize: '0.65rem',
                        fontWeight: 'bold',
                        color: isDone || isActive ? 'white' : '#8b949e'
                      }}>
                        {isDone ? "✓" : step.id}
                      </div>
                      <span style={{ 
                        fontSize: '0.75rem', 
                        color: isDone ? '#10b981' : isActive ? '#f3f4f6' : '#8b949e',
                        fontWeight: isActive ? '600' : 'normal'
                      }}>{step.name}</span>
                    </div>
                  );
                })}
              </div>
            </div>
          )}

          {/* IDLE state display */}
          {status === 'idle' && (
            <div className="glass-card" style={{ padding: '40px', textAlign: 'center', display: 'flex', flexDirection: 'column', alignItems: 'center', gap: '20px' }}>
              <div style={{
                background: 'linear-gradient(135deg, rgba(139, 92, 246, 0.1) 0%, rgba(59, 130, 246, 0.1) 100%)',
                border: '1px solid rgba(139, 92, 246, 0.2)',
                padding: '24px',
                borderRadius: '50%',
                boxShadow: '0 0 30px rgba(139, 92, 246, 0.15)',
                display: 'inline-flex',
                animation: 'float 3s infinite ease-in-out'
              }}>
                <Sparkles size={48} color="#a78bfa" />
              </div>
              <style dangerouslySetInnerHTML={{__html: `@keyframes float { 0%, 100% { transform: translateY(0); } 50% { transform: translateY(-8px); } }`}} />
              
              <div style={{ maxWidth: '550px' }}>
                <h3 className="font-outfit" style={{ fontSize: '1.6rem', margin: '0 0 10px 0', fontWeight: '700' }}>
                  Catalytic Viability Simulation Dashboard
                </h3>
                <p style={{ color: '#9ca3af', fontSize: '0.9rem', lineHeight: '1.5', margin: 0 }}>
                  Enter your candidate polymer monomer SMILES, select the target enzyme, and define the chain length. The pipeline will build coordinates, run pocket docking, grow coordinates, filter geometry, and validate stability via OpenMM.
                </p>
              </div>

              {/* Steps Info cards */}
              <div style={{ 
                display: 'grid', 
                gridTemplateColumns: 'repeat(auto-fit, minmax(120px, 1fr))', 
                gap: '12px', 
                width: '100%', 
                maxWidth: '650px', 
                marginTop: '10px' 
              }}>
                {[
                  { phase: "Phase 1", title: "Structure Build", desc: "Builds 3D coordinates & resolves linkages" },
                  { phase: "Phase 2", title: "Anchor Docking", desc: "Docks 3-mer anchor with GNINA/Vina" },
                  { phase: "Phase 3", title: "Geometry Filter", desc: "Verifies scissile attack coordinates" },
                  { phase: "Phase 4", title: "OpenMM MD", desc: "Equilibrates system with GAFF2/Amber forcefield" }
                ].map((s, idx) => (
                  <div key={idx} className="glass-card" style={{ padding: '12px', background: 'rgba(255,255,255,0.02)', textAlign: 'left', border: '1px solid rgba(255,255,255,0.03)' }}>
                    <span style={{ fontSize: '0.65rem', color: '#a78bfa', fontWeight: 'bold', display: 'block', marginBottom: '2px' }}>{s.phase}</span>
                    <span className="font-outfit" style={{ fontSize: '0.8rem', fontWeight: '600', display: 'block', marginBottom: '4px' }}>{s.title}</span>
                    <span style={{ fontSize: '0.7rem', color: '#8b949e', lineHeight: '1.2', display: 'block' }}>{s.desc}</span>
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* Running Console view */}
          {status === 'running' && (
            <div className="glass-card" style={{ display: 'flex', flexDirection: 'column', flex: 1, minHeight: '350px' }}>
              <div style={{ 
                display: 'flex', 
                alignItems: 'center', 
                gap: '8px', 
                padding: '12px 18px', 
                borderBottom: '1px solid rgba(255,255,255,0.06)'
              }}>
                <Terminal size={16} color="#8b5cf6" />
                <span className="font-outfit" style={{ fontSize: '0.85rem', fontWeight: '600' }}>Pipeline Execution Terminal Logs</span>
              </div>
              
              <div className="terminal-view" style={{ flex: 1, maxHeight: '420px', minHeight: '300px' }}>
                {logs.length === 0 ? (
                  <span style={{ color: '#8b949e' }}>Connecting to log stream handler...</span>
                ) : (
                  logs.map((log, index) => (
                    <div key={index} style={{ marginBottom: '4px' }}>{log}</div>
                  ))
                )}
                <div ref={terminalEndRef} />
              </div>
            </div>
          )}

          {/* COMPLETED State Displays 3D Viewer & Report card */}
          {status === 'completed' && results && (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
              
              {/* Summary Notification Bar */}
              <div className="glass-card" style={{ 
                padding: '14px 20px', 
                background: 'rgba(16, 185, 129, 0.1)', 
                border: '1px solid rgba(16, 185, 129, 0.2)', 
                color: '#10b981',
                display: 'flex',
                alignItems: 'center',
                gap: '12px',
                borderRadius: '12px'
              }}>
                <CheckCircle size={20} />
                <div style={{ flex: 1, fontSize: '0.88rem' }}>
                  <strong className="font-outfit">Simulation Succeeded!</strong> Best conformation pose found: <strong>Pose #{results.best_pose_num}</strong>. MD simulation verdict is <strong>{results.md_verdict}</strong>.
                </div>
              </div>

              {/* Main Outputs Grid */}
              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '20px' }}>
                {/* 3D viewer & downloads */}
                <div className="glass-card" style={{ padding: '20px', display: 'flex', flexDirection: 'column', gap: '16px' }}>
                  <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
                    <span className="font-outfit" style={{ fontWeight: '600', fontSize: '1rem' }}>3D Binding Conformation</span>
                    <span style={{ fontSize: '0.75rem', color: '#8b949e' }}>Highlighted magenta: {currentEnzymeData.nucleophile_res_name || 'SER'}{currentEnzymeData.nucleophile_res_num}</span>
                  </div>

                  <MoleculeViewer 
                    taskId={taskId} 
                    nucleophileResNum={currentEnzymeData.nucleophile_res_num}
                    nucleophileResName={currentEnzymeData.nucleophile_res_name}
                  />

                  {/* Downloads Section */}
                  <div style={{ borderTop: '1px solid rgba(255,255,255,0.06)', paddingTop: '16px' }}>
                    <span style={{ fontSize: '0.8rem', color: '#8b949e', fontWeight: '500', display: 'block', marginBottom: '10px' }}>
                      Downloads Center
                    </span>
                    <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '8px' }}>
                      <a href={`/api/files/${taskId}/complex`} download className="btn-primary" style={{ padding: '8px', fontSize: '0.75rem', textDecoration: 'none' }}>
                        <Download size={14} />
                        <span>Complex PDB</span>
                      </a>
                      <a href={`/api/files/${taskId}/trajectory`} download className="btn-primary" style={{ padding: '8px', fontSize: '0.75rem', textDecoration: 'none', background: 'linear-gradient(135deg, #4b5563 0%, #1f2937 100%)', boxShadow: 'none' }}>
                        <Download size={14} />
                        <span>MD Trajectory DCD</span>
                      </a>
                    </div>
                  </div>
                </div>

                {/* Score Summary Metrics & Traffic Lights */}
                <div className="glass-card" style={{ padding: '20px', display: 'flex', flexDirection: 'column', gap: '20px' }}>
                  <span className="font-outfit" style={{ fontWeight: '600', fontSize: '1.05rem', borderBottom: '1px solid rgba(255,255,255,0.06)', paddingBottom: '10px' }}>
                    Simulation Report Card
                  </span>

                  {/* Traffic Light Grid */}
                  <div style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
                    
                    {/* Catalytic Attack Viability */}
                    <div className="glass-card" style={{ padding: '12px 16px', display: 'flex', alignItems: 'center', justifyContent: 'space-between', background: 'rgba(255,255,255,0.02)' }}>
                      <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                        <div className={`traffic-indicator ${results.attack_distance < 4.5 ? 'indicator-green' : 'indicator-red'}`} />
                        <div>
                          <span style={{ fontSize: '0.82rem', fontWeight: '600', display: 'block' }}>Catalytic Proximity Filter</span>
                          <span style={{ fontSize: '0.7rem', color: '#8b949e' }}>Nucleophile distance to carbonyl &lt; 4.5 Å</span>
                        </div>
                      </div>
                      <span className="font-outfit" style={{ fontSize: '0.9rem', fontWeight: 'bold', color: results.attack_distance < 4.5 ? '#10b981' : '#ef4444' }}>
                        {results.attack_distance < 4.5 ? 'PASS' : 'FAIL'}
                      </span>
                    </div>

                    {/* MD Stability */}
                    <div className="glass-card" style={{ padding: '12px 16px', display: 'flex', alignItems: 'center', justifyContent: 'space-between', background: 'rgba(255,255,255,0.02)' }}>
                      <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                        <div className={`traffic-indicator ${results.md_verdict === 'STABLE' ? 'indicator-green' : 'indicator-yellow'}`} />
                        <div>
                          <span style={{ fontSize: '0.82rem', fontWeight: '600', display: 'block' }}>Molecular Dynamics Stability</span>
                          <span style={{ fontSize: '0.7rem', color: '#8b949e' }}>Trajectory equilibration stability check</span>
                        </div>
                      </div>
                      <span className="font-outfit" style={{ fontSize: '0.9rem', fontWeight: 'bold', color: results.md_verdict === 'STABLE' ? '#10b981' : '#f59e0b' }}>
                        {results.md_verdict}
                      </span>
                    </div>

                    {/* OpenMM status */}
                    <div className="glass-card" style={{ padding: '12px 16px', display: 'flex', alignItems: 'center', justifyContent: 'space-between', background: 'rgba(255,255,255,0.02)' }}>
                      <div style={{ display: 'flex', alignItems: 'center', gap: '10px' }}>
                        <div className={`traffic-indicator ${results.is_mock ? 'indicator-yellow' : 'indicator-green'}`} />
                        <div>
                          <span style={{ fontSize: '0.82rem', fontWeight: '600', display: 'block' }}>Physics Solver Backend</span>
                          <span style={{ fontSize: '0.7rem', color: '#8b949e' }}>
                            {results.is_mock ? 'Conda environment fallback solver' : 'GPU OpenMM MD engine'}
                          </span>
                        </div>
                      </div>
                      <div style={{ display: 'flex', flexDirection: 'column', alignItems: 'flex-end' }}>
                        <span className="font-outfit" style={{ fontSize: '0.85rem', fontWeight: 'bold', color: results.is_mock ? '#f59e0b' : '#10b981' }}>
                          {results.is_mock ? 'MOCK DATA' : 'REAL OPENMM'}
                        </span>
                        {results.is_mock && (
                          <div style={{ fontSize: '0.62rem', color: '#f59e0b', maxWidth: '140px', textAlign: 'right', display: 'flex', alignItems: 'center', gap: '2px' }}>
                            <Info size={8} />
                            <span>Python imports missed openff/openmm</span>
                          </div>
                        )}
                      </div>
                    </div>
                  </div>

                  {/* Quantitative metrics grid */}
                  <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '12px' }}>
                    <div className="glass-card" style={{ padding: '12px', textAlign: 'center', background: 'rgba(255,255,255,0.01)' }}>
                      <span style={{ fontSize: '0.72rem', color: '#8b949e', display: 'block', marginBottom: '2px' }}>MM-GBSA Energy</span>
                      <span className="font-outfit" style={{ fontSize: '1.2rem', fontWeight: '800', color: '#a78bfa' }}>
                        {results.score.toFixed(1)} <span style={{ fontSize: '0.75rem', fontWeight: 'normal' }}>kcal/mol</span>
                      </span>
                    </div>

                    <div className="glass-card" style={{ padding: '12px', textAlign: 'center', background: 'rgba(255,255,255,0.01)' }}>
                      <span style={{ fontSize: '0.72rem', color: '#8b949e', display: 'block', marginBottom: '2px' }}>Buried SASA</span>
                      <span className="font-outfit" style={{ fontSize: '1.2rem', fontWeight: '800', color: '#3b82f6' }}>
                        {results.sasa.toFixed(1)} <span style={{ fontSize: '0.75rem', fontWeight: 'normal' }}>Å²</span>
                      </span>
                    </div>

                    <div className="glass-card" style={{ padding: '12px', textAlign: 'center', background: 'rgba(255,255,255,0.01)' }}>
                      <span style={{ fontSize: '0.72rem', color: '#8b949e', display: 'block', marginBottom: '2px' }}>Avg Attack Distance</span>
                      <span className="font-outfit" style={{ fontSize: '1.2rem', fontWeight: '800', color: '#10b981' }}>
                        {results.avg_distance.toFixed(2)} <span style={{ fontSize: '0.75rem', fontWeight: 'normal' }}>Å</span>
                      </span>
                    </div>

                    <div className="glass-card" style={{ padding: '12px', textAlign: 'center', background: 'rgba(255,255,255,0.01)' }}>
                      <span style={{ fontSize: '0.72rem', color: '#8b949e', display: 'block', marginBottom: '2px' }}>Avg Ligand RMSD</span>
                      <span className="font-outfit" style={{ fontSize: '1.2rem', fontWeight: '800', color: '#f59e0b' }}>
                        {results.avg_rmsd.toFixed(2)} <span style={{ fontSize: '0.75rem', fontWeight: 'normal' }}>Å</span>
                      </span>
                    </div>
                  </div>

                  {/* Funnel yield metrics */}
                  <div style={{ fontSize: '0.75rem', color: '#8b949e', display: 'flex', justifyContent: 'space-between', borderTop: '1px solid rgba(255,255,255,0.06)', paddingTop: '12px' }}>
                    <span>Docked: {results.all_poses_count} poses</span>
                    <span>Passed Geometry: {results.passing_poses_count} poses</span>
                    <span>Validated: {n_to_validate} pose{n_to_validate > 1 ? 's' : ''}</span>
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* FAILED State details card */}
          {status === 'failed' && (
            <div style={{ display: 'flex', flexDirection: 'column', gap: '20px' }}>
              <div className="glass-card" style={{ 
                padding: '20px', 
                background: 'rgba(239, 68, 68, 0.1)', 
                border: '1px solid rgba(239, 68, 68, 0.2)', 
                color: '#ef4444',
                borderRadius: '12px'
              }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '10px', marginBottom: '10px' }}>
                  <XCircle size={24} />
                  <h3 className="font-outfit" style={{ margin: 0, fontSize: '1.1rem', fontWeight: '700' }}>Simulation Run Failed</h3>
                </div>
                <p style={{ margin: 0, fontSize: '0.85rem', color: '#f3f4f6', lineHeight: '1.4' }}>
                  {error}
                </p>
              </div>

              {/* Show logs to help trace the failure */}
              <div className="glass-card" style={{ display: 'flex', flexDirection: 'column', minHeight: '200px' }}>
                <div style={{ 
                  display: 'flex', 
                  alignItems: 'center', 
                  gap: '8px', 
                  padding: '12px 18px', 
                  borderBottom: '1px solid rgba(255,255,255,0.06)'
                }}>
                  <Terminal size={16} color="#ef4444" />
                  <span className="font-outfit" style={{ fontSize: '0.85rem', fontWeight: '600' }}>Traceback Terminal Logs</span>
                </div>
                
                <div className="terminal-view" style={{ maxHeight: '300px', flex: 1 }}>
                  {logs.map((log, index) => (
                    <div key={index} style={{ marginBottom: '4px' }}>{log}</div>
                  ))}
                  <div ref={terminalEndRef} />
                </div>
              </div>
            </div>
          )}
        </main>
      </div>
    </div>
  );
}
