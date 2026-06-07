import React, { useEffect, useRef, useState } from 'react';

export default function MoleculeViewer({ taskId, nucleophileResNum, nucleophileResName }) {
  const containerRef = useRef(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!taskId) return;

    let viewer = null;
    let isMounted = true;

    const initViewer = async () => {
      try {
        setLoading(true);
        setError(null);

        // Fetch PDB complex
        const response = await fetch(`/api/files/${taskId}/complex`);
        if (!response.ok) {
          throw new Error("Failed to load PDB structure file.");
        }
        const pdbData = await response.text();

        if (!isMounted) return;

        // Check if 3Dmol is available
        if (!window.$3Dmol) {
          throw new Error("3Dmol.js library is not loaded. Please check your internet connection.");
        }

        // Create viewer
        viewer = window.$3Dmol.createViewer(containerRef.current, {
          defaultcolors: window.$3Dmol.rasmolElementColors
        });

        // Add PDB model
        viewer.addModel(pdbData, "pdb");

        // Style protein
        viewer.setStyle({ protein: true }, { cartoon: { color: 'spectrum', opacity: 0.85 } });

        // Style ligand (UNL, LIG, POL, or HETATM)
        viewer.setStyle({ resname: 'UNL' }, { stick: { colorscheme: 'cyanCarbon', radius: 0.25 } });
        viewer.setStyle({ resname: 'LIG' }, { stick: { colorscheme: 'cyanCarbon', radius: 0.25 } });
        viewer.setStyle({ resname: 'POL' }, { stick: { colorscheme: 'cyanCarbon', radius: 0.25 } });

        // Highlight Catalytic Nucleophile Residue
        if (nucleophileResNum) {
          viewer.setStyle(
            { resseq: parseInt(nucleophileResNum) },
            { stick: { colorscheme: 'magentaCarbon', radius: 0.25 }, cartoon: { color: 'magenta' } }
          );
        }

        viewer.zoomTo();
        viewer.render();
        setLoading(false);
      } catch (err) {
        console.error("3Dmol rendering error:", err);
        if (isMounted) {
          setError(err.message);
          setLoading(false);
        }
      }
    };

    initViewer();

    return () => {
      isMounted = false;
      if (viewer) {
        viewer.clear();
      }
    };
  }, [taskId, nucleophileResNum, nucleophileResName]);

  return (
    <div style={{ position: 'relative', width: '100%', height: '450px' }}>
      {loading && (
        <div style={{
          position: 'absolute', top: 0, left: 0, right: 0, bottom: 0,
          background: 'rgba(9, 11, 15, 0.8)', display: 'flex', alignItems: 'center',
          justifyContent: 'center', borderRadius: '12px', zIndex: 10
        }}>
          <div className="font-outfit" style={{ display: 'flex', flexDirection: 'column', alignItems: 'center', gap: '12px' }}>
            <span style={{ fontSize: '1.5rem', fontWeight: 'bold' }}>🧬 Loading 3D structure...</span>
            <span style={{ color: '#8b949e', fontSize: '0.85rem' }}>Fetching PDB data from backend</span>
          </div>
        </div>
      )}
      {error && (
        <div style={{
          position: 'absolute', top: 0, left: 0, right: 0, bottom: 0,
          background: 'rgba(9, 11, 15, 0.95)', display: 'flex', alignItems: 'center',
          justifyContent: 'center', borderRadius: '12px', zIndex: 10, padding: '24px',
          border: '1px solid rgba(239, 68, 68, 0.2)', textAlign: 'center'
        }}>
          <div>
            <span style={{ fontSize: '2.5rem' }}>⚠️</span>
            <h4 className="font-outfit" style={{ color: '#ef4444', margin: '12px 0 6px 0' }}>3D Viewer Error</h4>
            <p style={{ color: '#9ca3af', fontSize: '0.85rem', margin: 0 }}>{error}</p>
          </div>
        </div>
      )}
      <div 
        ref={containerRef} 
        style={{ 
          width: '100%', 
          height: '100%', 
          borderRadius: '12px', 
          overflow: 'hidden', 
          border: '1px solid rgba(255, 255, 255, 0.05)'
        }}
      />
    </div>
  );
}
