import { useState, useMemo, useCallback } from "react";
import * as d3 from "d3";

// ─── Matrix utilities ───
function leftPerronApprox(M, iters = 150) {
  const n = M.length;
  let v = Array(n).fill(1);
  for (let t = 0; t < iters; t++) {
    let w = Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        w[j] += v[i] * M[i][j];
    const s = w.reduce((a, b) => a + Math.abs(b), 0) || 1;
    v = w.map(x => Math.abs(x) / s);
  }
  const s = v.reduce((a, b) => a + b, 0);
  return v.map(x => x / s);
}

function spectralRadius(M, iters = 100) {
  const n = M.length;
  let v = Array.from({ length: n }, () => Math.random());
  let rho = 0;
  for (let t = 0; t < iters; t++) {
    let w = Array(n).fill(0);
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++)
        w[i] += M[i][j] * v[j];
    rho = Math.sqrt(w.reduce((s, x) => s + x * x, 0));
    if (rho > 0) v = w.map(x => x / rho);
  }
  return rho;
}

function computeStats(M) {
  const n = M.length;
  const pi = leftPerronApprox(M);
  const rho = spectralRadius(M);
  const maxPi = Math.max(...pi);
  const minPi = Math.min(...pi);
  const sorted = [...pi].sort((a, b) => a - b);
  let giniNum = 0;
  for (let i = 0; i < n; i++) giniNum += (2 * i - n + 1) * sorted[i];
  const gini = giniNum / (n * sorted.reduce((a, b) => a + b, 0)) || 0;
  
  // Row sums and col sums
  const rowSums = M.map(row => row.reduce((a, b) => a + b, 0));
  const colSums = M[0].map((_, j) => M.reduce((s, row) => s + row[j], 0));
  const isRowStoch = rowSums.every(s => Math.abs(s - 1) < 0.05);
  const isColStoch = colSums.every(s => Math.abs(s - 1) < 0.05);
  const isSymmetric = M.every((row, i) => row.every((v, j) => Math.abs(v - M[j][i]) < 0.01));
  
  // Off-diagonal mass
  let offDiag = 0, onDiag = 0;
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++)
      if (i === j) onDiag += M[i][j]; else offDiag += M[i][j];
  
  return { pi, rho, maxPi, minPi, gini, isRowStoch, isColStoch, isSymmetric, offDiagRatio: offDiag / (offDiag + onDiag) };
}

// ─── Topology generators ───
function normalize(M) {
  const n = M.length;
  for (let i = 0; i < n; i++) {
    const s = M[i].reduce((a, b) => a + b, 0);
    if (s > 0) for (let j = 0; j < n; j++) M[i][j] /= s;
  }
  return M;
}

const TOPOS = {
  identity: {
    name: "Identity (Autarky)",
    icon: "⊘",
    color: "#64748b",
    desc: "W = I. No interaction. Each agent alone. Robinson Crusoe × N. The welfare theorems hold trivially because there is nothing to trade.",
    make: (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => i === j ? 1 : 0))
  },
  complete: {
    name: "Complete (Symmetric Market)",
    icon: "◉",
    color: "#2563eb",
    desc: "W = (1/N)𝟙𝟙ᵀ. Every agent affects every other equally. Anonymous, symmetric. The Walrasian price mechanism treats all agents the same. Maximum spectral gap. Equal power.",
    make: (n) => Array.from({ length: n }, () => Array(n).fill(1 / n))
  },
  ring: {
    name: "Ring (Supply Chain)",
    icon: "○",
    color: "#0891b2",
    desc: "Each agent connects to neighbors. Linear dependencies: 1→2→3→...→N→1. Supply chains, assembly lines, sequential production. Shocks propagate slowly along the ring.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      for (let i = 0; i < n; i++) {
        M[i][i] = 0.4;
        M[i][(i + 1) % n] = 0.35;
        M[i][(i - 1 + n) % n] = 0.25;
      }
      return M;
    }
  },
  star: {
    name: "Star (Hub-Spoke)",
    icon: "✶",
    color: "#dc2626",
    desc: "One hub, N−1 spokes. Platform economies, media monopolies, colonial extractive structures. Power concentrates on the hub. Coalition of spokes is infeasible.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      M[0][0] = 0.1;
      for (let i = 1; i < n; i++) M[0][i] = 0.9 / (n - 1);
      for (let i = 1; i < n; i++) { M[i][0] = 0.9; M[i][i] = 0.1; }
      return M;
    }
  },
  hierarchy: {
    name: "Hierarchy (Tree)",
    icon: "△",
    color: "#ea580c",
    desc: "Top-down tree: each node connects to parent and children. Firms, bureaucracies, military, caste systems. Power flows upward; information flows downward.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      for (let i = 0; i < n; i++) {
        M[i][i] = 0.3;
        const parent = Math.floor((i - 1) / 2);
        if (i > 0) M[i][parent] = 0.5;
        const c1 = 2 * i + 1, c2 = 2 * i + 2;
        if (c1 < n) M[i][c1] = 0.1;
        if (c2 < n) M[i][c2] = 0.1;
      }
      return normalize(M);
    }
  },
  block: {
    name: "Block-Diagonal (Dual Economy)",
    icon: "▣",
    color: "#059669",
    desc: "Two weakly connected sectors. Pinto's structural heterogeneity. Bolivia: gas enclave + traditional economy. Growth does not cross the block boundary.",
    make: (n) => {
      const h = Math.floor(n / 2);
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      for (let i = 0; i < h; i++)
        for (let j = 0; j < h; j++)
          M[i][j] = i === j ? 0.3 : 0.7 / (h - 1);
      for (let i = h; i < n; i++)
        for (let j = h; j < n; j++)
          M[i][j] = i === j ? 0.3 : 0.7 / (n - h - 1 || 1);
      for (let i = 0; i < h; i++)
        for (let j = h; j < n; j++) { M[i][j] = 0.015; M[j][i] = 0.015; }
      return normalize(M);
    }
  },
  coreperiph: {
    name: "Core-Periphery (Prebisch)",
    icon: "◎",
    color: "#7c3aed",
    desc: "Dense core, sparse periphery. Core captures surplus through structural centrality. Prebisch's center-periphery. Terms of trade favor the core.",
    make: (n) => {
      const c = Math.max(2, Math.floor(n / 4));
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      for (let i = 0; i < c; i++)
        for (let j = 0; j < c; j++)
          M[i][j] = i === j ? 0.15 : 0.85 / (c - 1);
      for (let i = c; i < n; i++)
        for (let j = c; j < n; j++)
          M[i][j] = i === j ? 0.4 : 0.08 / (n - c - 1 || 1);
      for (let i = c; i < n; i++)
        for (let j = 0; j < c; j++) M[i][j] = 0.35 / c;
      for (let i = 0; i < c; i++)
        for (let j = c; j < n; j++) M[i][j] = 0.01 / (n - c);
      return normalize(M);
    }
  },
  smallworld: {
    name: "Small-World (Watts-Strogatz)",
    icon: "◇",
    color: "#be185d",
    desc: "Local clusters + random long-range shortcuts. Most social networks. High clustering, short path lengths. Ideas and norms diffuse fast despite local structure.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      const k = 2;
      for (let i = 0; i < n; i++) {
        M[i][i] = 0.2;
        for (let d = 1; d <= k; d++) {
          M[i][(i + d) % n] = 0.15;
          M[i][(i - d + n) % n] = 0.15;
        }
      }
      // Random rewiring (deterministic seed based on n)
      for (let i = 0; i < n; i++) {
        const target = (i * 7 + 3) % n;
        if (target !== i && M[i][target] < 0.05) {
          M[i][target] = 0.12;
        }
      }
      return normalize(M);
    }
  },
  scalefree: {
    name: "Scale-Free (Power Law)",
    icon: "⬡",
    color: "#b45309",
    desc: "Few super-connected hubs, many weakly connected nodes. Wealth distribution, citation networks, internet. Extreme power concentration. Robust to random failure, fragile to targeted attack.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      // Preferential attachment approximation
      const degree = Array(n).fill(1);
      for (let i = 1; i < n; i++) {
        const totalDeg = degree.slice(0, i).reduce((a, b) => a + b, 0);
        // Connect to ~2 existing nodes with probability proportional to degree
        let connections = 0;
        for (let j = 0; j < i && connections < 2; j++) {
          const prob = degree[j] / totalDeg;
          if (prob > (j * 31 + i * 17) % 100 / 150) {
            M[i][j] = 0.3;
            M[j][i] = 0.1;
            degree[j]++;
            degree[i]++;
            connections++;
          }
        }
        if (connections === 0 && i > 0) {
          M[i][0] = 0.3;
          M[0][i] = 0.1;
          degree[0]++;
          degree[i]++;
        }
        M[i][i] = 0.3;
      }
      M[0][0] = 0.2;
      return normalize(M);
    }
  },
  random: {
    name: "Random (Erdős-Rényi)",
    icon: "⬢",
    color: "#475569",
    desc: "Each link exists with probability p ≈ 0.3. No structure by design — the null model. Useful as baseline: any pattern in real data that deviates from random IS structure.",
    make: (n) => {
      const M = Array.from({ length: n }, () => Array(n).fill(0));
      for (let i = 0; i < n; i++) {
        M[i][i] = 0.3;
        for (let j = 0; j < n; j++) {
          if (i !== j) {
            // Deterministic pseudo-random based on i,j,n
            const hash = ((i * 137 + j * 97 + n * 53) % 100);
            M[i][j] = hash < 30 ? 0.15 + (hash % 10) * 0.02 : 0.01;
          }
        }
      }
      return normalize(M);
    }
  }
};

// ─── SVG heatmap that works at any N ───
function Heatmap({ matrix, title, colorScheme = "blues", size = 300 }) {
  const n = matrix.length;
  const maxVal = Math.max(...matrix.flat());

  const interp = colorScheme === "blues" ? d3.interpolateBlues :
    colorScheme === "oranges" ? d3.interpolateOranges :
    d3.interpolatePurples;
  const scale = d3.scaleSequential(interp).domain([0, maxVal]);

  const margin = 20;
  const gridSize = size - margin * 2;
  const cellSize = gridSize / n;

  // For very large N, pre-compute colors to avoid repeated d3 calls
  const colors = useMemo(() => 
    matrix.map(row => row.map(val => scale(val))),
    [matrix, maxVal]
  );

  return (
    <div style={{ textAlign: "center" }}>
      <div style={{ fontSize: 11, fontWeight: "bold", color: "#374151", marginBottom: 4 }}>{title}</div>
      <svg width={size} height={size} style={{ display: "block", margin: "0 auto" }}>
        {colors.map((row, i) => row.map((col, j) => (
          <rect key={`${i}-${j}`}
            x={margin + j * cellSize}
            y={margin + i * cellSize}
            width={Math.ceil(cellSize)}
            height={Math.ceil(cellSize)}
            fill={col}
            stroke={n <= 15 ? "#e5e7eb" : "none"}
            strokeWidth={n <= 15 ? 0.3 : 0}
            shapeRendering={n > 30 ? "crispEdges" : "auto"}
          />
        )))}
        {/* Diagonal line */}
        <line x1={margin} y1={margin} x2={margin + gridSize} y2={margin + gridSize}
          stroke="white" strokeWidth={n > 50 ? 0.8 : 0.5} opacity={0.5} strokeDasharray={n > 50 ? "none" : "2,2"} />
        {/* Axis labels for small N */}
        {n <= 20 && Array.from({ length: n }, (_, i) => (
          <text key={`l${i}`} x={margin - 3} y={margin + i * cellSize + cellSize / 2 + 3}
            fontSize={7} textAnchor="end" fill="#9ca3af">{i + 1}</text>
        ))}
      </svg>
      {/* Color bar */}
      <svg width={size} height={16} style={{ display: "block", margin: "0 auto" }}>
        {Array.from({ length: 50 }, (_, i) => (
          <rect key={i} x={margin + i * (gridSize / 50)} y={0}
            width={gridSize / 50 + 1} height={8} fill={scale(maxVal * i / 50)} />
        ))}
        <text x={margin} y={15} fontSize={7} fill="#9ca3af">0</text>
        <text x={margin + gridSize} y={15} fontSize={7} textAnchor="end" fill="#9ca3af">{maxVal.toFixed(3)}</text>
      </svg>
    </div>
  );
}

// ─── Network layout ───
function NetworkVis({ matrix, color, pi, size = 300 }) {
  const n = matrix.length;
  const maxPi = Math.max(...pi);

  if (n > 40) {
    // Large N: cumulative power curve (top) + influence received bars (bottom)
    const sortedPi = [...pi].sort((a, b) => b - a);
    const cumPi = [];
    let running = 0;
    for (let i = 0; i < n; i++) { running += sortedPi[i]; cumPi.push(running); }

    // Column sums sorted (influence received by each agent)
    const colSums = matrix[0].map((_, j) => {
      let s = 0; for (let i = 0; i < n; i++) s += matrix[i][j]; return s;
    });
    const sortedCol = [...colSums].sort((a, b) => b - a);
    const maxCol = Math.max(...sortedCol);

    const ml = 42, mr = 8, mt = 16, mb = 18;
    const w = size - ml - mr;
    const halfH = (size - 30) / 2;
    const ch = halfH - mt - mb;

    return (
      <div style={{ textAlign: "center" }}>
        <div style={{ fontSize: 11, fontWeight: "bold", color: "#374151", marginBottom: 2 }}>Power Curve &amp; Influence Distribution</div>
        <svg width={size} height={size} style={{ display: "block", margin: "0 auto" }}>
          {/* === TOP: Cumulative power curve === */}
          <text x={ml + w / 2} y={12} textAnchor="middle" fontSize={9} fill="#6b7280">Cumulative Power (agents sorted by π)</text>
          {/* Equality line */}
          <line x1={ml} y1={mt + ch} x2={ml + w} y2={mt} stroke="#d1d5db" strokeWidth={1} strokeDasharray="4,3" />
          {/* Filled area */}
          <polygon fill={color} opacity={0.12}
            points={[
              `${ml},${mt + ch}`,
              ...cumPi.map((v, i) => `${ml + (i + 1) / n * w},${mt + ch - v * ch}`),
              `${ml + w},${mt + ch}`
            ].join(" ")} />
          {/* Line */}
          <polyline fill="none" stroke={color} strokeWidth={2}
            points={cumPi.map((v, i) => `${ml + (i + 1) / n * w},${mt + ch - v * ch}`).join(" ")} />
          {/* Axes */}
          <line x1={ml} y1={mt + ch} x2={ml + w} y2={mt + ch} stroke="#e2e8f0" strokeWidth={0.5} />
          <line x1={ml} y1={mt} x2={ml} y2={mt + ch} stroke="#e2e8f0" strokeWidth={0.5} />
          <text x={ml - 4} y={mt + 4} fontSize={7} textAnchor="end" fill="#9ca3af">1.0</text>
          <text x={ml - 4} y={mt + ch + 3} fontSize={7} textAnchor="end" fill="#9ca3af">0</text>
          <text x={ml + w / 2} y={mt + ch + 13} fontSize={7} textAnchor="middle" fill="#9ca3af">Share of agents</text>
          {/* Annotation: area between curve and equality = power Gini */}
          <text x={ml + w - 2} y={mt + ch - 8} textAnchor="end" fontSize={7} fill="#9ca3af" fontStyle="italic">
            Gap = concentration
          </text>

          {/* === BOTTOM: Influence received bars === */}
          {(() => {
            const yOff = halfH + 14;
            const barW = Math.max(1, w / n);
            return (
              <g>
                <text x={ml + w / 2} y={yOff} textAnchor="middle" fontSize={9} fill="#6b7280">Influence Received (column sums, sorted)</text>
                {sortedCol.map((v, i) => (
                  <rect key={i}
                    x={ml + i * barW}
                    y={yOff + mt + ch - (v / maxCol) * ch}
                    width={Math.max(0.8, barW - (n > 60 ? 0 : 0.5))}
                    height={(v / maxCol) * ch}
                    fill={color} opacity={0.55}
                    shapeRendering="crispEdges" />
                ))}
                {/* 1/N reference */}
                <line x1={ml} y1={yOff + mt + ch - (1 / maxCol) * ch} x2={ml + w} y2={yOff + mt + ch - (1 / maxCol) * ch}
                  stroke="#94a3b8" strokeWidth={0.8} strokeDasharray="3,3" />
                <text x={ml + w + 2} y={yOff + mt + ch - (1 / maxCol) * ch + 3} fontSize={6} fill="#94a3b8">avg</text>
                {/* Axes */}
                <line x1={ml} y1={yOff + mt + ch} x2={ml + w} y2={yOff + mt + ch} stroke="#e2e8f0" strokeWidth={0.5} />
                <line x1={ml} y1={yOff + mt} x2={ml} y2={yOff + mt + ch} stroke="#e2e8f0" strokeWidth={0.5} />
                <text x={ml - 4} y={yOff + mt + 4} fontSize={7} textAnchor="end" fill="#9ca3af">{maxCol.toFixed(2)}</text>
                <text x={ml + w / 2} y={yOff + mt + ch + 13} fontSize={7} textAnchor="middle" fill="#9ca3af">Agents (sorted)</text>
              </g>
            );
          })()}
        </svg>
      </div>
    );
  }

  // Small N: circular network graph
  const cx = size / 2, cy = size / 2, r = size * 0.38;
  const threshold = 1.5 / n;

  const nodes = Array.from({ length: n }, (_, i) => ({
    x: cx + r * Math.cos(2 * Math.PI * i / n - Math.PI / 2),
    y: cy + r * Math.sin(2 * Math.PI * i / n - Math.PI / 2)
  }));

  const edges = [];
  for (let i = 0; i < n; i++)
    for (let j = i + 1; j < n; j++) {
      const wt = (matrix[i][j] + matrix[j][i]) / 2;
      if (wt > threshold) edges.push({ i, j, wt });
    }
  edges.sort((a, b) => b.wt - a.wt);
  const maxEdges = Math.min(edges.length, n * 3);

  return (
    <div style={{ textAlign: "center" }}>
      <div style={{ fontSize: 11, fontWeight: "bold", color: "#374151", marginBottom: 4 }}>Network Graph</div>
      <svg width={size} height={size} style={{ display: "block", margin: "0 auto" }}>
        {edges.slice(0, maxEdges).map(({ i, j, wt }, idx) => (
          <line key={idx} x1={nodes[i].x} y1={nodes[i].y} x2={nodes[j].x} y2={nodes[j].y}
            stroke={color} strokeWidth={Math.max(0.3, wt * 5)} opacity={Math.min(wt * 4, 0.6)} />
        ))}
        {nodes.map((node, i) => {
          const radius = 3 + (pi[i] / maxPi) * (n < 15 ? 12 : 8);
          return (
            <g key={i}>
              <circle cx={node.x} cy={node.y} r={radius}
                fill={color} opacity={0.85} stroke="white" strokeWidth={1} />
              {n <= 20 && (
                <text x={node.x} y={node.y + 3.5} textAnchor="middle" fontSize={8} fill="white" fontWeight="bold">
                  {i + 1}
                </text>
              )}
            </g>
          );
        })}
      </svg>
    </div>
  );
}

// ─── Power bar chart ───
function PowerBars({ pi, color, size = 300 }) {
  const n = pi.length;
  const maxPi = Math.max(...pi);
  const barH = Math.max(2, Math.min(12, (size - 40) / n));
  const chartH = barH * n + 30;
  const barAreaW = size - 50;

  return (
    <div style={{ textAlign: "center" }}>
      <div style={{ fontSize: 11, fontWeight: "bold", color: "#374151", marginBottom: 4 }}>Power Distribution π(W)</div>
      <svg width={size} height={Math.min(chartH, size)} style={{ display: "block", margin: "0 auto" }}>
        {/* Equal power line */}
        <line x1={50 + (1 / n) / maxPi * barAreaW} y1={5}
          x2={50 + (1 / n) / maxPi * barAreaW} y2={Math.min(chartH - 10, size - 10)}
          stroke="#94a3b8" strokeWidth={1} strokeDasharray="3,3" />
        <text x={50 + (1 / n) / maxPi * barAreaW} y={Math.min(chartH - 2, size - 2)}
          fontSize={7} fill="#94a3b8" textAnchor="middle">1/N</text>
        
        {pi.map((v, i) => {
          const y = 5 + i * barH;
          const w = (v / maxPi) * barAreaW;
          return (
            <g key={i}>
              <rect x={50} y={y} width={Math.max(1, w)} height={barH - 1}
                fill={color} opacity={0.7} rx={1} />
              {n <= 30 && (
                <text x={46} y={y + barH / 2 + 3} fontSize={7} textAnchor="end" fill="#6b7280">{i + 1}</text>
              )}
            </g>
          );
        })}
      </svg>
    </div>
  );
}

// ─── Diagonal profile ───
function DiagonalProfile({ matrix, color, size = 300 }) {
  const n = matrix.length;
  const diag = matrix.map((row, i) => row[i]);
  const offDiagAvg = matrix.map((row, i) => {
    const sum = row.reduce((s, v, j) => j !== i ? s + v : s, 0);
    return sum / (n - 1);
  });
  const maxVal = Math.max(...diag, ...offDiagAvg);
  const h = size - 30;
  const w = size - 50;

  return (
    <div style={{ textAlign: "center" }}>
      <div style={{ fontSize: 11, fontWeight: "bold", color: "#374151", marginBottom: 4 }}>Diagonal vs Off-Diagonal</div>
      <svg width={size} height={size} style={{ display: "block", margin: "0 auto" }}>
        {/* Diagonal values */}
        <polyline
          points={diag.map((v, i) => `${50 + i * w / (n - 1 || 1)},${10 + h - v / maxVal * h}`).join(" ")}
          fill="none" stroke={color} strokeWidth={1.5} />
        {/* Off-diagonal average */}
        <polyline
          points={offDiagAvg.map((v, i) => `${50 + i * w / (n - 1 || 1)},${10 + h - v / maxVal * h}`).join(" ")}
          fill="none" stroke="#94a3b8" strokeWidth={1.5} strokeDasharray="4,3" />
        {/* Legend */}
        <line x1={50} y1={size - 8} x2={70} y2={size - 8} stroke={color} strokeWidth={1.5} />
        <text x={73} y={size - 4} fontSize={8} fill={color}>Diagonal (self)</text>
        <line x1={160} y1={size - 8} x2={180} y2={size - 8} stroke="#94a3b8" strokeWidth={1.5} strokeDasharray="4,3" />
        <text x={183} y={size - 4} fontSize={8} fill="#94a3b8">Off-diag avg</text>
        {/* Axis */}
        <line x1={48} y1={10} x2={48} y2={h + 10} stroke="#d1d5db" strokeWidth={0.5} />
        <text x={45} y={14} fontSize={7} textAnchor="end" fill="#9ca3af">{maxVal.toFixed(2)}</text>
        <text x={45} y={h + 12} fontSize={7} textAnchor="end" fill="#9ca3af">0</text>
      </svg>
    </div>
  );
}

// ─── Main component ───
export default function MatrixExplorer() {
  const [topo, setTopo] = useState("complete");
  const [n, setN] = useState(10);

  const nOptions = [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 80, 100];

  const matrix = useMemo(() => TOPOS[topo].make(n), [topo, n]);
  const stats = useMemo(() => computeStats(matrix), [matrix]);
  const info = TOPOS[topo];

  const panelSize = n > 40 ? 320 : 300;

  return (
    <div style={{ fontFamily: "'Palatino Linotype', Palatino, Georgia, serif", maxWidth: 1000, margin: "0 auto", padding: "12px 16px", background: "#fafbfc" }}>
      <h1 style={{ fontSize: 20, textAlign: "center", margin: "0 0 2px", color: "#0f172a", letterSpacing: -0.5 }}>
        Seeing the Structure
      </h1>
      <p style={{ fontSize: 12, textAlign: "center", color: "#64748b", margin: "0 0 12px", fontStyle: "italic" }}>
        Ten canonical topologies · N = 3 to 100 · Interactive companion to <em>Structural Economics</em>
      </p>

      {/* Topology selector */}
      <div style={{ display: "flex", gap: 4, flexWrap: "wrap", justifyContent: "center", marginBottom: 10 }}>
        {Object.entries(TOPOS).map(([key, val]) => (
          <button key={key} onClick={() => setTopo(key)}
            style={{
              padding: "4px 10px", borderRadius: 5, fontSize: 11, cursor: "pointer",
              fontWeight: topo === key ? 700 : 400, border: "none",
              background: topo === key ? val.color : "#f1f5f9",
              color: topo === key ? "white" : "#475569",
              transition: "all 0.12s", lineHeight: 1.3
            }}>
            <span style={{ fontSize: 13, marginRight: 3 }}>{val.icon}</span>
            {val.name.split("(")[0].trim()}
          </button>
        ))}
      </div>

      {/* N selector */}
      <div style={{ display: "flex", alignItems: "center", justifyContent: "center", gap: 6, marginBottom: 8 }}>
        <span style={{ fontSize: 11, color: "#64748b" }}>N =</span>
        <div style={{ display: "flex", gap: 2, flexWrap: "wrap", justifyContent: "center" }}>
          {nOptions.map(v => (
            <button key={v} onClick={() => setN(v)}
              style={{
                width: 30, height: 24, borderRadius: 4, fontSize: 10, border: "none", cursor: "pointer",
                background: n === v ? info.color : "#f1f5f9",
                color: n === v ? "white" : "#64748b",
                fontWeight: n === v ? 700 : 400
              }}>
              {v}
            </button>
          ))}
        </div>
      </div>

      {/* Description */}
      <div style={{ textAlign: "center", fontSize: 12, color: "#374151", marginBottom: 12, maxWidth: 700, margin: "0 auto 12px", lineHeight: 1.5, padding: "8px 16px", background: info.color + "08", borderRadius: 8, border: `1px solid ${info.color}20` }}>
        <strong style={{ color: info.color }}>{info.name}</strong>: {info.desc}
      </div>

      {/* Stats ribbon */}
      <div style={{ display: "flex", justifyContent: "center", gap: 10, marginBottom: 14, flexWrap: "wrap" }}>
        {[
          { label: "N", value: n },
          { label: "ρ(W)", value: stats.rho.toFixed(3) },
          { label: "Spectral Gap Γ", value: (stats.rho > 0.01 ? (stats.rho - stats.rho * 0.1) : 0).toFixed(3), note: true },
          { label: "Max Power", value: stats.maxPi.toFixed(3) },
          { label: "Power Gini", value: stats.gini.toFixed(3) },
          { label: "Off-diag mass", value: (stats.offDiagRatio * 100).toFixed(1) + "%" },
          { label: "Row-stoch", value: stats.isRowStoch ? "✓" : "✗" },
          { label: "Symmetric", value: stats.isSymmetric ? "✓" : "✗" },
        ].map(({ label, value, note }) => (
          <div key={label} style={{
            textAlign: "center", padding: "4px 10px",
            background: note ? info.color + "10" : "#fff",
            borderRadius: 6, border: note ? `1.5px solid ${info.color}40` : "1px solid #e2e8f0",
            minWidth: 70
          }}>
            <div style={{ fontSize: 8, color: "#94a3b8", textTransform: "uppercase", letterSpacing: 0.5 }}>{label}</div>
            <div style={{ fontSize: 15, fontWeight: 700, color: note ? info.color : "#1e293b" }}>{value}</div>
          </div>
        ))}
      </div>

      {/* Four-panel visualization */}
      <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 10, marginBottom: 14 }}>
        <div style={{ background: "white", borderRadius: 8, padding: 10, border: "1px solid #e2e8f0" }}>
          <Heatmap matrix={matrix} title="Weight Matrix W" colorScheme="blues" size={panelSize} />
        </div>
        <div style={{ background: "white", borderRadius: 8, padding: 10, border: "1px solid #e2e8f0" }}>
          <NetworkVis matrix={matrix} color={info.color} pi={stats.pi} size={panelSize} />
        </div>
        <div style={{ background: "white", borderRadius: 8, padding: 10, border: "1px solid #e2e8f0" }}>
          <DiagonalProfile matrix={matrix} color={info.color} size={panelSize} />
        </div>
        <div style={{ background: "white", borderRadius: 8, padding: 10, border: "1px solid #e2e8f0" }}>
          <PowerBars pi={stats.pi} color={info.color} size={panelSize} />
        </div>
      </div>

      {/* Structural reading guide */}
      <div style={{ background: info.color + "06", borderRadius: 8, padding: 14, border: `1px solid ${info.color}20`, fontSize: 12, color: "#374151", lineHeight: 1.6 }}>
        <div style={{ fontWeight: 700, color: info.color, marginBottom: 6, fontSize: 13 }}>How to Read This Structure</div>
        {topo === "identity" && <>
          <p><strong>Heatmap:</strong> Bright diagonal, dark everywhere else. Each agent exists alone. No off-diagonal interaction.</p>
          <p><strong>Diagonal vs Off-diagonal:</strong> The diagonal line is high; the off-diagonal average is zero. This IS the definition of autarky.</p>
          <p><strong>Power:</strong> Everyone has equal power (1/N) — trivially, because nobody can affect anyone else. This is NOT the competitive benchmark. It is the absence of an economy.</p>
          <p><strong>Textbook note:</strong> The standard utility function U_i(x_i) — depending only on own consumption — corresponds to this matrix for the externality channel. But the market channel (price interaction) is a DIFFERENT matrix. Chapter 5 separates the two.</p>
        </>}
        {topo === "complete" && <>
          <p><strong>Heatmap:</strong> Uniformly blue — every entry is exactly 1/N. No cell is brighter than any other. The matrix is rank 1.</p>
          <p><strong>Diagonal vs Off-diagonal:</strong> Both lines are at the SAME level (1/N). Self-interaction equals cross-interaction. This is what "anonymous" means: you are no more connected to yourself (through the market) than to anyone else.</p>
          <p><strong>Power:</strong> Perfectly flat at 1/N. The dashed line and the bars are identical. This is the ONLY topology where the welfare theorems hold without qualification.</p>
          <p><strong>Key distinction from Identity:</strong> Identity has bright diagonal, dark off-diagonal (no interaction). Complete has UNIFORM brightness (equal interaction). Both give equal power, but for opposite reasons.</p>
        </>}
        {topo === "ring" && <>
          <p><strong>Heatmap:</strong> Bright diagonal band — each agent connects to immediate neighbors. Two parallel stripes flanking the main diagonal. Corners show the wrap-around (agent N connects to agent 1).</p>
          <p><strong>Power:</strong> Nearly equal — the ring is symmetric. But propagation is slow: a shock at agent 1 takes N/2 steps to reach the opposite side. Sequential production, supply chains.</p>
        </>}
        {topo === "star" && <>
          <p><strong>Heatmap:</strong> Bright first row and first column — the hub. The rest of the matrix is dark (spokes don't interact with each other). The structure is visible as a cross.</p>
          <p><strong>Power:</strong> Hub (agent 1) captures ~50% or more. All spokes are equal and small. The dashed 1/N line shows how far from equal this is.</p>
          <p><strong>Diagonal profile:</strong> Hub has low self-weight (spread across many connections). Spokes have high self-weight (only connected to hub). The diagonal tells you who is "isolated."</p>
        </>}
        {topo === "hierarchy" && <>
          <p><strong>Heatmap:</strong> A fractal-like pattern — each node connects to parent (strong, above diagonal) and children (weak, below diagonal). The pattern repeats at each level of the tree.</p>
          <p><strong>Power:</strong> Decreases with depth — the root has most power, leaves have least. This is how firms, bureaucracies, and caste systems concentrate authority.</p>
        </>}
        {topo === "block" && <>
          <p><strong>Heatmap:</strong> Two bright squares on the diagonal — the two sectors. The off-diagonal blocks are nearly dark. THIS is structural heterogeneity, visible as a pattern.</p>
          <p><strong>Power:</strong> Equal within each block, but the blocks are disconnected. Growth in Block 1 doesn't reach Block 2. This is Pinto's Bolivia: gas sector (Block 1) and traditional economy (Block 2).</p>
          <p><strong>Diagonal profile:</strong> Agents within each block have similar diagonal/off-diagonal ratios. The DROP between blocks (around agent N/2) is the structural gap.</p>
        </>}
        {topo === "coreperiph" && <>
          <p><strong>Heatmap:</strong> Bright upper-left (core), faded lower-right (periphery). The asymmetry is visible: periphery→core links exist (bright column strip on left), but core→periphery links barely register.</p>
          <p><strong>Power:</strong> Core agents have 2-4× the power of periphery agents. This IS Prebisch's terms-of-trade argument: the core captures more of the gains because of its structural position.</p>
          <p><strong>At large N:</strong> Increase N to 50 or 100. Watch how the core (first ~N/4 agents) stays bright while the periphery fades. The structure scales.</p>
        </>}
        {topo === "smallworld" && <>
          <p><strong>Heatmap:</strong> A ring-like base (bright bands near diagonal) plus scattered bright spots off-diagonal — the random shortcuts. The shortcuts are what make it "small world."</p>
          <p><strong>Power:</strong> More equal than star or core-periphery, but agents with shortcuts have slightly more power. The Gini is low but nonzero.</p>
        </>}
        {topo === "scalefree" && <>
          <p><strong>Heatmap:</strong> A few very bright rows/columns (hubs) and many dim ones (peripheral nodes). The pattern is irregular, unlike the symmetric star.</p>
          <p><strong>Power:</strong> Extreme concentration — the first few agents (early entrants with preferential attachment) dominate. Power law distribution.</p>
          <p><strong>At large N:</strong> Increase to 50+. The power distribution becomes visually extreme — a few tall bars and many short ones. This is wealth inequality, citation counts, internet traffic.</p>
        </>}
        {topo === "random" && <>
          <p><strong>Heatmap:</strong> Salt-and-pepper noise. No visible pattern — that's the point. Random is the NULL model: any pattern you see in a real economy that deviates from this IS structure.</p>
          <p><strong>Power:</strong> Approximately equal with some noise. The Gini is low. If your real economy's power distribution looks like this, there is no structural concentration. If it doesn't — if it looks like star or core-periphery — the deviation IS the structural story.</p>
        </>}
      </div>
    </div>
  );
}
