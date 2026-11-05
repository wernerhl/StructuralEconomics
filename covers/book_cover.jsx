import { useState, useEffect, useRef } from "react";

export default function BookCover() {
  const canvasRef = useRef(null);
  const [hovered, setHovered] = useState(false);
  const [time, setTime] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => setTime(t => t + 0.02), 50);
    return () => clearInterval(interval);
  }, []);

  // Draw the animated matrix heatmap
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const w = canvas.width;
    const h = canvas.height;
    const n = 32;
    const cellW = w / n;
    const cellH = h / n;

    // Interpolate between competitive (uniform) and block-diagonal + core-periphery
    const blend = (Math.sin(time * 0.5) + 1) / 2; // 0 to 1

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        // Competitive: uniform 1/n
        const uniform = 1 / n;

        // Structural: block-diagonal with core-periphery
        const core = 8;
        let structural = 0;
        if (i < core && j < core) {
          structural = i === j ? 0.06 : 0.08;
        } else if (i >= core && i < 20 && j >= core && j < 20) {
          structural = i === j ? 0.04 : 0.03;
        } else if (i >= 20 && j >= 20) {
          structural = i === j ? 0.05 : 0.02;
        } else if (i >= core && j < core) {
          structural = 0.025;
        } else if (i < core && j >= core) {
          structural = 0.005;
        } else {
          structural = 0.003 + Math.sin(i * 0.3 + j * 0.2 + time) * 0.001;
        }

        const val = uniform * (1 - blend) + structural * blend;
        const normalized = Math.min(val / 0.09, 1);

        // Color: from deep indigo (low) to amber/gold (high)
        const r = Math.floor(10 + normalized * 220);
        const g = Math.floor(8 + normalized * 160);
        const b = Math.floor(40 + (1 - normalized) * 80 + normalized * 20);
        const alpha = 0.3 + normalized * 0.7;

        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${alpha})`;
        ctx.fillRect(j * cellW, i * cellH, cellW + 0.5, cellH + 0.5);
      }
    }

    // Draw subtle grid lines for the block boundaries
    ctx.strokeStyle = "rgba(255, 200, 100, 0.15)";
    ctx.lineWidth = 1;
    [8, 20].forEach(k => {
      const x = k * cellW;
      const y = k * cellH;
      ctx.globalAlpha = blend * 0.4;
      ctx.beginPath();
      ctx.moveTo(x, 0); ctx.lineTo(x, h);
      ctx.moveTo(0, y); ctx.lineTo(w, y);
      ctx.stroke();
    });
    ctx.globalAlpha = 1;

    // Diagonal glow
    ctx.strokeStyle = `rgba(255, 220, 140, ${0.1 + blend * 0.2})`;
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(w, h);
    ctx.stroke();

  }, [time]);

  // Draw network nodes
  const nodes = [];
  const n = 18;
  const blend = (Math.sin(time * 0.5) + 1) / 2;
  for (let i = 0; i < n; i++) {
    const angle = (2 * Math.PI * i) / n - Math.PI / 2;
    const baseR = 120;
    // In structural mode, core nodes pull inward, periphery nodes push outward
    const isCore = i < 5;
    const rAdj = isCore ? baseR * (1 - blend * 0.3) : baseR * (1 + blend * 0.15);
    nodes.push({
      x: 160 + rAdj * Math.cos(angle),
      y: 160 + rAdj * Math.sin(angle),
      size: isCore ? 3 + blend * 4 : 3 - blend * 1,
      core: isCore
    });
  }

  return (
    <div style={{
      width: "100%",
      minHeight: "100vh",
      display: "flex",
      alignItems: "center",
      justifyContent: "center",
      background: "#08060e",
      padding: 20,
      fontFamily: "'Palatino Linotype', Palatino, 'Book Antiqua', Georgia, serif"
    }}>
      {/* Book */}
      <div
        onMouseEnter={() => setHovered(true)}
        onMouseLeave={() => setHovered(false)}
        style={{
          width: 440,
          height: 660,
          position: "relative",
          background: "linear-gradient(170deg, #0a0812 0%, #12101e 30%, #0e0c18 70%, #080610 100%)",
          borderRadius: "0 8px 8px 0",
          boxShadow: hovered
            ? "12px 12px 60px rgba(0,0,0,0.8), -2px 0 20px rgba(180,140,60,0.1), 0 0 80px rgba(180,140,60,0.05)"
            : "8px 8px 40px rgba(0,0,0,0.7), -2px 0 12px rgba(100,80,40,0.08)",
          overflow: "hidden",
          transition: "box-shadow 0.5s ease",
          cursor: "default",
          borderLeft: "3px solid #1a1520"
        }}
      >
        {/* Spine edge detail */}
        <div style={{
          position: "absolute", left: 0, top: 0, bottom: 0, width: 6,
          background: "linear-gradient(180deg, #2a2030, #1a1520 20%, #2a2030 50%, #1a1520 80%, #2a2030)"
        }} />

        {/* Animated matrix heatmap background */}
        <canvas
          ref={canvasRef}
          width={320}
          height={320}
          style={{
            position: "absolute",
            top: 80,
            left: 60,
            width: 320,
            height: 320,
            opacity: 0.7,
            mixBlendMode: "screen",
            borderRadius: 2,
            filter: "blur(0.5px)"
          }}
        />

        {/* Subtle vignette over matrix */}
        <div style={{
          position: "absolute", top: 80, left: 60, width: 320, height: 320,
          background: "radial-gradient(ellipse at center, transparent 30%, #0a0812 100%)",
          pointerEvents: "none"
        }} />

        {/* Network overlay on the matrix */}
        <svg style={{
          position: "absolute", top: 80, left: 60, width: 320, height: 320,
          opacity: 0.5 + blend * 0.3
        }}>
          {/* Edges */}
          {nodes.map((a, i) => nodes.map((b, j) => {
            if (j <= i) return null;
            const bothCore = a.core && b.core;
            const crossLink = a.core !== b.core;
            const opacity = bothCore ? 0.25 + blend * 0.3 : crossLink ? 0.05 + (1 - blend) * 0.1 : 0.03;
            if (opacity < 0.04) return null;
            return (
              <line key={`${i}-${j}`}
                x1={a.x} y1={a.y} x2={b.x} y2={b.y}
                stroke={bothCore ? "#c9a84c" : "#4a5568"}
                strokeWidth={bothCore ? 0.8 : 0.4}
                opacity={opacity} />
            );
          }))}
          {/* Nodes */}
          {nodes.map((node, i) => (
            <circle key={i} cx={node.x} cy={node.y} r={node.size}
              fill={node.core ? "#d4a843" : "#5a6577"}
              opacity={0.8}
              stroke={node.core ? "#e8c55a" : "#3a4555"}
              strokeWidth={0.5} />
          ))}
        </svg>

        {/* Top line decoration */}
        <div style={{
          position: "absolute", top: 38, left: 50, right: 50,
          height: 1,
          background: "linear-gradient(90deg, transparent, #c9a84c40, #c9a84c80, #c9a84c40, transparent)"
        }} />

        {/* Title block */}
        <div style={{ position: "absolute", top: 44, left: 0, right: 0, textAlign: "center", zIndex: 2 }}>
          <div style={{
            fontSize: 11,
            letterSpacing: 6,
            color: "#c9a84c",
            textTransform: "uppercase",
            fontWeight: 400,
            marginBottom: 2,
            opacity: 0.8
          }}>
            A Textbook
          </div>
        </div>

        {/* Main title area - below the matrix */}
        <div style={{
          position: "absolute",
          bottom: 100,
          left: 30,
          right: 30,
          zIndex: 2
        }}>
          {/* The equation as visual element */}
          <div style={{
            fontSize: 13,
            color: "#8a7a5a",
            fontFamily: "'Palatino Linotype', serif",
            fontStyle: "italic",
            marginBottom: 16,
            letterSpacing: 0.5,
            textAlign: "center",
            opacity: 0.7
          }}>
            The neoclassical model is the special case
          </div>

          <div style={{
            textAlign: "center",
            marginBottom: 20,
            fontFamily: "'Palatino Linotype', serif",
          }}>
            <span style={{
              fontSize: 22,
              color: "#e8dcc8",
              fontStyle: "italic",
              letterSpacing: 1
            }}>
              W
            </span>
            <span style={{ fontSize: 18, color: "#8a7a5a", margin: "0 8px" }}>=</span>
            <span style={{
              fontSize: 16,
              color: "#8a7a5a",
              position: "relative",
              top: -2,
            }}>
              <sup style={{ fontSize: 12 }}>1</sup>
              ⁄
              <sub style={{ fontSize: 12 }}><em>N</em></sub>
            </span>
            <span style={{ fontSize: 22, color: "#e8dcc8", fontWeight: 700, letterSpacing: 1, marginLeft: 4 }}>
              𝟙𝟙
            </span>
            <span style={{ fontSize: 14, color: "#e8dcc8", position: "relative", top: -8 }}>⊤</span>
          </div>

          {/* Title */}
          <h1 style={{
            fontSize: 40,
            fontWeight: 400,
            color: "#f0e8d8",
            letterSpacing: 2,
            lineHeight: 1.1,
            textAlign: "center",
            margin: "0 0 6px",
            textShadow: "0 2px 20px rgba(0,0,0,0.8)"
          }}>
            STRUCTURAL
          </h1>
          <h1 style={{
            fontSize: 40,
            fontWeight: 400,
            color: "#f0e8d8",
            letterSpacing: 2,
            lineHeight: 1.1,
            textAlign: "center",
            margin: "0 0 14px",
            textShadow: "0 2px 20px rgba(0,0,0,0.8)"
          }}>
            ECONOMICS
          </h1>

          {/* Subtitle */}
          <div style={{
            fontSize: 13,
            color: "#c9a84c",
            fontStyle: "italic",
            textAlign: "center",
            letterSpacing: 0.5,
            opacity: 0.85,
            marginBottom: 4
          }}>
            This textbook is the general case.
          </div>
        </div>

        {/* Bottom line decoration */}
        <div style={{
          position: "absolute", bottom: 70, left: 50, right: 50,
          height: 1,
          background: "linear-gradient(90deg, transparent, #c9a84c40, #c9a84c80, #c9a84c40, transparent)"
        }} />

        {/* Author */}
        <div style={{
          position: "absolute",
          bottom: 32,
          left: 0,
          right: 0,
          textAlign: "center",
          zIndex: 2
        }}>
          <div style={{
            fontSize: 14,
            letterSpacing: 3,
            color: "#a09080",
            fontWeight: 400
          }}>
            Werner L. Hernani-Limarino
          </div>
        </div>

        {/* Subtle grain overlay */}
        <div style={{
          position: "absolute", inset: 0,
          backgroundImage: `url("data:image/svg+xml,%3Csvg viewBox='0 0 256 256' xmlns='http://www.w3.org/2000/svg'%3E%3Cfilter id='n'%3E%3CfeTurbulence type='fractalNoise' baseFrequency='0.9' numOctaves='4' stitchTiles='stitch'/%3E%3C/filter%3E%3Crect width='100%25' height='100%25' filter='url(%23n)' opacity='0.03'/%3E%3C/svg%3E")`,
          pointerEvents: "none",
          mixBlendMode: "overlay"
        }} />
      </div>
    </div>
  );
}
