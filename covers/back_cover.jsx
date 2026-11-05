import { useState, useEffect, useRef } from "react";

export default function BackCover() {
  const canvasRef = useRef(null);
  const [time, setTime] = useState(0);

  useEffect(() => {
    const interval = setInterval(() => setTime(t => t + 0.015), 50);
    return () => clearInterval(interval);
  }, []);

  // Draw a subtle faded matrix that's clearly block-diagonal — the "solved" state
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const w = canvas.width;
    const h = canvas.height;
    const n = 24;
    const cellW = w / n;
    const cellH = h / n;

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const core = 6;
        const mid = 14;
        let val = 0;
        // Three-block structure: center, semi-periphery, periphery
        if (i < core && j < core) val = 0.07 + Math.sin(time + i * 0.5) * 0.005;
        else if (i >= core && i < mid && j >= core && j < mid) val = 0.04;
        else if (i >= mid && j >= mid) val = 0.025;
        else if (j < core) val = 0.015; // periphery depends on core
        else val = 0.003;

        const normalized = Math.min(val / 0.08, 1);
        const r = Math.floor(10 + normalized * 180);
        const g = Math.floor(8 + normalized * 130);
        const b = Math.floor(35 + (1 - normalized) * 60);
        const alpha = 0.15 + normalized * 0.35;

        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${alpha})`;
        ctx.fillRect(j * cellW, i * cellH, cellW + 0.5, cellH + 0.5);
      }
    }
  }, [time]);

  const lineage = [
    { year: "1950", name: "Prebisch", work: "The Economic Development of Latin America", idea: "Center-periphery. Terms of trade." },
    { year: "1961", name: "Furtado", work: "Development and Underdevelopment", idea: "Underdevelopment as historical production." },
    { year: "1970", name: "Pinto", work: "Structural Heterogeneity", idea: "The dual economy within the periphery." },
    { year: "1970", name: "Sunkel", work: "Transnational Capitalism", idea: "Dependency as endogenous structure." },
    { year: "1983", name: "Taylor", work: "Structuralist Macroeconomics", idea: "The first formalization." },
    { year: "1990", name: "Fajnzylber", work: "Systemic Competitiveness", idea: "Truncated industrialization." },
    { year: "2009", name: "Ocampo", work: "Growth and Policy", idea: "The neo-structuralist synthesis." },
  ];

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
      <div style={{
        width: 440,
        height: 660,
        position: "relative",
        background: "linear-gradient(170deg, #0a0812 0%, #12101e 30%, #0e0c18 70%, #080610 100%)",
        borderRadius: "8px 0 0 8px",
        boxShadow: "-8px 8px 40px rgba(0,0,0,0.7), 2px 0 12px rgba(100,80,40,0.08)",
        overflow: "hidden",
        borderRight: "3px solid #1a1520"
      }}>
        {/* Spine edge */}
        <div style={{
          position: "absolute", right: 0, top: 0, bottom: 0, width: 6,
          background: "linear-gradient(180deg, #2a2030, #1a1520 20%, #2a2030 50%, #1a1520 80%, #2a2030)"
        }} />

        {/* Faded matrix in background — bottom right corner */}
        <canvas
          ref={canvasRef}
          width={240}
          height={240}
          style={{
            position: "absolute",
            bottom: 30,
            right: 30,
            width: 240,
            height: 240,
            opacity: 0.35,
            mixBlendMode: "screen",
            filter: "blur(0.5px)"
          }}
        />
        <div style={{
          position: "absolute", bottom: 30, right: 30, width: 240, height: 240,
          background: "radial-gradient(ellipse at center, transparent 20%, #0a0812 85%)",
          pointerEvents: "none"
        }} />

        {/* Top decoration */}
        <div style={{
          position: "absolute", top: 28, left: 40, right: 40, height: 1,
          background: "linear-gradient(90deg, transparent, #c9a84c40, #c9a84c60, #c9a84c40, transparent)"
        }} />

        {/* Main content */}
        <div style={{ position: "relative", zIndex: 2, padding: "42px 38px 40px" }}>

          {/* Opening quote */}
          <div style={{
            fontSize: 15,
            color: "#d8ceb8",
            lineHeight: 1.65,
            marginBottom: 20,
            fontStyle: "italic",
            textIndent: -8
          }}>
            <span style={{ fontSize: 28, color: "#c9a84c", fontWeight: 300, lineHeight: 0 }}>"</span>
            In 1950, Raúl Prebisch described a global economy organized by structure — 
            where the center captured gains and the periphery transferred them outward 
            through deteriorating terms of trade. He was right. But he wrote before the 
            mathematics of networks existed. This book gives the Latin American 
            structuralist tradition the tools it always deserved: weight matrices, 
            spectral theory, and the Perron-Frobenius theorem.
            <span style={{ fontSize: 28, color: "#c9a84c", fontWeight: 300, lineHeight: 0 }}>"</span>
          </div>

          {/* Separator */}
          <div style={{
            width: 40, height: 1, background: "#c9a84c50", margin: "0 0 18px"
          }} />

          {/* Book description */}
          <div style={{ fontSize: 12.5, color: "#a09888", lineHeight: 1.7, marginBottom: 20 }}>
            <span style={{ color: "#d4c8b0", fontStyle: "italic" }}>Structural Economics</span> shows 
            that the standard neoclassical model is a special case — the case where the 
            interaction matrix is symmetric, rank one, and anonymous. Replace that matrix 
            with the one the world actually has, and power becomes endogenous, preferences 
            become products of the structure, the welfare theorems fail, and the question 
            shifts from <span style={{ color: "#c9a84c", fontStyle: "italic" }}>is the equilibrium efficient?</span> to{" "}
            <span style={{ color: "#c9a84c", fontStyle: "italic" }}>what is the structure, who does it benefit, and can it be changed?</span>
          </div>

          {/* The lineage — the intellectual genealogy */}
          <div style={{ marginBottom: 20 }}>
            <div style={{
              fontSize: 9,
              letterSpacing: 4,
              color: "#c9a84c",
              textTransform: "uppercase",
              marginBottom: 10,
              opacity: 0.7
            }}>
              The Tradition
            </div>

            {lineage.map((entry, idx) => (
              <div key={idx} style={{
                display: "flex",
                alignItems: "baseline",
                marginBottom: 5,
                opacity: 0.65 + idx * 0.02
              }}>
                <span style={{
                  fontSize: 10,
                  color: "#8a7a5a",
                  width: 32,
                  flexShrink: 0,
                  fontVariantNumeric: "tabular-nums"
                }}>
                  {entry.year}
                </span>
                <span style={{
                  fontSize: 11.5,
                  color: "#c9b890",
                  width: 70,
                  flexShrink: 0,
                  fontWeight: 500
                }}>
                  {entry.name}
                </span>
                <span style={{
                  fontSize: 10,
                  color: "#7a6e5e",
                  fontStyle: "italic"
                }}>
                  {entry.idea}
                </span>
              </div>
            ))}

            {/* This book */}
            <div style={{
              display: "flex",
              alignItems: "baseline",
              marginTop: 8,
              paddingTop: 8,
              borderTop: "1px solid #c9a84c20"
            }}>
              <span style={{
                fontSize: 10,
                color: "#c9a84c",
                width: 32,
                flexShrink: 0,
                fontVariantNumeric: "tabular-nums"
              }}>
                2027
              </span>
              <span style={{
                fontSize: 11.5,
                color: "#e8d8b8",
                width: 70,
                flexShrink: 0,
                fontWeight: 600
              }}>
                This book
              </span>
              <span style={{
                fontSize: 10,
                color: "#c9a84c",
                fontStyle: "italic"
              }}>
                The eigenvalues they never had.
              </span>
            </div>
          </div>

          {/* What's inside */}
          <div style={{ fontSize: 11, color: "#8a8070", lineHeight: 1.7, marginBottom: 16 }}>
            Twenty-six chapters unify network economics, general equilibrium, 
            endogenous preference formation, bargaining theory, and the political 
            economy of structural change through a single mathematical object: the 
            weight matrix <span style={{ color: "#d4c8b0", fontStyle: "italic" }}>W</span>. 
            With applications to gender, caste, global trade, and development — 
            and Bolivia as a running case study throughout.
          </div>

          {/* Bottom equation echo */}
          <div style={{
            textAlign: "center",
            marginTop: 12,
            paddingTop: 12,
          }}>
            <span style={{
              fontSize: 11,
              color: "#6a5e4e",
              fontStyle: "italic",
              letterSpacing: 0.5
            }}>
              Newton is not refuted by Einstein. He is nested as the low-velocity limit.
            </span>
          </div>
        </div>

        {/* Bottom decoration */}
        <div style={{
          position: "absolute", bottom: 28, left: 40, right: 40, height: 1,
          background: "linear-gradient(90deg, transparent, #c9a84c40, #c9a84c60, #c9a84c40, transparent)"
        }} />

        {/* ISBN / Publisher placeholder */}
        <div style={{
          position: "absolute",
          bottom: 8,
          left: 38,
          fontSize: 8,
          color: "#4a4238",
          letterSpacing: 0.5
        }}>
          Draft Manuscript — 2027
        </div>

        {/* Barcode placeholder area */}
        <div style={{
          position: "absolute",
          bottom: 38,
          right: 40,
          width: 70,
          height: 42,
          border: "1px solid #2a2520",
          borderRadius: 2,
          display: "flex",
          alignItems: "center",
          justifyContent: "center"
        }}>
          <div style={{ fontSize: 7, color: "#3a342e", textAlign: "center" }}>
            ISBN<br />forthcoming
          </div>
        </div>

        {/* Grain overlay */}
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
