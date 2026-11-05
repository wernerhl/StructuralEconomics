# STRUCTURAL ECONOMICS: CONSTRUCTION MANUAL
## Version 3.0 — March 2026

This file governs every chapter. Read it before writing. Follow it without being asked.

---

## THE IRON RULES

### Rule 1: NEVER CUT EXISTING MATH
When upgrading a chapter, start from the existing .tex file. Add to it. Never rewrite from scratch. Every equation, proof, example, appendix, and exercise that exists stays. The math was tested. It works.

### Rule 2: THE TEMPLATE IS ADDITIVE
Three elements are added to every chapter. They do not replace anything:
1. **Visualization panels** (TikZ heatmaps) inserted at the mathematical moment they illustrate
2. **Structuralist case study** (1.5–2 pages, one thinker, developed argument) before Further Reading
3. **Epigraph** in original Spanish/Portuguese with English translation

### Rule 3: NO OPENING STORIES
The chapter opens with the epigraph, then goes straight into math. No salteña vendors, no minibus rides, no altitude metaphors in the opening paragraphs. Stories belong in exercises, where they force formalization. The reader is a PhD student. Respect that.

### Rule 4: STRUCTURALIST CASE STUDY IS AN ARGUMENT, NOT A LIST
The case study develops ONE thinker's argument in full:
- State their thesis in 2–3 sentences
- Reconstruct their argument in 2–3 numbered parts
- Map each part to a specific equation/proposition in the chapter
- Close with a mapping table: | Thinker's claim | Formal object | This chapter |
- Final sentence connects the prose to the proof

**Bad:** "Prebisch said X. Pinto said Y. Marini said Z."
**Good:** "Pinto's 1970 essay has three parts. Part 1 corresponds to Proposition 3.5..."

### Rule 5: VISUALIZATION IS MATHEMATICAL
Heatmap panels are mathematical objects rendered visually. They go at the exact point in the argument where the student needs to SEE the matrix property being discussed. They are not decoration. Each panel has:
- TikZ-rendered heatmap(s) with numerical values in cells where appropriate
- A `\vizpanel{Title}` box environment
- Caption connecting the visual pattern to the theorem/proposition
- Typical placement: after a key worked example or after a key insight box

---

## CHAPTER TEMPLATE

Every chapter follows this skeleton:

```
1. EPIGRAPH (original language + translation)
   - One structuralist thinker per chapter
   - Choose the thinker whose work most directly states the chapter's main result
   - Optional: second epigraph (Dijkstra-style, from math/CS/philosophy)

2. OPENING PARAGRAPH
   - 3–5 sentences. States what the chapter does mathematically.
   - References the previous chapter's result and says what comes next.
   - No narrative, no scene-setting. "This chapter derives X from Y."

3. MAIN BODY (sections and subsections)
   - Definitions → Propositions → Proofs → Examples → Intuition boxes
   - Key Structural Insight boxes at major punchlines
   - Prebisch Connection boxes at results that formalize CEPAL concepts
   - VISUALIZATION PANELS at the mathematical moments they illustrate
   - Worked numerical examples with full computations shown

4. SUMMARY (5 numbered results)

5. EXERCISES (7–9 problems)
   - Mix: computation, proof, data, interpretation
   - Bolivian stories appear HERE as formalization exercises
   - At least one exercise uses real data (INE, ICIO, CEPAL)
   - At least one exercise previews a concept from a later chapter
   - Exercise numbering: 1.1, 1.2, ... (chapter.problem)

6. MATHEMATICAL APPENDICES
   - Self-contained. Skippable. Teach exactly the linear algebra needed.
   - Cumulative across chapters (Ch.1: vectors/matrices, Ch.5: graphs, Ch.6: PF, Ch.10: difference eqs)
   - Include Python code where relevant

7. STRUCTURALIST CASE STUDY (1.5–2 pages)
   - Developed argument, not a list
   - One thinker per chapter (see assignment below)
   - Mapping table at the end

8. FURTHER READING (0.5 page)
   - Four categories: Standard micro/macro, Network economics, Structuralist tradition, Empirical
```

---

## TCOLORBOX ENVIRONMENTS

Every chapter uses these five box types:

```latex
% Blue: mathematical intuition
\newtcolorbox{intuition}[1][]{colback=blue!3, colframe=blue!40!black,
  fonttitle=\bfseries, title={Intuition: #1}, sharp corners, boxrule=0.5pt}

% Orange: key structural insight (chapter punchlines)
\newtcolorbox{keyinsight}[1][]{colback=orange!5, colframe=orange!60!black,
  fonttitle=\bfseries, title={Key Structural Insight}, sharp corners, boxrule=0.5pt}

% Green: Prebisch Connection (results that formalize CEPAL concepts)
\newtcolorbox{prebisch}[1][]{colback=green!3, colframe=green!40!black,
  fonttitle=\bfseries, title={Prebisch Connection: #1}, sharp corners, boxrule=0.5pt}

% Gray: Visualization panels (TikZ heatmaps)
\newtcolorbox{vizpanel}[1][]{colback=gray!4, colframe=gray!50!black,
  fonttitle=\bfseries, title={Visualization: #1}, sharp corners, boxrule=0.5pt}

% Green variant: Case study (used once per chapter, at the end)
% Not a tcolorbox — it's a subsection in the appendices area
```

---

## STRUCTURALIST CASE STUDY ASSIGNMENTS

| Chapter | Thinker | Work | Core mapping |
|---------|---------|------|--------------|
| 1 | Prebisch | 1949 Manifesto | 3-step argument → $W \neq \frac{1}{N}\mathbf{11}'$ |
| 2 | Furtado | *Desenvolvimento e Subdesenvolvimento* 1961 | Underdevelopment as fixed point of $W$ |
| 3 | Pinto | "Heterogeneidad estructural" 1970 | Tripartite structure → block-diagonal $W^P$, small $\Gamma$ |
| 4 | Sunkel & Paz | *El subdesarrollo latinoamericano* 1970 | Catalog of assumptions as properties of $W$ |
| 5 | Dos Santos | "Structure of Dependence" 1970 | Three regimes of $W$: colonial/financial/technological |
| 6 | Prebisch (late) | *Capitalismo Periférico* 1981 | Power as appropriation of surplus → Perron eigenvector |
| 7 | Marini | *Dialéctica de la dependencia* 1973 | Superexploitation → $\pi_{\text{worker}} \ll 1/N$ |
| 8 | Quijano | "Colonialidad del poder" 2000 | Coloniality as $W^I$ structure → epistemic Perron dominance |
| 9 | Furtado (late) | *O Mito do Desenvolvimento* 1974 | Cultural dependency → endogenous preferences |
| 10 | Cardoso & Faletto | *Dependência e Desenvolvimento* 1969 | Situations of dependency as local network topology |
| 11 | — | Empirical chapter, no case study | — |
| 12 | Sen | *Development as Freedom* 1999 | Capabilities under endogenous preferences |
| 13 | Sunkel | "Capitalismo transnacional" 1971 | Integration/disintegration as network rewiring |
| 14 | Cardoso & Faletto | Class alliances as coalitions | Coalition payoffs from $\Delta\boldsymbol{\pi}$ |
| 15 | Tavares | *Da substituição de importações* 1972 | ISI as operator design on $W^P$ |
| 16 | Prebisch + Hirschman | 1949 + 1958 | Global $W^P$: linkages as Katz centrality |
| 17 | Pinto + Fajnzylber | 1970 + 1990 | Domestic $\Gamma$ + casillero vacío |
| 18 | Quijano + Marini | 2000 + 1973 | Race/gender as node attributes constraining $W$ |
| 19 | Sunkel | "Five processes, one process" | Three-scale invariance of $\boldsymbol{\pi}$ |
| 20–23 | Taylor | *Structuralist Macroeconomics* 1983 | SAMs as weight matrices |
| 24 | Gramsci + Skocpol | Hegemony + state crisis | Counter-hegemony as alternative $W^I$ |
| 25 | Fajnzylber | *Industrialización trunca* 1983 | Missing dense subgraph = truncated $W^P$ |
| 26 | Ocampo | Neo-structuralist synthesis 2009 | High $\Gamma$ + high ECI = growth with equity |
| 27 | — | AI chapter, contemporary | Algorithmic $W^I$ |
| 28 | Furtado | Return to 1961 thesis | The general case restated |
| Appendix | All 15 | Full translation table | Prose → equation mapping |

---

## EPIGRAPH LANGUAGE RULES

- Always in original language first (Spanish or Portuguese)
- English translation in \textnormal{[brackets]} below
- Attribution: author, *title in italics*, (year)
- If the thinker wrote in both languages, use the language of the specific work cited
- Non-LAC thinkers (Gramsci, Sen, etc.): English is fine

---

## VISUALIZATION GUIDELINES

### When to insert a heatmap panel:
- After a worked example that computes a matrix (show the matrix as colors)
- After a key insight box that makes a claim about matrix topology (show the topology)
- When comparing two structures (side-by-side heatmaps with "vs." between them)
- When the spectral gap, Perron vector, or Leontief inverse is the punchline

### How to build heatmap panels:
- Use TikZ \fill[blue!XX] rectangles where XX = percentage brightness
- Cell size: \s = 0.3–0.4 for N=8–10, \s = 1.0–1.3 for N=2–3
- White gridlines: draw=white, line width=0.1pt (large N) or 0.3pt (small N)
- Numerical values in cells: \node[font=\tiny, white] for small N
- Labels below each heatmap in \scriptsize or \small \bfseries
- Diagnostic values ($\Gamma$, $L$, $\Var(\pi_i)$) below labels in \scriptsize

### Standard heatmap color mapping:
- Value 0 → blue!0 (black/dark)
- Value 0.5 → blue!50 (medium)
- Value 1.0 → blue!100 (bright)
- Scale linearly. Normalize to [0,1] range of the specific matrix.

---

## EXERCISE DESIGN RULES

### Bolivian stories as exercises:
- The Calle Comercio story → Exercise "From story to matrix" (Ch.1)
- The minibus from El Alto → Exercise "El Alto vs. Zona Sur" (Ch.2)
- Gas well to quinoa bowl → Exercise "Bolivia I-O table" (Ch.1, Ch.3)
- Luisa in Potosí → Exercise "The phone as structural change" (Ch.5)
- The landlord and farmer → Exercise "Bargaining with asymmetric power" (Ch.7)

### Structure of each exercise:
- Bold title. 3–5 sub-parts labeled (a), (b), (c)...
- Part (a): setup/computation
- Part (b): extend or vary a parameter
- Part (c): interpret economically
- Part (d) or (e): connect to a structuralist concept or preview a later chapter
- Data exercises: cite the specific source (INE, ICIO, CEPAL repositorio)

---

## EXISTING INVENTORY

| Component | Status | Pages | File |
|-----------|--------|-------|------|
| Preface | Written | 9 | preface.pdf |
| Visual Prologue | Written | 9 | visual_prologue.pdf |
| Ch. 1 v3 (template applied) | Written | 29 | chapter1_v3_final.pdf |
| Ch. 2 v2 (template applied) | Written | 17 | chapter2_v2.pdf |
| Ch. 3 v2 (template applied) | Written | 24 | chapter3_v2.pdf |
| Ch. 4 | Written (template) | 17 | chapter4.pdf |
| Ch. 5 v2 | Written (template applied) | 28 | chapter5_v2.pdf |
| Ch. 6 v2 | Written (template applied) | 21 | chapter6_v2.pdf |
| Ch. 7 v2 | Written (template applied) | 21 | chapter7_v2.pdf |
| Ch. 8–28 | NOT WRITTEN | — | — |
| Essay series (6) | Written | 82 | various .pdf |
| Structuralist atlas | Written | ~30K words | from deep research |
| Matrix visualizer | Built | — | matrix_explorer.jsx |
| Book covers | Built | — | book_cover.jsx, back_cover.jsx |
| Architecture v2 | Written | 14 | architecture_v2.pdf |

### All template upgrades complete for existing chapters.

### Next to write (in order):
1. Ch. 8 (Belief Manipulation and Epistemic Power — new)
2. Ch. 9–12 (Formation Model — Part III, the intellectual heart)

---

## LATEX CONVENTIONS

- Document class: article, 12pt, a4paper
- Font: mathpazo (Palatino)
- Margins: 1.15in
- Spacing: onehalfspacing
- Header: chapter title (left), page number (right)
- Theorem numbering: Definition/Proposition/Theorem share a counter, reset per section
- Equation numbering: continuous within chapter
- No affiliations in author block
- Disclaimer: "Opinions expressed do not reflect the views of affiliated institutions; all errors are the author's own."
- Email: wernerhl@gmail.com

---

## THE MASTER EQUATION

$$W^{CE} = \frac{1}{N}\mathbf{1}\mathbf{1}'$$

The neoclassical model is this matrix. Structural economics is every other matrix. The spectral gap $\Gamma$ measures how far the world is from the textbook. The Perron eigenvector $\boldsymbol{\pi}$ measures who benefits from the distance.

The structuralists built the cathedral in words. This book provides the engineering drawings. But the blueprint was always theirs.
