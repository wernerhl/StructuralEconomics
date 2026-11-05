# STRUCTURAL ECONOMICS: CONSTRUCTION MANUAL
## Version 4.0 — November 2026

**Read `docs/conventions.md` first.** It is authoritative for notation, the meaning of π and Γ, the allocation rule, and the editorial rules. Where this file and the conventions disagree, the conventions win.

This file governs every chapter. Read it before writing. Follow it without being asked.

---

## THE IRON RULES

### Rule 0: THEOREM MEANS PROVED; NUMBERS ARE REPRODUCIBLE; Γ IS NOT POWER
- Theorem/Proposition/Lemma/Corollary only for statements proved in the text or an appendix. Otherwise: Conjecture, Stylized Fact, Empirical Result, or Definition.
- Every number in the text is printed by `code/chNN.py`. Unsourced numbers carry the `\illustrative` tag and must be internally consistent.
- Γ = 1 − |λ₂| measures integration/mixing speed. Concentration of power is measured by π (N_eff, Var(π), power Gini). Never use one for the other.
- Never invent a source, a dataset, a citation, or a quotation.

### Rule 1: NEVER CUT MATH SILENTLY — AND NEVER KEEP WRONG MATH
When upgrading a chapter, start from the existing .tex file. Do not rewrite from scratch. But every result and every number must be *checked*: worked examples are reproduced by `code/chNN.py`, and a statement that is false is corrected or removed (never kept because "it was there"). Record the reason for any removal in the commit message.

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
- A `\begin{vizpanel}{Title}` box (the title is a mandatory argument)
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

Every chapter uses these box types. They are defined once, in `structecon.sty`; chapters never redefine them.

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

% Gray: Visualization panels (TikZ heatmaps) — MANDATORY title argument
\newtcolorbox{vizpanel}[1]{colback=gray!4, colframe=gray!50!black,
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
| 26 | Ocampo | Neo-structuralist synthesis (2007; Ocampo, Rada & Taylor 2009) | Conjecture: integration (Γ) + capabilities (ECI) + information diversity (N_eff) for growth with equity |
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
- Value 0 → blue!0 (**white** in xcolor, not dark)
- Value 0.5 → blue!50 (medium)
- Value 1.0 → blue!100 (saturated)
- Scale linearly. Normalize to [0,1] range of the specific matrix.
- Cell text: black on cells below blue!50, white above.
- Cell values must be the values of the matrix in the text (generate them with the chapter's `code/chNN.py`). Schematic heatmaps are captioned "schematic", never as data.

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
| Ch. 8–28 | Written (first draft; under revision per docs/review.md) | — | chapters/chapterNN.tex |
| Essay series (6) | Written | 82 | various .pdf |
| Structuralist atlas | Written | ~30K words | from deep research |
| Matrix visualizer | Built | — | matrix_explorer.jsx |
| Book covers | Built | — | book_cover.jsx, back_cover.jsx |
| Architecture v2 | Written | 14 | architecture_v2.pdf |

### Status (November 2026)
All 28 chapters are drafted. The revision follows `docs/review.md` and `docs/conventions.md`; the plan of work is in `docs/review.md` §9.

---

## LATEX CONVENTIONS

- The book is `structural_economics_book.tex` (book class, 12pt, a4paper). Each chapter is a `subfiles` child: it starts with `\documentclass[../structural_economics_book.tex]{subfiles}` and compiles standalone with `make chNN`.
- All packages, boxes, theorem environments and macros live in `structecon.sty`. Chapters have no preamble.
- Labels are prefixed by chapter: `\label{ch06:prop:star}`. Cross-chapter references use `\ref`, never hard-coded numbers.
- Wide tables use `tabularx` with `L` (wrapping) columns.
- Font: mathpazo (Palatino)
- Margins: 1.15in
- Spacing: onehalfspacing
- Header: chapter title (left), page number (right)
- Theorem numbering: all theorem-like environments share one counter, numbered per chapter (Proposition 6.4)
- Equation numbering: continuous within chapter
- No affiliations in author block
- Disclaimer: "Opinions expressed do not reflect the views of affiliated institutions; all errors are the author's own."
- Email: wernerhl@gmail.com

---

## THE MASTER EQUATION

$$W^{CE} = \frac{1}{N}\mathbf{1}\mathbf{1}'$$

The equal-endowment competitive benchmark is this matrix: anonymous (rank one) and symmetric. Walrasian markets are anonymous in general ($W^{mkt} = \mathbf{1}s'$, $s$ = wealth shares); $W^{CE}$ is the equal-wealth case. Structural economics studies everything else, along two independent dimensions:
- **Power** — the Perron vector $\boldsymbol{\pi}$ and its concentration ($N_{\text{eff}}$, $\Var(\pi)$): who benefits.
- **Integration** — the spectral gap $\Gamma = 1 - |\lambda_2|$: how fast the structure mixes, propagates, and reaches consensus.
Neither determines the other (Proposition `ch06:prop:gap-power-independent`).

The structuralists built the cathedral in words. This book provides the engineering drawings. But the blueprint was always theirs.
