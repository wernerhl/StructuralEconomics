# Conventions and Foundational Fixes (authoritative)

This file fixes the book's mathematical conventions. Every chapter follows it.
- Where a chapter contradicts it, the chapter changes.
- Chapter 1 gets a "Conventions" section (`\label{ch01:sec:conventions}`) that states items 1–6 for the reader.
- Chapter 6 states and proves items 3 and 7 (labels given below).

Every formula here has been checked numerically.

---

## 1. Matrices and normalization

**Interaction matrices are row-stochastic.**
- This covers W^I, W^U, any composite, and any generic "W" whose Perron vector is called power.
- They are nonnegative N×N matrices with `W·1 = 1`; `W_ij` is the weight agent i places on agent j.
- If a raw matrix K (counts, flows, overlaps) is not stochastic, define `W = D⁻¹K` with `D = diag(K·1)`, and say so.

**Production uses two different objects. Never confuse them.**
- **Technical coefficients A** (Chapter 1 convention, keep it):
  - `a_ij` = units of good j needed per unit of good i, so row i is i's recipe.
  - Balance: `y = A'y + d`.
  - Leontief inverse: `L = (I − A')⁻¹`, which requires `ρ(A) < 1` (Hawkins–Simon).
  - Multipliers and backward/forward linkages are computed **only** from A (or from a SAM coefficient matrix), never from a row-normalized matrix.
  - For any row-stochastic W, `(I − γW)⁻¹1 = 1/(1−γ)·1`, so all row-sum linkages coincide.
- **Production weight matrix** `W^P := D⁻¹A`: the input-cost shares of each sector. It is row-stochastic. This is the input-share matrix of Acemoglu, Carvalho, Ozdaglar & Tahbaz-Salehi (2012).

**Composites.**
- `W = Σ_k α_k W^k`, with `α_k ≥ 0`, `Σα_k = 1`, and every layer row-stochastic **on the same node set**.
- If the node sets differ (sectors vs. media outlets vs. groups vs. individuals), there is no composite. Say so, and do not compare N_eff or Γ across different node sets without normalizing by N.

## 2. Power: the Perron vector π, and how to measure its concentration

**Definition.** For an irreducible row-stochastic W, π is the unique vector with `π'W = π'`, `π > 0` and `Σπ_i = 1`.

**What π means. State which reading is in force in each chapter.**
1. **Information (DeGroot).** π is exactly the vector of consensus weights: `b* = π'b⁽⁰⁾`. If W is aperiodic as well as irreducible, this is the limit of the dynamics.
2. **Production.** π is the γ→1 limit of the normalized Domar/influence vector:
   - `π' = lim_{γ→1} ((1−γ)/N)·1'(I − γW^P)⁻¹`.
   - For γ < 1, the Acemoglu et al. (2012) influence vector is `(α/N)·1'(I − (1−α)W^P)⁻¹`. That is a Katz/Leontief object, **not** the Perron vector.
3. **Bargaining.** Using π as Nash bargaining weights is a **modelling assumption**. State it as an Assumption, not a result. Note that in strategic network bargaining (Manea 2011; Abreu & Manea 2012) payoffs are generally *not* eigenvector centrality.

**Concentration measures.** These, not Γ, are what "unequal power" means.
- Effective number of agents: `N_eff = 1/Σπ_i²`, which lies in `[1, N]` (macro `\Neff`).
- Variance of power: `Var(π) = (1/N)Σ(π_i − 1/N)² = (1/N)(1/N_eff − 1/N)`.
  - Its maximum is `(N−1)/N²`, at `π = e_1`.
  - That is about 1/N for large N: 0.0475 at N = 20, 0.1875 at N = 4, and about 9×10⁻⁸ at N = 11 million.
- Power Gini `G(π)`.
- Power distance from the benchmark: `d_π := ½‖π − 1/N‖₁`, which lies in `[0, 1 − 1/N]`.

**When is π uniform?** For row-stochastic irreducible W: **π = 1/N ⇔ W is doubly stochastic** (all column sums equal 1). For non-stochastic nonnegative W the condition is "all column sums equal", not "doubly stochastic". Avoid non-stochastic W when π is read as power.

## 3. The spectral gap Γ measures integration (mixing speed), never concentration

**Definition.** `Γ(W) := 1 − |λ₂(W)|`, where λ₂ is the eigenvalue of second-largest modulus (macro `\Gap`). Use this one definition everywhere. Delete `|λ₁| − |λ₂|`, `ρ − |λ₂|`, and Γ of sub-stochastic matrices.

**What it means.**
- Asymptotic speed of consensus or mixing: deviations decay like |λ₂|^t, with half-life `ln 2 / (−ln|λ₂|)`.
- **Γ = 0** if W is reducible (several unit eigenvalues) **or** periodic (an eigenvalue of modulus 1 other than 1, e.g. bipartite stars or even rings). A ring or star *without self-loops* has Γ = 0 for that reason.
- For a *lazy* walk, use `½(I + W)`.

**Proposition: Γ and π are independent.** Label: `ch06:prop:gap-power-independent`.
- **Statement.** For every π in the open simplex and every γ ∈ (0, 1], the matrix `W = γ·1π' + (1−γ)I` is irreducible, aperiodic and row-stochastic, with `π(W) = π` and `Γ(W) = γ`.
- **Proof.** `π'W = γπ' + (1−γ)π' = π'`. The spectrum is {1, (1−γ) with multiplicity N−1}.
- **Consequences.**
  - Any degree of power concentration is compatible with any spectral gap.
  - Examples: the near-dictatorship `1π'` with π₁ = 0.91 has Γ = 1; a doubly stochastic weakly coupled block matrix has uniform π and Γ ≈ 0.02.
  - The path `W_t = (1−t)W + t·1π'` holds π fixed while Γ rises to 1. So "maximize Γ" never costs anyone power.

**Correct uses of Γ.**
- Speed of consensus, convergence and mixing.
- Fragmentation. For reversible W with conductance Φ (Cheeger): `Φ²/2 ≤ 1 − λ₂ ≤ 2Φ`. Weak cross-block coupling (small Φ) therefore implies small Γ.
- **The converse fails:** small Γ can come from periodicity, not fragmentation.

**Never say**
- "small Γ ⇒ concentrated power" or "large Γ ⇒ π near uniform";
- "Γ measures distance from the textbook" or "Γ bounds Var(π)" (that bound is false);
- "star / hub ⇒ small Γ". A hub everyone listens to gives Γ near 1. A pure bipartite star gives Γ = 0 only because it is periodic.

**Time-varying W** (products `W_t ⋯ W_1`):
- Per-period Γ does not control the product. Counterexample: Γ = 0.85 and 0.776, but the product has Γ = 0.042.
- Use the Dobrushin coefficient `τ₁(W) = ½ max_{i,j} Σ_k |W_ik − W_jk|`. It satisfies `|λ₂(W)| ≤ τ₁(W)` and `τ₁(AB) ≤ τ₁(A)τ₁(B)`.
- Weak ergodicity holds if `Σ_t (1 − τ₁(W_t)) = ∞` (Hajnal; Seneta 2006).

## 4. Δ(W) and the benchmark

**Δ(W).** Define `Δ(W) := W − W^CE` (a matrix), with `W^CE = (1/N)11'` (macro `\Wce`). Summarize it by **two** coordinates, which Prop `ch06:prop:gap-power-independent` shows are independent:
- **power distance** `d_π` (or N_eff / N);
- **integration deficit** `1 − Γ`.

**Facts.**
- W^CE is the **unique** row-stochastic W that is rank one with uniform π. Every rank-one row-stochastic W equals `1π'` for some π.
- Do not claim that a scalar "‖Δ(W)‖" determines welfare. Outcomes depend on W *and* on primitives (κ, ρ, μ, endowments, the game) that are not functions of Δ(W).

## 5. W^CE and the welfare theorems: the corrected statement

- **Walrasian price formation is anonymous.**
  - Each agent interacts with the market aggregate.
  - Aggregate demand weights agents by expenditure/wealth shares s.
  - So the implied interaction matrix is `W^mkt = 1s'`: rank one, Γ = 1, π = s.
  - **W^CE = (1/N)11' is the equal-endowment special case.**
- **The First Welfare Theorem** holds for *any* wealth distribution s. It requires local nonsatiation, complete markets, price-taking and **no externalities**; in this book's language, **W^U = 0**.
- **What breaks efficiency vs. what is distributional.**
  - Utility interdependence (`W^U ≠ 0`, i.e. externalities) and missing markets break efficiency.
  - Non-anonymous interaction with unequal π (power) is **distributional**: weighted Nash bargaining is Pareto efficient. The associated loss is a *utilitarian (distributional) loss*, never "Pareto inefficiency".
- **The book's thesis becomes:**
  - The textbook benchmark is anonymous interaction with no externalities and, in its equal-endowment form, `W = W^CE`.
  - Structural economics studies non-anonymous W (structure: Γ, blocks) and the distribution of power π, and how both are produced and changed.
- **Ch2's "doubly stochastic" result** is correct only with `W_ii = 0`: then equal externality weights `c_j` ⇔ equal column sums.
- **Autarky** (W = I) does *not* make the welfare theorems "trivially hold".

## 6. The allocation rule and utility calibration

**Bargaining (transferable utility).**
- Joint surplus S, disagreement payoffs `d_i`, asymmetric Nash with weights π:
  - `x_i = d_i + π_i·(S − Σ_j d_j)`.
- With equal `d_i = d̄` and `Σx = Ω`, define the bargaining surplus `Φ := Ω − N·d̄`. Then:
  - **`x_i = Ω/N + (π_i − 1/N)·Φ`**.
- This is the book's allocation rule. Derive it (Ch7) and cite it (`ch07:prop:allocation-rule`).
- Φ is the surplus over disagreement, **not** 1/κ. One good with endowments as the disagreement point gives Φ = 0, i.e. no bargaining.

**Welfare evaluation.**
- Use concave U. With `U = x − (κ/2)x²`, **every allocation must satisfy κ·x_i < 1** (below satiation). Recalibrate examples that violate this, or use log utility.
- For quadratic U the distributional loss is exact:
  - **`L_power := Σ_i [U(Ω/N) − U(x_i)] = (κN/2)·Φ²·Var(π)`**.

## 7. Perturbation lemma

State it once, in a Chapter 6 appendix, with label `ch06:lem:perturbation`. Cite it everywhere else (Ch7, 8, 13, 14, 15, 23, 24).

**Hypotheses.** W is irreducible and row-stochastic with stationary π. E is a perturbation with `E·1 = 0` that keeps W + tE nonnegative.

1. **Perron vector (no diagonalizability needed):**
   - `dπ'/dt = π'·E·Z`, with `Z := (I − W + 1π')⁻¹` (the fundamental matrix; Kemeny–Snell; Meyer 1980).
   - **The formula `Σ_k u_k'·E·v_1 …` used in the old Ch8/Ch14 appendices is identically zero**, because `E·v_1 = E·1 = 0`.
   - Equivalently, `Σ_{k≥2} [π'E v_k/(1−λ_k)] u_k'`, with biorthogonal normalization `u_k'v_l = δ_kl`.
   - π is most sensitive (power is cheapest to move) when an eigenvalue is close to **+1** (weak coupling). Eigenvalues near −1 (periodicity) do not create sensitivity. So manipulation or escape cost scales *with* `1 − λ₂`, not with 1/Γ.
2. **Eigenvalue.** For a simple eigenvalue λ with right vector v, left vector u and `u'v = 1`:
   - `dλ = u'Ev`;
   - `d|λ| = Re(conj(λ)·dλ)/|λ|`. The old `sign(λ)·dλ` is wrong for complex λ.
   - Γ is differentiable only if λ₂ is simple and no other eigenvalue (other than its conjugate) shares its modulus.
   - For a directional derivative, state the compensation rule that keeps the rows stochastic.
3. **Gradient intuition.** The largest entries of `|u₂|` and `|v₂|` lie deep inside blocks, not on bridge nodes. Raising Γ means "connect the cores of the blocks".

## 8. Results, numbers, quotations (editorial rules)

**Labelling results.**
- **Theorem / Proposition / Lemma / Corollary** only for statements proved in the text or an appendix (a correct sketch is OK if it is labelled "Proof sketch").
- If a statement is false, fix it (add hypotheses, weaken it, or replace it with the correct statement) or delete it.
- If it is plausible but unproved, use **Conjecture**.
- Data claims use **Stylized Fact** / **Empirical Result**, with a source.
- Definitions that were labelled theorems become **Definition**.

**Numbers.**
- Every number in a worked example, table or exercise solution must be reproduced by `code/chNN.py` (numpy only; run as `python3 code/chNN.py`; print every value the text uses).
- A number with no data source gets the `\illustrative` tag in its caption or example title, and must still be internally consistent: π sums to 1, linkages ≥ 1, mean BL = mean FL, Var(π) ≤ (N−1)/N², and so on.
- **Never invent a data source, sample, citation or quotation.** If the text attributes numbers to a dataset or paper that does not contain them, either correct the attribution (only if certain) or relabel the numbers as illustrative and remove the attribution.

**Quotations.**
- Fix an epigraph only when the correct text and source are known with high confidence (see docs/review.md §4.2).
- Otherwise keep it, and add a LaTeX comment `% TODO-VERIFY: <what to check, against which edition/page>`.
- Never fabricate replacement wording.
- A paraphrase must be presented as a paraphrase, not in quotation marks.

**Tone.** Remove straw-man phrasing ("the standard model assumes W^I does not exist") and hyperbole ("this is not a metaphor", "the single most important result"). Use "the benchmark Arrow–Debreu model abstracts from …".

**Rule 3 (no opening stories)** stands. Move opening stories into exercises.

**Summaries** have 5 numbered items. Keep the existing exercises; fix broken ones.

## 9. Notation decisions made during the revision

- **Formation model** (Ch9–10).
  - The formation map is **T**, with anchor θ° and Jacobian J = ρI + μφW^I C P.
  - Ψ_i(ΔW, T) is defined once, as `ch09:def:psi`.
  - Ch11's variance-ratio statistic is **Ψ̂_var**, which is explicitly *not* Ψ.
  - ℳ = μ/(1−ρ) is the "formation multiplier"; (1+Ψ) is the "formation amplifier".
- **Network evolution map** (Ch13–14, 24): **ℰ**. It was renamed from Ψ and then from 𝒯 to avoid clashing with Ψ and T.
- **Welfare totals:** 𝒲 (calligraphic), to avoid clashing with the weight matrix W.
- **Welfare losses** use Ch7's names: L_power, L_externality, L_interaction (= −Cov(c, π)·…, which *offsets* when the powerful are also the most cared-about). L_formation (Ch18) and L_block (Ch2) are separate objects.
- **Production:** technology in growth models is Z (Ch3), so A always means the I-O coefficient matrix.
- The full table is in `frontmatter/notation.tex`, printed after the table of contents.
