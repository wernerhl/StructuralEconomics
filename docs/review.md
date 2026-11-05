# *Structural Economics: The General Case*: Manuscript Review

**Scope.** I reviewed every file in `structural_economics_source.zip`: `BOOK.md`, `structural_economics_book.tex`, the 28 chapter files, and the 3 cover/explorer `.jsx` files. All chapters were compiled with pdflatex (about 490 pages in total).

**How the numbers were checked.** Every worked numerical example I could reconstruct was recomputed in numpy. Line numbers refer to the `.tex` source as uploaded.

**On the quotation checks.** Some judgements about quotations and attributions come from reviewers' recall plus web checks. Each carries a confidence level and should be verified against the printed source before you change anything.

---

## 0. Executive summary

The book's ambition is real and the core idea is good. Re-expressing allocation, production, learning and preferences as operations on interaction matrices gives students one vocabulary for Prebisch, Hirschman, Pinto, Furtado and modern network economics. The exposition in the early chapters (Ch1's allocation and Leontief material, Ch6's centrality progression, Ch8's DeGroot chapter) is clear and well paced. The structuralist case studies are a distinctive asset.

As it stands, however, the manuscript would not survive a referee. The problems are structural, not cosmetic, and fall into seven groups:

1. **The central diagnostic is wrong.** The spectral gap Γ does not measure "distance from the textbook" or "symmetry of power".
   - Any rank-one matrix `W = 1π'` has Γ = 1 for *every* π, including a near-dictatorship.
   - A doubly stochastic, weakly coupled block matrix has perfectly uniform π and Γ ≈ 0.
   - High Γ is therefore neither necessary nor sufficient for equal power. The claim is load-bearing in 20+ chapters, the master equation in `BOOK.md`, the back cover, and the Ch26 and Ch28 "theorems".
2. **W^CE = (1/N)11' is asserted, not derived.**
   - Ch1's own Cobb-Douglas example implies an interaction matrix `1s'`, where s is the vector of wealth shares. That matrix has Γ = 1 but unequal π, and the welfare theorems still hold.
   - The welfare theorems do not "require" W^CE, double stochasticity or π = 1/N. They require local nonsatiation, complete markets, price-taking and no externalities.
   - The book's statements of what they require also change from chapter to chapter (Ch1 vs Ch2 vs Ch4 vs Ch28).
3. **The allocation engine is broken.**
   - The rule `xᵢ = Ω/N + (πᵢ − 1/N)Φ` is used everywhere but never derived, and it contradicts the book's own Nash-bargaining model.
   - The quadratic utility `U = x − (κ/2)x²` is calibrated past satiation in almost every example (κΩ/N ≥ 1). The "powerful" agent who gets more is often made *worse off*, and several bargaining examples are degenerate.
4. **Many worked examples have wrong numbers.** Recomputation found errors in roughly two-thirds of the chapters, and several reverse their own conclusion.
5. **Most "theorems" are unproved, and many are false.**
   - Chapters 20–28 contain two proof sketches, both invalid.
   - Several stated results have explicit counterexamples, given below.
   - Definitions, narratives and arithmetic evaluations are labelled "Theorem".
6. **The empirical numbers lack provenance.**
   - Nearly every Bolivia, Latin America or global figure is unsourced and presented as an estimate (Γ, π, N_eff, Ψ, λ, the ICIO linkage tables, survey tables).
   - Several are mathematically impossible: forward linkages below 1, π summing to 0.65, Var(π) above its maximum possible value.
   - Values drift between chapters (Bolivia's N_eff is 3.2 / 3.4 / 2.8; the thresholds Γ* and N* conflict).
7. **The scholarship needs tightening.**
   - There is no bibliography and no `\cite` anywhere.
   - Several epigraphs are misquoted, paraphrased or misattributed: Gramsci, Quijano (×2), Sunkel & Paz, Prebisch, Ocampo, Marini, and John Paul II (presented as "social movements").
   - Several thinkers are misread: Skocpol inverted, Marini's superexploitation recast as marginal product, Fajnzylber's *casillero vacío* placed in the wrong cell, Crenshaw's intersectionality reversed.

**The good news.** Most of these trace back to a handful of root causes. Fixing them at the source (§1) repairs large parts of the book at once. My recommended order of work is in §9.

---

## 1. Foundational issues (fix these first; they propagate everywhere)

### 1.1 Γ is not a measure of power symmetry or of distance from W^CE

Verified numerically (N = 10):

| Matrix | Γ = 1 − \|λ₂\| | max πᵢ | What the book's reading predicts |
|---|---|---|---|
| `11'/N` (competitive) | 1.000 | 0.100 | equal power ✓ |
| `1π'`, π₁ = 0.91 (near-dictator) | **1.000** | **0.910** | should be "unequal" ✗ |
| Ring, no self-loops (periodic) | **0.000** | 0.100 | book says ≈0.19 |
| Ring with self-loops | 0.127 | 0.100 | — |
| Star, no self-loops (bipartite) | **0.000** | 0.500 | book says ≈0.05 |
| Star with self-loops | **0.500** | 0.357 | large Γ *and* concentrated π |
| Doubly stochastic weak blocks (ε = 0.01) | **0.020** | 0.100 (uniform) | small Γ *and* equal power |

**What the counterexamples show.**
- *Ch15:* the path `W_t = (1−t)W + t·1π'` holds π **exactly fixed** while Γ rises monotonically to 1. "Maximize Γ" therefore never costs anyone power, so the political-cost results in Ch15 are empty.
- *Ch16:* the chapter's own core–periphery schematic has **Γ ≈ 0.90**, not 0.20. "Everyone imports from the core" is close to rank one, which raises Γ. Frictionless gravity trade is exactly rank one: Γ = 1, with π equal to world export shares.
- *Ch24:* the hegemonic W^I drawn in its own viz panel has Γ = 0.72, larger than the counter-hegemonic one (0.67).
- *Ch27:* a viral hub taking 60% of attention *raises* Γ from 0.44 to 0.78 while N_eff collapses from 6 to 2. Winner-take-all concentration registers as improved "structural health".
- *Ch28:* "neoclassical = Γ = 1" is wrong, since `Γ(1e₁') = 1` with N_eff = 1 (a dictatorship).

**Where the false claim originates.**
- `BOOK.md` ("The spectral gap Γ measures how far the world is from the textbook").
- Ch2 Prop "Spectral Gap Bounds the Power Loss" (`chapter02.tex:464–479`; the bound Var(π) ≤ (1−Γ)/(NΓ) is false, and its Courant–Fischer step needs symmetry).
- Ch6 §"Definition and Meaning" (`chapter06.tex:429–460`, item 3: "Large Γ implies π close to uniform").

**Where it carries weight later.**
- Ch1 App. 1B; Ch3 (Pinto as a spectral property).
- Ch5 case study (colonial regime "Γ ≈ 0"; that W ≈ `1e₁'` has Γ = 1).
- Ch6 (star proposition, viz caption "Γ predicts the shape").
- Ch7 (escape cost ∝ 1/Γ); Ch8 (manipulation cost); Ch13 (equilibrium classification); Ch15 (objective, political cost).
- Ch16 (casillero vacío bound); Ch17 (the dashboard's "Γ alone captures most information").
- Ch18 (assortative mixing ⇒ unequal power); Ch19 ("development = high Γ"); Ch24; Ch25 (truncation theorem); Ch26 (proof of the growth-equity theorem); Ch27; Ch28.

**Recommended repair** (one conceptual change, applied everywhere):

- **Separate two concepts the book currently fuses:**
  - **Integration / mixing / propagation speed.** Use Γ, or better the Dobrushin coefficient τ₁ (which composes well over time-varying products) or conductance/Cheeger bounds (which connect to "bridges" and "fragmentation").
  - **Concentration of power.** Use Var(π), N_eff = 1/Σπᵢ², the power Gini, or ‖π − 1/N‖. This is the true "distance from W^CE" on the power dimension.
- **Redefine Δ(W).** Make it a norm of W − 11'/N, or the pair (‖π − 1/N‖, 1 − Γ), and state plainly that neither coordinate determines the other.
- **Use one definition of Γ.** The book currently has four: `|λ₁| − |λ₂|` (Ch1), `1 − |λ₂|` (most chapters), `ρ − |λ₂|` (Ch6), and `ρ − |λ₂|` on sub-stochastic matrices (Ch3, where it just scales with overall linkage intensity).
- **Account for periodicity.** Either assume aperiodicity (lazy walks) or use 1 − Re λ₂ / 1 − λ₂ where relevant. Ch8 App. 8B shows that perturbation fragility is governed by 1 − λ₂ (λ₂ near +1), not by |λ₂|. Eigenvalues near −1, as in bipartite graphs, do not create fragility.

### 1.2 What is π, and which "power" is meant?

The book uses at least three different objects as one "power vector π":

| Object | Matrix it comes from | Chapters |
|---|---|---|
| Eigenvector centrality | Non-stochastic W^P / composite W | Ch5–7 |
| DeGroot influence weights | Row-stochastic W^I | Ch5, Ch8 |
| Nash bargaining weights | Postulated | Ch7, App. 7B ("π as Nash weights" is asserted, not derived) |

Consequences:
- **Ch6 Prop 3.2 is false as stated.** For general nonnegative irreducible W, "π uniform ⇔ competitive" holds iff all **column sums are equal**. Counterexample: `[[0.2, 0.6], [0.4, 0]]`. Yet most Ch5–7 examples are not row-stochastic.
- **Ch5's composite W has no invariance.** For `α_U W^U + α_P W^P + α_I W^I`, π changes if any layer is rescaled. The α's are unidentified, and allowing W^U < 0 (Ch5 L185) breaks Perron–Frobenius.
- **For a row-normalized balanced SAM, π = y/Σy exactly** (account size). Ch20's "distributional power" G(π) is just the Gini of account totals.
- **For symmetric homophily π ∝ degree.** Ch18's assortative mixing therefore cannot create a power gap with equal group sizes. A gap needs unequal group sizes or an asymmetric W^max.
- **The network-bargaining literature disagrees with π-as-bargaining-weight.** Manea (2011) and Abreu–Manea (2012) show that network bargaining payoffs are *not* eigenvector centrality.

**Fix.**
- Pick one normalization (row-stochastic throughout is cleanest) and state which interpretation of π is in force in each chapter.
- Present "π as bargaining weight" as a modelling postulate, with a paragraph on what microfoundation could deliver it.
- Cite Acemoglu et al. (2012) correctly. Their influence vector is `(α/n)[I − (1−α)W']⁻¹1`, a Leontief/Katz object, not the Perron vector (Ch6 L778).

### 1.3 W^CE and the welfare theorems

- **Ch1 (L603–612, L688–694).** Price effects are not symmetric, and aggregate demand weights agents by wealth share, not 1/N. Derive the implicit W from the model (you get `1s'`) and confront the consequence. Otherwise present W^CE as a *normative benchmark*, not "the matrix the standard model assumes".
- **Ch2 "FWT — Matrix Version" (L139–149) is false.**
  - In a one-good economy, the competitive equilibrium (CE) is the endowment.
  - The proof compares to a *utilitarian* planner, not to Pareto efficiency.
  - Restate it as "the V-utilitarian optimum coincides with the U-utilitarian optimum iff αⱼ are equal", or build a multi-good economy with consumption externalities.
- **Ch2 "Doubly Stochastic Condition" (L151–157) is false** because cⱼ excludes the diagonal. Counterexample: `[[.8,.1,.1],[.1,.5,.4],[.4,.4,.2]]`, which has equal cⱼ and column sums (1.3, 1.0, 0.7). Assume Wᵢᵢ = 0 and it holds.
- **Ch4 "A5 is a consequence of A1–A4" is false** under the chapter's own composite, which does not contain W^market.
- **Ch1's table row "Autarky — welfare theorems trivially hold" (L683) is false.**
- **Settle one precise, proven statement** and use it in Ch1, Ch2, Ch4, Ch6 and Ch28.

### 1.4 The allocation rule and the utility calibration

**Satiation.** With `U = x − (κ/2)x²`, satiation is at x = 1/κ. Examples with κΩ/N ≥ 1 put agents at or past it:
- Ch1 (κ = 1, Ω = N).
- Ch2 and Ch7 (κ = 0.5, Ω/N = 5, which is 2.5× the satiation point).

**Symptoms.**
- **Ch1 "Bargaining with Asymmetric Power" (L370–400).** The only Pareto-efficient allocation is (1, 1, 1) for any λ. On the book's own numbers the "powerful" agent gets U(1.3) = 0.455 < U(1) = 0.5.
  - Suggested recalibration: κ = 1, Ω = 2.4, ω = (1.8, 0.3, 0.3), λ = (0.6, 0.2, 0.2).
  - This gives x ≈ (0.88, 0.76, 0.76), W^B ≈ 1.435, W^FB = 1.44.
- **Ch1 Ex. 1.4** is degenerate for every λ. Use Ω = 1.6 with ω = (1.3, 0.3).
- **Ch7 LQ "fundamental allocation equation" (L222–250, App. 7A) has no valid derivation.** With one good and endowments as the disagreement point, x = ω is the only individually rational allocation. The appendix's μ < 0 then makes some Sᵢ < 0 inside a log.

**Fix.**
- **Use the transferable-utility (TU) landlord–farmer structure the book already has.** With joint surplus S and disagreement payoffs d, asymmetric Nash gives `xᵢ = dᵢ + πᵢS`. With equal d this is exactly `xᵢ = (d + S/N) + (πᵢ − 1/N)S`, i.e. Φ = S.
- **Keep concave U only for welfare evaluation.**
- **Relabel losses.** Weighted Nash is always Pareto efficient, so the "loss" is a *utilitarian* (distributional) loss, not inefficiency. This affects Ch1, Ch2, Ch7 Prop "direction" and the Marini case study.
- **Make Φ explicit.** Φ = 1/κ under the rule, giving `L_power = N·Var(π)/(2κ)`. Then the "for some positive function Φ, Ψ, Ξ" results in Ch2 become falsifiable.

### 1.5 Input-output conventions

The convention drifts:
- Ch1: `y = A'y + d`, with aᵢⱼ meaning "i needs j".
- Ch3: `Y = W^P Y + f`, so W^P plays the role of A′.
- Ch4: W^P = A.
- Ch5 (L249–259): the index is transposed relative to its own definition.
- Ch16: row = buyer, but App. 16B uses row = user, the reverse of the OECD convention.
- Ch17: row-normalized technical coefficients mix different denominators.

**Consequences.**
- **Rasmussen/Hirschman linkages degenerate under row-stochastic W.** `(I − γW)⁻¹1 = 1/(1−γ)·1` for every row-stochastic W, so all forward linkages are identical. Ch15's Ω₃ is constant (8.0) for every feasible ΔW, and its multiplier column cannot come from the stated object.
- **Demand and supply shocks are conflated** (Ch3). Demand shocks propagate upstream via the quantity inverse; productivity shocks propagate downstream via the price dual.

**Fix.**
- Fix one convention in Ch1 (OECD: row = supplier) and use explicit transposes elsewhere.
- Keep the **shares matrix W** (for π and Γ) separate from the **coefficient matrix A** (for multipliers and linkages).
- Never compute linkages from a row-normalized matrix.

### 1.6 One shared perturbation lemma

**The current formulas are wrong.**
- Ch8 App. 8B (L498) and Ch14 App. 14B project ΔW onto v₁ = 1. Row-sum-preserving perturbations satisfy ΔW·1 = 0, so **these formulas are identically zero**. Verified numerically: the formula returns 0 while the finite difference is non-zero.
- Ch15 App. 15A's ∂Γ/∂Wᵢⱼ is wrong for complex λ₂. The correct form is `∂|λ₂| = Re(conj(λ₂)·∂λ₂)/|λ₂|`.
  - Worked case: λ₂ = 0.453 + 0.088i. The book gives 0.1369; the finite difference and the corrected formula both give 0.1450.
  - It also needs biorthogonal normalization `u_k'v_l = δ_kl`, a simple λ₂, and a stated compensation direction.

**The fix.** Write one lemma, put it in Ch6's appendix, and cite it from Ch7, 8, 13, 14, 15, 23 and 24:
- **The Perron vector:** `dπ' = π'·dW·Z` with `Z = (I − W + 1π')⁻¹` (the fundamental matrix; no diagonalizability needed). Equivalently, `Σ_{k≥2} [π'dW v_k/(1−λ_k)] u_k'`.
- **Γ:** the complex-safe derivative of |λ₂|, with hypotheses.

**Consequences once the lemma is correct.**
- **Ch7's "structural trap" has the wrong direction.** Ch7 says ε* ∝ π_h/Γ, i.e. escape is *hard* when Γ is small. In fact π is most sensitive (cheapest to move) when 1 − λ₂ is small, which is what Ch8 says. The trap narrative needs rethinking.
- **Ch8's manipulation-cost formula** `C ∝ 1/Γ^I` has the wrong sign relative to its own intuition and summary. Ch24 inherits it, and its maintenance-cost summary contradicts the formula.
- **Ch15's "bridge-builder" intuition is backwards.** |u₂| and |v₂| are largest *deep inside* blocks and near zero on bridge nodes. The gradient says "connect the cores of the two blocks".

### 1.7 Time-varying W (Ch13, Ch23)

**Per-period Γ(W_t) does not control products of matrices.** Verified counterexample:
- X = [[1,0,0],[0,0,1],[.85,0,.15]] has Γ = 0.850.
- Y = [[0,.95,.05],[0,1,0],[1,0,0]] has Γ = 0.776.
- Their product XY has **Γ = 0.042**.

**Fix.**
- Add a section on inhomogeneous products: Dobrushin τ₁ (submultiplicative), Hajnal/Seneta weak and strong ergodicity, and the absolute-probability sequence as time-varying "power".
- Report the realized-product Γ alongside per-period values.

---

## 2. Proof status and the use of "Theorem"

**What the counts show.**

| Chapter | Results (Prop/Thm/Lem/Cor) | `proof` environments | Notes |
|---|---|---|---|
| 4 | 8 | 0 | |
| 5 | 4 | 0 | |
| 11 | 11 | 0 | |
| 15 | 7 | 0 | |
| 16 | 9 | 0 | |
| 17 | 2 | 0 | |
| 20–24 | 22 | 0 | |
| 25–26 | 9 | 2 | Both are "sketches", and both are invalid. |
| 27–28 | 10 | 0 | |

**Results with explicit counterexamples (false as stated).**
- **Ch2:** FWT-matrix version; doubly-stochastic condition; spectral-gap bound; Transfer Trap part 4.
- **Ch3:** growth-propagation bound. `[[.5,.9],[0,.5]]` has entry 0.9 > ρ = 0.5. A productive 3×3 example gives an average off-diagonal of 56.5 against a "bound" of 12.7.
- **Ch3:** "block-diagonal ⇒ small Γ" under ρ − |λ₂|. diag(0.6, 0.05) has Γ = 0.55.
- **Ch6:**
  - Prop 3.2, which needs row-stochasticity.
  - Prop welfare_revisited.
  - PF Part II item 5, which needs primitivity, not irreducibility.
  - Prop star_power. The hub weight `π_h = a/(a+1−s)` is independent of N, so "Γ → 0" is false.
  - Prop block_power.
- **Ch7:**
  - Prop "nesting" (symmetric Nash ≠ Walras). Counterexample: with Cobb-Douglas utility, endowments (0.9, 0.3) and (0.1, 0.7), Walras gives 0.600 and symmetric Nash gives 0.631.
  - Prop "direction" (weighted Nash is Pareto efficient).
- **Ch8:** wisdom bound (holds only with σ̄ = min σⱼ); perturbation formula ≡ 0; manipulation-cost sign.
- **Ch13:**
  - The integration "iff" fails. Adding within-block weight 0.005 → 0.02 *raises* Γ from 0.05 to 0.08.
  - The Brouwer existence argument fails: π is discontinuous at reducible W, and the map does not land in the row-stochastic set.
  - Path dependence: |λ₂| is only Hölder at defective eigenvalues.
- **Ch14:**
  - The game is constant-sum, because Σπ = 1 makes ΣΔπ ≡ 0.
  - The core is generically empty: v(N) = 0 while complementary coalitions each gain.
  - The grand coalition is indifferent among all ΔW.
  - The equilibrium characterization does not cover the chapter's own baseline.
  - Bernheim–Peleg–Whinston does *not* guarantee existence.
- **Ch15:** optimal edge selection. Under a Frobenius budget the optimum spreads over all entries: gradient step ΔΓ = 0.069 vs best single edge 0.038. The single-edge answer holds for an ℓ₁ budget.
- **Ch16:** casillero-vacío bound; ECI variance claim (ECI is z-scored, so its variance ≡ 1).
- **Ch18:** Prop assortative in general. With within-group 0.05 and cross-group 1, suppressing cross-group edges *raises* Γ. Thm gender_gap is unprovable with symmetric φ.
- **Ch19:** App. 19B eigenvalues. The correct distinct values are 1, 0.763, 0.368, 0; the book gives 0.167 and 0.30 for the last two.
- **Ch20:** "Γ_SAM ≤ Γ_IO". The chapter's own data give 0.40–0.50 vs 0.39.
- **Ch22:** alignment theorem. With rank-one layers, composite π = Σαₖπₖ, so Var ≤ max Varₖ, never above. The example's 0.25 < min(0.46, 0.35, 0.30) is impossible.
- **Ch23:** SCI ∈ [0,1], which fails for Argentina and Venezuela on the chapter's own numbers.
- **Ch24:**
  - The "Hegemonic Cycle" does not cycle: stable fixed points with monotone N_eff give no oscillation.
  - "Pitchfork" should be saddle-node.
  - The counter-hegemony threshold is ill-posed under linear dynamics.
- **Ch25:** Truncation-Γ theorem.
  - Interlacing needs symmetric principal submatrices.
  - The appendix bound runs in the wrong direction.
  - A *dense module* lowers Γ, unlike a *bridge*. Counterexample: removing a dense D raises Γ from 0.22 to 0.80.
- **Ch26:** Growth-Equity theorem (depends on the false Var-Γ bound; see §1.1).
- **Ch27:**
  - The echo-chamber theorem's proof is circular.
  - The "diverse equilibrium" is impossible: a primitive W(b*) forces consensus.
- **Ch28:**
  - "General Case" is a tautology.
  - "Var(π) > 0 iff Δ ≠ 0" is false: any doubly stochastic W ≠ 11'/N has uniform π.

**Recommended convention.**
- Reserve **Theorem / Proposition** for statements proved in the text or in an appendix.
- Use **Definition** for definitional "theorems" (Ch13 co-evolution, Ch19 "five processes", Ch23 decomposition).
- Use **Conjecture** or **Hypothesis** (with a testable prediction) for Ch16 Prebisch–Singer, Ch24 cycle, Ch26 growth-equity.
- Use **Illustration** or **Stylized fact** for Ch14 alliance cycling and Ch15 sequencing.
- Add a short "Status of results" table to each Part opener. Referees will respect honesty here far more than a theorem label they can break in five minutes.

---

## 3. Empirical provenance

**Internally impossible numbers.**

| Location | Problem |
|---|---|
| Ch16 L189–210 (ICIO table) | Forward linkage below 1; mean backward linkage ≠ mean forward linkage; Bolivia and Ethiopia are not separate ICIO economies; "67 economies, 36 industries" matches no release. |
| Ch16 L111 | Column sums: US ≈ 0.15 and Bolivia ≈ 0.001 cannot both hold, since columns of a row-stochastic matrix average 1. |
| Ch17 L116–153 | Forward linkages 0.95 and 0.90; mean backward linkage ≠ mean forward linkage; "Informal" and "Household production" as I-O sectors. |
| Ch17 L266–291 | π sums to 0.65 (stated "≈1.00"); the power Gini of the listed π is 0.24, not 0.52. |
| Ch18 L409–417 | Var(π) "≈0.003 higher" at N = 11M; the maximum is about 9×10⁻⁸. |

**Unsourced numbers presented as estimates.** These cover essentially every Bolivia, Latin America or global value of Γ, π, N_eff, Ψ and λ, plus:
- Ch21's media-overlap results and validation regressions (r = 0.65 "significant at 1%" from 6 groups).
- Ch22's transfer, experiment, referendum and Latinobarómetro tables.
- Ch23's country trajectories.
- Ch25–26's thresholds, PCI and φ̄ values, and Γ* and N*.

**Where the chapter's own table contradicts its conclusion.** Recomputing Ch21's media table gives Γ = 0.74 and N_eff = 5.99 (out of a maximum of 6), i.e. near-information democracy. The text reports Γ = 0.14 and N_eff = 3.4.

**Values that drift between chapters.**
- Bolivia N_eff^I: 3.2 (Ch17/24/26), 3.4 (Ch21/27/28), 2.8 → 3.2 (Ch23).
- Γ*: 0.50 (Ch23 SCI), ≥0.40 (Ch23 viz), 0.30–0.35 (Ch26). N*: 2.5–3.5 (Ch24) vs 8–12 (Ch26).
- China and Chile Γ differ between Ch23 and Ch26.
- Bolivia's structure is "two blocks" in Ch17 but "three disconnected blocks" in Ch25, which would mean Γ = 0.
- λ_ethnicity: 0.30 (Ch18), 0.28 (Ch21), 0.72 (recomputed from Ch21's own table).
- G(π): 0.52 (Ch17) vs 0.58 (Ch20) for the same years.

**Errors that carry across chapters.** Ch18's 0.80 should be 0.57 per capita (or 0.69 as power/population). That makes the intersectional 0.46 wrong (it should be 0.33). Ch19 then uses 0.46 as its within-sector factor, and Ch22 reuses it as "π/(1/N)".

**Bolivia facts that cut against the narrative and need engaging head-on.**
- **Gini.** Bolivia's Gini fell from about 0.60 to about 0.42 during the 2006–14 gas boom, one of the largest declines in the region. Ch3 L262 cites Bolivia as "grew and inequality rose". Ch16 says transfers cannot reduce inequality.
- **Gas share of GDP.** "Gas ≈ 30% of GDP" (Ch3 L181) is likely far above INE's extraction share, roughly 5–8% at peak; verify.
- **Hydrocarbons take.** The "18% → 82%" take (Ch23 L162) applied to two mega-fields under a transitional decree. Law 3058 had already set 50%.
- **Lithium.** Bolivia has the largest lithium *resources*, not reserves (Ch25).
- **Remittances.** Bolivia is a net *recipient* of remittances. Ch22 inverts the flow.

**Recommendation.** Adopt a strict two-tier convention:
- **Data tables** carry a source line: dataset, release, year, aggregation, normalization, code repo.
- **Illustrative tables** are tagged *Illustrative* in the caption, and their numbers must at least be internally consistent.
- Build one canonical `bolivia_parameters.tex` (macros like `\BolNeffI`) so a value can only change in one place.
- Ship a companion code repository that computes every number from public sources (INE, Eora/ICIO, Atlas of Economic Complexity, EH).

---

## 4. Scholarship: references, quotations, attributions

### 4.1 Bibliography
- There is no `\cite` and no bibliography in ~490 pages. Further Reading sections are prose with inconsistent (and sometimes wrong) venues.
- **Fix:** move to `biblatex` (authoryear) with one shared `.bib`, and cite inline wherever a result or fact is borrowed.
- **Venue and date errors found:**
  - Kern & Hainmueller 2009 is *Political Analysis*.
  - Hausmann & Hidalgo 2011 is *J. Econ. Growth*.
  - Bursztyn & Cantoni 2016 is *REStat*.
  - Banerjee et al. 2019 ("Using Gossips") is *REStud*.
  - Lang & Lehmann 2012 is *JEL*.
  - Abreu–Manea 2012 is *GEB*.
  - Boyd, Diaconis & Xiao 2004 is *SIAM Review*, not "Boyd 2006".
  - Ocampo's 2014 chapter is in Stiglitz & Heymann (eds.), *Life After Debt*.
  - Ocampo's bonanza article is *Revista CEPAL* 93 (2007), not 2009.

### 4.2 Epigraphs that need correction (verify each against the printed source and add page numbers)

| Ch | Issue | Confidence |
|---|---|---|
| 1 | Prebisch translation truncated; add the full clause. Dijkstra source: "The Humble Programmer" (1972). | high |
| 2 | Furtado: drops "em si"; check the edition (ch. 4, ~p. 191 1st ed.). | medium-high |
| 3 | The Prebisch Spanish looks like a back-translation of the 1950 UN English. **Delete "Che boludo, it's obvious." — What Prebisch would have said**: fabricated attribution, vulgar, inappropriate for a PhD text. | medium-high |
| 4 | Sunkel & Paz appears to be a secondary-source paraphrase. | medium |
| 5 | Dos Santos: the English is his 1970 AER wording; the Spanish looks back-translated, and the title "La estructura de la dependencia" appears invented. | medium |
| 6 | Prebisch (1981): reads like a paraphrase; give a page. | low-medium |
| 7 | Marini: concerns *fuerza de trabajo* (labour-power) paid below its value. A Furtado English paraphrase violates the BOOK.md rule. | medium |
| 8 | Quijano: the text is the opening of "Colonialidad del poder y clasificación social" (*JWSR* 6(2), 2000), not the Lander essay, and the wording is altered. | medium-high |
| 15 | Tavares: first published 1964 in Spanish (CEPAL *Boletín* IX(1)); cite "(1964; in Tavares 1972)". | medium-high |
| 16 | Prebisch epigraph has no title (BOOK.md requires one). | medium-low |
| 18 | The Marini "principio fundamental" line is from "Subdesarrollo y revolución" (1968/69), not *Dialéctica*. The Quijano line is unverified. | medium |
| 19 | Sunkel & Paz "five processes" sentence not found; probably not verbatim. | medium |
| 20 | Taylor "three pillars" unverifiable; likely a paraphrase. | medium |
| 21 | "Attributed to various Latin American media critics". Replace with a sourced quote, e.g. Luis Ramiro Beltrán, *Comunicación dominada* (1980). | high |
| 22 | **This is John Paul II, *Sollicitudo Rei Socialis* §38 (1987)**, not "Latin American social movements", and text is omitted. | high |
| 24 | **Gramsci "monsters" is a popular misquote.** Q3 §34: "La crisi consiste appunto nel fatto che il vecchio muore e il nuovo non può nascere: in questo interregno si verificano i fenomeni morbosi più svariati." Also recurs at L394. | high |
| 25 | Fajnzylber: likely not verbatim (his term is "núcleo endógeno de dinamización tecnológica"). | medium |
| 26 | Ocampo: wrong year (2007), and the sentence reads like CEPAL 1990. Consider Ocampo, Rada & Taylor (2009), *Growth and Policy in Developing Countries*. | medium-high |
| 27 | A self-quote in English with no thinker. | — |
| 28 | Furtado: possibly *Teoria e Política* (1967) rather than 1961. | medium |

### 4.3 Thinkers and literature misrepresented

- **Prebisch.**
  - "Unequal exchange" is Emmanuel (1969). Prebisch's mechanism is terms-of-trade deterioration via asymmetric wage and price behaviour over the cycle plus income elasticities, not weak linkages (Ch3, Ch7, Ch16).
  - The 1949 study was written at CEPAL in Santiago.
- **Hirschman.**
  - He was the *critic* of the big push (Rosenstein-Rodan 1943), not its author (Ch15 L84).
  - He identified mining enclaves as *low*-linkage, so "Hirschman says invest in gas" is backwards.
  - His measures were Chenery–Watanabe direct coefficients; Leontief-inverse indices are Rasmussen (1956).
- **Marini.** Superexploitation is labour-power paid below its *value* (intensified work, longer days), not below marginal product (Ch7, Ch18). His main critics were Cardoso & Serra (1978).
- **Skocpol** is explicitly structural and non-voluntarist (international pressure breaking agrarian-bureaucratic states, plus agrarian structures enabling peasant revolt). An "alliance of elite defectors and mass movements driven by information diversity" inverts her (Ch24). Bolivia 1952 fits her scope; 2003 does not.
- **Gramsci.** Missing: consent "armoured by coercion", the historic bloc, and passive revolution (the most relevant Gramscian concept for Latin America). Organic intellectuals ≠ "high-π^I nodes".
- **Fajnzylber.**
  - The *casillero vacío* is the growth-*and*-equity cell, which no Latin American country occupied in 1965–84. Ch16 and Ch17 each put it somewhere different (three different regions in total).
  - "Spurious vs genuine competitiveness" is from 1988, not 1983.
  - Imitative consumption and transnational-corporation leadership are missing.
  - He appears at length in Ch16, Ch17 and Ch25; keep him in one chapter.
- **Ocampo.** His pillars (structural dynamics, countercyclical macro and the balance-of-payments constraint, equality) do not include "information democracy". Ch26 has three different triples all labelled "the three pillars".
- **Sunkel.** "Five processes, one process" and transnational integration / national disintegration are Sunkel 1971, not Sunkel & Paz 1970 (Ch4). Ch19 says he wrote 1970 "in Chilean exile after 1973".
- **Furtado.** He held a Sorbonne chair from 1965 to 1985. "Could not give him tenure" (Ch28 L352) is false. *O Longo Amanhecer* is not his autobiography. SUDENE dates from 1959.
- **Cardoso & Faletto.** The four-alliance typology is the author's. Map their actual categories (enclave vs nationally controlled; associated-dependent development). Cardoso launched the Plano Real as Finance Minister, not President.
- **Taylor.** "Neoclassical closure ≡ doubly stochastic W" has no basis; closures are about which variable is exogenous. Stone, Pyatt–Thorbecke and Pyatt–Round preceded his SAM work. A structuralist CGE *has* been built (Taylor 1990), contrary to Ch28 L272.
- **Crenshaw.** A product of independent λₖ is additive in logs, which is exactly the *no-interaction* case. Intersectionality is the interaction term `φ − Πφₖ`, which is what Ch18 App. 18B describes. Move it into the main text.
- **Golub & Jackson (2010)** do not show that DeGroot is a rational-learning limit. Their result is wisdom iff max πᵢ → 0. The foundations are DeMarzo–Vayanos–Zwiebel (2003) and Molavi–Tahbaz-Salehi–Jadbabaie (2018).
- **Banerjee–Chandrasekhar–Duflo–Jackson (2013).**
  - Villages average about 200 households, not "10–30".
  - The λ_caste and N_eff figures attributed to the paper are not in it.
  - The finding concerns injection-point centrality.
- **Henrich et al. / Cardenas (2003).** Neither varied recipient ethnicity. Cardenas is about rural Colombia, and the Tsimane are Amazonian, not ayllu.

### 4.4 Literature the book must engage (it is your strongest ally, not your opponent)

The dichotomy "neoclassical = 11'/N vs structural = general W" will be read as a straw man. Production-network macro and network learning already run on arbitrary W.

**Production networks.**
- Long–Plosser (1983).
- Acemoglu, Carvalho, Ozdaglar & Tahbaz-Salehi (2012).
- Baqaee–Farhi (2019).
- Oberfield (2018); Acemoglu–Azar (2020), on endogenous input-output structure.

**Industrial policy.**
- Liu (2019, *QJE*), "Industrial Policies in Production Networks". This is essentially the rigorous version of Ch15.
- Juhász, Lane & Rodrik (2024).
- Lane (2025, Korea's heavy-and-chemical-industry drive).
- Hausmann–Hwang–Rodrik (2007); Rodrik (2008); Hausmann–Rodrik–Sabel (2008).

**Optimizing the spectral gap.** Boyd, Diaconis & Xiao (2004). Fastest-mixing Markov chain is convex (a semidefinite program) for symmetric W, and nonconvex otherwise. This gives Ch15 a correct optimization core.

**Network learning.**
- DeGroot (1974); DeMarzo–Vayanos–Zwiebel (2003).
- Golub–Jackson (2010, 2012). The 2012 *QJE* paper links homophily to learning speed via the spectral gap, which is directly Ch18's topic.
- Chandrasekhar–Larreguy–Xandri (2020).

**Network bargaining.** Manea (2011); Abreu–Manea (2012).

**Measurement and identification.**
- Manski (1993).
- Bramoullé–Djebbari–Fortin (2009).
- Chandrasekhar & Lewis (2016).
- Breza et al. (2020, aggregated relational data).
- Comola & Fafchamps (2014).
- de Paula (2017).
- Shalizi & Thomas (2011).
- Griffith (2022) on censored name generators.

**Platforms and polarization.**
- Bakshy, Messing & Adamic (2015).
- The 2023 Meta election studies (Guess et al.; Nyhan et al.; González-Bailón et al.).
- Huszár et al. (2022).
- Global South: WhatsApp in Brazil (Machado et al.; Resende et al.); Badrinathan (2021); the UN Myanmar mission; Bolivia 2019 (DFRLab and the CLS Strategies takedown).
- WhatsApp is not a recommender system. This is load-bearing for Ch27.

**Time-varying Markov chains.** Seneta (2006); Hajnal (1958); Touri–Nedić.

**Opinion dynamics.** Hegselmann–Krause (2002); Deffuant et al. (2000), on bounded confidence (Ch27 App. B).

**Perturbation theory for π.** Meyer (1980); Cho–Meyer (2001); Schweitzer (1968).

**Inequality decomposition.** Vargas & Garriga (2015, IMF) attribute Bolivia's Gini fall mostly to labour income.

---

## 5. Production, LaTeX and build

1. **The master file does not include the chapters.**
   - `structural_economics_book.tex` contains only placeholders (`\textit{29 pages...}`), and each chapter is a standalone `article` with its own ~70-line duplicated preamble. The preambles have already drifted: `structroots` is defined only in Ch1 and never used, and the box definitions are formatted differently.
   - **Fix:** create `preamble.sty` (packages, theorem environments, boxes, notation macros). Convert each chapter to a body file starting with `\chapter{...}`, demote `\section → \section` (they're already sections under a book `\chapter`), drop `\maketitle`/`\title`, and `\include{chapters/chapterNN}` from the master.
   - For standalone drafts, use the `subfiles` package.
   - Number theorems per chapter (`\newtheorem{definition}{Definition}[chapter]`) so "Proposition 6.3.2" is unambiguous.
2. **Every visualization panel loses its title (all 34).**
   - `vizpanel` is defined as `\newtcolorbox{vizpanel}[1][]{...title={Visualization: #1}}`, with an optional `[...]` argument, but every panel is called as `\begin{vizpanel}{Title}`. The title bar shows "Visualization:" and the intended title appears as body text.
   - The bug originates in BOOK.md's template.
   - **Fix:** `\newtcolorbox{vizpanel}[1]{...}` (mandatory argument).
3. **Hard errors** (the PDF still builds in nonstopmode): `\boldsymbol` inside `\section`/`\subsection*` titles with hyperref (Ch4 L237, Ch14 L528, Ch21 L231). Use `\texorpdfstring{$\boldsymbol{\pi}$}{pi}`.
4. **Undefined references:** Ch7 (`eq:fundamental`, `prop:LQ_allocation`, which should be `eq:LQ_equal_omega` / `ex:LQ_bargaining`), Ch12 (`eq:gap_bound`).
5. **Hard-coded cross-chapter numbers.** There are 376 "Chapter N" mentions, plus many hard-coded "Proposition 6.6" style references. Many are wrong because numbering resets per section:
   - Ch7's "Prop 6.6" renders as 3.6.
   - Ch4's references to Ch2 Props 2.4/2.3/2.2 should be 2.3/2.1/1.4.
   - Ch3's "Ch1 Prop 1B.1" does not exist.
   - Use `\ref` once the chapters live in one document. Duplicate labels (`sec:summary`, `sec:exercises`, `sec:properties`, `sec:bolivia`, `eq:consensus`) will collide on merge, so prefix labels with the chapter (`ch06:prop:star`).
6. **280 overfull boxes.**
   - The worst are the case-study mapping tables, which run 60–150pt (2–5 cm) into the margin.
   - **Ch28's "Book in One Page" table overflows the page by 259pt (≈9 cm), so its bottom rows don't print.**
   - **Fix:** `tabularx` with `X` columns, or `p{}` columns, for every mapping table; `longtable` for Ch28.
7. **Draft text left in published examples:** "Wait — check…" (Ch1 L829–838; Ch6 L563).
8. **Heatmap colours.**
   - In xcolor, `blue!0` is **white**, not dark as BOOK.md says. Low-value cells therefore render white, the opposite of the intended dark-to-bright mapping.
   - White text on `blue!3–20` cells is unreadable.
   - Several heatmaps print values that do not match the matrices in the text (Ch1, Ch5, Ch6, Ch16, Ch17, Ch23, Ch24).
   - Some are formula-generated but labelled as data ("Bolivia W^P (2016)", "South Korea W^P (2016)"), and their spectra contradict the stated Γ: the Bolivia pattern gives 0.36 (stated 0.13) and the Korea pattern 0.98 (stated 0.45).
   - Ch16 uses `rnd`, so the figure changes on every compile.
   - **Fix:** generate every heatmap from the actual matrix with a script (Python → TikZ/pgfplots `matrix plot`). That guarantees the picture equals the math.
9. **Front matter.**
   - The title page says "DRAFT — September 2027" with © 2027.
   - The dedication is buried mid-chapter in Ch28 (L305–308); move it to the front matter.
   - The "Answers to the Prologue Challenge" section sits at the end of Ch1 (L1048); move it to the Prologue.
   - Add a **notation table**. Symbols currently collide:
     - W is welfare *and* the weight matrix.
     - Δ is depreciation *and* the deviation.
     - L is labour, the Leontief inverse *and* the loss.
     - Ψ is the evolution map *and* the formation ratio.
     - Ω is world output *and* the objective; ω is endowments *and* relatedness density.
     - ρ is spectral radius *and* persistence; γ, φ, δ and B are also reused.

---

## 6. `BOOK.md` (the construction manual) itself

- **Update the inventory.** It says Chapters 8–28 are "NOT WRITTEN", but they all exist.
- **Fix the template bugs:** the `\vizpanel{Title}` argument (§5.2) and the colour mapping ("0 → blue!0 (black/dark)" is wrong).
- **"Rule 1: NEVER CUT EXISTING MATH" should go.** The review shows a lot of existing math is wrong. Replace it with: "never cut math *without recording why*; every change to a result must be checked numerically."
- **Add rules the book needs:**
  1. Every number is either sourced or tagged *Illustrative*, and recomputed by script.
  2. Theorem means proved.
  3. State the normalization of W and the interpretation of π in every chapter.
  4. Γ is never used as a measure of concentration.
  5. Exercises with numbers ship with a solutions script.
- **Template compliance is uneven.** It is not urgent, but it is easy to fix:
  - Summaries have 6–7 items in most chapters (the rule is 5).
  - Many chapters lack the real-data exercise.
  - Opening stories violate Rule 3 (Ch3 §2.1, Ch5 "Before the Math: Three Stories", Ch1 §4.1, Ch24, Ch18's Potosí story).
  - Ch21 and Ch22 abandon the Taylor assignment and the single-thinker rule.
  - Ch21, Ch27 and Ch28 have no viz panel; Ch15, Ch19 and Ch22 panels are not heatmaps.
  - Ch28 has no summary, exercises or appendices.
  - Further Reading never uses the four required categories.
  - "Prebisch Connection" boxes in Ch24 and Ch27 have no Prebisch content.
  - Revisit the "one thinker per chapter" rule where it forces fit. Several case studies bend a thinker to the apparatus (Skocpol, Fajnzylber, Ocampo, Marini). It is better to say explicitly "this is my formal reinterpretation" and mark where the thinker would disagree.

---

## 7. Framing and tone

- **Retitle Ch28's "General Case Theorem" (L83–94).** Every W equals W^CE + Δ(W), so the content lies in the maps from W to outcomes. Those maps need extra primitives (κ, ρ, μ, the coalition game, the product space) that are not functions of Δ(W). Call it a **Framing Principle**. The back-cover claim that the book "proves that every result in development economics … is a corollary of Δ(W)" should go.
- **Recast the conclusion as a research programme.** The best material in Ch28 is "What the book does not do" (L233–243). Build the ending on that honesty.
- **Drop the straw man.** Replace "the standard model assumes W^I doesn't exist / A is technological fact / preferences are unexplained" with "the benchmark Arrow–Debreu model abstracts from…". Asymmetric information, Bayesian network learning, endogenous input-output structure and Becker–Stigler consumption capital are mainstream, and several of them are cited in your own Further Reading.
- **Cut the hyperbole.** For example: "This is not a metaphor. It is a mathematical statement"; "This book is the proof"; "the single most important result"; "remarkably similar" (for 0.13 vs 0.20); "Che boludo".
- **Rework the sensitive material for care and accuracy.**
  - Ch18 contrasts gender ("biological component") with race, contradicting its own social treatment of gender.
  - Ch22 calls the Santa Cruz elite "W^U ≈ 0", which is an unsupported label.
  - Ch22 frames ayllu solidarity as "fragmentation that prevents broader coalitions", which ignores CONAMAQ and the Katarista federations. Consider John Murra's verticality as the thinker for that case study.
  - Ch3's OLG "men vs women patience" framing risks reading as essentialist; use neutral group labels.

---

## 8. Chapter-by-chapter findings

These are the most important items per chapter, beyond the cross-cutting issues above. "→" gives the corrected value from numerical recomputation.

### Ch1: Economics in Matrix Notation
*Works:* the allocation-matrix and Leontief exposition is clear and well paced.
- **L370–400:** asymmetric bargaining example is degenerate (satiation); recalibrate as in §1.4.
- **L512–533:** Leontief inverse is wrong.
  - Correct L = (I−A′)⁻¹ = [[1.138, 0.246, 0.118], [0.088, 1.084, 0.246], [0.131, 0.088, 1.138]].
  - Multipliers → (1.357, 1.417, 1.501).
  - Output for d = (0, 0, 100) → (11.78, 24.58, 113.79).
- **L542–567:** heatmap values match neither A nor L, and the label says (I−A)⁻¹ while the text uses (I−A′)⁻¹.
- **L978–994:** Neumann table is wrong; the exact inverse is (1.304, 0.435, 0.145, 1.159).
- **L485:** convergence reasoning is wrong. aᵢⱼ < 1 doesn't imply convergence; you need row sums < 1 or ρ(A) < 1.
- **L960:** "entries of Aᵏ ≤ ρᵏ" is false; use the Gelfand bound.
- **L1006–1028:** the Prebisch case study maps to objects not defined in Ch1 (π, Var(π), terms of trade). Map the three steps to Ch1 objects and flag the rest as forward-looking.
- **Ex. 1.7 / 3.5:** Bolivia is not in OECD ICIO. Use INE supply-use tables, Eora, or CEPAL regional I-O tables.

### Ch2: The Welfare Theorems in Matrix Form
*Works:* the FWT proof and the block-loss algebra are cleanly presented.
- FWT-matrix, doubly-stochastic and spectral-gap propositions are false (see §1.3, §1.1).
- **L204–234:** corrected values.
  - Var(π) → 1.70×10⁻³ / 30.0×10⁻³.
  - L → 0.0068 / 0.120; ratio 17.6×, not 16×.
- **L390–416:** "Transfer vs Structural" example is inconsistent.
  - The allocations imply κ = 0.125, not 0.5.
  - Uniform cross-block entries give block masses exactly ½/½, so the stated π cannot arise.
  - Scaling cross-block entries leaves block masses unchanged (Simon–Ando), contradicting Policy B.
- **L372–388:** Transfer Trap part 4 is false within the model. Nash with disagreement point (Ω/N)1 returns the first-best for any π.
- **L656–684:** Furtado case study. BOOK.md's core mapping (underdevelopment as a stable fixed point) appears only in an exercise with an undefined f. Give a 2×2 example in the main text.
- **Ex. 2.6:** mixes a 2×2 matrix with 8 agents, and κ = 0.001 makes U ≈ −5625.

### Ch3: Growth and Distribution in the Weight Matrix
*Works:* the Bolivia enclave Leontief example and its heatmap are the clearest pedagogy in the opening chapters.
- **L291–307:** corrected inverses.
  - L_A → [[1.184, 0.329], [0.263, 1.184]].
  - L_B → [[1.112, 0.025], [0.025, 1.112]].
- **L376–381:** the symmetric W must have a symmetric inverse.
  - L_conn → [[1.217, 0.228], [0.228, 1.293]].
  - Multiplier 1.445; spill 15.8%; 7.8× the enclave.
- **L642–656:** eigenvalues → 0.21 and 0.14 (gap 0.070). Connected case → 0.327 and 0.023 (gap 0.304); ratio 4.3×, not 5×.
- **L309–329:** the growth-propagation bound and its proof step are false. Use `γ‖W‖/(1−γ‖W‖)` with any submultiplicative norm.
- **L426–449:** "Pinto heterogeneity is a spectral property" fails (§1.1). Measure decomposability by cross-block mass instead.
- **L181, L262:** gas share, export share, quinoa and Bolivia's inequality path (§3).
- **L484–501:** "Transfer Perversity" has no model (no wage premium or labour mobility). Add a Lewis or Harris–Todaro block, or downgrade it.
- **L97–99:** delete the "Che boludo" pseudo-epigraph.

### Ch4: What the Standard Model Assumes
*Works:* the A1–A6 catalogue is a useful organizing device.
- **L237:** compile error (`\boldsymbol` in `\section`).
- **L338–345:** "A5 follows from A1–A4" is false (§1.3).
- **L394–412:** "Empirical Departures" is unsourced (§3). Tag it illustrative or source it. *(Correction: an earlier version of this review called Var(π^I) ≈ 0.03 for N = 20 impossible. It is not: the maximum of (1/N)Σ(π_i − 1/N)² is (N−1)/N² = 0.0475.)*
- **L161–167:** "(Jacobian terms)" is not mathematics.
  - Correct derivative: dL = L·γdW·L.
  - Comparative static: ΔK* = (Δ − SLF′)⁻¹SLγΔA·Lf.
- **L575–579:** π is analytic in W for irreducible W, so there are no "phase transitions" except at the reducible boundary.
- **Case study:** conflates Sunkel 1971 with Sunkel & Paz 1970, and the mapping rows are mismatched.

### Ch5: Three Networks, One Economy
*Works:* the three-network framing and the W^P vs W^U distinction.
- **L627–705:** La Paz example.
  - The composite entry W₂₁ → 0.31.
  - π → (0.311, 0.365, 0.325). The Executive leads and the Vendor is *not* marginal.
  - π(W^I) → (0.28, 0.27, **0.45**): the Vendor has the *most* epistemic power, so the caption is reversed.
- **L609–615:** eigen-identity is missing 1/ρ (ρ = 0.566). π is not invariant to layer rescaling.
- **L249–259:** Leontief equation transposed relative to the definition.
- **L526, L715–723:** welfare-theorem straw man.
- **L119–153:** opening stories violate Rule 3.
- **L898–928:** case study's Γ column fails (colonial W ≈ `1e₁'` has Γ = 1). It also has no "This chapter" column and imports Sunkel and Pinto.

### Ch6: Power and the Perron Eigenvector
*Works:* degree → Katz → Bonacich → eigenvector is a good arc; the rich-farmer 2×2 example is correct (π = (9/11, 2/11)).
- **L241–263:** five-agent example; correct values (ρ = 0.810):

  | Measure | Agent 1 | Agent 2 | Agent 3 | Agent 4 | Agent 5 |
  |---|---|---|---|---|---|
  | Eigenvector π | 0.139 | 0.205 | 0.230 | 0.214 | 0.212 |
  | Katz (δ = 0.5) | 0.498 | 0.735 | 0.761 | 0.682 | 0.676 |

  Agent 2 ranks 4th by eigenvector, not 2nd.
- **L184–204:** Bonacich as defined is exactly Katz/β (verified), so "Bonacich makes neighbours' centrality matter, Katz doesn't" is wrong.
- **L328–332:** PF item 5 needs primitivity.
- **L355–367 / App. 6B:** Hilbert-metric proof has gaps (variables mismatched, complex-λ step, simplicity not shown).
- **L431–515:** topology table and viz (§1.1). The constructions must be specified; the viz also says N = 8 while the table says N = 10.
- **L563:** remove "Wait — check…".
- **L778:** Acemoglu et al.'s influence vector is not the Perron vector.

### Ch7: The Bargaining Equilibrium
*Works:* the TU landlord–farmer split `xᵢ = dᵢ + πᵢS` is correct and intuitive. Build the chapter on it (§1.4).
- **L222–250 / App. 7A:** LQ allocation equation is underivable (§1.4).
- **L274–295:** "nesting" and "direction" propositions are false (§2).
- **L400–455:** worked economies.
  - Economy A: π → (0.244, 0.277, 0.243, 0.237).
  - Economy B: π → (0.346, 0.386, 0.134, 0.134).
  - Welfare loss B → 0.362, not 0.18. Ratio 45×, not 18×.
  - The B inequality comes from coefficient scale, not block structure: row-normalizing B gives nearly uniform π.
- **L121–136:** landlord–farmer inconsistencies.
  - K^0.5 L^0.5 with K = 10, L = 2 gives 4.47, not 20.
  - Cobb-Douglas gives zero output from land alone, contradicting the Landlord's autarky payoff.
  - (14, 6) *is* the equal split.
- **L533–547:** escape cost ∝ π_h/Γ has the wrong direction (§1.6).
- **L350:** "unequal exchange" is Emmanuel. Only 6 exercises (minimum 7).

### Ch8: Belief Manipulation and Epistemic Power
*Works:* the cleanest chapter. The DeGroot exposition, N_eff, and the Quijano case study follow the template.
- **L154–165:** star example; all non-unit eigenvalues are 0.1.
  - Γ → **0.9**; half-life → 0.30 rounds.
  - π → (2/3, 1/9, 1/9, 1/9); b* → 0.722.
  - This is itself a counterexample to "star ⇒ small Γ".
- **L297–313 / App. 8B:** manipulation-cost sign and the zero perturbation formula (§1.6).
  - Reachability constraint: b_T must lie in [min b⁽⁰⁾, max b⁽⁰⁾].
- **L110, L232–246:** Golub–Jackson misattributed; the "wisdom" bound needs σ̄ = min σⱼ.
- **L272–276:** "π_IMF + π_WB + π_UST ≈ 0.6" is invented; three agents at 0.2 each give N_eff ≈ 8.3, not 3.
- **Ex. 8.4(d):** ill-posed; b* = 0.5 for all β, and only the half-life diverges.

### Part III overview (Ch9–12): The Formation Model

**Cross-cutting problems.**
- **Ψ has four incompatible definitions:**
  - a derivative ratio F/E (Ch9);
  - a sum of absolute effects (Ch10);
  - a *variance* ratio (Ch11; in squared units, so Ψ̂ = 2.08 corresponds to an effect ratio of about 1.44);
  - an R² ratio (App. 11B).
- **Ψ values drift:** 0.5–2.1, "0.5–4.0", "2–5", 1.2, 3.5, 1.6.
- **ρ and the time unit are inconsistent:** ≈0.9 (Ch9), 0.85 (Ch10), 0.62 "per generation" but estimated on decade cohorts (Ch11), ≈0.2 (the Dohmen citation).
- **Fix:** define one estimand Ψᵢ(ΔW, T) and derive the rest from it.

**Preferences vs beliefs is never addressed.** Most of the cited evidence is about *beliefs*:
- Fernández (2013) and Fogli–Veldkamp model learning.
- Bursztyn et al. (2020) is about misperceived norms: men's beliefs about others' support were corrected; their own support did not "increase dramatically".
- Malmendier–Nagel is personal experience, not W^I.

This duplicates Ch8's DeGroot learning. Add a section separating taste from belief (Bursztyn, Ederer, Ferman & Yuchtman 2014; Dietrich & List 2013), or reinterpret θ explicitly.

**The model cannot sign the welfare effect of W^I edges.** Edges liberate in the Ch9 intuition, harm in Ch9's Furtado box and in Ch8, and are "always good" in Ch12. Edges need content, i.e. the relevance or attainability of the observed returns.

**Literature to engage:**
- Becker (1996); Becker & Murphy (2000).
- Mookherjee, Napel & Ray (2010); Genicot & Ray (2017); Dalton, Ghosal & Mani (2016).
- Durlauf (1996); Bénabou (1996); Brock & Durlauf (2001); Panebianco & Verdier (2017).
- On welfare under endogenous tastes: von Weizsäcker (1971); Hammond (1976); Pollak (1978); Bernheim–Rangel (2009, properly, not only in Further Reading); Sen's intersection quasi-ordering; Nussbaum (2001); Khader (2011).

### Ch9: Where Preferences Come From
*Works:* the disciplining criteria for Becker–Stigler (L122–127).
- **L270–282, L570: "Two Steady States" contradicts its own theory.**
  - It uses a *reducible* block-diagonal W^I, while the proposition assumes irreducibility.
  - Each map is a contraction, so each agent has a *unique* steady state. The two θ* come from heterogeneous returns, not multiplicity or the network: with W^I = I the gap is unchanged.
  - The numbers themselves (1.0882, 0.9227) are right.
- **L333–338:** the connected 0.6/0.4 gap → 0.031, not 0.04. For Ex. 9.2(b), θ → (1.033, 0.971).
- **L294–305:** phase diagram.
  - True slopes are 0.83 and 0.806; the drawn lines have 0.35 and 0.15.
  - Neither marked fixed point lies on its line.
- **L349–362: "Hysteresis" theorem.**
  - The decay rate is ρ(ρI + μW^I·Dg), not ρ: 0.83 vs 0.8 in your own example, giving a half-life of 3.7 periods, not 3.1.
  - The result is local only.
  - With a unique stable steady state this is persistence, not hysteresis.
- **L140–199: citations misdescribed.**
  - Falk et al. (2018) is a survey; "factor of 10" on a z-score is meaningless.
  - Kosse et al. (2020) is about prosociality.
  - The immigrant LFP result is Fernández & Fogli (2009) and Fernández (2007). Fernández (2013) is in the *AER*.
  - Bronnenberg et al. find about 60% of the gap closes quickly, which is evidence *against* ρ ≈ 0.9.
  - Atkin (2013) finds that tastes *persist*.
- **L255–267:** Brouwer needs bounded g and a compact invariant Θ, but every example uses an unbounded linear g.
- **Ex. 9.3:** ill-posed. There is a unique steady state (0.542), so the "war shock" has no permanent effect, contradicting part (d). The WWII result is Fernández, Fogli & Olivetti (2004).
- **Furtado (*O Mito*):** income concentration lets a minority imitate centre consumption, and that demand shapes technology. "Low mass aspirations" is the book's gloss. "Demonstration effect" is Duesenberry (1949) and Nurkse (1953), not Furtado. Verify the epigraph.

### Ch10: The Formal Formation Model
*Works:* the formation Leontief inverse (I−J)⁻¹ is a clean, teachable bridge to Ch1/Ch3.
- **L349–357: the poverty trap is impossible.**
  - Φ(θ) − θ is a quadratic, so there are at most two fixed points: **θ = 1.258 (stable) and 8.742 (unstable)**.
  - The stated 0.78, 1.02 and 1.35 are not fixed points (Φ − θ = +0.048, +0.023, −0.008). Every start below 8.74 converges to 1.258.
  - Use an S-shaped (logistic) g to get three fixed points, and redo Ex. 10.2 and 10.6.
- **L214–235: the Bolivia three-block table is inconsistent.**
  - Returns don't match eq. (returns): it gives 0.833, 0.367, 0.123, not 0.182, 0.112, 0.063.
  - With the *stated* W^I, θ* → (1.090, 1.018, 0.968). The table's (1.109, 1.016, 0.951) is what pure block-diagonal W^I gives.
  - W^P, γ, κ and Φ are never given, so Ex. 10.4 ("verify") is impossible.
  - The block π values sum to 1, so they are totals, not means.
- **L464–502: the policy example.**
  - The cross weight is set "0 → 0.15", but the baseline was 0.025.
  - The direct effect → 0.024, not 0.06.
  - The information lever has E = 0 by Prop "levers", so Ψ = F/E is undefined, not 3.5.
- **L456, L477: transfers.**
  - r is linear in x, so transfers have a first-order effect, not second-order.
  - In eq. (bargaining), x doesn't depend on ω at all.
  - Put ω into h, and engage Macours & Vakis (2014).
- **L139–142: Jacobian.** ∂h/∂θ = (Φ/Nκ)(I − 11′/N), so J·1 = ρ·1. The *aggregate* preference level is never amplified, which rules out economy-wide "virtuous circles". That is a substantive result; state it.
- **L302–313: stability bound.** ρ(·) is not submultiplicative, so use norms. The condition is sufficient, not "required".
- **L412–424:** the Ψ LQ formula turns *negative* as ρ → 1, rather than diverging; it blows up at ρ = 1 − μβη·ρ(W^IW^P).
- **L340–347:** "Basin Structure" is false; your own convex example has one steady state.
- **L643–650:** the bifurcation threshold is the contraction bound, not a saddle-node condition.
- **L666–692: Cardoso & Faletto.**
  - They were explicitly anti-formalist (historical-structural method; Palma 1978; Cardoso 1977 "Consumption of Dependency Theory"). Engage the objection rather than assert "no longer narratives; coordinates".
  - Their ideologies are class-alliance political projects, not a scalar θ.
  - "Associated-dependent development" is Cardoso (1973).
  - "Rejected Dos Santos" is anachronistic.
  - Cut "Cardoso's presidency as an exercise for the reader".
  - The mapping table overflows by 149pt.

### Ch11: The Ψ Ratio, Measuring Formation
*Works:* the epidemiological design is the right anchor; second-generation immigrants are handled correctly.
- **L121–136: age–period–cohort non-identification makes Ψ̂ arbitrary.**
  - On Ex. 11.1's own table, Ψ̂ = 8.05.
  - Shifting a linear trend between cohort and year effects, which leaves the fit unchanged, gives **Ψ̂ ∈ {0.03, 0.62, 8.0, 507}**.
  - "Cohort = formation, year = enforcement" fails both ways. The pill operated by age at access; formation can be contemporaneous.
  - Cite Hall (1971), Heckman–Robb (1985), Deaton–Paxson (1994), McKenzie (2006), and Fosse–Winship (2019); state the restriction and show sensitivity.
- **W^I disappears from the econometrics.**
  - The steady state is a spatial-autoregressive model with an endogenous effect, and cohort-mean regressions are Manski's non-identified case.
  - Bramoullé–Djebbari–Fortin requires I, W, W² to be independent. The block-diagonal W^I used in Ch10 fails this (W² = W).
  - Add a section on measuring and identifying W^I: Lee (2007); Calvó-Armengol et al. (2009); Blume et al. (2011); Goldsmith-Pinkham–Imbens (2013); de Paula–Rasul–Souza; partial-population designs.
- **L158–178: the epidemiological derivation.**
  - It skips the daughter's θ.
  - Origin LFP includes origin enforcement.
  - App. 11B's β ≈ ρ^g contradicts ρ ≈ 0.625 from the 0.25/0.40 ratio.
  - Controlling for education, marriage and children is a bad-control problem.
- **L173–257: unsourced estimates presented as findings.**
  - Fernández (2013) is a calibrated model in the *AER*. The immigrant evidence is Fernández & Fogli (2009): 1970 Census, about 25 countries.
  - The variance numbers, DHS-28, Bolivia EH estimates with SEs, savings β̂ and PSID split have no sources.
  - Source them or tag them illustrative, and drop the SEs.
- **L391–407: bounds table.** Combining the stated reductions gives 0.83 multiplicatively or 0.41 additively, not 1.05. So "Ψ > 1 in the worst case" fails on your own numbers, and 0.8 × 2.08 = 1.66, not 1.60.
- **L294–356: structural estimation.**
  - Var(θ) is missing from the correlation formula.
  - One moment is used for three parameters.
  - θ̂_c is a generated regressor that is identified only up to the age–period–cohort (APC) trend.
  - A lagged dependent variable in a short panel suffers Nickell bias.
  - ρ is per decade, not per generation; 0.62 per decade is about 0.30 per 25 years.
  - SE(M̂) → 0.165.
- **L281, L455:** "Ψ = 0 rejected everywhere" is vacuous, since a variance ratio is ≥ 0 and sampling noise is positive. Use F-tests and bias-correct.
- **L389, L499–504: twin and adoption studies misread.**
  - Cesarini et al. (2009) find common environment ≈ 0.
  - Sacerdote (2007) never observed biological parents, so "β̂_bio = 0.04" is impossible by design. Rewrite Ex. 11.3.
- **Labels:** 9 of the 11 "Propositions" are empirical claims; relabel them. The Bragg epigraph needs "attributed to".

### Ch12: Welfare Under Endogenous Preferences
*Works:* posing the circularity explicitly and ordering the approaches from conservative to radical is good pedagogy.
- **L129–198: the worked example is wrong, and the correction reverses its conclusion.**
  - With U = θx − (κ/2)x², κ = 0.5: **U → 0.228 (gas), 0.605 (urban), 0.832 (rural)**. The table has 2.24 for gas, a factor-of-10 error.
  - At θ = 1: −0.232, 0.551, 0.952. The *rich* block has negative welfare because all allocations exceed the bliss point θ/κ.
  - Raising rural θ *raises* rural utility, so the LQ form encodes the *opposite* of adaptive preferences.
  - Use reference-dependent utility u(x − aᵢ) with aspiration aᵢ = θᵢ (Genicot–Ray 2017; Clark, Frijters & Shields 2008), and recalibrate below satiation.
- **L105–116: "Satisfaction Paradox" theorem.** The proof derives U = θ²/(2κ), which is increasing in θ. That contradicts the theorem.
- **L321–339, L475–487: the "Structural Welfare Theorem" is false.**
  - (b) and (c) rest on Γ ⇒ equal π. Counterexample (N = 3): W_A = [[.43,.51,.06],[.01,.96,.03],[.39,.02,.59]] and W_B = [[.02,.45,.53],[0,.12,.88],[.11,.88,.01]].
    - Γ_A = 0.38 > Γ_B = 0.18, and Δ_A = 0.61 < Δ_B = 1.07.
    - Yet Var(π_A) = 0.140 > Var(π_B) = 0.040.
  - The algebra in the proof of (a) is wrong. The identity is U = θ²/(2κ) − d²/(2κ); the book's "+θd/κ" is wrong (the cross-term cancels).
  - Connecting the blocks *lowers* mean θ in Ch9's own example.
  - (d) holds "by assumption", which is a tautology.
  - **Verdict:** the circularity is renamed, not solved. Choosing Γ, N_eff or W^I_just *is* the value judgment, so say so.
- **L236–265: the capability formalization is not Sen's.**
  - O_i is defined as the image of preferences, but Sen's capability set is the feasible functioning set given entitlements and conversion factors.
  - Σᵢ Cᵢ is constant in π, so Prop (3) is false.
  - "Jensen, h concave in π" is wrong: h is linear in π.
  - Engage Pattanaik & Xu (1990), Sen (1991) and Sugden (2004, *AER*).
- **L292:** "adding W^I edges is robustly welfare-improving" contradicts Ch8 and Ch9.
- **L389–393:** the Pareto sets are not nested (two-agent counterexample).
- **L395–401:** Prebisch (1981) predicted a *crisis*; he was not "puzzled by stability".
- **L483:** undefined reference to the false Ch2 bound.

### Ch13: How Networks Evolve
*Works:* the three-mechanism taxonomy.
- **L218–242:** both "equilibria" are symmetric and non-stochastic.
  - Γ_A → 0.101; π_A → nearly uniform, with the rural block *highest*.
  - Γ_B → 0.649.
  - θ* is underivable because Φ, Ψ and h are unspecified, so Ex. 13.2(a) is impossible.
- **L194–200:** Brouwer existence fails (§2). Restrict to {W ≥ ε} and renormalize inside Ψ.
- **L335–344:** integration "iff" is false. **L305–317:** path dependence is a heuristic.
- **L357–359 / Ex. 13.5:** Sunkel's prediction holds in only 57% of random starts. Specify the starting W and how rows are rebalanced.
- **L379–410:** phase diagram is linear (one fixed point), contradicting the multiplicity claim; Γ also runs to 8.

### Ch14: Coalitions and Class Alliances
*Works:* the Acemoglu commitment-problem motivation.
- **L360–370:** constant-sum game, generically empty core, indifferent grand coalition (§2). Payoffs need a surplus term Ω(W).
- **L168–179:** existence is not guaranteed by Bernheim–Peleg–Whinston. Use Bogomolnaia–Jackson (2002) conditions, or downgrade to an example.
- **L377–385:** the characterization does not cover its own baseline; the regions in the figure don't match the inequalities.
- **App. 14B:** perturbation formula ≡ 0 (§1.6).
- **L570:** Cardoso launched the Plano Real as Finance Minister (1993–94).
- **L528:** compile error.

### Ch15: Operator Design
*Works:* the ISI-A vs ISI-B teaching device.
- **L376–407:** the political-cost results are empty (π-preserving path to Γ = 1; §1.1). Rebuild the objective around Var(π) or the power Gini.
- **L128–139:** edge-selection theorem is wrong under a Frobenius budget, and eq. L136 contains a literal "correction" placeholder.
- **App. 15A:** derivative fixes (§1.6).
- **L115–119, L212:** multipliers degenerate under row-stochastic W (§1.5).
- **L174–193:** ISI example can't be replicated.
  - With a symmetric calibration, ISI-A gives ΔΓ = 0.000 and ISI-B +0.047 (text +0.12).
  - ΔVar(π) is positive, not negative.
  - Print the W.
- **L283–306:** greedy sequencing isn't optimal without submodularity; the information operator has ΔΓ(W^P) = 0.
- **Tavares:** fix the date (1964), and present the mechanism as external strangulation (*estrangulamento externo*).
- **Engage Liu (2019) and Boyd–Diaconis–Xiao (2004).**

### Ch16: The Global Weight Matrix
*Works:* the ECI, RCA and proximity definitions.
- **L113–172:** own schematic gives Γ ≈ 0.90 (§1.1).
- **L98–104:** the W^P definition is inconsistent about the diagonal (domestic absorption). Specify source, year, aggregation and diagonal treatment, and report sensitivity.
- **L177–210:** wrong orientation, and the linkage table is impossible (§3).
- **L243–306:** ECI is not a function of W^P and is z-scored. The casillero bound is false. 6 of 14 countries plot *below* the "lower bound" and none sits in the box the caption describes.
  - Useful correct facts to add: M̃ has real eigenvalues in [0,1], and its left Perron vector ∝ diversity.
- **L393–405:** Prebisch–Singer "theorem" has no ToT link, and USD-valued W^P *does* move with commodity prices.
- **Viz:** titled 20×20 but N = 10; uses `rnd`; white labels are invisible.

### Ch17: Domestic Structural Diagnostics
*Works:* the dashboard idea; σ_B arithmetic is correct.
- **L108–114:** normalization mixes denominators. Report Γ under two normalizations and two aggregations.
- **L116–207:** Bolivia I-O example and heatmaps look invented and are impossible (§3, §5.8).
- **L266–291:** π sums to 0.65; Gini → 0.24. π divided by employment share mostly measures capital intensity, not "structural privilege".
- **L307–343:** N_eff depends on how many outlets are enumerated. The Ψ ratio faces the age–period–cohort identification problem.
- **L412–418, L507:** "Γ alone captures most information" is false.
- **L486–492:** the inequality-fall narrative conflicts with decomposition evidence and with Ch16.

### Ch18: Gender, Race, and Structure
*Works:* the node-attribute / mixing-function formalism.
- **L120–131:** the assortative proposition holds only for an equitable two-block partition; prove it via the quotient matrix [[a,b],[b,a]].
- **L127–176:** symmetric homophily gives π ∝ degree, so there is *no* gender gap at 50/50. A minority gap arises from group size: `d_min/d_maj = (p + λ(1−p))/((1−p) + λp)`, which is 0.84 for Bolivia (p = 0.41, λ = 0.35).
- **Ex. 18.3:** asks students to verify something false (π is uniform). Weight W^max by population.
- **L299–312:** arithmetic.
  - Power share → 0.284, not 0.33.
  - Power/population → 0.69, not 0.80.
  - Intersection → 0.33, not 0.46.
- **L409–417:** welfare example is impossible (Var bound).
- **L316–341:** intersectionality is reversed (§4.3). Marini is misrepresented (§4.3).
- **L284:** the 2012 Census had no "mestizo/white" category.
- **L176:** contradicts Ch17's sector π ranking.
- **Ex. 18.4:** degenerate (Δθ* = 0 for every λ).
- **L558:** half-life arithmetic. 500 years is 2.5–4.7 half-lives, not "2–3".

### Ch19: Scale Invariance
*Works:* the nesting question is the right one.
- **L148–167:** "scale invariance" is undefined, and its conclusions are near-vacuous.
  - A real result is available: under exact lumpability, Γ_coarse ≥ Γ_fine.
  - That shows the scale ordering is partly an aggregation artefact.
- **L124–139:** the individual-scale column has no possible data source.
- **L200–239:** the Kronecker algebra is correct, but the justification is Simon–Ando near-decomposability, not Kronecker structure. Per-capita factors omit size: the sector factor is about 32×, not 1.8×.
- **L287–310:** "Theorem five" is a dictionary; label it a Definition.
- **L463–472:** App. 19B eigenvalues are wrong.
- **L505:** "Chilean exile" contradicts the 1970 date.

### Ch20: SAMs as Weight Matrices
*Works:* choosing the SAM as the data object is right.
- **L167–191:** the SAM does not balance. Only Gas, Lab and Cap balance; there is no household consumption and no savings account. A balancing fill is in the reviewer notes.
- **L257–279:** the reported eigenvalues can't come from this table. The leading pair is complex (circular flow) and Γ = 0.40–0.50.
- **L284:** Γ_SAM ≤ Γ_IO is refuted by the chapter's own data.
- **L126, L286:** π = account size for any balanced SAM.
- **L306–317, L489–506:** Pyatt–Round decomposition is misstated.
  - Correct form: M = M₃M₂M₁, with M₂ = I + A* + A*².
  - The block pattern is transposed.
  - Column sums of a product ≠ product of column sums: 5.70 vs a true 5.18. Use Stone's additive form.
  - The "M₃" row equals Ch17's Rasmussen indices.
- **L396–414:** Labour → poor-household share → 0.40. ΔΓ from the reallocation ≈ +0.0001, not +15%.
- **Taylor representation** (§4.3).

### Ch21: Measuring the Information Network
*Works:* the four-layer scheme; the limitations appendix.
- **L99–343:** recomputed media table contradicts the chapter (§3).
- **L167–177:** misdescribes Banerjee et al. (2013).
- **L289–298:** validation results appear invented.
- **Missing:** the network-measurement and identification literature (§4.4). App. 21B misstates Manski.
- **L69–72:** unattributed epigraph. **L231:** compile error.

### Ch22: Measuring the Utility Network
*Works:* the five-strategy menu.
- **L390–395:** the altruism FOC is inverted; it should be `Wᵢⱼ = (cⱼ/cᵢ)^σ`. Altruism is confounded with exchange and insurance motives (Cox 1987; Foster & Rosenzweig 2001).
- **L304–328:** alignment theorem is false (§2).
- **L113–131:** the transfer matrix is not EH data, inverts the remittance flow, and is 4×6.
- **L144–217:** experiments misattributed; referendum and Latinobarómetro tables invented. The dictator share s is not W: W = s/(1−s).
- **L67–70:** the epigraph is papal (§4.2).
- **L259–296:** N_eff compared across different node sets.
- **L414–442:** ayllu case study (§7).

### Ch23: Tracking Structural Change
*Works:* the RAS, HP and Bai–Perron toolkit; SCI arithmetic is correct.
- **No theory of products of time-varying W** (§1.7).
- **L297–316:** the decomposition "theorem" is a definition, and its taxonomy is inconsistent.
- **L128–134:** π is badly conditioned at Γ ≈ 0.13, so bootstrap the break tests.
- **L431:** HP cutoffs. λ = 400–1600 gives 28–40-year cycles; 10–20 years needs λ ≈ 6.25–100.
- **L162 (gas take), L324 (Norway's fund invests only abroad).** The case study uses CEPAL rather than the assigned Taylor.

### Ch24: Hegemony and Counter-Hegemony
*Works:* the Gramsci/Skocpol pairing and the "three loops" device.
- **L127–137:** cost sign error, inherited from Ch8. Its own viz falsifies "concentration = small Γ".
- **L168–176, L283–294:** threshold is ill-posed; the "cycle" doesn't cycle. **L379:** saddle-node, not pitchfork.
- **Skocpol and Gramsci** (§4.2–4.3).
- **L179–181:** the miners' radios were union-founded; the Water War was led by the Coordinadora.

### Ch25: Truncated Industrialization
*Works:* τ_c is simple and testable; Ex. 25.7 (lithium) is well designed.
- **L89–95:** "dense subgraph" is ill-defined. Pick Goldberg/Charikar densest subgraph with |D| ≥ k, a weighted k-core, or modularity with stated thresholds.
- **L118–130:** a dense module is not a bridge, so the theorem fails (§2).
- **L97–105:** Hidalgo et al. used about 775 SITC-4 products, not "774 HS-4"; HS 90 is listed twice.
- **L226–238:** the "three-quarters" truncation penalty depends on an arbitrary Γ(0). Report the 55–86% range.
- **Fajnzylber fidelity; Brazil's II PND** (1974–79).

### Ch26: The Neo-Structuralist Synthesis
*Works:* the plan arithmetic is internally consistent; the pathology table teaches well.
- **The theorem is not well-defined.** ECI is not a function of W and is z-scored (positional). Γ* and N* are "estimated" with no data. Ex. 26.1 estimates Γ from ECI, which makes the conditions collinear.
- **Proof (1)** rests on the false bound (§1.1).
- **Falsified by Korea.** Korea grew with equity under authoritarian media and with Γ below Γ* before 1980.
- **Viz:** Chile is plotted inside the green zone but the caption says it is outside; Venezuela's ECI is inconsistent across the table, plot and exercise; the casillero is mislabelled.
- **Ocampo** (§4.2–4.3).

### Ch27: Algorithmic W^I
*Works:* the operator framing and the state-dependent DeGroot model.
- **L124–132 / App.:** echo-chamber proof is circular, and "for all k" is false. Use biased assimilation (Dandekar–Goel–Lee 2013).
- **App. 27B L316:** the diverse equilibrium is impossible; recast with bounded confidence.
- **Γ vs N_eff move in opposite directions** (§1.1).
- **Empirical literature absent; WhatsApp is not a recommender system** (§4.4).
- **N_eff is whole-network.** Per-user claims need row entropy. N_eff is not additive.
- **Ch24 says social media raised N_eff; Ch27 says it lowers it.** Reconcile.

### Ch28: The General Case
*Works:* "What the book does not do" (L233–243) is honest and well pitched.
- **L104–141:** the summary table overflows by 259pt and its bottom rows don't print.
- **L83–96:** the theorem is a tautology, and several "corollaries" are false (§2).
- **L180–184:** "neoclassical = Γ = 1" (§1.1).
- **Furtado biography errors; a structuralist CGE exists (Taylor 1990)** (§4.3).
- **No summary, exercises, appendices or viz panel.** The dedication sits mid-chapter.

### Covers and the matrix explorer (`covers/*.jsx`)
- **`matrix_explorer.jsx` never computes Γ,** the book's headline diagnostic.
- **`leftPerronApprox` is unreliable.** It uses 150 steps of plain power iteration with `Math.abs`, which does not converge on periodic matrices (the star, even rings). Iterate on the lazy chain `½(I + W)` (same π, guaranteed convergence) or use a proper eigen-solver.
- **`spectralRadius` uses a random start,** so results are non-deterministic.
- **Preset descriptions repeat the false claims:** "Complete … Maximum spectral gap. Equal power"; "Scale-free … extreme power concentration".
- **Add a Γ readout next to N_eff.** The rank-one "everyone listens to one hub" preset becomes a live demonstration that Γ and concentration are different things, which is a great teaching tool once §1.1 is fixed.
- **`book_cover.jsx` / the TikZ cover:** low-value cells render white on black because of the `blue!0` issue. The cover heatmap is ad-hoc rather than a real matrix; consider using Bolivia's actual (sourced) W^P.

---

## 9. Recommended order of work

1. **Decide the foundations** (§1.1–1.5). Write a 3–4 page "Conventions" section in Ch1 covering:
   - the normalization of W and the I-O orientation;
   - what π means and when;
   - Γ (mixing) vs N_eff / Var(π) (concentration), and Δ(W) defined as the pair;
   - the allocation rule derived from TU Nash;
   - a non-satiated utility calibration.
2. **Fix the source propositions** (Ch2 FWT and spectral bound, Ch6 Prop 3.2, star_power and Γ section, Ch8 App. 8B). Write the shared perturbation lemma. Then propagate the changes chapter by chapter; the Γ-dependency list in §1.1 is the checklist.
3. **Build a numerics repository.** It should hold one script per chapter that produces every matrix, number and heatmap in the text, plus the exercise solutions. This would have caught most of the arithmetic errors listed in §8, and it prevents regressions.
4. **Proof pass.** For each result: prove it, downgrade it, or delete it (§2).
5. **Data pass.** Source or tag every number, centralize Bolivia parameters in macros, and reconcile the cross-chapter drift (§3).
6. **Scholarship pass.** Bibliography, epigraph verification with page numbers, and thinker fidelity (§4).
7. **Production pass.** Single preamble, `\include` the chapters, fix the viz-panel title, `tabularx` the tables, use `\ref` everywhere, add a notation table (§5).
8. **Reframe the ending and the back cover** around "a unifying language and research programme" (§7).

Steps 1–2 are where the intellectual work is. Steps 3–8 are mostly mechanical once 1–2 are settled.
