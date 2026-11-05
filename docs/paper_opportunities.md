# Paper opportunities from the book

> Prior-art assessment of which parts of the revised manuscript could become standalone papers. It was produced with targeted literature searches (November 2026). Citations marked "from memory" or "unverified" were not re-checked. The working log follows the ranked summary.

# Paper-extraction assessment: "Structural Economics: The General Case"

Working log (appended per candidate as completed; final ranked report at end).

---
## C1. What eigenvectors and spectral gaps of IO tables / SAMs do and do not measure

**Material (file:line)**
- Ch20 `chapters/chapter20.tex:61-81`: SAM row-normalization; Prop `ch20:prop:pi_size` (π = y/1'y for a balanced SAM, one-line proof); `:83-93` (SAM ⊃ W^P); `:192-233` (example: full-SAM Γ 0.402 > I-O block 0.391/0.248; complex λ₂ from the circular flow; Prop `ch20:prop:more_structure`, no ordering between Γ_SAM and Γ_IO); `:246-270` Pyatt–Round decomposition (textbook; correct proof).
- Ch19 `chapters/chapter19.tex:118-141`: Prop `ch19:prop:lumpability` (exact lumpability ⇒ spec(W̄) ⊂ spec(W) ⇒ Γ_coarse ≥ Γ_fine; π aggregates by summation).
- Ch3 `chapters/chapter03.tex:317-353`: ρ(A) − |λ₂(A)| is a scale artefact; Prop `ch03:prop:block_spectral` (2×2: Γ = a+b, π = (b,a)/(a+b); Cheeger for reversible W^P).
- Ch17 `chapters/chapter17.tex:47-59, 84-89`: domestic Γ must report normalization (input shares 0.13 vs sales shares 0.075, illustrative).
- Ch16 `chapters/chapter16.tex:37-62` (diagonal treatment; frictionless gravity ⇒ W = 1s′, Γ = 1, π = s), `:132-144` (linkages from A only; for row-stochastic W, (I−γW)⁻¹1 = 1/(1−γ)), `:200-218` (ECI facts: M̃ eigenvalues in [0,1], π(M̃) ∝ diversity, ECI variance ≡ 1).
- Ch6 `chapters/chapter06.tex:383-462`: Γ definition, periodicity caveat, Prop `ch06:prop:gap-power-independent` (W = γ1π′ + (1−γ)I).
- Code: `code/ch20.py`, `code/ch19.py`, `code/ch16.py`, `code/ch17.py`, `code/ch03.py`.

**What is mathematically there.** Every statement is correct, and every one is a short textbook fact:
- π of a row-normalized *balanced* flow matrix equals throughput. This is the stationary distribution of a random walk on a circulation (Eulerian) network.
- Exact lumpability preserves eigenvalues (Kemeny–Snell 1960).
- (I − γW)⁻¹1 = 1/(1−γ)·1 for row-stochastic W.
- Rank-one W has a single non-zero eigenvalue.
- W = γ1π′ + (1−γ)I is a one-line construction.
- Cheeger's inequality.
- The ECI facts are in Mealy, Farmer & Teytelboym (2019).

None of this is a theorem a referee would credit as new.

**Prior art found**
- Wirkierman, Bianchi & Torriero (2021), "Leontief Meets Markov: Sectoral Vulnerabilities Through Circular Connectivity", *Networks and Spatial Economics*, doi:10.1007/s11067-021-09551-8 (OA: https://eprints-gro.gold.ac.uk/30489/).
  - They map a *closed* IO system into an ergodic regular Markov chain, following Bródy (1970), *Proportions, Prices and Planning*.
  - They state explicitly that "stationary probabilities would coincide with (normalised) observed activity levels".
  - They find that their mean-first-passage betweenness ranking is almost the same as a ranking by gross output.
  - **This pre-empts the book's headline SAM result (π = account size) and its "centrality ≈ size" warning.**
- Blöchl, Theis, Vega-Redondo & Fisher (2011), "Vertex centralities in input-output networks reveal the structure of modern economies", *Phys. Rev. E* 83, 046127 (WP: CESifo 3175, https://www.cesifo.org/DocDL/cesifo1_wp3175.pdf). Random-walk centralities on IO networks, with an explicit treatment of self-loops and open vs. closed systems.
- Cerina, Zhu, Chessa & Riccaboni (2015), "World Input-Output Network", *PLOS ONE* 10(7): e0134025 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4519177/). PageRank on WIOD.
- Kostoska, Mitikj, Jovanovski & Kocarev (2020), absorbing-Markov-chain approach to the world economy (https://arxiv.org/pdf/2003.05204).
- **Second-eigenvalue / aggregation literature in ESR itself.** This is the most important finding: the "Γ depends on aggregation" point is an active ESR debate.
  - Bródy (1997), "The Second Eigenvalue of the Leontief Matrix", *ESR* 9(3):253-258 (https://ideas.repec.org/a/taf/ecsysr/v9y1997i3p253-258.html).
  - Białas & Gurgul (1998), *ESR* 10(3):285-290.
  - Molnár & Simonovits (1998), *ESR* 10(1):79-82.
  - Bidard & Schatteman (2001), *ESR* 13(3):289-298.
  - Mariolis & Tsoulfidis (2014), "On Bródy's conjecture: theory, facts and figures about instability of the US economy", *ESR* 26(2):209-223.
  - Gurgul & Wójtowicz (2015), "On the Economic Interpretation of the Bródy Conjecture", *ESR* 27(1):122-131 (https://econpapers.repec.org/RePEc:taf:ecsysr:v:27:y:2015:i:1:p:122-131). They show by simulation, and on EU tables, that *aggregation lowers* |λ₂/λ₁|.
  - Shaikh, Nassif-Pires & Coronado (2023), "A new empirical contribution to an old theoretical puzzle…", *ESR* 35(1):118-135 (https://ideas.repec.org/a/taf/ecsysr/v35y2023i1p118-135.html). On 307 US matrices of different sizes, subdominant moduli *rise* with matrix size.
  - The book's lumpability proposition is the exact-aggregation special case of what this literature documents empirically.
- McNerney, Fath & Silverberg (2013), "Network structure of inter-industry flows", *Physica A* 392:6427-6441 (arXiv:1204.4122): community structure of 20 national IO networks.
- McNerney, Savoie, Caravelli, Carvalho & Farmer (2022), "How production networks amplify economic growth", *PNAS* 119(1) e2106031118 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8740592). The economically meaningful centrality is the output multiplier, i.e. a Leontief/Katz object.
- Mealy, Farmer & Teytelboym (2019), "Interpreting economic complexity", *Science Advances* 5(1) eaau1705 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6326748). ECI = spectral clustering / second eigenvector of the same M̃. The Ch16 ECI facts are known.
- Acemoglu, Carvalho, Ozdaglar & Tahbaz-Salehi (2012, *Econometrica*): the influence vector is Katz/Leontief, not Perron. Ch6/Ch16 already concede this.

**What is genuinely new (small).** No source I found collects these points as a checklist of pitfalls for IO/SAM practitioners. The pieces are scattered across Markov-chain IO, ECI and Bródy-conjecture papers. Potentially new or under-stated points:
- (a) **Normalization.** For a *balanced closed* SAM T, the receipt-share matrix D⁻¹T and the expenditure-share matrix TD⁻¹ are similar (TD⁻¹ = D·(D⁻¹T)·D⁻¹). So Γ does not depend on the choice in a closed SAM, but it does in an open IO block. The book does not state this. I could not find it stated in the IO literature, though Bródy (1970) may have it. Verify.
- (b) **Γ vs π.** In a near-rank-one (gravity / hub) trade or IO network, a high Γ coexists with concentrated π. So any "integration index" built on Γ will score hub dependence as integration. Trade-network papers that use spectral measures do not seem to make this point.
- (c) **Aggregation beyond exact lumpability.** Exact lumpability ⇒ Γ_coarse ≥ Γ_fine. For reversible chains the stationary-weighted collapsed chain also has a gap at least as large as the original. That is standard Markov-chain theory (Aldous–Fill; Levin–Peres–Wilmer, ch. 13); verify the exact lemma before citing. For *non-reversible* IO share matrices, no general ordering holds. A counterexample plus an empirical frequency on real tables would be a modest new contribution, and it ties directly to the Gurgul–Wójtowicz / Shaikh et al. findings.

**What it would take**
1. Reframe the note as a contribution to the ESR eigenvalue/aggregation debate, not as "structuralism meets networks". Cite Bródy, Gurgul–Wójtowicz, Shaikh et al. and Wirkierman et al. up front, and concede that π = size and λ₂-falls-with-aggregation are known.
2. Add one real theorem beyond the textbook. Candidates:
   - aggregation bounds for approximately lumpable, non-reversible share matrices, e.g. via Dobrushin τ₁ or a perturbation bound around the exactly lumpable case;
   - the normalization-similarity result for closed systems, with a clean characterization of when the input-share and sales-share Γ coincide (balanced/closed) and how far apart they can be when the system is open (bounds in terms of value-added and final-demand leak shares).
3. Do the empirics on public data: OECD ICIO 2023 (76 countries, 45 industries), WIOD 2016, and US BEA at 15/71/405 sectors for aggregation paths. Report, for each country:
   - Γ under both normalizations;
   - Γ across aggregation levels;
   - corr(π, gross-output share), which should be ≈1 when the system is closed;
   - corr(PageRank, gross output);
   - how much "integration" rankings change with the diagonal and normalization treatment.

   Ship replication code. The repository's `code/` scripts are numpy toy checks, not data pipelines.
4. Drop the Bolivia narrative and the "power" vocabulary. ESR referees will read "π = power" as a category error unless it is called "stationary throughput share".

**Effort.** 2-3 months for one person with IO data experience: about one month for the theory (aggregation/normalization bounds), 1-1.5 months for the data pipeline, plus writing. Length: a 15-20 page article, or a 6-8 page note.

**Venues and odds**
- *Economic Systems Research*: **moderate**, if framed as above with real data and at least one non-textbook result. ESR has a tradition of short notes on exactly these eigenvalue properties (Bródy, Białas–Gurgul, Molnár–Simonovits, Soklis 2008). As a pure pitfalls note with only the current propositions: **low**. Referees will cite Wirkierman et al. and Gurgul–Wójtowicz.
- *Structural Change and Economic Dynamics*: **low-moderate**. It needs an economic finding (e.g. that measured "integration" rankings of LatAm vs East Asia flip with normalization), not just measurement hygiene.
- *Journal of Economic Structures* (Springer OA, PAPAIOS) or *Networks and Spatial Economics*: **moderate-good** for a careful empirical/methodological version. Lower prestige.

**Verdict.** A viable but modest paper. Its best asset is that it plugs into an existing ESR conversation. Its main risk is that most of its content is already known.

**Numerical verification (my script: /tmp/claude-0/-home-user/f7a0cabe-4458-593f-968b-7c82f8f4f63e/scratchpad/agg_check.py)**
- (a) Confirmed. On a random balanced 6-account circulation, receipt-share and expenditure-share Γ agree to machine precision (0.5489).
- (c) Confirmed. With stationary-weighted 2-block aggregation of random *non-reversible* 6×6 row-stochastic matrices:
  - aggregation *lowered* Γ in 494 of 20,000 draws (about 2.5%), with a maximum drop of 0.21;
  - for reversible matrices it never did (0 of 5,000).

  So "aggregation raises Γ" holds under exact lumpability or reversibility, but not in general. This is a small, correct and publishable-as-a-note clarification of the Gurgul–Wójtowicz / Shaikh et al. aggregation finding.

---
## C2. Network-mediated aspiration/preference formation and welfare (Ch9–12)

**Material (file:line)**
- Ch9 `chapters/chapter09.tex:175-200`: formation equation (Def `ch09:def:formation`).
  - `:220-244`: steady state (Brouwer / affine).
  - `:322-345`: persistence and half-life (Prop `ch09:thm:hysteresis`).
  - `:384-417`: E/F decomposition and Ψ (Def `ch09:def:psi`, Remark `ch09:rem:psi`).
  - `:425-437`: block-wise steady states.
- Ch10 `chapters/chapter10.tex:40-75`: model (TU rule + affine returns + formation).
  - `:83-104`: J = ρI + μφW^I C P; J1 = ρ1; "only relative preferences are amplified" (Prop `ch10:prop:jacobian`).
  - `:364-385`: formation Leontief inverse (I − J)⁻¹ (Thm `ch10:thm:formation_effect`).
  - `:400-417`: one-loop Ψ(T) = κ(1−λ^T)/(1−λ) (Prop `ch10:prop:psi_lq`).
- Ch11 `chapters/chapter11.tex:42-104`: omitted-variable bias; APC non-identification; when Ψ̂_var maps to Ψ.
  - `:124-140`: epidemiological β = ρλ.
  - `:265-280`: reflection / Bramoullé–Djebbari–Fortin.
  - `:291-316`: multigenerational moments.
- Ch12 `chapters/chapter12.tex:52-66`: satisfaction paradox, gap ratio = 1 − ℳ.
  - `:273-287`: fixed-evaluation ranking by Var(π) / N_eff, not Γ.
  - `:289-300`: endogenous welfare Ω − ℳ c′x.
  - `:335-346`: satisfaction-bias corollary (Easterlin).
  - `:371-383`: endogenous-Pareto ⊂ laundered-Pareto (Prop `ch12:prop:consent`).

**Assessment of the mathematics.** Everything is correct after the revision, and everything is shallow. Each headline result is a one-line consequence of a modelling assumption:
- **J1 = ρ1.** This holds because the bargaining rule was *assumed* to depend on θ_i − θ̄. "Networks amplify only relative preferences" is that assumption restated, not a finding.
- **Ψ(T) closed form.** A geometric series in a scalar AR(1).
- **(I − J)⁻¹ as a "formation Leontief inverse".** The standard Neumann-series reading of any linear recursion.
- **Satisfaction gap ratio = 1 − ℳ.** Linear aspirations a = ℳW^I x, two blocks, quadratic v: algebra.
- **Fixed-evaluation ranking by Var(π).** Quadratic utility is mean–variance. Var(x) = Φ²Var(π) under the linear allocation rule.
- **Consent proposition.** A Neumann-series positivity argument. The economic point (others' gains are experienced as losses, so fewer changes command unanimity) is the standard relative-income externality.
- **Ch11.** A competent survey of textbook identification results (Manski 1993; Bramoullé et al. 2009; APC non-identification; the epidemiological approach) applied to the book's model. Useful pedagogy, not research.

**Prior art**
- Genicot & Ray (2017), "Aspirations and Inequality", *Econometrica* 85(2):489-519 (https://www.econometricsociety.org/doi/10.3982/ECTA13865). Aspirations formed from the social income distribution, with a feedback to investment and inequality. Much richer: endogenous investment, frustration region, steady-state polarization.
- Mookherjee, Napel & Ray (2010), "Aspirations, Segregation, and Occupational Choice", *JEEA* 8(1):139-168 (https://ideas.repec.org/a/tpr/jeurec/v8y2010i1p139-168.html). Neighbourhood aspirations plus segregation steady states: exactly the book's block-wise steady-state story, done with equilibrium wages.
- Dalton, Ghosal & Mani (2016), "Poverty and Aspirations Failure", *Economic Journal* 126(590):165-188 (from memory; not re-checked).
- Ghiglino & Goyal (2010), "Keeping Up with the Neighbors: Social Interaction in a Market Economy", *JEEA* 8(1):90-119 (https://ideas.repec.org/a/tpr/jeurec/v8y2010i1p90-119.html).
  - Network reference consumption in general equilibrium; prices and consumption depend on Bonacich centrality.
  - "When segregated communities become integrated, the poor lose utility while the rich gain."
  - This is the rigorous version of the book's satisfaction-paradox and integration-welfare point, with prices.
- Bisin & Verdier (2001, *JET*); Panebianco & Verdier (2017), "Paternalism, homophily and cultural transmission in random networks", *GEB* 105:155-176 (https://econpapers.repec.org/RePEc:eee:gamebe:v:105:y:2017:i:c:p:155-176). Network cultural transmission with homophily.
- Friedkin & Johnsen (1990, 1999). The anchored-DeGroot opinion model x_{t+1} = ΛWx_t + (I−Λ)x₀ is structurally the book's formation equation with ρ, the anchor and μW.
- Golub & Jackson (2010, 2012). Network learning, homophily and the speed of convergence.
- Endogenous-preference welfare:
  - Pollak (1978), *AER P&P* 68(2):374-379 (https://profiles.wustl.edu/en/publications/endogenous-tastes-in-demand-and-welfare-analysis/);
  - Hammond (1976), *JET* 13(2):329-340 (https://ideas.repec.org/a/eee/jetheo/v13y1976i2p329-340.html);
  - Gintis (1974);
  - von Weizsäcker (2005), "The welfare economics of adaptive preferences" (MPI Bonn, https://www.coll.mpg.de/51271/2005_11online.pdf);
  - Bernheim & Rangel (2009, *QJE*).
- Relative income and satisfaction: Luttmer (2005, *QJE*); Clark, Frijters & Shields (2008, *JEL*). The local-reference-group explanation of the Easterlin pattern is standard.

**What is new.** Very little. The cleanest new-looking item is the "laundered vs endogenous Pareto" nesting combined with the "fixed-evaluation ranks by N_eff, not Γ" remark. That is a footnote in a theory paper, not a paper.

**What it would take.** An actual theory contribution, for example:
- nonlinear returns with investment (as in Genicot–Ray), so that W^I shapes the *selection* among steady states;
- welfare criteria compared across those steady states, with W^I-dependent comparative statics that Ghiglino–Goyal and Mookherjee–Napel–Ray do not already deliver.

Or an empirical paper that estimates network-mediated aspiration formation with randomized ΔW^I. That needs field data the author does not have.

**Effort.** 6-12 months for the theory route, with uncertain novelty. More for the empirical route.

**Venues and odds**
- *Journal of Development Economics* / *JEBO* / *Journal of Economic Inequality* with the current content: **low**. Referees will cite Genicot–Ray and Ghiglino–Goyal and ask what is new.
- *Economics & Philosophy* or *Journal of Economic Methodology*, as a conceptual piece on welfare under network-formed aspirations and the "rubber ruler": **low-moderate**. It would need a much stronger philosophy-of-welfare literature review (Sen, Elster's sour grapes, Nussbaum, Bernheim–Rangel, Sugden).

**Verdict.** Reject as a standalone research paper. At most, the satisfaction-paradox and consent material could be a section inside C6.

---
## C3. "Spectral industrial policy": which links to build, and its political economy (Ch15, with Ch14)

**Material (file:line)**
- Ch15 `chapters/chapter15.tex:30-58`: operator design problem and three objectives (Γ, Var(π), Leontief output; the output objective is constant on row-stochastic W).
- `:64-79`: Prop `ch15:prop:pi_path` (W_t = (1−t)W + t1π′ keeps π and raises Γ).
- `:84-106`: gradient G_ij = −Re(λ̄₂u₂ᵢv₂ⱼ)/|λ₂|; Prop `ch15:thm:optimal_edge` (ℓ₁ ⇒ single within-row move; Frobenius ⇒ row-centred gradient).
- `:108-118`: "connect the cores".
- `:125-176`: ISI first-order condition and the ISI-A vs ISI-B example.
- `:188-193`: road placement.
- `:355-388`: political constraint Δπ_i ≥ 0 for i ∈ S*; Prop `ch15:prop:political_cost` (integration has zero political cost).
- Ch14 `chapters/chapter14.tex:60-95`: payoff = ΔΩ/N + ΦΔπ_i (exact), formation amplifier (assumed).
- `:315-330`: constant-sum power, empty core when ΔΩ ≡ 0.
- Ch6 `chapters/chapter06.tex:729ff`: perturbation lemma dπ′ = π′EZ.

**Assessment**
- **π-preserving path.** A two-line construction. Economically it is also nearly vacuous: W = 1π′ means every sector buys inputs in the *same* proportions π. It is a degenerate technology, reachable only with a budget ‖1π′ − W‖ of order one. So "integration costs no one power" is an artefact of measuring power as π and integration as Γ. It is not a political-economy result. A referee will say so in one paragraph.
- **ℓ₁ vs Frobenius edge selection.** Linear optimization over a cross-polytope vs a ball. Textbook.
- **The gradient formula.** Standard simple-eigenvalue perturbation (Kato; Stewart–Sun; Schweitzer 1968; Meyer 1980; Cho–Meyer 2001). The perturbation lemma is a correct textbook statement.
- **"Connect the cores".** This is the Fiedler-vector heuristic of Ghosh & Boyd (2006). Add the edge maximising (v_i − v_j)², i.e. join nodes with large opposite-signed Fiedler entries, which sit deep in the two clusters. Known.
- **Ch14.** Constant-sum Δπ and an empty core when ΔΩ ≡ 0 are immediate from Σπ = 1.
- **ISI-A vs ISI-B.** A 4×4 illustration. The one genuinely interesting qualitative point is in it: Γ, power and output objectives disagree, and the Leontief objective must use A, not W. But it is an example, not a result.

**Prior art**
- Liu (2019), "Industrial Policies in Production Networks", *QJE* 134(4):1883-1948, doi:10.1093/qje/qjz024 (https://ideas.repec.org/a/oup/qjecon/v134y2019i4p1883-1948..html).
  - The rigorous welfare-based version: distortion centrality (influence / Domar weight) is a sufficient statistic for which sector to promote.
  - Tested on 1970s Korea and on China.
  - There is an IPEA application to Brazil (https://repositorio.ipea.gov.br/handle/11058/18525), which is relevant given the author's IPEA background.
- Boyd, Diaconis & Xiao (2004), "Fastest Mixing Markov Chain on a Graph", *SIAM Review* 46(4):667-689. Convex/SDP for symmetric W.
- Ghosh & Boyd (2006), "Growing Well-Connected Graphs", *IEEE CDC* (https://web.stanford.edu/~boyd/papers/growing_fast_graphs.html). Greedy Fiedler-vector edge addition plus an SDP bound.
- PageRank/stationary-distribution manipulation:
  - de Kerchove, Ninove & Van Dooren (2008), "Maximizing PageRank via outlinks", *Linear Algebra Appl.* 429 (arXiv:0711.2867);
  - Csáji, Jungers & Blondel (2014), "PageRank optimization by edge selection", *Discrete Applied Math.* 169:73-87 (https://ar5iv.arxiv.org/html/0911.2280);
  - Fercoq, Akian, Bouhtou & Gaubert (2013, *IEEE TAC*).

  These already solve the "which edges move π" problem (Ω₂), with MDP formulations and complexity results.
- Production-network macro: Acemoglu, Carvalho, Ozdaglar & Tahbaz-Salehi (2012); Baqaee & Farhi (2019, *Econometrica*); Acemoglu & Azar (2020, *Econometrica*, endogenous production networks).
- Industrial policy: Juhász, Lane & Rodrik (2024, *Annual Review of Economics*); Lane (2025, *QJE*, Korea HCI).

**What is new.** Essentially nothing at the level of a research contribution. An *economically grounded* version would ask which input-substitution policies change measured mixing and measured centrality, using real IO tables and a real policy episode such as Korea HCI or Brazil's II PND. That would collide head-on with Liu (2019), which has the welfare foundation the spectral objective lacks. The spectral objective has no welfare justification: why should a planner care about the mixing time of input shares?

**What it would take.** A welfare foundation for Γ, e.g. Γ governs the speed at which sectoral productivity shocks are diffused or prices converge, inside a dynamic model with adjustment costs. Then show that it gives policy rankings different from Liu's distortion centrality on data. That is a new theory paper; the book supplies little of it.

**Effort.** 9-18 months and a co-author with production-network macro skills.

**Venues and odds**
- *JDE* / *SCED* / *Journal of Economic Structures* with the current content: **low**. "Known heuristics relabelled; π-preserving path is degenerate; no welfare criterion."
- An OR/control venue (e.g. *IEEE Transactions on Network Science and Engineering*): **low**. The mathematics is already there.

**Verdict.** Reject. The single defensible lesson is that integration and power objectives are distinct, and that multipliers must come from A, not W. It belongs in C1 or C6 as a remark.

---
## C4. Power-weighted TU bargaining with the L_power + L_externality + L_interaction decomposition (Ch7, Ch2)

**Material (file:line)**
- Ch7 `chapters/chapter07.tex:107-117`: Assumption 7.A (π as Nash weights is a postulate; Manea / Abreu–Manea cited against it).
- `:167-185`: Prop `ch07:prop:allocation-rule` (x_i = d_i + π_i(Ω − Σd)).
- `:191-198`: L_power = (κN/2)Φ²Var(π).
- `:226-238`: weighted Nash is Pareto efficient; the loss is distributional.
- `:259-283`: Thm `ch07:thm:decomposition`, with L_interaction = −γ_U N U′(Ω/N)ΦCov(c, π) + quadratic term.
- `:323-336`: comparative statics via dπ′ = π′EZ.
- Ch2 `chapters/chapter02.tex:141-230` (power, externality and block losses); `:387ff` (structural welfare loss).

**Assessment**
- The allocation rule is the textbook closed form of asymmetric Nash bargaining with transferable utility (Kalai 1977; Roth 1979).
- The decomposition is an exact second-order expansion of quadratic utility around equal division, weighted by 1 + γ_U c_j.
- The "offsetting interaction term" is −Cov(c, π) from that expansion. The chapter itself says it is an argument *against* evaluating with V. The interaction term is therefore a curiosity, not a result.
- The load-bearing economic content is Assumption 7.A, which is a postulate.

**Prior art**
- Manea (2011, *AER*), "Bargaining in Stationary Networks": payoffs are set by mutually estranged sets and shortage ratios, not by eigenvector centrality.
- Abreu & Manea (2012, *GEB* 75(1):1-16), "Markov equilibria in a model of bargaining in networks" (https://econpapers.repec.org/RePEc:eee:gamebe:v:75:y:2012:i:1:p:1-16).
- Elliott & Golub (2019, *JPE*), "A Network Approach to Public Goods": Lindahl contributions are proportional to eigenvector centrality of the benefits matrix. This is the closest existing *derivation* of an eigenvector-centrality share.
- C. Emerson Melo Sanchez, Caltech PhD thesis (2013, adviser Echenique): a networked-market bargaining chapter in which equilibrium payoffs coincide with Bonacich centrality (https://mathgenealogy.org/id.php?id=227663; summary at https://independent.academia.edu/EmersonSanchez33). I could not verify the exact model or whether it was published.
- Ballester, Calvó-Armengol & Zenou (2006, *Econometrica*): Bonacich centrality in network games.

**What would make it a paper.** Only a *microfoundation* of Assumption 7.A: a non-cooperative bargaining or exchange game on a row-stochastic W whose limit payoffs are the Perron vector π, with comparative statics under ΔW. That would be a genuine theory contribution, positioned between Manea and Elliott–Golub. The book has none of it, and it is a hard, uncertain theory project (6-12+ months for a network-theory specialist).

**Venues and odds**
- As is: **low** anywhere. Algebra on a postulate.
- With a microfoundation: *GEB* / *JET* / *TE*, **moderate**. That would be a different paper by a different skill set.

**Verdict.** Reject as a standalone paper.

---
## C5. Truncated industrialization as missing bridges, not missing dense modules (Ch25)

**Material (file:line)**
- Ch25 `chapters/chapter25.tex:29-41`: densest-subgraph definition (Goldberg / Charikar).
- `:44-53`: product-space core (Stylized Fact; numbers illustrative).
- `:58-65`: truncation index τ_c (illustrative values).
- `:67-68`: product-space proxy W_c = D_c⁻¹φ_c.
- `:71-82`: Prop `ch25:prop:bridge_bound` (Γ ≤ 1 − λ₂ ≤ 2Φ_c(S): easy half of Cheeger).
- `:84-89`: Remark: bridge vs module (Γ 0.548 → 0.009 when a bridge core is removed; 0.218 → 0.800 when a dense module is removed; densest-subgraph algorithms cannot tell them apart).
- `:91-94`: Conjecture `ch25:thm:truncation` + test design (participation coefficient, conductance with/without D, regress Γ_c on τ_c).
- `:176-196`: capital-goods-as-bridge claims (betweenness / participation).
- `:231-243`: Fajnzylber mechanisms.
- Related: Ch16 `chapters/chapter16.tex:231-242` (ECI→growth conjecture mediated by domestic Γ; testable with IO tables).
- Code: `code/ch25.py`.

**Assessment**
- **Theory.** The proposition is the easy direction of Cheeger's inequality (Lawler–Sokal 1988; Sinclair–Jerrum 1989): textbook. The bridge/module remark is correct and is a useful corrective to "dense core ⇒ integration" readings of the product space. It is still an observation, not a theorem.
- **The real content is an empirical hypothesis:**

  > the capital-goods/machinery/chemicals core of the product space (and of IO networks) is a *bridge set* whose absence leaves commodity exporters' activity clusters connected only through low-conductance cuts.

  That hypothesis is testable with public data.
- **Weaknesses.**
  - (i) The product-space proxy treats co-export proximity as if it were a production network. Hidalgo et al. are explicit that φ is capability similarity, and the book itself says so (`chapters/chapter16.tex:177`). A referee will insist the test be run on national IO tables, with the product space as secondary.
  - (ii) Γ of a country-specific subgraph of φ depends mechanically on how many products the country exports (diversity), so τ_c and Γ_c are confounded by size and diversity. This needs null models, e.g. random product subsets of equal size, degree-preserving rewiring, and controls for diversity/ECI.
  - (iii) The link from "low Γ" to growth or inequality outcomes is itself a conjecture (Ch26). The paper would have to stop at the structural fact, or add a reduced-form growth regression with all the usual identification caveats.

**Prior art**
- Hidalgo, Klinger, Barabási & Hausmann (2007), "The Product Space Conditions the Development of Nations", *Science* 317:482-487 (arXiv:0708.2090). Core–periphery, and the "hard to reach the core" diffusion argument. The core-as-hub idea is there qualitatively; *bridge vs module* and *conductance* are not.
- Hausmann & Klinger (2006/2007) on product-space structure and the evolution of comparative advantage. Hausmann, Hidalgo et al. (2014), *The Atlas of Economic Complexity*. Mealy, Farmer & Teytelboym (2019).
- Fadinger, Ghiglino & Teteryatnikova (2022), "Income Differences, Productivity, and Input-Output Networks", *AEJ: Macro* 14(2):367-415 (https://www.aeaweb.org/doi/10.1257/mac.20180342). IO structure and cross-country income on 38 countries; highly connected sectors are relatively more productive in poor countries.
- McNerney et al. (2022, *PNAS*): longer production chains, i.e. higher output multipliers, predict growth on WIOD.
- Bartelme & Gorodnichenko (2015), "Linkages and Economic Development", NBER WP 21251 (from memory; not re-checked). IO linkage strength correlates with income.
- Mutreja, Ravikumar & Sposi (2018), "Capital goods trade, relative prices, and economic development", *Review of Economic Dynamics* 27:101-122 (https://econpapers.repec.org/RePEc:red:issued:17-216). Frictionless capital-goods trade closes about 40% of the rich-poor income gap. This is a competing, *trade*-based account of why a missing domestic capital-goods sector matters less than structuralists think.
- Hypothetical-extraction literature: Dietzenbacher, van der Linden & Steenge (1993, *ESR*); Miller & Lahr (2001). This is the standard IO way to "remove the capital-goods core" and measure the consequences. The conductance/Γ version would be a spectral variant of it.
- Historical/structuralist: Fajnzylber (1983, 1988, 1990 "casillero vacío"); Tavares (1964); Amsden (2001) *The Rise of "the Rest"*; Lane (2025) on Korea HCI.

**What is new.** If the data cooperate:
- (a) Show that the product-space/IO "core" is a bridge rather than a module: high participation coefficient, and a conductance collapse on extraction relative to a matched null.
- (b) Show that countries lacking it have lower domestic IO mixing *conditional on diversity and size*.

Point (a) has value even alone. The densest-subgraph framing used informally in the complexity literature does not distinguish the two cases, and the 2-example remark in Ch25 shows why it matters.

**What it would take**
1. **Data.**
   - Atlas/BACI HS6 or SITC4 export data → RCA → φ (public, Growth Lab Dataverse).
   - OECD ICIO 2023 domestic tables (76 countries × 45 industries, 1995-2020) and/or Eora26 for Bolivia and other small LatAm economies. Eora's quality for Bolivia is poor; say so.
   - National supply-use tables at a finer level (e.g. INE Bolivia, IBGE Brazil, INEGI Mexico) for robustness. 45 industries is coarse for a "capital-goods core" claim.
2. **Methods.**
   - Participation coefficient and betweenness of core products.
   - Conductance and Γ (lazy/reversible symmetrization for the product space, where Cheeger applies; Dobrushin τ₁ for non-reversible IO shares).
   - Hypothetical extraction of capital-goods industries in IO tables.
   - Matched null models of equal size and degree sequence.
   - Diversity/ECI controls.
   - LatAm vs East Asia event studies: Korea 1973-90 vs Brazil 1974-90 (II PND) vs Mexico (maquila).
3. **Positioning.** Present it as an empirical test of a Fajnzylber/Tavares hypothesis with modern network tools. Put the Cheeger bound in a short section and the null-model design up front. Drop the "truncation explains 55-86% of Bolivia's integration deficit" calculation (Ch25 `:205-224`). It rests on an unobserved counterfactual Γ(0) and illustrative inputs.

**Effort.** About 4-6 months with a research assistant; 6-9 months solo. The pipeline (RCA, proximity, ICIO parsing, nulls) is standard but fiddly. There is a risk of a null or ambiguous result. Mitigation: pre-specify the tests, and write the paper so that "the core is a module, not a bridge" is also a publishable finding, because it would cut against a common reading of the product space.

**Venues and odds**
- *Structural Change and Economic Dynamics*: **moderate**. Good topical fit (structural change, product space, IO), and receptive to structuralist framing with network methods.
- *CEPAL Review* / *Revista CEPAL*: **good**, if written for that audience (Fajnzylber connection, LatAm focus, Spanish/English). Lower impact.
- *World Development* or *Journal of Development Studies*: **low-moderate**. They would need a clearer development-outcome link.
- *Research Policy* / *Industrial and Corporate Change*: **low-moderate**. Possible if the capability/technology angle is strengthened.
- *Economic Systems Research*: **low-moderate**, as an IO-only version (hypothetical extraction + conductance).

**Verdict.** The best *empirical* paper in the book: a clear, falsifiable hypothesis, a tight link to the author's comparative advantage (LatAm structuralism plus development-bank data experience), and public data. The current chapter supplies the question and the counterexample logic, not the evidence.

---
## C6. A formal translation of Latin American structuralism into network economics (book as a whole; Ch1, Ch28)

**Material (file:line)**
- Ch1 `chapters/chapter01.tex:513-560`: the implicit market interaction matrix W^mkt = 1s′; W^CE as the equal-endowment case.
- `:660-755`: conventions; π uniform ⇔ doubly stochastic; allocation rule.
- Ch28 `chapters/chapter28.tex:20-41`: framing principle (two independent coordinates d_π and 1 − Γ).
- `:43-91`: "Book in one page" table (result + primitives for each chapter). A ready-made skeleton for a mapping table.
- `:93-135`: Furtado's three claims → formal content, with "reinterpretation, not proof".
- `:177-197`: "What the book does not do". Unusually candid; a referee asset.
- Mapping material across chapters:
  - Pinto → cross-block mass / conductance (Ch3 `:237-353`; Ch17 `:134-188`).
  - Prebisch centre–periphery → concentration of π, not Γ (Ch16 `:57-62`).
  - Furtado / Nurkse demonstration effect → W^I edges to the centre (Ch9 `:448-457`; Ch12 `:300-301`).
  - Fajnzylber → bridge vs module (Ch25).
  - Hirschman linkages → must come from A, not W (Ch15 `:194-202`).
  - Sunkel / Cardoso–Faletto → coalitions (Ch14).
  - Gramsci → hegemony as formation-aligned π^I (Ch24).
  - Taylor → SAM multipliers (Ch20).
- Bibliography: `references.bib` has 365 entries; chapters now use natbib. Quotations still flagged in `docs/verification_checklist.md` (the TODO-VERIFY items).

**Assessment.** This is the most *distinctive* output of the book. I found no existing paper that maps structuralist propositions systematically onto network objects. The closest work formalizes structuralism with *macro* models (BoP-constrained growth, North–South, structuralist CGE), or links it to economic complexity empirically.

The essay's strongest thesis comes out of the corrections themselves:

> Latin American structuralism makes two logically independent claims about social structure. Centre–periphery dependence is about the **concentration** of interaction weight (π, N_eff). Pinto/Furtado heterogeneity is about weak **coupling** (conductance, Γ). A network formalization that collapses them into one index (as the original draft did with Γ) mis-scores hub dependence as integration. Frictionless gravity trade is fully "integrated" (Γ = 1) and maximally centre-dependent.

That is a clean, non-trivial conceptual point, supported by simple but correct propositions (Ch6 gap-power independence; Ch16 gravity rank-one; Ch3 2×2 case; Ch25 bridge vs module).

**Prior art**
- Formal structuralism (macro):
  - Taylor (1983), *Structuralist Macroeconomics*; Taylor (1990), *Socially Relevant Policy Analysis*;
  - Ocampo, Rada & Taylor (2009), *Growth and Policy in Developing Countries: A Structuralist Approach*;
  - Ros (2000; 2013);
  - Cimoli & Porcile (2014), "Technology, structural change and BOP-constrained growth: a structuralist toolbox", *CJE* 38(1):215-237, doi:10.1093/cje/bet020 (https://ideas.repec.org/a/oup/cambje/v38y2014i1p215-237..html);
  - Findlay (1980) and Dutt (1990) North–South models.
- Interpretive and history of thought:
  - Rodríguez (2006), *El estructuralismo latinoamericano*; an English UN paper is at https://www.un.org/esa/usg_ocampo/articles/pdf/lcg2150i_Rodriguez.pdf;
  - Bielschowsky (1998), "Cincuenta años del pensamiento de la CEPAL";
  - Love (1980, *LARR*) on Prebisch's origins.
- Complexity ↔ structuralism: Gala, Rocha & Magacho (2018), "The structuralist revenge: economic complexity as an important dimension to evaluate growth and development", *Brazilian Journal of Political Economy* 38(2):219-236, doi:10.1590/0101-31572018v38n02a01 (https://econpapers.repec.org/RePEc:ekm:repojs:v:38:y:2018:i:2:p:219-236:id:63). Also Hartmann et al. (2017, *World Development*) on ECI and inequality.
- Network allies the essay must credit: Acemoglu et al. (2012); Golub & Jackson (2010, 2012); Liu (2019); Ghiglino & Goyal (2010); Elliott & Golub (2019).

**What it would take**
1. Cut to about 9,000-12,000 words. Use one mapping table: claim → formal object → what is proved → what is conjecture → fidelity caveat ("where the thinker would disagree").
2. Lead with the two-coordinate thesis and the negative results; they are what make it more than relabelling.
3. Fix fidelity before submission. `docs/review.md` §4.3 lists serious misreadings: Prebisch vs Emmanuel, Marini, Skocpol, Fajnzylber's *casillero vacío*, Cardoso–Faletto, Ocampo's pillars. Check that each chapter's corrected text is right, and clear the TODO-VERIFY quotations, or the essay will be rejected by a structuralist referee on scholarship alone.
4. Drop every Bolivia number tagged *Illustrative*. An essay can be number-free.
5. Decide the audience:
   - for CEPAL/BJPE: Spanish or Portuguese is an option, the tone is policy-relevant, and Ocampo/Bárcena-era neo-structuralism should be engaged;
   - for CJE/ROKE: the post-Keynesian structuralist toolbox (Cimoli–Porcile, Taylor) should be engaged, and the essay must say what the network view explains that BoP-constrained models do not.

**Effort.** 6-10 weeks of writing and checking, given the book. Most of the cost is verifying the thinkers and quotations.

**Venues and odds**
- *CEPAL Review* / *Revista CEPAL*: **good**. Natural home; receptive to formal reinterpretations of Prebisch and Furtado; the author's IDB/World Bank/IPEA profile fits.
- *Brazilian Journal of Political Economy (Revista de Economia Política)*, *El Trimestre Económico*, *Investigación Económica*: **good** / **moderate-good**.
- *PSL Quarterly Review*: **moderate**. It publishes history-of-thought plus formal-reinterpretation essays.
- *Review of Keynesian Economics*, *Cambridge Journal of Economics*: **low-moderate**. Referees will push on "what new explanation?" and on fidelity.
- *Structural Change and Economic Dynamics*: **low** for a pure essay.
- *Journal of Economic Surveys*: **low**. It is not a literature survey; it would need comprehensive coverage of formal structuralist models.
- *History of Political Economy*: **low**. HOPE wants archival or textual history-of-thought scholarship, not formal translation, and the book's quotation record is a liability there.

**Verdict.** A good, realistic, low-effort publication. Lower prestige than C1 or C5, but the most "on brand" for the author.

---
## C7. Quick checks (rejected as standalone papers)

- **Ch23: time-varying W** (`chapters/chapter23.tex:66-114`).
  - The counterexample (Γ(X) = 0.85, Γ(Y) = 0.776, Γ(XY) = 0.042) is a nice teaching device.
  - The Dobrushin properties, the Hajnal condition and absolute-probability sequences are textbook: Seneta (2006) *Non-negative Matrices and Markov Chains*; Hajnal (1958); Kolmogorov (1936); Touri & Nedić (2014).
  - "Realized-product gap" is just |λ₂| of a product.
  - Nothing new. Reject. At most a paragraph in C1 ("report τ₁ or the realized-product gap for IO time series").
- **Ch24: hegemony with conformity and a saddle-node** (`chapters/chapter24.tex:111-163`).
  - A logistic best-response map with a saddle-node. This is Granovetter (1978) thresholds, Brock & Durlauf (2001) binary choice with social interactions, and Kuran (1989/1995) preference falsification, in a different vocabulary.
  - The proof is fine but standard.
  - The Bolivian N_eff^I series is illustrative.
  - Reject. It could be a vignette in C6.
- **Ch18: homophily** (`chapters/chapter18.tex:66-150`).
  - π ∝ degree for symmetric K: the textbook stationary distribution of a reversible random walk.
  - The minority power ratio (p + λ(1−p))/((1−p) + λp) is just a degree ratio. Network-science prior art shows homophily plus minority size ⇒ ranking disadvantage:
    - Avin et al. (2015, ITCS), "Homophily and the glass ceiling effect in social networks";
    - Karimi, Génois, Wagner, Singer & Strohmaier (2018, *Scientific Reports*), "Homophily influences ranking of minorities in social networks";
    - Golub & Jackson (2012, *QJE*) on homophily and convergence speed.
  - "Intersectionality as an interaction term" is a modelling remark.
  - Reject.
- **Ch11: formation vs enforcement identification** (`chapters/chapter11.tex:42-140, 265-316`).
  - A correct, well-referenced synthesis of APC non-identification, omitted-variable bias, the epidemiological approach (Fernández), reflection (Manski 1993) and network identification (Bramoullé et al. 2009).
  - The multigenerational-moment result is the latent-AR(1) logic of multigenerational mobility (Clark 2014; Braun & Stuhler 2018, *EJ*).
  - The variance-ratio statistic Ψ̂_var is (correctly) shown to be not Ψ.
  - Reject as research. Good teaching material.
- **Ch19: power nesting** π_I = π_G·π_{N|c}·π_within. This is stochastic complementation / aggregation–disaggregation (Meyer 1989, *SIAM Review*; Simon & Ando 1961). Textbook.
- **Ch27: engagement ranking lowers mixing** (`chapters/chapter27.tex:70-95`). A corollary of the Ch18 two-group formula. It would need real platform data, and the review notes that WhatsApp is not a recommender. Reject.

---
---
# FINAL RANKED REPORT

**Bottom line.** The corrected book contains almost no research-grade theorems. After the revision, every proved result is a short textbook fact (Perron–Frobenius, lumpability, Cheeger, the asymmetric-Nash closed form, geometric series, quadratic-utility algebra), a modelling assumption restated, or an illustrative calculation. The publishable value lies in three things:
- **(i)** measurement clarifications for IO/SAM practitioners, which plug into an existing ESR debate;
- **(ii)** one falsifiable empirical hypothesis (Ch25);
- **(iii)** the conceptual translation of structuralism, which is distinctive.

The "theory" candidates (C2 formation/welfare, C3 spectral industrial policy, C4 bargaining decomposition) should not be submitted as standalone papers.

## 1. ESR methodological note on spectral measures of IO tables and SAMs (from C1)

**Working title.** "What Perron Vectors and Spectral Gaps of Input-Output Tables and SAMs Measure: Size, Normalization and Aggregation"

**Contribution.** Eigenvector centrality, PageRank and spectral-gap "integration" indices are increasingly computed on IO tables, SAMs and trade-share matrices. The paper shows what they mechanically measure:
- **Size.** In any balanced (closed) flow matrix the Perron vector of the row-normalized table is exactly the account-throughput share. Centrality rankings of SAMs are therefore rankings by size.
- **Normalization.** In closed systems the receipt-share and expenditure-share normalizations are similar matrices and give the same spectral gap. In open IO blocks they do not, and the gap can differ substantially; the paper gives bounds in terms of leakage shares.
- **Aggregation.** Under exact lumpability, or for reversible chains, aggregation weakly raises the gap. For general non-reversible share matrices it need not: in my simulation about 2.5% of random 6×6 cases fall, by up to 0.21. This qualifies the Bródy-conjecture literature's finding that subdominant moduli rise with disaggregation.
- **Concentration vs mixing.** Mixing speed and concentration are independent (W = γ1π′ + (1−γ)I). Near-rank-one, hub-dependent structures such as frictionless gravity trade (W = 1s′) have Γ = 1 with π equal to export shares. So spectral "integration" scores reward hub dependence.
- **Linkages.** Rasmussen/Hirschman linkages are degenerate on row-stochastic matrices and must come from A.

Everything is demonstrated on OECD ICIO 2023, WIOD 2016 and multi-size US BEA tables, with replication code.

**Closest prior art**
- Wirkierman, A.L., Bianchi, M. & Torriero, A. (2021), "Leontief Meets Markov: Sectoral Vulnerabilities Through Circular Connectivity", *Networks and Spatial Economics*, doi:10.1007/s11067-021-09551-8, https://eprints-gro.gold.ac.uk/30489/. Closed IO → regular Markov chain; stationary probabilities = normalized activity levels, following Bródy (1970). **Pre-empts the "π = size" result. Must be cited and conceded.**
- Gurgul, H. & Wójtowicz, T. (2015), "On the Economic Interpretation of the Bródy Conjecture", *ESR* 27(1):122-131, https://econpapers.repec.org/RePEc:taf:ecsysr:v:27:y:2015:i:1:p:122-131. Aggregation lowers |λ₂/λ₁|.
- Shaikh, A., Nassif-Pires, L. & Coronado, J.A. (2023), "A new empirical contribution to an old theoretical puzzle…", *ESR* 35(1):118-135, https://ideas.repec.org/a/taf/ecsysr/v35y2023i1p118-135.html. 307 US matrices; subdominant moduli rise with size.
- Bródy, A. (1997), "The Second Eigenvalue of the Leontief Matrix", *ESR* 9(3):253-258, https://ideas.repec.org/a/taf/ecsysr/v9y1997i3p253-258.html.
- Also: Białas & Gurgul (1998, *ESR*); Molnár & Simonovits (1998, *ESR*); Mariolis & Tsoulfidis (2014, *ESR* 26(2)).
- Blöchl, F., Theis, F.J., Vega-Redondo, F. & Fisher, E.O'N. (2011), "Vertex centralities in input-output networks reveal the structure of modern economies", *Phys. Rev. E* 83:046127, https://www.cesifo.org/DocDL/cesifo1_wp3175.pdf.
- Cerina, F., Zhu, Z., Chessa, A. & Riccaboni, M. (2015), "World Input-Output Network", *PLOS ONE* 10(7):e0134025, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4519177/.
- McNerney, J., Fath, B.D. & Silverberg, G. (2013), *Physica A* 392:6427-6441.
- McNerney, J., Savoie, C., Caravelli, F., Carvalho, V.M. & Farmer, J.D. (2022), *PNAS* 119(1):e2106031118, https://pmc.ncbi.nlm.nih.gov/articles/PMC8740592.
- Mealy, P., Farmer, J.D. & Teytelboym, A. (2019), "Interpreting economic complexity", *Science Advances* 5(1):eaau1705, https://pmc.ncbi.nlm.nih.gov/articles/PMC6326748.

**What is new.**
- The closed-system normalization-invariance statement and its open-system failure.
- The reversible vs non-reversible aggregation result, with a counterexample and its empirical incidence.
- The concentration/mixing independence applied to trade and IO "integration" indices.
- A single practitioner checklist with real-data magnitudes.

Each piece is small. The package is a legitimate ESR note.

**Gaps and effort.**
- One genuinely non-textbook proposition: open-system normalization bounds, or an approximate-lumpability bound.
- The full data pipeline.
- Drop all "power"/structuralist vocabulary: call π the "stationary throughput share".
- Verify the reversible-collapse lemma (Levin–Peres–Wilmer ch. 13 / Aldous–Fill) before citing it.
- About 2-3 months.

**Venues and odds**
- *Economic Systems Research*: **moderate**. It has a tradition of short eigenvalue notes; with only the book's current propositions, **low**.
- *Journal of Economic Structures* or *Networks and Spatial Economics*: **moderate-good**.
- *SCED*: **low-moderate**, unless an economic finding emerges, e.g. that LatAm vs East Asia integration rankings flip with normalization.

**Chapters supplying material.**
- Ch20 (`chapter20.tex:61-81, 192-233`)
- Ch19 (`chapter19.tex:118-141`)
- Ch3 (`chapter03.tex:317-353`)
- Ch17 (`chapter17.tex:47-89`)
- Ch16 (`chapter16.tex:37-62, 132-144, 200-218`)
- Ch6 (`chapter06.tex:383-462`)
- Ch23 (`chapter23.tex:66-114`), for time series

## 2. Empirical test: is the capital-goods core a bridge or a module? (from C5)

**Working title.** "Truncated Industrialization as Missing Bridges: Testing Fajnzylber's Capital-Goods Hypothesis with Product-Space and Input-Output Networks"

**Contribution.** Fajnzylber and Tavares argued that Latin American industrialization was "truncated" because it never built a domestic capital-goods core. The product-space literature locates such a core: machinery, metals, chemicals. But a *dense* module and a *bridge* set have opposite implications for integration. Removing a bridging core collapses the conductance between activity clusters. Removing a dense but peripheral module can *raise* mixing. Densest-subgraph methods cannot tell the two apart.

The paper:
- measures, on HS/SITC product-space data and OECD ICIO / Eora domestic IO tables, whether the core is a bridge (participation coefficient, betweenness, conductance collapse under hypothetical extraction, against size- and degree-matched null models);
- tests whether countries lacking the core have lower domestic IO mixing, conditional on diversity and size;
- contrasts Korea's HCI drive with Brazil's II PND and Mexico's maquila path.

**Closest prior art**
- Hidalgo, C.A., Klinger, B., Barabási, A.-L. & Hausmann, R. (2007), "The Product Space Conditions the Development of Nations", *Science* 317:482-487, https://arxiv.org/abs/0708.2090.
- Mealy, Farmer & Teytelboym (2019) (above).
- Fadinger, H., Ghiglino, C. & Teteryatnikova, M. (2022), "Income Differences, Productivity, and Input-Output Networks", *AEJ: Macro* 14(2):367-415, https://www.aeaweb.org/doi/10.1257/mac.20180342.
- McNerney et al. (2022, *PNAS*) (above).
- Mutreja, P., Ravikumar, B. & Sposi, M. (2018), "Capital goods trade, relative prices, and economic development", *RED* 27:101-122, https://econpapers.repec.org/RePEc:red:issued:17-216. A competing trade-based view.
- Gala, Rocha & Magacho (2018, *BJPE*) (below).
- Hypothetical extraction: Dietzenbacher, van der Linden & Steenge (1993, *ESR*); Miller & Lahr (2001). Not re-checked online.
- Lane (2025, *QJE*) on Korea HCI. Not re-checked.

**What is new.** The bridge-vs-module distinction (Ch25's counterexample logic) turned into an empirical test with null models, linking a classic structuralist hypothesis to measurable network properties. Either answer is informative. "The core is a module, not a bridge" would cut against a common reading of the product space.

**Gaps and effort.**
- Drop the product-space-as-production-network proxy as the main object, or justify it. Use IO tables as primary.
- Control mechanically for diversity and size.
- Replace the Ch25 "55-86% of Bolivia's deficit" calculation, which rests on an unobserved counterfactual.
- Eora quality for Bolivia is weak. ICIO at 45 industries is coarse; national supply-use tables (IBGE, INEGI, INE) are needed for robustness.
- About 4-6 months with a research assistant (6-9 months solo). Real risk of an ambiguous result. Pre-specify the tests.

**Venues and odds**
- *Structural Change and Economic Dynamics*: **moderate**.
- *CEPAL Review*: **good**.
- *World Development* / *Journal of Development Studies*: **low-moderate**, needing an outcome link.
- *Research Policy* / *Industrial and Corporate Change*: **low-moderate**.

**Chapters supplying material.**
- Ch25 (`chapter25.tex:29-94, 176-243`)
- Ch16 (`chapter16.tex:157-242`)
- Ch3 (`chapter03.tex:331-353`)
- Ch17 (`chapter17.tex:43-89`)
- `code/ch25.py` (toy bridge/module examples only)

**Synergy.** Build one ICIO data pipeline and use it for both #1 and #2.

## 3. Essay: structuralism's two dimensions in network form (from C6)

**Working title.** "Concentration and Coupling: A Network Formalization of Latin American Structuralism, and What It Gets Wrong"

**Contribution.** The essay maps core structuralist propositions to network objects, with explicit fidelity caveats:
- Prebisch's centre–periphery;
- Pinto's structural heterogeneity;
- Furtado's demonstration effect and self-reproducing underdevelopment;
- Hirschman's linkages;
- Fajnzylber's truncation;
- Taylor's SAM multipliers.

Its central thesis is that structuralism makes two logically independent claims:
- **dependence = concentration** of interaction weight (Perron vector, N_eff);
- **heterogeneity = weak coupling** (conductance, spectral gap).

Network indices that fuse them misread hub dependence as integration. Frictionless gravity trade is fully "integrated" and maximally centre-dependent. The essay states precisely which translations are theorems (few, simple), which are modelling postulates (π as bargaining power), and which are open conjectures (truncation, growth-with-equity thresholds). It closes with a research agenda.

**Closest prior art**
- Cimoli, M. & Porcile, G. (2014), "Technology, structural change and BOP-constrained growth: a structuralist toolbox", *CJE* 38(1):215-237, https://ideas.repec.org/a/oup/cambje/v38y2014i1p215-237..html.
- Ocampo, J.A., Rada, C. & Taylor, L. (2009), *Growth and Policy in Developing Countries: A Structuralist Approach* (Columbia UP).
- Taylor, L. (1983), *Structuralist Macroeconomics*.
- Rodríguez, O. (2006), *El estructuralismo latinoamericano*; English UN paper at https://www.un.org/esa/usg_ocampo/articles/pdf/lcg2150i_Rodriguez.pdf.
- Gala, P., Rocha, I. & Magacho, G. (2018), "The structuralist revenge…", *BJPE* 38(2):219-236, https://econpapers.repec.org/RePEc:ekm:repojs:v:38:y:2018:i:2:p:219-236:id:63.
- Network allies: Acemoglu et al. (2012); Golub & Jackson (2010, 2012); Liu (2019); Ghiglino & Goyal (2010).

I found no prior systematic *network* translation.

**Gaps and effort.**
- Fidelity: the `docs/review.md` §4.3 misreadings and the TODO-VERIFY quotations must be resolved.
- Cut to 9-12k words.
- Remove the illustrative Bolivia numbers.
- Engage BoP-constrained structuralism, and say what the network view adds.
- About 6-10 weeks.

**Venues and odds**
- *CEPAL Review*: **good**.
- *Brazilian Journal of Political Economy* / *El Trimestre Económico*: **good / moderate-good**.
- *PSL Quarterly Review*: **moderate**.
- *ROKE* / *CJE*: **low-moderate**.
- *JES*, *HOPE*: **low**.

**Chapters supplying material.**
- Ch1 (`chapter01.tex:513-560, 660-755`)
- Ch28 (`chapter28.tex:20-135, 177-227`)
- Ch6 (`chapter06.tex:424-462`)
- Ch16 (`chapter16.tex:57-62`)
- Ch3 (`chapter03.tex:237-353`)
- Ch9 (`chapter09.tex:448-457`)
- Ch12 (`chapter12.tex:52-66, 371-390`)
- Ch25 (`chapter25.tex:84-89`)

## Rejected candidates

- **C2 (formation/aspirations and welfare, Ch9–12).** Correct but shallow.
  - J1 = ρ1 restates the assumption that allocations depend on θ_i − θ̄.
  - Ψ(T) is a geometric series.
  - The satisfaction-gap factor (1 − ℳ) and the Var(π) ranking are quadratic-utility algebra.
  - The consent nesting is a Neumann-series argument.
  - Genicot & Ray (2017, *Econometrica*), Mookherjee, Napel & Ray (2010, *JEEA*) and Ghiglino & Goyal (2010, *JEEA*; integration hurts the poor and helps the rich under network reference consumption) already do richer versions. Ch11's identification material is a good textbook synthesis (APC, Manski, Bramoullé et al.), not research.
  - **Low** odds anywhere without a new nonlinear theory or randomized ΔW^I data.
- **C3 (spectral industrial policy, Ch15/14).**
  - The π-preserving path to W = 1π′ is a two-line construction toward a degenerate technology in which every sector buys the same input mix. "Integration has zero political cost" is an artefact of the definitions.
  - The ℓ₁/Frobenius results are LP-vs-ball facts. "Connect the cores" is the Ghosh–Boyd (2006) Fiedler heuristic.
  - Stationary-distribution manipulation is solved in the PageRank-optimization literature (de Kerchove, Ninove & Van Dooren 2008; Csáji, Jungers & Blondel 2014).
  - Liu (2019, *QJE*) provides the welfare-founded version that the spectral objective lacks.
  - **Low**.
- **C4 (power-weighted TU bargaining and loss decomposition, Ch7/2).**
  - The allocation rule is the textbook asymmetric-Nash TU solution.
  - The three-term decomposition is a second-order expansion. The "offsetting" −Cov(c, π) term is, by the chapter's own account, a reason not to use V.
  - Everything rests on Assumption 7.A (π = Nash weights), which Manea (2011) and Abreu–Manea (2012) contradict for strategic network bargaining.
  - It would become a paper only with a microfoundation delivering π as limit payoffs, in the spirit of Elliott & Golub (2019, *JPE*) or the Melo Sanchez (Caltech 2013) thesis. That is a different project for a different skill set.
  - **Low** as is.
- **C7 quick checks.**
  - Ch23: textbook Dobrushin/Hajnal (Seneta 2006). The product counterexample is a good teaching device only.
  - Ch24: a logistic conformity saddle-node, i.e. Granovetter, Brock–Durlauf and Kuran in a different vocabulary.
  - Ch18: π ∝ degree is the textbook reversible-walk fact. The minority ratio is a degree ratio, already in Avin et al. (2015) and Karimi et al. (2018, *Sci. Rep.*).
  - Ch11: textbook identification.
  - Ch19: power nesting is stochastic complementation (Meyer 1989).
  - Ch27: a corollary of Ch18.
  - None is a standalone paper. Ch23 can feed #1, and Ch24 can be a vignette in #3.

## Suggested sequence

1. Write #3 first (CEPAL Review), about 2 months. It forces the fidelity clean-up the book needs anyway.
2. Build the ICIO / product-space pipeline once.
3. Write #1 (ESR note), about 2-3 months, mostly theory plus descriptive statistics.
4. Then attempt #2 (SCED), about 4-6 months, the highest-upside and highest-risk paper.
5. Do not invest in C2–C4 as standalone theory papers unless a co-author brings a genuine new theorem (e.g. a microfoundation for π as bargaining power).
