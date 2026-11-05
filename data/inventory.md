# Inventory of illustrative / unsourced quantitative claims

Built 2026-11-05 from `grep -n 'illustrative\|Illustrative' chapters/*.tex` (222 hits), the tables
next to those hits, and a second grep for untagged numbers (Gini, poverty, remittances, GDP,
ECI, Bs amounts, `TODO-VERIFY`).

**Line numbers refer to commit `f8e594c`** (chapters were being edited concurrently;
re-grep the quoted text if lines have moved).

Classes:
- **A**: obtainable from public data fetchable now. "A ✔" = fetched/computed in `data/` (see `README.md`); "A ✗" = fetchable but not done in this pass (reason given).
- **B**: public but not fetchable here (registration, paywall, microdata portal, blocked host, or no machine-readable source).
- **C**: inherently illustrative (hypothetical example, toy matrix, model parameter, or a number derived only from other illustrative numbers).

A row can cover several adjacent hits in one example. If a row mixes classes, it is split.

| # | chapter:line | claim | value used in book | class | sourced value / note |
|---|---|---|---|---|---|
| 1 | ch01:649-650 | stylized dual-economy heatmap | Γ≈0.16, N_eff≈7.9 | C | toy matrix |
| 2 | ch02:343, 405 | 4-agent examples | N=4, Ω=20, κ=0.1, Φ=8 | C | toy |
| 3 | ch02:545 | El Alto vs Zona Sur 8×8 block matrix (5 and 3 agents) | ε coupling | C | toy |
| 4 | ch02:620 | one-parameter linkage dynamics | ε_t map | C | toy |
| 5 | ch03:93 | gas economy: "large share of export revenue in boom years and a small share of employment" | qualitative | A ✔ | fuel exports 54.5% of merchandise exports (2013), 32.1% (2019), 18.9% (2024) [WDI]; mining & quarrying incl. hydrocarbons 1.3% of employment (2019) [ILOSTAT] |
| 6 | ch03:260-275 | A_Bolivia 2×2 (0.15, 0.03; 0.02, 0.20), χ=0.125 | toy | C | a real 2-block aggregation of the INE table is possible (`processed/bolivia_Z/`) |
| 7 | ch04:340-357 | departure-measure ranges A1–A6 | 0.1–0.5, 5–20%, ... | C | placeholders |
| 8 | ch05:525 | Minister / Gas Executive / Vendor | 3 agents | C | toy |
| 9 | ch06:45 | 5-agent production network | W^P | C | toy |
| 10 | ch07:354 | two 4-agent economies | Ω=20, d̄=3 | C | toy |
| 11 | ch08:219 | π_IMF = π_WB = π_UST = 0.2 | 0.2 each | C | toy |
| 12 | ch08:286 | consumer priors | (0.8, 0.5, 0.2) | C | toy |
| 13 | ch09:518 | documentary raised aspirations 0.1 SD | 0.1 SD | C | exercise premise |
| 14 | ch10:172-173, 329, 446, 507 | Bolivia three blocks / poverty trap / three policies | parameters | C | toy |
| 15 | ch11:21 | blanket tag "all magnitudes illustrative" | – | C | – |
| 16 | ch11:165-166 | epidemiological coefficients β2=0.40, β3=0.25 | 0.625 | C | stated as not from the papers |
| 17 | ch11:169-170 | US female LFP cohort–year variance decomposition | 0.042, 0.013, 0.027 → Ψ̂=2.08 | C | a real one needs IPUMS microdata (B) |
| 18 | ch11:183-218, 235-256 | five illustrative variance ratios | 1.3, 0.8, 0.5, 1.6, 2.08 | C | – |
| 19 | ch11:225-233 | Bolivia self-employment 30–65% across departments; Ψ̂=1.6 | hypothetical | B | departmental rates need EH microdata (INE ANDA portal; not fetched). National level A ✔: self-employed 68.1% of employment (2019, ILO modelled, WDI SL.EMP.SELF.ZS) |
| 20 | ch11:309-317, 394-395, 464, 527 | (ρ̂, μ̂) = (0.62, 0.18), sensitivity | – | C | – |
| 21 | ch11:503 | adoption gradients 0.10 vs 0.30 | – | C | "suppose" |
| 22 | ch11:519 | US female LFP ~28% (1940), ~36% (1945), ~33% (1950) | approx. | B | Census/BLS historical tables; FRED LNS11300002 starts 1948; not fetched |
| 23 | ch12:73, 303, 362 | satisfaction examples | – | C | toy |
| 24 | ch14:133-134, 262 | class-alliance π = (0.35, 0.30, 0.20, 0.15), Δπ | – | C | stylization |
| 25 | ch14:332-372 | alliance-region thresholds; Bolivia 1985/2006 coordinates | – | C | – |
| 26 | ch15:228, 264-322, 382 | Ψ ranges (info 3–5, institutions 1–2, infrastructure 0.5–1); schematic Γ_t | – | C | – |
| 27 | ch15:339-342 | conventional ISI ΔΓ≈0.02 | 0.02 | C | – |
| 28 | ch16:53-58, 377 | global W^P π: US 0.18, China 0.15, Germany 0.08; top-10 ≈65%; Bolivia 0.0005 | as listed | A ✔ | 2019 trade-based W^P (imports-only, 230 economies, Growth Lab HS92): π_US 0.089, π_CHN 0.130, π_DEU 0.082, top-10 0.519, π_BOL 0.00053 |
| 29 | ch16:56 | core–periphery: column sum of US ≫ 1, Bolivia ≪ 1 | qualitative | A ✔ | 2019 column sums of W: USA 24.2, CHN 30.7, DEU 11.1, BOL 0.10 (`processed/global_wp_pi.csv`) |
| 30 | ch16:116 | schematic global heatmap | – | C | schematic |
| 31 | ch16:132-147 | country-level BL/FL (US 2.45/1.98; China 2.68/2.31; Germany 2.51/2.12; Japan 2.33/1.85) | as listed | A ✗ | needs the full ICIO aggregated to countries (zips are in `raw/oecd_icio/`); not computed in this pass |
| 32 | ch16:140-141 | Bolivia 1.28/1.12, Ethiopia 1.15/1.08 | as listed | B | not separate in OECD ICIO; Eora needs registration |
| 33 | ch16:270-280 | complexity–inequality scatter | anonymous points | C | schematic |
| 34 | ch16:305-320 | relatedness density ω Bolivia vs Korea (0.04/0.52, 0.03/0.41, 0.02/0.68, 0.08/0.45; avg 0.06/0.55) | as listed | A ✗ | needs Growth Lab country×product RCA (451 MB) + product-space proximities; not done |
| 35 | ch16:343 | Bolivia ECI −0.85 (2006), −0.80 (2018) | −0.85, −0.80 | A ✔ | Growth Lab HS92: −1.21 (2006), −1.46 (2018); SITC: −1.08, −1.50. ECI **fell** |
| 36 | ch16:358 | N_eff of global economic policy institutions 3–5 | 3–5 | C | – |
| 37 | ch17:39 | blanket tag for Bolivian numbers | – | C | – |
| 38 | ch17:57-90 | 12-sector Bolivia: output shares, BL, FL; λ2=0.873, λ3=0.817, Γ=0.127 (sales shares 0.075) | as listed | A ✔ | INE COU 2019 (base 2017), domestic, basic prices, book12: Γ=0.467 (λ2=0.533, λ3=0.313); sales shares Γ=0.694; linkages in `processed/bolivia_io_linkages.csv` |
| 39 | ch17:92 | Korea Γ≈0.45, Germany 0.52, USA 0.48 | as listed | A ✔ | OECD ICIO 2023, 2019, same 12-sector aggregation: KOR 0.511, DEU 0.544, USA 0.463 (Bolivia like-for-like 0.475) |
| 40 | ch17:107-129 | Bolivia and Korea W^P heatmaps | illustrative matrices | A ✔ | real 12×12 domestic Z in `processed/bolivia_Z/`, Korea in `processed/oecd_icio_domestic/KOR.csv.gz` |
| 41 | ch17:162-178 | between-block inequality: block shares 0.08/0.35/0.57, mean incomes 18,500/5,200/1,800 Bs; V_total 32M; σ_B=0.62; Korea 0.20, Brazil 0.45, US 0.25 | as listed | B | needs EH 2019 microdata (INE ANDA portal, not fetched) and foreign microdata |
| 42 | ch17:182-185 | Pinto modern/primitive income ratio 10.3 | derived | C | derived from row 41 |
| 43 | ch17:198-215 | sector employment shares (hydrocarbons 0.01, mining 0.04, public adm 0.08, manuf 0.10, constr 0.07, commerce 0.15, transport 0.06, finance 0.02, agric 0.27, other 0.08, personal 0.09, domestic 0.03) | as listed | A ✔ | ILOSTAT 2019 (ISIC sections): agric 27.4%, mining+hydrocarbons (B) 1.3%, manuf 11.2%, construction 7.9%, commerce 17.1%, transport 7.1%, public adm (O) 2.6%, finance 0.7%, domestic (T) 1.9% |
| 44 | ch17:198-215 | sectoral π_s (hydrocarbons 0.195, public adm 0.308, agriculture 0.037 …) | as listed | A ✔ | INE 2019 book12: π_hydro 0.039, π_public adm 0.005, π_manuf 0.159, π_agric 0.091, π_other services 0.270, π_finance 0.162 |
| 45 | ch17:224 | power Gini 0.52, N_eff 5.82 of 12; Korea ≈0.28, Germany ≈0.22 | as listed | A ✔ | Bolivia 2019: G(π)=0.445, N_eff=6.48 of 11 (domestic services has no intermediate flows); like-for-like (book12_isic): BOL 0.435, KOR 0.658, DEU 0.611 |
| 46 | ch17:228-235 | π_gas − π_agr = 0.158 | 0.158 | A ✔ | 2019: 0.039 − 0.091 = −0.052 |
| 47 | ch17:249-264 | media ecosystems π^I = (0.34, 0.37, 0.12, 0.17), N_eff 3.4; US 15–20, Scandinavia 8–12 | as listed | B | audience measurement (Kantar IBOPE Bolivia: paywall; Latinobarómetro: registration) |
| 48 | ch17:262 | population 11 million | 11 M | A ✔ | 11.67 M (2019), 12.41 M (2024) [WDI SP.POP.TOTL] |
| 49 | ch17:281 | Ψ for self-employment 1.56 | 1.56 | C | model-derived |
| 50 | ch17:328-343 | dashboard 2019: Γ 0.13 vs Korea 0.45; G(π) 0.52 vs 0.28 | as listed | A ✔ | Γ 0.475 vs 0.511; G(π) 0.435 vs 0.658 (book12_isic) |
| 51 | ch17:328-343 | dashboard: σ_B 0.62/0.20, N_eff(W^I) 3.4/≈12 | as listed | B | see rows 41, 47 |
| 52 | ch17:328-343 | dashboard: Ψ 1.56/≈0.6 | as listed | C | – |
| 53 | ch17:400-415 | trajectory Γ: 0.10 (1990), 0.11 (2000), 0.12 (2006), 0.13 (2010), 0.13 (2019); G(π): 0.58, 0.55, 0.54, 0.52, 0.52 | as listed | A ✔ | INE base 1990 book12: Γ 0.567 (1990), 0.494 (2000), 0.518 (2006), 0.478 (2011), 0.484 (2013); base 2017: 0.467 (2019). G(π): 0.420, 0.420, 0.424, 0.369, 0.445. Base break between 2014 and 2017 |
| 54 | ch17:400-415 | trajectory σ_B and N_eff | as listed | B | microdata / media data |
| 55 | ch17:417 | state take 50% (Law 3058, 2005); 82% (DS 28701, 2006) | 50%, 82% | B | legal texts, not a data series; not checked here |
| 56 | ch17:421 | extreme poverty ~38% → ~15%; Gini ~0.60 → ~0.42; growth ~5% (2006–2019) | as listed | A ✔ | INE extreme poverty 38.2% (2005), 37.7% (2006) → 15.2% (2018, old lines); new lines 13.2% (2019). INE Gini 0.600 (2005), 0.592 (2006) → 0.426 (2019). Mean GDP growth 2006–2019 4.52%; 2006–2014 5.06% [WDI] |
| 57 | ch17:438-440, 578 | summary repeats Γ 0.13, σ_B 0.62, G(π) 0.52, N_eff 3.4, Ψ 1.56 | – | A ✔/B/C | as rows 50–52 |
| 58 | ch17:498, 573 | ECI ≈ −0.8 (2019) | −0.8 | A ✔ | −1.30 (HS92, 2019); −1.31 (SITC) |
| 59 | ch18:59 | gender openness λ_g = 0.3 | 0.3 | C | – |
| 60 | ch18:129-140 | female share of employment: gas/mining 5%, manuf 38%, construction 4%, commerce 65%, agric 45%, domestic 95%, government 42% | as listed | A ✔ | ILOSTAT 2019: mining+hydrocarbons 6.5%, manuf 44.1%, construction 3.1%, commerce 61.4%, agric 46.0%, domestic 93.9%, public adm 34.3% |
| 61 | ch18:245 | correlated λ_r ranges | ranges | C | – |
| 62 | ch18:255-256 | Census 2012: 41% aged 15+ declare indigenous belonging | 41% | A ✗ | INE Censo 2012 tables not fetched in this pass (flagged TODO-VERIFY in text) |
| 63 | ch18:258-275 | λ_r^U 0.25, λ_r^P 0.35, λ_r^I 0.30; per-capita π 0.008/0.014 | as listed | C | stated as not estimated |
| 64 | ch18:302, 369, 395-406, 446, 450, 488, 561 | intersectional products, Ψ=2.3, ratios 0.84/0.58 | – | C | – |
| 65 | ch18:352-363 | Findex: account ownership women 42%, men 58% | 42/58 | A ✔ | Findex via WDI: 63.3% / 74.0% (2021); 53.5% / 60.3% (2024) |
| 66 | ch18:352-363 | Findex: saved for business 12/22%, borrowed for business 8/18% | as listed | A ✗ | in the full Findex database (not in WDI); not fetched |
| 67 | ch19:68-80 | three scales: global Γ 0.20; max π/min π 360 (US/Laos) | as listed | A ✔ | 2019 trade W^P: Γ_global 0.452; π_US/π_LAO = 193 |
| 68 | ch19:68-80 | national Γ 0.13; max/min π 88 (public adm / domestic services) | as listed | A ✔ | Γ 0.467; max/min 54 (other services 0.270 / public adm 0.005; domestic services has π undefined) |
| 69 | ch19:68-80 | σ_B global ≈0.55, national 0.62; individual-scale column | as listed | B / C | σ_B needs microdata (B); individual column hypothetical (C) |
| 70 | ch19:196-207 | Germany π 0.08 / Bolivia 0.0005 = 160; workforces 45 and 6 million | as listed | A ✔ | 2019: 0.0824 / 0.000526 = 157; labour force DEU 43.98 M, BOL 5.95 M [WDI] |
| 71 | ch19:203, 243, 310-320, 365 | compound per-capita power ratios (0.33, 141, 427×, 1/100) | derived | C | derived from illustrative inputs |
| 72 | ch19:338-342 | Korea Γ_N 0.45, Γ_I 0.35; Uruguay Γ_N ≈0.30, π_G 0.001 | as listed | A ✔ / B | Korea Γ_N 0.511 (2019, ICIO); Uruguay not in OECD ICIO 2023 (B); π_URY (trade W^P, 2019) in `processed/global_wp_pi.csv` |
| 73 | ch19:361, 455 | Bolivia Γ in [0.10, 0.20] at three scales | – | A ✔ | national 0.467, global 0.452 |
| 74 | ch20:122-147, 179 | 10×10 Bolivia SAM (billions of Bs) | toy | C | a real SAM is not published by INE; the SUT shows hydrocarbons gross output 4.7% of total output (2019) |
| 75 | ch20:194, 288-317, 331, 363, 392 | SAM spectral results and multipliers (gas 5.18, manuf 5.82, agric 6.20) | derived | C | from toy SAM |
| 76 | ch20:317, 396 | "gas buys few domestic inputs; income goes mainly to capital" | qualitative | A ✔ | INE 2019: hydrocarbons domestic intermediate inputs 0.225 per unit of output (manufacturing 0.558); compensation 6.0% of hydrocarbons VA; 70.1% of output exported |
| 77 | ch20:340 | G(π) of SAM account totals 0.18 | 0.18 | C | toy |
| 78 | ch21:43-62, 193-215 | 6-group × 5-channel audience shares; population shares p = (0.17, 0.18, 0.25, 0.10, 0.15, 0.15) | as listed | B | survey data (no public module); population shares by urban/rural × indigenous are A ✗ (Censo 2012, not fetched) |
| 79 | ch21:133 | smartphone / social media (earlier 72% / 65% removed; TODO-VERIFY) | – | A ✔ (partial) | internet users 47.5% (2019), 77.3% (2023), 79.7% (2024) [ITU via WDI]; smartphone/social-media shares (DataReportal) not fetched |
| 80 | ch21:156-215, 281-300, 314, 349, 409-411 | composite weights, Γ^I 0.74, N_eff 5.59/4.85, λ values | derived | C | from illustrative audience table |
| 81 | ch22:56-101, 289, 304 | 4-group transfer-share matrix; λ^U 0.20; Γ 0.30; N_eff 3.85 | as listed | C | EH records transfers received, not senders: no data source exists |
| 82 | ch22:79, 201 | private transfers ≈10 bn Bs, "mostly remittances" | 10 bn Bs | A ✔ | remittances received US$1.336 bn = 9.30 bn Bs (2019); 10.08 bn Bs (2023) [WDI; 6.96 Bs/US$] |
| 83 | ch22:119 | directed dictator game | hypothetical | C | – |
| 84 | ch22:163-178 | group-specific trust table (1–4 scale) | hypothetical | C | no such battery in Latinobarómetro (and that needs registration: B) |
| 85 | ch22:191-197 | IDH ≈8 bn Bs; Bono Juancito Pinto ≈0.5 bn Bs; Renta Dignidad ≈3.5 bn Bs | as listed | A ✗ | MEFP / UDAPE fiscal tables (udape.gob.bo responds 200); not fetched in this pass |
| 86 | ch22:199 | fiscal incidence 12% / 8% / 25% | as listed | B | needs a fiscal-incidence study (e.g. CEQ), not fetched |
| 87 | ch22:201 | government transfers ≈15 bn Bs | 15 bn Bs | A ✗ | MEFP; not fetched |
| 88 | ch22:218, 250-251, 364 | alignment example, population shares (0.35, 0.25, 0.15, 0.25), correlations 0.96–0.99 | – | C | – |
| 89 | ch23:148-180 | five epochs: Γ 0.11→0.10 (1985–90), 0.10→0.12 (1994–2005), 0.12–0.13 (2006–09), 0.13 (2010–19) | as listed | A ✔ | INE base 1990 book12: 1988 0.563, 1990 0.567, 1994 0.567, 2005 0.502, 2006 0.518, 2009 0.504, 2011 0.478, 2013 0.484 (2014 0.361, an outlier); base 2017: 2017 0.506, 2019 0.467 |
| 90 | ch23:149 | "Bolivia's I-O tables (1990, 2000, 2006, 2012, 2016)" | list of years | A ✔ | INE publishes annual tables 1988–2009, 2011–2014 (base 1990) and 2017–2024 (base 2017); there is no 2016 table and no 2010 table |
| 91 | ch23:150-176 | σ_B by epoch (0.65→0.68→0.63→0.61→0.62) | as listed | B | microdata |
| 92 | ch23:176-180 | ΔECI by epoch: +0.02 (1994–2005), −0.01 (2006–09), +0.05 (2010–19) | as listed | A ✔ | HS92: −0.24 (1995→2005), −0.14 (2006→2009), −0.02 (2010→2019) |
| 93 | ch23:159 | royalties 18% + IDH 32% = 50% | legal | B | legal text; not checked |
| 94 | ch23:186-232 | Korea path Γ 0.15 → 0.45 (1960–2000), σ_B 0.50 → 0.20 | as listed | B | pre-1995 Bank of Korea tables not fetched; ICIO gives Korea from 1995 (2021 edition) and 2011–2020 (2023 edition) |
| 95 | ch23:234 | Bolivia growth average 4% per year | 4% | A ✔ | 4.04% mean (1990–2019); 3.76% (1986–2019) [WDI] |
| 96 | ch23:255-270 | structural velocities: Korea +0.075, China +0.058, Chile +0.033, Bolivia +0.010, Argentina −0.005, Venezuela −0.030 per decade | as listed | A ✔ (partial) | Bolivia 1990–2013 (base 1990): −0.036/decade. within-edition endpoints per decade: Korea −0.055 (1995–2014, ICIO 2021) / +0.055 (2011–2020, ICIO 2023); Chile +0.091 / +0.107; Argentina +0.037 / −0.050; China −0.038 (2011–2020). Pre-1995 and Venezuela not available (B) |
| 97 | ch23:280 | ρ ≈ 0.85 for Bolivia | 0.85 | C | – |
| 98 | ch23:303-312 | π_gas 0.07 (2003) → 0.11 (2014) → 0.09 (2019) | as listed | A ✔ | INE book12 π_hydrocarbons: 0.047 (2003) → 0.117 (2014) (base 1990); 0.039 (2019, base 2017; not comparable across bases) |
| 99 | ch23:330-366 | SCI: Korea 0.86, China 0.61, Chile 0.31, Bolivia 0.075; Γ_ref = 0.50 ("Germany 0.52") | as listed | A ✔ (partial) | Germany 0.544 (2019). Bolivia's Γ (0.47–0.57) is already near Γ_ref, so the SCI as defined is ill-conditioned for Bolivia |
| 100 | ch24:83-100, 225-235 | N_eff^I by hegemonic period (1.5, 2.0, 1.8, 2.5, 3.0–3.4, 3.5; 2.0→2.8) | as listed | C | no historical audience data |
| 101 | ch25:45-50 | product space: ~775 SITC-4 products; core proximity 0.35 vs 0.08; ratio 4.4; PCI_D ≈1.8 | as listed | A ✗ | Growth Lab product-space files (`top_edges_hs92.csv`) plus PCI; not computed in this pass |
| 102 | ch25:62, 99-108, 306 | truncation τ: Bolivia 0.98 (2 of ~100), Korea 0.18 (82 of ~100); LA 0.70–0.95, East Asia 0.10–0.40 | as listed | A ✗ | needs country×product RCA (451 MB) |
| 103 | ch25:107 | ECI Bolivia −0.80, Korea +2.10 | as listed | A ✔ | 2019 HS92: −1.30 / +1.81 |
| 104 | ch25:108 | Γ(W^P) 0.13 / 0.45 | as listed | A ✔ | 0.475 / 0.511 (2019, like-for-like) |
| 105 | ch25:109 | σ_B 0.62 / 0.20 | as listed | B | microdata |
| 106 | ch25:110 | GDP per capita PPP 2019: \$8,800 / \$44,000 | as listed | A ✔ | \$11,005 / \$46,511 (current international \$, 2021 ICP vintage) [WDI NY.GDP.PCAP.PP.CD] |
| 107 | ch25:171, 195 | schematic product space; relaxation-time penalty | – | C | – |
| 108 | ch25:247-258 | Korea τ 0.70 → 0.20 and Γ 0.25 → 0.40 (1973–1993) | as listed | B | pre-1995 data |
| 109 | ch25:273-290 | Bolivia incremental path ω ≈ 0.08, 0.06, 0.04 | as listed | C / A ✗ | relatedness computable (row 34) |
| 110 | ch25:340 | PCI of lithium batteries ≈1.5 | 1.5 | A ✔ | HS92 8507 (electric accumulators incl. Li-ion): PCI 0.549 (2019), 0.811 (2023); lithium carbonate (2836) −0.57 (2019) |
| 111 | ch26:43, 47 | thresholds Γ* ≈ 0.30, N*_eff ≈ 8 | – | C | – |
| 112 | ch26:53-62 | pathology table (Brazil, Mexico, China, Cuba, Venezuela, Korea) | qualitative | C | – |
| 113 | ch26:73 | Korea Γ below 0.30 until ≈1980 | – | B | pre-1995 |
| 114 | ch26:76-113 | (Γ, ECI) scatter: BOL (0.13, −0.8), BRA (0.25, 0.8), COL (0.20, −0.3), MEX (0.28, 1.0), KOR (0.45, 2.1), JPN (0.56, 2.3), DEU (0.52, 1.8), CZE (0.40, 1.5), CHL (0.28, 0.3), ARG (0.15, 0.5), VEN (0.10, −0.5), MYS (0.42, 0.8), POL (0.38, 1.2), CHN (0.35, 1.0) | as listed | A ✔ | 2019, Γ like-for-like (book12_isic, domestic, input shares) and HS92 ECI: BOL (0.475, −1.30), BRA (0.496, 0.48), COL (0.460, 0.14), MEX (0.649, 0.95), KOR (0.511, 1.81), JPN (0.510, 2.05), DEU (0.544, 1.54), CZE (0.486, 1.48), CHL (0.557, −0.06), ARG (0.522, 0.06), MYS (0.381, 1.03), POL (0.552, 1.06), CHN (0.582, 1.30); VEN ECI −1.02, Γ not available (not in ICIO: B) |
| 115 | ch26:165-206, 222, 270 | 30-year plan: costs \$17B, phases, targets | – | C | plan |
| 116 | ch26:167 | current position: Γ 0.13, ECI −0.80, N_eff 3.4, σ_B 0.62 | as listed | A ✔ / B | Γ 0.467–0.475, ECI −1.30; N_eff and σ_B B |
| 117 | ch26:208 | Bolivia GDP ≈ \$40B | \$40B | A ✔ | \$49.1B (2019), \$54.9B (2024), current US\$ [WDI; reflects INE's base-2017 rebasing] |
| 118 | ch27:85-104 | λ ≈ 0.22; community radio 14%; N_eff^I 3.4 → 2.5–3.0 | as listed | C | derived from illustrative tables |
| 119 | ch27:116 | echo chamber vs viral hub | – | C | – |
| 120 | ch28:72 | multiple structural equilibria Γ = 0.072 and 0.643 | – | C | – |
| 121 | ch28:143-160 | portrait: Γ 0.13; G(π) 0.52; ECI −0.80; π_G 0.0005; v_Γ +0.01/decade; SCI 0.075 | as listed | A ✔ | Γ 0.467; G(π) 0.445; ECI −1.30; π_G 0.00053; v_Γ −0.036/decade (1990–2013) |
| 122 | ch28:143-160 | portrait: σ_B 0.62; N_eff(W^I) 3.4 | as listed | B | – |
| 123 | ch28:143-160 | portrait: N_eff(W^U) 3.85; Ψ̂ 1.56; τ 0.98; power ratio 1:467; truncation share 55–86% | as listed | C / A ✗ (τ) | – |

## Counts

Each row is counted once, in its first-listed class ("A ✔ / B" counts as A ✔). Counted by script from the table above:

| class | rows |
|---|---|
| A ✔ (fetched / computed) | 42 |
| A ✗ (fetchable, not done) | 8 |
| B | 18 |
| C | 55 |
| **total** | **123** |

Summed over classes: A = 50 (42 done + 8 not done), B = 18, C = 55.
