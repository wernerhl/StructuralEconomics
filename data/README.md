# data/: sourced figures for the data pass

This directory replaces the book's `\illustrative` Bolivia / Latin America / world numbers with
sourced values where public data could be fetched. Nothing in `chapters/` was edited.

- **Retrieval date for everything below: 2026-11-05.** Each download, with its URL, HTTP status, bytes
  and sha256 prefix, is logged in `raw/download_log.csv`.
- `inventory.md` lists every illustrative or unsourced claim (chapter:line at commit `f8e594c`),
  classified A/B/C. The status table at the end of this file follows the same row numbers.
- **Never-fabricate rule:** every value below is printed by a script in `fetch/` or `compute/`.
  Where something failed, it says so (§9).

## 0. How to reproduce

```
pip install pandas numpy openpyxl requests
python3 data/fetch/fetch_wdi.py              # World Bank WDI API
python3 data/fetch/fetch_pip.py              # World Bank PIP API
python3 data/fetch/fetch_ine_poverty.py      # INE poverty / Gini tables
python3 data/fetch/fetch_ine_sut.py          # INE supply-use tables, base 2017 + base 1990
python3 data/fetch/fetch_ilostat.py          # ILOSTAT SDMX
python3 data/fetch/fetch_atlas_eci.py        # Growth Lab ECI/PCI (Harvard Dataverse)
python3 data/fetch/fetch_atlas_bilateral.py  # Growth Lab bilateral trade (31 MB, git-ignored)
python3 data/fetch/fetch_oecd_icio.py                 # OECD ICIO 2023 (2011-2020; ~180 MB, git-ignored)
python3 data/fetch/fetch_oecd_icio.py --edition 2021  # OECD ICIO 2021 (1995-2014; ~460 MB, git-ignored)
python3 data/compute/compute_bolivia_io.py
python3 data/compute/compute_oecd_io.py
python3 data/compute/compute_global_wp.py
python3 data/compute/compute_eci_facts.py
python3 data/compute/compute_velocity.py
python3 data/compute/compute_bolivia_facts.py
```

Files over about 20 MB are listed in `.gitignore` (`raw/oecd_icio/`, `raw/atlas/hs92_country_country_year.csv`)
and are re-downloaded by their fetch scripts. Every file in `processed/` is under 14 MB (largest: `ine_sut_long.csv`, 13.2 MB).

## 1. Conventions used for Γ, π, N_eff and the linkages (read before quoting a number)

These follow `docs/conventions.md` and are implemented in `compute/iolib.py`.

- **Flows.** Published tables have rows = suppliers. With Z[s,s′] = sales of s to s′ and x = output:
  - data convention: A_data = Z / x (by column), L = (I − A_data)⁻¹;
  - book convention: A_book = A_dataᵀ (row i is i's recipe), so the book's L = (I − A_book′)⁻¹ is **the same matrix**.
- **Linkages** are Rasmussen raw sums of L: BL_j = column sum, FL_i = row sum (as `code/ch17.py`: `L.sum(0)`, `L.sum(1)`).
  - They are computed from A, never from W^P. mean(BL) = mean(FL).
  - Normalised indices (BL/mean, FL/mean) and Ghosh forward linkages are also in the CSVs.
- **W^P = D⁻¹A_book**: W^P[s,s′] = Z[s′,s] / Σ Z[·,s], the input-cost shares.
  - The sales-share alternative is B[s,s′] = Z[s,s′] / Σ Z[s,·].
  - Both are reported with the diagonal kept and with the diagonal zeroed and the rows renormalised.
- **Other measures.**
  - Γ = 1 − |λ₂|; π = left Perron vector; N_eff = 1/Σπ²; power Gini as in `code/ch06.py`.
  - A sector with no intermediate purchases or no intermediate sales has an undefined W^P row and is dropped. The CSV column `dropped_sectors` names it (in Bolivia: "Domestic services").
- **Γ depends on aggregation and normalisation.** Bolivia 2019 ranges from 0.33 to 0.69 across the eight variants (§2.2). Compare countries **only** like-for-like:
  - `book12_isic` aggregation;
  - domestic flows, basic prices;
  - input shares, diagonal kept.
  - That is the headline everywhere below.
  - `detail` Γ (39 INE industries vs 45 ICIO industries) is not comparable across sources.

## 2. Bolivia's input-output structure (INE supply-use tables)

**Source.** INE Bolivia, Cuentas Nacionales, current prices.
- **(a) Base 2017** "Cuadros de Oferta y Utilización" (2-digit, 42 products × 39 industries, millions of Bs), 2017–2024. 2023–2024 are marked (P) preliminary.
  - Page: https://www.ine.gob.bo/referencia2017/cuadros_oferta_utilizacion.html
  - Files: `https://www.ine.gob.bo/referencia2017/CUADROS/pagina_web/4.2.N.xlsx`. Cuadro 4.2.(2k−1) is supply and 4.2.(2k) is use for year 2016+k; 2019 = 4.2.5 / 4.2.6.
- **(b) Base 1990** "Matrices de Insumo Producto" (35 × 35, thousands of Bs), 1988–2009 and 2011–2014.
  - Page: https://www.ine.gob.bo/index.php/estadisticas-economicas/pib-y-cuentas-nacionales/matrices/matrices-de-insumo-producto/
  - Files: `https://www.ine.gob.bo/index.php/descarga/406/matrices-de-insumo-producto-corriente/<id>/matriz-de-insumo-producto-en-miles-de-bolivianos<year>.xlsx` (ids 44576–44601).

**Scripts.** `fetch/fetch_ine_sut.py` writes `processed/ine_sut_long.csv`. `compute/compute_bolivia_io.py` writes:
- `processed/bolivia_io_summary.csv`, `_linkages.csv`, `_pi.csv`;
- `processed/bolivia_io_sector_facts.csv`, `_industry_facts.csv`, `_checks.csv`;
- `processed/bolivia_Z/Z_<base>_<year>_<agg>.csv`.

**Construction.** The SUTs do not split domestic from imported use, so the I-O table is derived:
- **Valuation.** Purchasers' prices → basic prices. Each product's net taxes and trade margins are removed in proportion to every use, and the margins are reallocated to commerce (the only product with negative margins).
- **Imports.** Import proportionality across non-export uses. Exports are domestic up to domestic output; any excess is treated as re-exports.
- **Industry × industry.** Fixed product-sales-structure (Eurostat model D): Z = D·U_dom.
- **Base 1990 specifics.**
  - Make matrices are published only for 2002, 2003, 2006 and 2011–2014. Other years assume a diagonal make matrix (the 2006 make matrix is 98.6% diagonal).
  - The FISIM dummy industry is allocated to industries in proportion to intermediate consumption.
- **Checks.** Supply = use and output = intermediate + final sales hold to < 10⁻⁴ in every year except two inconsistencies in the published tables: 2011 product 2 (13%) and 2012 industry 21. See `bolivia_io_checks.csv`.
- **12-sector mapping** (`book12`, the 12 sectors of ch17:57). The mapping dicts in the script give the details.
  - Utilities, communication, real estate, business services, market education/health → Other services.
  - Hotels/restaurants and community/personal services → Personal services.
  - Refining → Manufacturing.
  - Government non-market education/health → Public administration in `book12`, and → Other services in `book12_isic` (to match ICIO).

### 2.1 Headline: Γ, π, N_eff, power Gini, 2019 (base 2017, domestic, input shares)

| measure | book (ch17:57-92, 198-224, 328-343) | INE 2019, `book12` | INE 2019, `book12_isic` | INE 2019, `detail39` |
|---|---|---|---|---|
| Γ = 1−\|λ₂\| | 0.127 (≈0.13) | **0.467** (λ₂ = 0.533, real) | 0.475 | 0.410 |
| \|λ₃\| | 0.817 | 0.313 | – | 0.440 |
| Γ, sales-share normalisation | 0.075 | 0.694 | 0.687 | 0.331 |
| N_eff (of N) | 5.82 of 12 | 6.48 of 11 | 6.67 of 11 | 11.86 of 38 |
| power Gini G(π) | 0.52 | 0.445 | 0.435 | 0.670 |
| π_max | Public adm 0.308 | Other services 0.270 | – | Financial services 0.184 |

### 2.2 Γ under every normalisation and diagonal treatment, 2019 (`processed/bolivia_io_summary.csv`)

| aggregation / flows | input shares | input, zero diag | sales shares | sales, zero diag |
|---|---|---|---|---|
| book12, domestic | 0.467 | 0.614 | 0.694 | 0.653 |
| book12, total (incl. imports) | 0.443 | 0.576 | 0.664 | 0.691 |
| detail39, domestic | 0.410 | 0.476 | 0.331 | 0.544 |
| detail39, total | 0.354 | 0.417 | 0.324 | 0.542 |

### 2.3 The 12-sector table of ch17:57-90 and ch17:198-215 with real 2019 values

Output shares, BL, FL and π are from INE 2019 (domestic, book12). Employment shares are from ILOSTAT 2019 (§6), by ISIC section; the mapping to the book's sectors is approximate.

| sector | output share (book) | output share | BL (book) | BL | FL (book) | FL | π_s (book) | π_s | employment share (book) | employment share (ILOSTAT ISIC) |
|---|---|---|---|---|---|---|---|---|---|---|
| Hydrocarbons | 0.08 | 0.047 | 1.34 | 1.340 | 1.54 | 1.195 | 0.195 | 0.039 | 0.01 | 0.013 (B, incl. mining) |
| Mining | 0.06 | 0.045 | 1.47 | 1.364 | 1.27 | 1.296 | 0.123 | 0.037 | 0.04 | (in B) |
| Public administration | 0.09 | 0.107 | 1.27 | 1.245 | 1.69 | 1.023 | 0.308 | 0.005 | 0.08 | 0.026 (O) |
| Manufacturing | 0.15 | 0.208 | 1.93 | 1.855 | 2.15 | 2.207 | 0.108 | 0.159 | 0.10 | 0.112 (C) |
| Construction | 0.07 | 0.093 | 1.85 | 1.574 | 1.38 | 1.128 | 0.036 | 0.029 | 0.07 | 0.079 (F) |
| Commerce | 0.10 | 0.078 | 1.51 | 1.416 | 1.76 | 1.576 | 0.070 | 0.095 | 0.15 | 0.171 (G) |
| Transport | 0.08 | 0.063 | 1.76 | 1.542 | 1.57 | 1.522 | 0.052 | 0.092 | 0.06 | 0.071 (H) |
| Financial services | 0.04 | 0.041 | 1.42 | 1.433 | 1.38 | 1.570 | 0.036 | 0.162 | 0.02 | 0.007 (K) |
| Agriculture | 0.12 | 0.084 | 1.50 | 1.549 | 1.52 | 1.704 | 0.037 | 0.091 | 0.27 | 0.274 (A) |
| Other services | 0.12 | 0.144 | 1.42 | 1.356 | 1.26 | 2.018 | 0.018 | 0.270 | 0.08 | – |
| Personal services | 0.06 | 0.087 | 1.21 | 1.661 | 1.18 | 1.097 | 0.012 | 0.021 | 0.09 | – |
| Domestic services | 0.03 | 0.004 | 1.07 | 1.000 | 1.05 | 1.000 | 0.003 | undefined (no intermediate flows) | 0.03 | 0.019 (T) |
| mean | | | 1.48 | 1.445 | 1.48 | 1.445 | | | | |

- **Hydrocarbons in the 2019 SUT** (`processed/bolivia_facts.csv`; used by ch17, ch20:317 and 396, and ch3:93):
  - 70.1% of output is exported;
  - domestic intermediate inputs are 0.225 per unit of output (manufacturing: 0.558);
  - compensation of employees is 6.0% of the sector's value added;
  - the sector is 6.14% of total gross value added and 4.7% of gross output.
- **Hydrocarbons over time** (VA share): 3.5% (2000), 7.9% (2006), 10.0% (2013) on base 1990; 6.3% (2017), 4.4% (2023), 3.3% (2024) on base 2017.

### 2.4 Trajectory (ch17:400-415, ch23:148-180, ch23:303-312); `book12`, domestic, input shares

| year | 1988 | 1990 | 1995 | 2000 | 2003 | 2005 | 2006 | 2009 | 2011 | 2013 | 2014 | ‖ | 2017 | 2019 | 2022 | 2024 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Γ (book: 0.10 in 1990 → 0.13 in 2019) | 0.563 | 0.567 | 0.578 | 0.494 | 0.514 | 0.502 | 0.518 | 0.504 | 0.478 | 0.484 | 0.361 | ‖ | 0.506 | 0.467 | 0.400 | 0.419 |
| G(π) (book: 0.58 → 0.52) | 0.442 | 0.420 | 0.411 | 0.420 | 0.444 | 0.366 | 0.424 | 0.353 | 0.369 | 0.371 | 0.352 | ‖ | 0.446 | 0.445 | 0.475 | 0.470 |
| π_hydrocarbons (book: 0.07 in 2003 → 0.11 in 2014 → 0.09 in 2019) | 0.025 | 0.037 | 0.025 | 0.036 | 0.047 | 0.074 | 0.064 | 0.077 | 0.088 | 0.085 | 0.117 | ‖ | 0.045 | 0.039 | 0.026 | 0.021 |

- ‖ marks the change of base (1990 → 2017): different classification and SNA vintage. Do not compute changes across it.
- In base-1990 diagonal-make years N = 10, because public administration has no intermediate sales. 2014 looks like an outlier.
- **Structural velocity** (`processed/gamma_velocity.csv`; replaces ch23:255-270 "+0.010/decade"):
  - base 1990, 1990–2013: −0.036 per decade (endpoints); −0.046 (OLS);
  - base 2017, 2017–2024: −0.12 per decade.

## 3. Comparators (OECD ICIO domestic blocks)

**Source.** OECD (2023), *OECD Inter-Country Input-Output Database*, http://oe.cd/icio. SML version; 76 economies + ROW; 45 ISIC Rev.4 industries; USD million at basic prices.
- Files used:
  - `https://stats.oecd.org/wbos/fileview2.aspx?IDFile=d1ab2315-298c-4e93-9a81-c6f2273139fe` (2016-2020_SML.zip);
  - `…IDFile=1f791bc6-befb-45c5-8b34-668d08a1702a` (2011-2015_SML.zip).
- The 2021 edition (1995–2014 used) came from `…IDFile=91d8e84b-7406-46b9-af5f-ec096242755c`, `8adf89dd-18b4-40fe-bc7f-c822052eb961`, `fe218690-0a3b-44aa-a82c-b3e3da6d24db` and `2c2f499f-5703-4034-9457-2f7518e8f2fc`.
- The IDs come from the iotr R package and pymrio 0.7.0, and each was verified by its Content-Disposition file name.

**Scripts.** `fetch/fetch_oecd_icio.py` writes `processed/oecd_icio{,2021}_domestic/<ISO3>.csv.gz` and `processed/oecd_icio{,2021}_output.csv.gz`. `compute/compute_oecd_io.py` writes `processed/oecd_io_summary.csv`, `oecd_io_linkages.csv`, `oecd_io_pi.csv` and `gamma_comparison.csv`. Bolivia is not in ICIO (it is in ROW).

### 3.1 Γ like-for-like, 2019 (book12_isic, domestic, input shares, diagonal kept), ICIO 2023

| country | Γ | N_eff (of 11–12) | G(π) | book value (ch17:92, ch19:341, ch25:108) |
|---|---|---|---|---|
| **Bolivia (INE)** | **0.475** | 6.67 | 0.435 | Γ 0.13, G(π) 0.52 |
| Korea | 0.511 | 3.62 | 0.658 | Γ 0.45, G(π) 0.28 |
| Germany | 0.544 | 3.89 | 0.611 | Γ 0.52, G(π) 0.22 |
| United States | 0.463 | 3.89 | 0.614 | Γ 0.48 |
| Chile | 0.557 | 4.22 | 0.568 | – |
| China | 0.582 | 3.44 | 0.622 | – |
| Mexico | 0.649 | 4.30 | 0.608 | – |
| Brazil | 0.496 | 4.35 | 0.597 | – |
| Argentina | 0.522 | 5.40 | 0.526 | – |
| Colombia | 0.460 | 4.78 | 0.537 | – |
| Peru | 0.407 | 5.15 | 0.558 | – |
| Japan | 0.510 | 3.77 | 0.664 | – |

- **The book's premise fails here.** It says Bolivia is a fragmented outlier (Γ 3.5× smaller than Korea's; power Gini twice as high). Like-for-like, Bolivia's Γ sits inside the comparator range, and its sectoral power concentration is lower than that of every comparator above.
- At 45-industry detail, Γ (2019, input shares) is: KOR 0.455, DEU 0.372, USA 0.465, MEX 0.560, ARG 0.373, BRA 0.436, COL 0.327, PER 0.363; CHL 0.483 at 43 and CHN 0.377 at 44 industries. These are not comparable with Bolivia's detail39.

### 3.2 Trajectories and velocities (ch23:255-270, ch23:333-366, ch26:73)

Γ, book12_isic, input shares (`processed/oecd_io_summary.csv`):

| | 1995 (ed.2021) | 2000 | 2005 | 2010 | 2014 | ‖ 2011 (ed.2023) | 2015 | 2019 | 2020 |
|---|---|---|---|---|---|---|---|---|---|
| Korea | 0.626 | 0.517 | 0.411 | 0.514 | 0.522 | 0.491 | 0.521 | 0.511 | 0.541 |
| Chile | 0.467 | 0.482 | 0.515 | 0.578 | 0.640 | 0.449 | 0.635 | 0.557 | 0.544 |
| Argentina | 0.407 | 0.423 | 0.505 | 0.477 | 0.476 | 0.538 | 0.484 | 0.522 | 0.492 |
| Germany | 0.554 | 0.516 | 0.490 | 0.534 | 0.549 | 0.502 | 0.545 | 0.544 | 0.541 |
| USA | 0.475 | 0.478 | 0.437 | 0.462 | 0.457 | 0.461 | 0.455 | 0.463 | 0.446 |
| China | n/a* | | | | | 0.613 | 0.607 | 0.582 | 0.579 |

\*China and Mexico are split (CN1/CN2, MX1/MX2) in the 2021 edition, and the split-block extraction returned empty: not available for 1995–2010 in this pass.

- **Velocities per decade** (endpoints, within edition):
  - Korea: −0.055 (1995–2014), +0.055 (2011–2020);
  - Chile: +0.091 (1995–2014), +0.107 (2011–2020);
  - Argentina: +0.037 (1995–2014), −0.050 (2011–2020);
  - China: −0.038 (2011–2020).
- Nothing before 1995 is available here, so the book's Korea 1960–2000 path (Γ 0.15 → 0.45) remains unsourced (B).
- The SCI with Γ_ref = 0.50 is undefined or ill-conditioned for most countries, because Γ_t0 ≥ 0.40 already.

## 4. Inequality and poverty

### 4.1 INE official series (Encuesta de Hogares)

**Source.** INE, "Pobreza y Desigualdad – Encuestas de Hogares", https://www.ine.gob.bo/index.php/estadisticas-economicas/encuestas-de-hogares/
- Cuadros 3.06.01.01 (poverty), 3.06.01.02 (extreme poverty), 3.06.01.03 (Gini of monthly per-capita household income).
- Files: `https://nube.ine.gob.bo/index.php/s/<id>/download` (ids in the script).
- The "-rev25" files use the 2025 population projections.

**Scripts.** `fetch/fetch_ine_poverty.py` writes `processed/ine_poverty_gini_long.csv`, with footnotes verbatim in `processed/ine_poverty_notes.csv`. Scalars are in `compute/compute_bolivia_facts.py` → `processed/bolivia_facts.csv`.

**Two line vintages; do not splice them:** old lines (2005–2018) and EPF 2015–16 lines (2016–2025). There was no survey in 2010.

| indicator | 2005 | 2006 | 2018 | 2019 | 2024 | 2025 | replaces |
|---|---|---|---|---|---|---|---|
| extreme poverty, old lines (%) | 38.2 | 37.7 | 15.2 | – | – | – | ch17:421 "≈38% → ≈15%" |
| extreme poverty, EPF lines, rev25 (%) | – | – | 15.8 | 13.2 | 13.3 | 16.7 | ch17:421 |
| moderate poverty, old lines (%) | 60.6 | 59.9 | 34.6 | – | – | – | – |
| moderate poverty, EPF lines, rev25 (%) | – | – | 40.5 | 37.7 | 38.1 | 44.7 | – |
| Gini, per-capita household income (rev25) | 0.600 | 0.592 | 0.435 | 0.426 | 0.418 | 0.426 | ch17:421 "≈0.60 → ≈0.42" |

The rev25 rural Gini was 0.659 in 2005 and 0.471 in 2019; the urban Gini was 0.539 and 0.386.

### 4.2 World Bank PIP (income-based, comparable spell 5 from 2005)

- **Source:** https://api.worldbank.org/pip/v1/pip?country=BOL&year=all&povline=3.0&format=json (also povline 4.2 and 8.3).
- **Script:** `fetch/fetch_pip.py` → `processed/pip_bolivia.csv`.
- **Gini:** 0.616 (2000), 0.585 (2005), 0.567 (2006), 0.416 (2019), 0.409 (2024).
- **Poverty at \$3.00 (2021 PPP):** 28.2% (2000), 15.9% (2006), 3.1% (2019), 2.1% (2024).
- The 1990 and 1992 surveys are urban-only.
- WDI `SI.POV.GINI` gives the same numbers (e.g. 41.6 in 2019).

## 5. National accounts, trade, remittances, digital (World Bank WDI API)

**Source.** `https://api.worldbank.org/v2/country/<ISO3>/indicator/<CODE>?format=json&per_page=20000&date=1960:2025`. The WDI database was last updated 2026-10-08.
- **Scripts:** `fetch/fetch_wdi.py` → `processed/wdi_long.csv`, `processed/wdi_bolivia_wide.csv`. Derived scalars: `compute/compute_bolivia_facts.py` → `processed/bolivia_facts.csv`.

| figure | value | code | replaces |
|---|---|---|---|
| GDP growth, mean 2006–2019 | 4.52%/yr (2006–2014: 5.06%) | NY.GDP.MKTP.KD.ZG | ch17:421 "close to 5%" |
| GDP growth, mean 1990–2019 | 4.04%/yr (1986–2019: 3.76%) | NY.GDP.MKTP.KD.ZG | ch23:234 "average 4%" |
| GDP, current US\$ | 49.1 bn (2019), 54.9 bn (2024) | NY.GDP.MKTP.CD | ch26:208 "≈\$40B" (WDI now reflects INE's base-2017 rebasing) |
| GDP per capita PPP, 2019 | Bolivia \$11,005; Korea \$46,511 (current intl \$) | NY.GDP.PCAP.PP.CD | ch25:110 "\$8,800 / \$44,000" |
| fuel exports, % of merchandise exports | 48.8 (2006), 54.5 (2013), 32.1 (2019), 18.9 (2024) | TX.VAL.FUEL.ZS.UN | ch3:93, ch16:340 |
| ores and metals exports, % | 22.0 (2006), 28.7 (2019), 43.6 (2024) | TX.VAL.MMTL.ZS.UN | ch16:340 |
| natural gas rents, % GDP | 3.56 (2006), 5.56 (2013), 2.07 (2019) | NY.GDP.NGAS.RT.ZS | ch17, ch23 |
| remittances received | US\$1.336 bn = 9.30 bn Bs (2019); 2.72% of GDP; peak 8.0% of GDP (2007) | BX.TRF.PWKR.CD.DT, BX.TRF.PWKR.DT.GD.ZS | ch22:201 "≈10 bn Bs" (Bs at 6.96/US\$) |
| population | 10.52 M (2012), 11.67 M (2019), 12.41 M (2024) | SP.POP.TOTL | ch17:262 "11 million" |
| labour force 2019 | Bolivia 5.95 M; Germany 43.98 M | SL.TLF.TOTL.IN | ch19:207 "6 and 45 million" |
| self-employed, % of employment | 68.1 (2019) | SL.EMP.SELF.ZS | ch11:225 (national level only) |
| account ownership, women / men (Findex) | 63.3 / 74.0 (2021); 53.5 / 60.3 (2024) | FX.OWN.TOTL.FE.ZS / .MA.ZS | ch18:357 "42 / 58" |
| internet users, % population | 47.5 (2019), 77.3 (2023), 79.7 (2024) | IT.NET.USER.ZS | ch21:133 |

Comparators are in `wdi_long.csv` for KOR, CHL, CHN, MEX, DEU, USA, BRA, ARG, COL, PER, URY, VEN, JPN, LAC and the World.

## 6. Employment by sector and sex (ILOSTAT)

**Source.** ILOSTAT SDMX, `https://sdmx.ilo.org/rest/data/ILO,DF_EMP_TEMP_SEX_ECO_NB,1.0/BOL.A...?format=csvfile`, dataflow DF_EMP_TEMP_SEX_ECO_NB. Values are thousands of persons from the household surveys (Encuesta Continua de Empleo / EH); status in employment is DF_EMP_TEMP_SEX_STE_NB.
- **Script:** `fetch/fetch_ilostat.py` → `processed/ilostat_bolivia_emp_eco_sex.csv`, `processed/ilostat_bolivia_emp_status.csv`.

Female share of employment, 2019 (replaces ch18:129-140):

| sector | book | ILOSTAT 2019 (ISIC section) |
|---|---|---|
| Gas / Mining | 5% | 6.5% (B, mining & quarrying incl. hydrocarbons) |
| Manufacturing | 38% | 44.1% (C) |
| Construction | 4% | 3.1% (F) |
| Commerce | 65% | 61.4% (G) |
| Agriculture | 45% | 46.0% (A) |
| Domestic service | 95% | 93.9% (T) |
| Government | 42% | 34.3% (O, public administration; education and health excluded) |

Employment shares, 2019: see §2.3 (agriculture 27.4%, commerce 17.1%, manufacturing 11.2%, mining & quarrying 1.3%).

## 7. Economic complexity (Growth Lab, Harvard Dataverse)

**Source.** Growth Lab, "Growth Projections and Complexity Rankings", doi:10.7910/DVN/XTAQMC, v12 (2026-02-08).
- File `growth_proj_eci_rankings.csv`: https://dataverse.harvard.edu/api/access/datafile/13439575
- "International Trade Data (HS, 92)", doi:10.7910/DVN/T4CHWJ, v18 (2026-04-22): `hs92_product_year_4.csv` (datafile 13685116) for PCI and `hs92_country_year.csv` (13685109).
- SITC: doi:10.7910/DVN/H8SFD2 v18, `sitc_country_year.csv` (13685143).

**Scripts.** `fetch/fetch_atlas_eci.py` and `compute/compute_eci_facts.py` write `processed/eci_rankings_selected.csv`, `eci_selected.csv`, `eci_facts.csv`, `eci_epochs_bolivia.csv` and `pci_selected_products.csv`.
- ECI is z-scored across countries each year and is revised with every release: quote the release.

ECI (HS92) and rank:

| | 1995 | 2000 | 2006 | 2010 | 2014 | 2018 | 2019 | 2024 |
|---|---|---|---|---|---|---|---|---|
| Bolivia | −0.85 (111) | −0.78 (111) | −1.21 (121) | −1.28 (127) | −1.34 (130) | −1.46 (131) | −1.30 (131) | −1.10 (123) |
| Korea | 1.10 | 1.29 | 1.51 | 1.62 | 1.70 | 1.80 | 1.81 (2) | 1.60 |
| Chile | 0.10 | 0.19 | 0.07 | −0.05 | 0.01 | −0.08 | −0.06 | −0.04 |
| China | 0.63 | 0.79 | 0.98 | 1.15 | 1.26 | 1.31 | 1.30 | 1.27 |
| Mexico | 0.80 | 0.90 | 0.93 | 0.97 | 0.96 | 0.92 | 0.95 | 0.94 |

- **Replaces:**
  - ch16:343 (−0.85 in 2006, −0.80 in 2018) → −1.21 and −1.46;
  - ch17:498, ch17:573, ch19:313, ch25:107, ch26:167, ch28:155 (−0.80) → −1.30 (2019);
  - ch25:107 Korea +2.10 → +1.81.
- **Epoch changes for ch23:176-180** (book: +0.02, −0.01, +0.05): −0.24 (1995→2005), −0.14 (2006→2009), −0.02 (2010→2019).
- **PCI** (replaces ch25:340 "lithium batteries ≈1.5"): HS 8507 electric accumulators (incl. Li-ion) 0.549 (2019) and 0.811 (2023); lithium carbonate (2836) −0.573; natural gas (2711) −2.226.

## 8. Global production weight matrix from trade (Chapter 16, ch19:68-80)

**Source.** Growth Lab bilateral goods trade, HS92, `hs92_country_country_year.csv`, doi:10.7910/DVN/T4CHWJ v18, datafile 13685107 (31 MB, git-ignored).

**Scripts.** `fetch/fetch_atlas_bilateral.py` → `processed/atlas_bilateral_selected_years.csv.gz`; `compute/compute_global_wp.py` → `processed/global_wp_summary.csv`, `processed/global_wp_pi.csv`.

**Construction.**
- W[i,j] = share of importer i's goods imports coming from j (imports only, zero diagonal, 230 economies).
- Countries with zero imports or exports inside the set are dropped iteratively.

| figure | book | 2019 | replaces |
|---|---|---|---|
| π_USA | 0.18 | 0.089 | ch16:58, ch16:377 |
| π_CHN | 0.15 | 0.130 | ch16:58 |
| π_DEU | 0.08 | 0.082 | ch16:58, ch19:207 |
| π_BOL | 0.0005 | 0.00053 | ch16:58, ch19:313, ch28:157 |
| top-10 share | ≈65% | 51.9% (CHN, USA, DEU, JPN, FRA, KOR, NLD, ITA, GBR, RUS) | ch16:58 |
| π_DEU / π_BOL | 160 | 157 | ch19:207 |
| π_USA / π_LAO | 360 | 193 | ch19:77 |
| Γ_global | 0.20 | 0.452 | ch19:76 |
| column sum of W, USA / BOL | "well above 1" / "well below 1" | 24.2 / 0.10 | ch16:56 |

π_BOL by year: 0.00018 (1995), 0.00047 (2006), 0.00071 (2014), 0.00053 (2019), 0.00055 (2023). π tracks world export shares closely (gravity), as ch16:59 anticipates.

## 9. What failed, what was not fetched, and why

- **www.oecd.org** (IOTS landing page): HTTP 403 Cloudflare browser challenge from this machine, and WebFetch also gets 403.
  - The OECD harmonised national IOTs ("IOTS") were therefore not fetched.
  - ICIO domestic blocks are used instead (§3).
- **OECD ICIO 2023 year blocks 1995–2000, 2001–2005 and 2006–2010:** HTTP 404 on stats.oecd.org (the files were withdrawn).
  - Replaced by the 2021 edition for 1995–2014, flagged as a separate vintage.
  - China and Mexico could not be extracted from the 2021 edition (split CN1/CN2, MX1/MX2 blocks).
- **ICIO coverage gaps:** Bolivia, Venezuela, Uruguay and Ethiopia are not separate economies in ICIO. Eora needs registration (B).
- **Pre-1995 I-O tables** (Korea 1960–1995, Bolivia before 1988): not fetched (B).
- **INE Encuesta de Hogares microdata** (σ_B, departmental self-employment, block mean incomes): on INE's ANDA portal; not fetched (B).
- **Media audience shares and survey-based trust:**
  - Kantar IBOPE ratings are paywalled; Latinobarómetro requires registration (B);
  - a group-specific trust battery and dyadic transfer data do not exist (C).
- **Not done in this pass, but fetchable (A ✗):**
  - Growth Lab country×product RCA (451 MB) for relatedness density, truncation τ and product-space core proximities (ch16:305, ch25:45-110);
  - full ICIO aggregated to countries, for the country BL/FL of ch16:132;
  - INE Censo 2012 ethnicity tables (ch18:255, ch21:62);
  - full Findex database ("saved/borrowed for business", ch18:357);
  - MEFP/UDAPE fiscal tables (IDH, Bono Juancito Pinto, Renta Dignidad; ch22:194-201). udape.gob.bo responds, but the tables were not located in time.
- **INE poverty-lines table** (`lineas_pobreza_rev25.xlsx`): downloaded but not parsed (different layout); it is not needed for any claim.
- **api.imf.org** (502) was not needed.
- **Approximate values from memory are deliberately not given anywhere in this directory.**

## 10. Status summary (same rows as `inventory.md`)

Status values: **sourced** (A ✔: a value from a script in this directory is given), **fetchable, not done** (A ✗), **not fetchable here** (B), **remains illustrative** (C).

| # | chapter:line | claim | value used in book | class | status | sourced value / note |
|---|---|---|---|---|---|---|
| 1 | ch01:649-650 | stylized dual-economy heatmap | Γ≈0.16, N_eff≈7.9 | C | remains illustrative | toy matrix |
| 2 | ch02:343, 405 | 4-agent examples | N=4, Ω=20, κ=0.1, Φ=8 | C | remains illustrative | toy |
| 3 | ch02:545 | El Alto vs Zona Sur 8×8 block matrix (5 and 3 agents) | ε coupling | C | remains illustrative | toy |
| 4 | ch02:620 | one-parameter linkage dynamics | ε_t map | C | remains illustrative | toy |
| 5 | ch03:93 | gas economy: "large share of export revenue in boom years and a small share of employment" | qualitative | A ✔ | sourced | fuel exports 54.5% of merchandise exports (2013), 32.1% (2019), 18.9% (2024) [WDI]; mining & quarrying incl. hydrocarbons 1.3% of employment (2019) [ILOSTAT] |
| 6 | ch03:260-275 | A_Bolivia 2×2 (0.15, 0.03; 0.02, 0.20), χ=0.125 | toy | C | remains illustrative | a real 2-block aggregation of the INE table is possible (`processed/bolivia_Z/`) |
| 7 | ch04:340-357 | departure-measure ranges A1–A6 | 0.1–0.5, 5–20%, ... | C | remains illustrative | placeholders |
| 8 | ch05:525 | Minister / Gas Executive / Vendor | 3 agents | C | remains illustrative | toy |
| 9 | ch06:45 | 5-agent production network | W^P | C | remains illustrative | toy |
| 10 | ch07:354 | two 4-agent economies | Ω=20, d̄=3 | C | remains illustrative | toy |
| 11 | ch08:219 | π_IMF = π_WB = π_UST = 0.2 | 0.2 each | C | remains illustrative | toy |
| 12 | ch08:286 | consumer priors | (0.8, 0.5, 0.2) | C | remains illustrative | toy |
| 13 | ch09:518 | documentary raised aspirations 0.1 SD | 0.1 SD | C | remains illustrative | exercise premise |
| 14 | ch10:172-173, 329, 446, 507 | Bolivia three blocks / poverty trap / three policies | parameters | C | remains illustrative | toy |
| 15 | ch11:21 | blanket tag "all magnitudes illustrative" | – | C | remains illustrative | – |
| 16 | ch11:165-166 | epidemiological coefficients β2=0.40, β3=0.25 | 0.625 | C | remains illustrative | stated as not from the papers |
| 17 | ch11:169-170 | US female LFP cohort–year variance decomposition | 0.042, 0.013, 0.027 → Ψ̂=2.08 | C | remains illustrative | a real one needs IPUMS microdata (B) |
| 18 | ch11:183-218, 235-256 | five illustrative variance ratios | 1.3, 0.8, 0.5, 1.6, 2.08 | C | remains illustrative | – |
| 19 | ch11:225-233 | Bolivia self-employment 30–65% across departments; Ψ̂=1.6 | hypothetical | B | not fetchable here | departmental rates need EH microdata (INE ANDA portal; not fetched). National level A ✔: self-employed 68.1% of employment (2019, ILO modelled, WDI SL.EMP.SELF.ZS) |
| 20 | ch11:309-317, 394-395, 464, 527 | (ρ̂, μ̂) = (0.62, 0.18), sensitivity | – | C | remains illustrative | – |
| 21 | ch11:503 | adoption gradients 0.10 vs 0.30 | – | C | remains illustrative | "suppose" |
| 22 | ch11:519 | US female LFP ~28% (1940), ~36% (1945), ~33% (1950) | approx. | B | not fetchable here | Census/BLS historical tables; FRED LNS11300002 starts 1948; not fetched |
| 23 | ch12:73, 303, 362 | satisfaction examples | – | C | remains illustrative | toy |
| 24 | ch14:133-134, 262 | class-alliance π = (0.35, 0.30, 0.20, 0.15), Δπ | – | C | remains illustrative | stylization |
| 25 | ch14:332-372 | alliance-region thresholds; Bolivia 1985/2006 coordinates | – | C | remains illustrative | – |
| 26 | ch15:228, 264-322, 382 | Ψ ranges (info 3–5, institutions 1–2, infrastructure 0.5–1); schematic Γ_t | – | C | remains illustrative | – |
| 27 | ch15:339-342 | conventional ISI ΔΓ≈0.02 | 0.02 | C | remains illustrative | – |
| 28 | ch16:53-58, 377 | global W^P π: US 0.18, China 0.15, Germany 0.08; top-10 ≈65%; Bolivia 0.0005 | as listed | A ✔ | sourced | 2019 trade-based W^P (imports-only, 230 economies, Growth Lab HS92): π_US 0.089, π_CHN 0.130, π_DEU 0.082, top-10 0.519, π_BOL 0.00053 |
| 29 | ch16:56 | core–periphery: column sum of US ≫ 1, Bolivia ≪ 1 | qualitative | A ✔ | sourced | 2019 column sums of W: USA 24.2, CHN 30.7, DEU 11.1, BOL 0.10 (`processed/global_wp_pi.csv`) |
| 30 | ch16:116 | schematic global heatmap | – | C | remains illustrative | schematic |
| 31 | ch16:132-147 | country-level BL/FL (US 2.45/1.98; China 2.68/2.31; Germany 2.51/2.12; Japan 2.33/1.85) | as listed | A ✗ | fetchable, not done | needs the full ICIO aggregated to countries (zips are in `raw/oecd_icio/`); not computed in this pass |
| 32 | ch16:140-141 | Bolivia 1.28/1.12, Ethiopia 1.15/1.08 | as listed | B | not fetchable here | not separate in OECD ICIO; Eora needs registration |
| 33 | ch16:270-280 | complexity–inequality scatter | anonymous points | C | remains illustrative | schematic |
| 34 | ch16:305-320 | relatedness density ω Bolivia vs Korea (0.04/0.52, 0.03/0.41, 0.02/0.68, 0.08/0.45; avg 0.06/0.55) | as listed | A ✗ | fetchable, not done | needs Growth Lab country×product RCA (451 MB) + product-space proximities; not done |
| 35 | ch16:343 | Bolivia ECI −0.85 (2006), −0.80 (2018) | −0.85, −0.80 | A ✔ | sourced | Growth Lab HS92: −1.21 (2006), −1.46 (2018); SITC: −1.08, −1.50. ECI **fell** |
| 36 | ch16:358 | N_eff of global economic policy institutions 3–5 | 3–5 | C | remains illustrative | – |
| 37 | ch17:39 | blanket tag for Bolivian numbers | – | C | remains illustrative | – |
| 38 | ch17:57-90 | 12-sector Bolivia: output shares, BL, FL; λ2=0.873, λ3=0.817, Γ=0.127 (sales shares 0.075) | as listed | A ✔ | sourced | INE COU 2019 (base 2017), domestic, basic prices, book12: Γ=0.467 (λ2=0.533, λ3=0.313); sales shares Γ=0.694; linkages in `processed/bolivia_io_linkages.csv` |
| 39 | ch17:92 | Korea Γ≈0.45, Germany 0.52, USA 0.48 | as listed | A ✔ | sourced | OECD ICIO 2023, 2019, same 12-sector aggregation: KOR 0.511, DEU 0.544, USA 0.463 (Bolivia like-for-like 0.475) |
| 40 | ch17:107-129 | Bolivia and Korea W^P heatmaps | illustrative matrices | A ✔ | sourced | real 12×12 domestic Z in `processed/bolivia_Z/`, Korea in `processed/oecd_icio_domestic/KOR.csv.gz` |
| 41 | ch17:162-178 | between-block inequality: block shares 0.08/0.35/0.57, mean incomes 18,500/5,200/1,800 Bs; V_total 32M; σ_B=0.62; Korea 0.20, Brazil 0.45, US 0.25 | as listed | B | not fetchable here | needs EH 2019 microdata (INE ANDA portal, not fetched) and foreign microdata |
| 42 | ch17:182-185 | Pinto modern/primitive income ratio 10.3 | derived | C | remains illustrative | derived from row 41 |
| 43 | ch17:198-215 | sector employment shares (hydrocarbons 0.01, mining 0.04, public adm 0.08, manuf 0.10, constr 0.07, commerce 0.15, transport 0.06, finance 0.02, agric 0.27, other 0.08, personal 0.09, domestic 0.03) | as listed | A ✔ | sourced | ILOSTAT 2019 (ISIC sections): agric 27.4%, mining+hydrocarbons (B) 1.3%, manuf 11.2%, construction 7.9%, commerce 17.1%, transport 7.1%, public adm (O) 2.6%, finance 0.7%, domestic (T) 1.9% |
| 44 | ch17:198-215 | sectoral π_s (hydrocarbons 0.195, public adm 0.308, agriculture 0.037 …) | as listed | A ✔ | sourced | INE 2019 book12: π_hydro 0.039, π_public adm 0.005, π_manuf 0.159, π_agric 0.091, π_other services 0.270, π_finance 0.162 |
| 45 | ch17:224 | power Gini 0.52, N_eff 5.82 of 12; Korea ≈0.28, Germany ≈0.22 | as listed | A ✔ | sourced | Bolivia 2019: G(π)=0.445, N_eff=6.48 of 11 (domestic services has no intermediate flows); like-for-like (book12_isic): BOL 0.435, KOR 0.658, DEU 0.611 |
| 46 | ch17:228-235 | π_gas − π_agr = 0.158 | 0.158 | A ✔ | sourced | 2019: 0.039 − 0.091 = −0.052 |
| 47 | ch17:249-264 | media ecosystems π^I = (0.34, 0.37, 0.12, 0.17), N_eff 3.4; US 15–20, Scandinavia 8–12 | as listed | B | not fetchable here | audience measurement (Kantar IBOPE Bolivia: paywall; Latinobarómetro: registration) |
| 48 | ch17:262 | population 11 million | 11 M | A ✔ | sourced | 11.67 M (2019), 12.41 M (2024) [WDI SP.POP.TOTL] |
| 49 | ch17:281 | Ψ for self-employment 1.56 | 1.56 | C | remains illustrative | model-derived |
| 50 | ch17:328-343 | dashboard 2019: Γ 0.13 vs Korea 0.45; G(π) 0.52 vs 0.28 | as listed | A ✔ | sourced | Γ 0.475 vs 0.511; G(π) 0.435 vs 0.658 (book12_isic) |
| 51 | ch17:328-343 | dashboard: σ_B 0.62/0.20, N_eff(W^I) 3.4/≈12 | as listed | B | not fetchable here | see rows 41, 47 |
| 52 | ch17:328-343 | dashboard: Ψ 1.56/≈0.6 | as listed | C | remains illustrative | – |
| 53 | ch17:400-415 | trajectory Γ: 0.10 (1990), 0.11 (2000), 0.12 (2006), 0.13 (2010), 0.13 (2019); G(π): 0.58, 0.55, 0.54, 0.52, 0.52 | as listed | A ✔ | sourced | INE base 1990 book12: Γ 0.567 (1990), 0.494 (2000), 0.518 (2006), 0.478 (2011), 0.484 (2013); base 2017: 0.467 (2019). G(π): 0.420, 0.420, 0.424, 0.369, 0.445. Base break between 2014 and 2017 |
| 54 | ch17:400-415 | trajectory σ_B and N_eff | as listed | B | not fetchable here | microdata / media data |
| 55 | ch17:417 | state take 50% (Law 3058, 2005); 82% (DS 28701, 2006) | 50%, 82% | B | not fetchable here | legal texts, not a data series; not checked here |
| 56 | ch17:421 | extreme poverty ~38% → ~15%; Gini ~0.60 → ~0.42; growth ~5% (2006–2019) | as listed | A ✔ | sourced | INE extreme poverty 38.2% (2005), 37.7% (2006) → 15.2% (2018, old lines); new lines 13.2% (2019). INE Gini 0.600 (2005), 0.592 (2006) → 0.426 (2019). Mean GDP growth 2006–2019 4.52%; 2006–2014 5.06% [WDI] |
| 57 | ch17:438-440, 578 | summary repeats Γ 0.13, σ_B 0.62, G(π) 0.52, N_eff 3.4, Ψ 1.56 | – | A ✔/B/C | sourced | as rows 50–52 |
| 58 | ch17:498, 573 | ECI ≈ −0.8 (2019) | −0.8 | A ✔ | sourced | −1.30 (HS92, 2019); −1.31 (SITC) |
| 59 | ch18:59 | gender openness λ_g = 0.3 | 0.3 | C | remains illustrative | – |
| 60 | ch18:129-140 | female share of employment: gas/mining 5%, manuf 38%, construction 4%, commerce 65%, agric 45%, domestic 95%, government 42% | as listed | A ✔ | sourced | ILOSTAT 2019: mining+hydrocarbons 6.5%, manuf 44.1%, construction 3.1%, commerce 61.4%, agric 46.0%, domestic 93.9%, public adm 34.3% |
| 61 | ch18:245 | correlated λ_r ranges | ranges | C | remains illustrative | – |
| 62 | ch18:255-256 | Census 2012: 41% aged 15+ declare indigenous belonging | 41% | A ✗ | fetchable, not done | INE Censo 2012 tables not fetched in this pass (flagged TODO-VERIFY in text) |
| 63 | ch18:258-275 | λ_r^U 0.25, λ_r^P 0.35, λ_r^I 0.30; per-capita π 0.008/0.014 | as listed | C | remains illustrative | stated as not estimated |
| 64 | ch18:302, 369, 395-406, 446, 450, 488, 561 | intersectional products, Ψ=2.3, ratios 0.84/0.58 | – | C | remains illustrative | – |
| 65 | ch18:352-363 | Findex: account ownership women 42%, men 58% | 42/58 | A ✔ | sourced | Findex via WDI: 63.3% / 74.0% (2021); 53.5% / 60.3% (2024) |
| 66 | ch18:352-363 | Findex: saved for business 12/22%, borrowed for business 8/18% | as listed | A ✗ | fetchable, not done | in the full Findex database (not in WDI); not fetched |
| 67 | ch19:68-80 | three scales: global Γ 0.20; max π/min π 360 (US/Laos) | as listed | A ✔ | sourced | 2019 trade W^P: Γ_global 0.452; π_US/π_LAO = 193 |
| 68 | ch19:68-80 | national Γ 0.13; max/min π 88 (public adm / domestic services) | as listed | A ✔ | sourced | Γ 0.467; max/min 54 (other services 0.270 / public adm 0.005; domestic services has π undefined) |
| 69 | ch19:68-80 | σ_B global ≈0.55, national 0.62; individual-scale column | as listed | B / C | not fetchable here | σ_B needs microdata (B); individual column hypothetical (C) |
| 70 | ch19:196-207 | Germany π 0.08 / Bolivia 0.0005 = 160; workforces 45 and 6 million | as listed | A ✔ | sourced | 2019: 0.0824 / 0.000526 = 157; labour force DEU 43.98 M, BOL 5.95 M [WDI] |
| 71 | ch19:203, 243, 310-320, 365 | compound per-capita power ratios (0.33, 141, 427×, 1/100) | derived | C | remains illustrative | derived from illustrative inputs |
| 72 | ch19:338-342 | Korea Γ_N 0.45, Γ_I 0.35; Uruguay Γ_N ≈0.30, π_G 0.001 | as listed | A ✔ / B | sourced | Korea Γ_N 0.511 (2019, ICIO); Uruguay not in OECD ICIO 2023 (B); π_URY (trade W^P, 2019) in `processed/global_wp_pi.csv` |
| 73 | ch19:361, 455 | Bolivia Γ in [0.10, 0.20] at three scales | – | A ✔ | sourced | national 0.467, global 0.452 |
| 74 | ch20:122-147, 179 | 10×10 Bolivia SAM (billions of Bs) | toy | C | remains illustrative | a real SAM is not published by INE; the SUT shows hydrocarbons gross output 4.7% of total output (2019) |
| 75 | ch20:194, 288-317, 331, 363, 392 | SAM spectral results and multipliers (gas 5.18, manuf 5.82, agric 6.20) | derived | C | remains illustrative | from toy SAM |
| 76 | ch20:317, 396 | "gas buys few domestic inputs; income goes mainly to capital" | qualitative | A ✔ | sourced | INE 2019: hydrocarbons domestic intermediate inputs 0.225 per unit of output (manufacturing 0.558); compensation 6.0% of hydrocarbons VA; 70.1% of output exported |
| 77 | ch20:340 | G(π) of SAM account totals 0.18 | 0.18 | C | remains illustrative | toy |
| 78 | ch21:43-62, 193-215 | 6-group × 5-channel audience shares; population shares p = (0.17, 0.18, 0.25, 0.10, 0.15, 0.15) | as listed | B | not fetchable here | survey data (no public module); population shares by urban/rural × indigenous are A ✗ (Censo 2012, not fetched) |
| 79 | ch21:133 | smartphone / social media (earlier 72% / 65% removed; TODO-VERIFY) | – | A ✔ (partial) | sourced | internet users 47.5% (2019), 77.3% (2023), 79.7% (2024) [ITU via WDI]; smartphone/social-media shares (DataReportal) not fetched |
| 80 | ch21:156-215, 281-300, 314, 349, 409-411 | composite weights, Γ^I 0.74, N_eff 5.59/4.85, λ values | derived | C | remains illustrative | from illustrative audience table |
| 81 | ch22:56-101, 289, 304 | 4-group transfer-share matrix; λ^U 0.20; Γ 0.30; N_eff 3.85 | as listed | C | remains illustrative | EH records transfers received, not senders: no data source exists |
| 82 | ch22:79, 201 | private transfers ≈10 bn Bs, "mostly remittances" | 10 bn Bs | A ✔ | sourced | remittances received US$1.336 bn = 9.30 bn Bs (2019); 10.08 bn Bs (2023) [WDI; 6.96 Bs/US$] |
| 83 | ch22:119 | directed dictator game | hypothetical | C | remains illustrative | – |
| 84 | ch22:163-178 | group-specific trust table (1–4 scale) | hypothetical | C | remains illustrative | no such battery in Latinobarómetro (and that needs registration: B) |
| 85 | ch22:191-197 | IDH ≈8 bn Bs; Bono Juancito Pinto ≈0.5 bn Bs; Renta Dignidad ≈3.5 bn Bs | as listed | A ✗ | fetchable, not done | MEFP / UDAPE fiscal tables (udape.gob.bo responds 200); not fetched in this pass |
| 86 | ch22:199 | fiscal incidence 12% / 8% / 25% | as listed | B | not fetchable here | needs a fiscal-incidence study (e.g. CEQ), not fetched |
| 87 | ch22:201 | government transfers ≈15 bn Bs | 15 bn Bs | A ✗ | fetchable, not done | MEFP; not fetched |
| 88 | ch22:218, 250-251, 364 | alignment example, population shares (0.35, 0.25, 0.15, 0.25), correlations 0.96–0.99 | – | C | remains illustrative | – |
| 89 | ch23:148-180 | five epochs: Γ 0.11→0.10 (1985–90), 0.10→0.12 (1994–2005), 0.12–0.13 (2006–09), 0.13 (2010–19) | as listed | A ✔ | sourced | INE base 1990 book12: 1988 0.563, 1990 0.567, 1994 0.567, 2005 0.502, 2006 0.518, 2009 0.504, 2011 0.478, 2013 0.484 (2014 0.361, an outlier); base 2017: 2017 0.506, 2019 0.467 |
| 90 | ch23:149 | "Bolivia's I-O tables (1990, 2000, 2006, 2012, 2016)" | list of years | A ✔ | sourced | INE publishes annual tables 1988–2009, 2011–2014 (base 1990) and 2017–2024 (base 2017); there is no 2016 table and no 2010 table |
| 91 | ch23:150-176 | σ_B by epoch (0.65→0.68→0.63→0.61→0.62) | as listed | B | not fetchable here | microdata |
| 92 | ch23:176-180 | ΔECI by epoch: +0.02 (1994–2005), −0.01 (2006–09), +0.05 (2010–19) | as listed | A ✔ | sourced | HS92: −0.24 (1995→2005), −0.14 (2006→2009), −0.02 (2010→2019) |
| 93 | ch23:159 | royalties 18% + IDH 32% = 50% | legal | B | not fetchable here | legal text; not checked |
| 94 | ch23:186-232 | Korea path Γ 0.15 → 0.45 (1960–2000), σ_B 0.50 → 0.20 | as listed | B | not fetchable here | pre-1995 Bank of Korea tables not fetched; ICIO gives Korea from 1995 (2021 edition) and 2011–2020 (2023 edition) |
| 95 | ch23:234 | Bolivia growth average 4% per year | 4% | A ✔ | sourced | 4.04% mean (1990–2019); 3.76% (1986–2019) [WDI] |
| 96 | ch23:255-270 | structural velocities: Korea +0.075, China +0.058, Chile +0.033, Bolivia +0.010, Argentina −0.005, Venezuela −0.030 per decade | as listed | A ✔ (partial) | sourced | Bolivia 1990–2013 (base 1990): −0.036/decade. within-edition endpoints per decade: Korea −0.055 (1995–2014, ICIO 2021) / +0.055 (2011–2020, ICIO 2023); Chile +0.091 / +0.107; Argentina +0.037 / −0.050; China −0.038 (2011–2020). Pre-1995 and Venezuela not available (B) |
| 97 | ch23:280 | ρ ≈ 0.85 for Bolivia | 0.85 | C | remains illustrative | – |
| 98 | ch23:303-312 | π_gas 0.07 (2003) → 0.11 (2014) → 0.09 (2019) | as listed | A ✔ | sourced | INE book12 π_hydrocarbons: 0.047 (2003) → 0.117 (2014) (base 1990); 0.039 (2019, base 2017; not comparable across bases) |
| 99 | ch23:330-366 | SCI: Korea 0.86, China 0.61, Chile 0.31, Bolivia 0.075; Γ_ref = 0.50 ("Germany 0.52") | as listed | A ✔ (partial) | sourced | Germany 0.544 (2019). Bolivia's Γ (0.47–0.57) is already near Γ_ref, so the SCI as defined is ill-conditioned for Bolivia |
| 100 | ch24:83-100, 225-235 | N_eff^I by hegemonic period (1.5, 2.0, 1.8, 2.5, 3.0–3.4, 3.5; 2.0→2.8) | as listed | C | remains illustrative | no historical audience data |
| 101 | ch25:45-50 | product space: ~775 SITC-4 products; core proximity 0.35 vs 0.08; ratio 4.4; PCI_D ≈1.8 | as listed | A ✗ | fetchable, not done | Growth Lab product-space files (`top_edges_hs92.csv`) plus PCI; not computed in this pass |
| 102 | ch25:62, 99-108, 306 | truncation τ: Bolivia 0.98 (2 of ~100), Korea 0.18 (82 of ~100); LA 0.70–0.95, East Asia 0.10–0.40 | as listed | A ✗ | fetchable, not done | needs country×product RCA (451 MB) |
| 103 | ch25:107 | ECI Bolivia −0.80, Korea +2.10 | as listed | A ✔ | sourced | 2019 HS92: −1.30 / +1.81 |
| 104 | ch25:108 | Γ(W^P) 0.13 / 0.45 | as listed | A ✔ | sourced | 0.475 / 0.511 (2019, like-for-like) |
| 105 | ch25:109 | σ_B 0.62 / 0.20 | as listed | B | not fetchable here | microdata |
| 106 | ch25:110 | GDP per capita PPP 2019: \$8,800 / \$44,000 | as listed | A ✔ | sourced | \$11,005 / \$46,511 (current international \$, 2021 ICP vintage) [WDI NY.GDP.PCAP.PP.CD] |
| 107 | ch25:171, 195 | schematic product space; relaxation-time penalty | – | C | remains illustrative | – |
| 108 | ch25:247-258 | Korea τ 0.70 → 0.20 and Γ 0.25 → 0.40 (1973–1993) | as listed | B | not fetchable here | pre-1995 data |
| 109 | ch25:273-290 | Bolivia incremental path ω ≈ 0.08, 0.06, 0.04 | as listed | C / A ✗ | remains illustrative | relatedness computable (row 34) |
| 110 | ch25:340 | PCI of lithium batteries ≈1.5 | 1.5 | A ✔ | sourced | HS92 8507 (electric accumulators incl. Li-ion): PCI 0.549 (2019), 0.811 (2023); lithium carbonate (2836) −0.57 (2019) |
| 111 | ch26:43, 47 | thresholds Γ* ≈ 0.30, N*_eff ≈ 8 | – | C | remains illustrative | – |
| 112 | ch26:53-62 | pathology table (Brazil, Mexico, China, Cuba, Venezuela, Korea) | qualitative | C | remains illustrative | – |
| 113 | ch26:73 | Korea Γ below 0.30 until ≈1980 | – | B | not fetchable here | pre-1995 |
| 114 | ch26:76-113 | (Γ, ECI) scatter: BOL (0.13, −0.8), BRA (0.25, 0.8), COL (0.20, −0.3), MEX (0.28, 1.0), KOR (0.45, 2.1), JPN (0.56, 2.3), DEU (0.52, 1.8), CZE (0.40, 1.5), CHL (0.28, 0.3), ARG (0.15, 0.5), VEN (0.10, −0.5), MYS (0.42, 0.8), POL (0.38, 1.2), CHN (0.35, 1.0) | as listed | A ✔ | sourced | 2019, Γ like-for-like (book12_isic, domestic, input shares) and HS92 ECI: BOL (0.475, −1.30), BRA (0.496, 0.48), COL (0.460, 0.14), MEX (0.649, 0.95), KOR (0.511, 1.81), JPN (0.510, 2.05), DEU (0.544, 1.54), CZE (0.486, 1.48), CHL (0.557, −0.06), ARG (0.522, 0.06), MYS (0.381, 1.03), POL (0.552, 1.06), CHN (0.582, 1.30); VEN ECI −1.02, Γ not available (not in ICIO: B) |
| 115 | ch26:165-206, 222, 270 | 30-year plan: costs \$17B, phases, targets | – | C | remains illustrative | plan |
| 116 | ch26:167 | current position: Γ 0.13, ECI −0.80, N_eff 3.4, σ_B 0.62 | as listed | A ✔ / B | sourced | Γ 0.467–0.475, ECI −1.30; N_eff and σ_B B |
| 117 | ch26:208 | Bolivia GDP ≈ \$40B | \$40B | A ✔ | sourced | \$49.1B (2019), \$54.9B (2024), current US\$ [WDI; reflects INE's base-2017 rebasing] |
| 118 | ch27:85-104 | λ ≈ 0.22; community radio 14%; N_eff^I 3.4 → 2.5–3.0 | as listed | C | remains illustrative | derived from illustrative tables |
| 119 | ch27:116 | echo chamber vs viral hub | – | C | remains illustrative | – |
| 120 | ch28:72 | multiple structural equilibria Γ = 0.072 and 0.643 | – | C | remains illustrative | – |
| 121 | ch28:143-160 | portrait: Γ 0.13; G(π) 0.52; ECI −0.80; π_G 0.0005; v_Γ +0.01/decade; SCI 0.075 | as listed | A ✔ | sourced | Γ 0.467; G(π) 0.445; ECI −1.30; π_G 0.00053; v_Γ −0.036/decade (1990–2013) |
| 122 | ch28:143-160 | portrait: σ_B 0.62; N_eff(W^I) 3.4 | as listed | B | not fetchable here | – |
| 123 | ch28:143-160 | portrait: N_eff(W^U) 3.85; Ψ̂ 1.56; τ 0.98; power ratio 1:467; truncation share 55–86% | as listed | C / A ✗ (τ) | remains illustrative | – |

Totals: sourced 42, fetchable but not done 8, not fetchable here 18, remains illustrative 55 (of 123 rows).
