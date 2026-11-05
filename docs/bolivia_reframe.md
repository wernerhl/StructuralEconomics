# Reframing the Bolivia case around the data (authoritative for the data pass)

The author has decided to **reframe the book around the sourced data**. Every number below is produced by a script in `data/` (retrieved 2026-11-05). `data/README.md` gives the source, URL and series for each, and `data/processed/*.csv` holds the values.

Cite data in the text as, for example, "INE supply-use tables, base 2017; author's calculations (data/compute/compute_bolivia_io.py)". Add bibliography entries for datasets as `@misc` (see `docs/bibliography.md`).

## Canonical values

### Domestic production network W^P
W^P = D⁻¹A, computed on domestic flows with input shares. Γ depends on aggregation and normalization, so always state both.

| Quantity | Value | Basis / source |
|---|---|---|
| Bolivia 2019, Γ | **0.467** | book's 12-sector classification (`book12`) |
| Bolivia 2019, N_eff | **6.48 of 11** | same |
| Bolivia 2019, power Gini | **0.445** | same |
| Bolivia 2019, Γ, cross-country basis | **0.475** | ISIC-based 12 sectors (`book12_isic`); use for cross-country comparison |
| Bolivia 2019, power Gini, cross-country basis | **0.435** | same |

Comparators, ICIO 2023, 2019, same `book12_isic` basis:

| | KOR | DEU | USA | CHL | CHN | MEX |
|---|---|---|---|---|---|---|
| Γ | 0.511 | 0.544 | 0.463 | 0.557 | 0.582 | 0.649 |
| Power Gini | 0.658 | 0.611 | — | — | — | — |

For the other countries' power Gini, see `data/processed/oecd_io_summary.csv`.

**Bolivia Γ over time** (INE base 1990): 0.567 (1990), 0.494 (2000), 0.518 (2006), 0.484 (2013). The 2019 value of 0.467 is on the new base 2017 and is not strictly comparable. Korea fell from 0.626 (1995) to 0.511 (2019), from ICIO.

**Normalization sensitivity, Bolivia 2019:** Γ ranges over 0.33–0.69 across input/sales shares and with or without the diagonal. Report this range wherever a single Γ is used to compare economies.

**Hydrocarbons, 2019** (INE):
- domestic power share π = **0.039**;
- **70%** of output exported;
- domestic intermediate inputs **0.225** per unit of output;
- compensation of employees **6%** of the sector's value added.

### Complexity
Growth Lab ECI, Dataverse:
- Bolivia: −1.21 (2006), −1.46 (2018), −1.30 (2019).
- Korea: +1.81 (2019).
- ΔECI by epoch: see `data/processed/eci_epochs_bolivia.csv`.

### Distribution
- Gini (INE): 0.600 (2005) → 0.426 (2019).
- Extreme poverty: 38.2% (2005) → 15.2% (2018) on the old lines; 13.2% (2019) on the updated lines.
- GDP growth averaged 4.52% a year over 2006–19 (WDI).

### Global trade W^P
Growth Lab bilateral trade, 230 economies, 2019:
- π: US 0.089, China 0.130, Germany 0.082, Bolivia 0.00053;
- top-10 share of π: 51.9%;
- Γ_global = 0.452.

### Other figures
Labour forces, GDP, remittances, export shares and the ILOSTAT female employment shares are in `data/processed/bolivia_facts.csv` and the `ilostat_*` CSVs.

## How the narrative changes
1. **Bolivia is not a fragmented outlier in integration.** Its domestic Γ is in the comparator range (the USA's is lower). Γ fell over 1990–2013 in Bolivia and over 1995–2019 in Korea. **Aggregate Γ does not separate development paths.** State this plainly wherever the book used "Bolivia Γ = 0.13" as evidence of dualism or fragmentation.
2. **The enclave is real, but it is measured by π and by export orientation, not by Γ.** The export sector sits at the periphery of the domestic network: hydrocarbons' π is 0.039 while 70% of its output is exported. It has a low domestic input intensity and a low labour share; rents are captured fiscally (IDH) rather than through linkages. This is Pinto's and Prebisch's enclave, measured correctly.
3. **Power concentration in W^P is not automatically bad.** Bolivia's W^P power Gini (0.435) is *lower* than Korea's (0.658) or Germany's (0.611). In industrial economies, π concentrates on manufacturing and business-service hubs, i.e. on industrial depth. Do not read high π concentration in the production network as "oppression". The distributional reading of π belongs to bargaining (W^U/W^I) contexts, under Assumption 7.A.
4. **Bolivia is distinctive in complexity:** ECI −1.30, declining since 2006. **And in distribution:** large falls in inequality and poverty in 2005–2019, while Γ was roughly flat. This is evidence against any claim that equity gains *require* high Γ or ECI > 0. Ch26's conjecture must confront it.
5. **Ch25's truncation conjecture** (missing bridge ⇒ low Γ) gets no support at the 12-sector level. Say so honestly. Keep it as a conjecture, and say what finer-grained test (39-sector INE vs comparable detail, or product-level data) would be needed.
6. **Toy and illustrative examples** (inventory class C) stay. They must stay labelled `\illustrative`, and must not be presented as describing Bolivia's actual structure. Where an illustrative example previously claimed to *be* Bolivia (e.g. "Bolivia Γ = 0.13"), either swap in the real values, or relabel it as a stylized economy ("a dual economy") and add the real Bolivia values next to it.

## Mechanics
- **Chapter scripts:** `code/chNN.py` should read sourced values from `data/processed/*.csv` with the csv module (pandas only if already used). Assert that the text's values match, to the rounding used in the text.
- **Sources:** each table using sourced data gets a source note under it, e.g. `\par\footnotesize Source: …`.
- **Inventory:** `data/inventory.md` lists every claim by chapter and line (as of commit f8e594c; line numbers have shifted, so search by content). Class A rows marked done are the swaps. Class B rows: keep the value as illustrative and add "(not available here)" only if useful.
