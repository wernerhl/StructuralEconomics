"""Growth Lab (Harvard), Atlas of Economic Complexity data on Harvard Dataverse.

Datasets (Dataverse API, file access by numeric id):
  * "Growth Projections and Complexity Rankings", doi:10.7910/DVN/XTAQMC (v12,
    released 2026-02-08): growth_proj_eci_rankings.csv + data dictionary.
  * "International Trade Data (HS, 92)", doi:10.7910/DVN/T4CHWJ (v18,
    2026-04-22): hs92_country_year.csv (ECI, COI, diversity by country-year),
    hs92_product_year_4.csv (PCI by product-year), data dictionary.
  * "International Trade Data (SITC, Rev. 2)", doi:10.7910/DVN/H8SFD2 (v18,
    2026-04-22): sitc_country_year.csv (long ECI series from 1962), dictionary.

Outputs: data/raw/atlas/*.csv (all < 3 MB) and
  data/processed/eci_selected.csv : ECI for selected countries, both
  classifications, long format (classification, country, year, eci, eci_rank,
  diversity, coi, gdp... whatever the file provides), and
  data/processed/eci_rankings_selected.csv (rankings file, selected countries).

Run: python3 data/fetch/fetch_atlas_eci.py
"""
from __future__ import annotations

import pandas as pd

from common import PROC, RAW, download

API = "https://dataverse.harvard.edu/api/access/datafile/{}"
FILES = {
    "growth_proj_eci_rankings.csv": 13439575,
    "rankings_data_dictionary.csv": 13439574,
    "hs92_country_year.csv": 13685109,
    "hs92_data_dictionary.csv": 13685113,
    "hs92_product_year_4.csv": 13685116,
    "sitc_country_year.csv": 13685143,
    "sitc_data_dictionary.csv": 13685147,
}
SEL = ["BOL", "KOR", "CHL", "CHN", "MEX", "DEU", "USA", "BRA", "ARG", "COL", "PER",
       "JPN", "URY", "VEN", "ECU", "PRY", "NOR", "MYS", "POL", "CZE", "ETH", "LAO"]


def main():
    d = RAW / "atlas"
    paths = {n: download(API.format(i), d / n, timeout=300) for n, i in FILES.items()}
    out = []
    for cls in ("hs92", "sitc"):
        df = pd.read_csv(paths[f"{cls}_country_year.csv"])
        iso = next(c for c in df.columns if c.lower() in ("country_iso3_code", "iso3_code", "location_code", "country_code"))
        df = df[df[iso].isin(SEL)].copy()
        df.insert(0, "classification", cls)
        df = df.rename(columns={iso: "country_iso3"})
        out.append(df)
    E = pd.concat(out, ignore_index=True)
    E.to_csv(PROC / "eci_selected.csv", index=False)
    R = pd.read_csv(paths["growth_proj_eci_rankings.csv"])
    iso = next((c for c in R.columns if "iso3" in c.lower() or c.lower() in ("country_code", "code")), None)
    (R[R[iso].isin(SEL)] if iso else R).to_csv(PROC / "eci_rankings_selected.csv", index=False)
    print(E.columns.tolist())
    print(R.columns.tolist())


if __name__ == "__main__":
    main()
