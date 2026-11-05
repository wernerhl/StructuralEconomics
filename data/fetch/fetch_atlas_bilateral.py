"""Growth Lab bilateral goods trade totals (country x partner x year), HS92,
Harvard Dataverse doi:10.7910/DVN/T4CHWJ (v18, 2026-04-22), file
hs92_country_country_year.csv (id 13685107, ~31 MB -> git-ignored, re-downloaded
by this script).  Values are the Growth Lab's reconciled ("mirror-cleaned")
export/import values in current USD.

Writes data/processed/atlas_bilateral_selected_years.csv.gz : rows for the
years used by compute/compute_global_wp.py (exporter, importer, year,
trade_value), where trade_value = the exporter's reported/reconciled exports to
the importer.

Run: python3 data/fetch/fetch_atlas_bilateral.py
"""
from __future__ import annotations

import pandas as pd

from common import PROC, RAW, download

URL = "https://dataverse.harvard.edu/api/access/datafile/13685107"
YEARS = [1995, 2000, 2006, 2010, 2014, 2018, 2019, 2022, 2023]


def main():
    p = download(URL, RAW / "atlas" / "hs92_country_country_year.csv", timeout=900)
    df = pd.read_csv(p)
    print(df.columns.tolist(), df.year.min(), df.year.max())
    df = df[df.year.isin(YEARS)]
    out = pd.DataFrame({
        "exporter": df["country_iso3_code"], "importer": df["partner_iso3_code"],
        "year": df["year"], "trade_value": df["export_value"],
    })
    out = out[(out.trade_value > 0) & (out.exporter != out.importer)]
    out.to_csv(PROC / "atlas_bilateral_selected_years.csv.gz", index=False)
    print(len(out), "rows")


if __name__ == "__main__":
    main()
