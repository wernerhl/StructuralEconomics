"""ECI / PCI facts used by Chapters 16, 17, 23, 25, 26, 28.

Inputs : data/processed/eci_rankings_selected.csv (Growth Lab rankings file,
         doi:10.7910/DVN/XTAQMC v12) and data/raw/atlas/hs92_product_year_4.csv
         (doi:10.7910/DVN/T4CHWJ v18).  Run fetch/fetch_atlas_eci.py first.
Outputs: data/processed/eci_facts.csv   (country, year, eci_hs92, rank, eci_sitc, rank)
         data/processed/eci_epochs_bolivia.csv  (Delta ECI over the Ch.23 epochs)
         data/processed/pci_selected_products.csv
ECI is standardised across countries each year (mean 0, sd 1) and is revised
with every Growth Lab release: always quote the release.
"""
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parents[1]
PROC = DATA / "processed"

PRODUCTS = {"8507": "Electric accumulators (incl. lithium-ion batteries)",
            "2836": "Carbonates (incl. lithium carbonate)",
            "2825": "Hydrazine, metal oxides/hydroxides (incl. lithium hydroxide)",
            "2711": "Petroleum gases (natural gas)", "2603": "Copper ores", "2616": "Precious-metal ores (silver)",
            "2608": "Zinc ores", "8001": "Unwrought tin", "1201": "Soya beans", "2304": "Soybean oil-cake",
            "7108": "Gold, unwrought"}


def main():
    R = pd.read_csv(PROC / "eci_rankings_selected.csv")
    R = R.rename(columns={"country_iso3_code": "country"})
    cols = ["country", "year", "eci_hs92", "eci_rank_hs92", "eci_sitc", "eci_rank_sitc", "eci_hs12", "eci_rank_hs12"]
    R[cols].sort_values(["country", "year"]).to_csv(PROC / "eci_facts.csv", index=False)
    b = R[R.country == "BOL"].set_index("year")
    ep = []
    for name, (y0, y1) in {"2. Gas discovery 1994-2005 (1995->2005)": (1995, 2005),
                           "3. Nationalization 2006-2009": (2006, 2009),
                           "4. Redistribution 2010-2019": (2010, 2019),
                           "Full 1995-2019": (1995, 2019), "2006->2018 (ch16 claim)": (2006, 2018),
                           "2019->2024": (2019, 2024)}.items():
        ep.append(dict(epoch=name, y0=y0, y1=y1, eci_hs92_y0=b.eci_hs92[y0], eci_hs92_y1=b.eci_hs92[y1],
                       d_eci_hs92=b.eci_hs92[y1] - b.eci_hs92[y0], eci_sitc_y0=b.eci_sitc[y0],
                       eci_sitc_y1=b.eci_sitc[y1], d_eci_sitc=b.eci_sitc[y1] - b.eci_sitc[y0],
                       rank_hs92_y0=b.eci_rank_hs92[y0], rank_hs92_y1=b.eci_rank_hs92[y1]))
    pd.DataFrame(ep).to_csv(PROC / "eci_epochs_bolivia.csv", index=False)
    P = pd.read_csv(DATA / "raw" / "atlas" / "hs92_product_year_4.csv", dtype={"product_hs92_code": str})
    P = P[P.product_hs92_code.isin(PRODUCTS)].copy()
    P["description"] = P.product_hs92_code.map(PRODUCTS)
    P[["product_hs92_code", "description", "year", "pci", "export_value"]].sort_values(
        ["product_hs92_code", "year"]).to_csv(PROC / "pci_selected_products.csv", index=False)
    print(pd.DataFrame(ep).round(2).to_string())


if __name__ == "__main__":
    main()
