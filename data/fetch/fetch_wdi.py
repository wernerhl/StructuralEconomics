"""World Bank World Development Indicators (WDI) via the v2 API.

Downloads selected indicators for Bolivia and comparators, saves the raw JSON
in data/raw/wdi/ and writes a tidy long CSV data/processed/wdi_long.csv with
columns: country_iso3, country, indicator, indicator_name, year, value,
source_note, last_updated.

Run: python3 data/fetch/fetch_wdi.py
API docs: https://datahelpdesk.worldbank.org/knowledgebase/articles/889392
"""
from __future__ import annotations

import json

import pandas as pd

from common import PROC, RAW, get_json

COUNTRIES = ["BOL", "KOR", "CHL", "CHN", "MEX", "DEU", "USA", "BRA", "ARG",
             "COL", "PER", "URY", "VEN", "JPN", "LCN", "WLD"]

INDICATORS = {
    # inequality and poverty
    "SI.POV.GINI": "Gini index (World Bank estimate, PIP)",
    "SI.POV.DDAY": "Poverty headcount at international poverty line (% pop.)",
    "SI.POV.LMIC": "Poverty headcount at lower-middle-income line (% pop.)",
    "SI.POV.UMIC": "Poverty headcount at upper-middle-income line (% pop.)",
    "SI.POV.NAHC": "Poverty headcount at national poverty lines (% pop.)",
    "SI.DST.10TH.10": "Income share held by highest 10%",
    "SI.DST.FRST.20": "Income share held by lowest 20%",
    # output, growth, composition
    "NY.GDP.MKTP.CD": "GDP (current US$)",
    "NY.GDP.MKTP.CN": "GDP (current LCU)",
    "NY.GDP.MKTP.KD.ZG": "GDP growth (annual %)",
    "NY.GDP.PCAP.PP.CD": "GDP per capita, PPP (current international $)",
    "NY.GDP.PCAP.PP.KD": "GDP per capita, PPP (constant 2021 international $)",
    "NV.AGR.TOTL.ZS": "Agriculture, forestry, fishing value added (% GDP)",
    "NV.IND.TOTL.ZS": "Industry (incl. construction) value added (% GDP)",
    "NV.IND.MANF.ZS": "Manufacturing value added (% GDP)",
    "NV.SRV.TOTL.ZS": "Services value added (% GDP)",
    "NY.GDP.NGAS.RT.ZS": "Natural gas rents (% GDP)",
    "NY.GDP.PETR.RT.ZS": "Oil rents (% GDP)",
    "NY.GDP.MINR.RT.ZS": "Mineral rents (% GDP)",
    "NY.GDP.TOTL.RT.ZS": "Total natural resources rents (% GDP)",
    # exports
    "TX.VAL.FUEL.ZS.UN": "Fuel exports (% of merchandise exports)",
    "TX.VAL.MMTL.ZS.UN": "Ores and metals exports (% of merchandise exports)",
    "TX.VAL.FOOD.ZS.UN": "Food exports (% of merchandise exports)",
    "TX.VAL.MANF.ZS.UN": "Manufactures exports (% of merchandise exports)",
    "NE.EXP.GNFS.ZS": "Exports of goods and services (% GDP)",
    # remittances
    "BX.TRF.PWKR.CD.DT": "Personal remittances, received (current US$)",
    "BX.TRF.PWKR.DT.GD.ZS": "Personal remittances, received (% GDP)",
    # employment (ILO modelled estimates)
    "SL.AGR.EMPL.ZS": "Employment in agriculture (% of total employment)",
    "SL.IND.EMPL.ZS": "Employment in industry (% of total employment)",
    "SL.SRV.EMPL.ZS": "Employment in services (% of total employment)",
    "SL.EMP.SELF.ZS": "Self-employed (% of total employment)",
    "SL.EMP.VULN.ZS": "Vulnerable employment (% of total employment)",
    "SL.TLF.TOTL.IN": "Labor force, total",
    "SL.TLF.CACT.FE.ZS": "Labor force participation rate, female (% 15+, ILO)",
    "SP.POP.TOTL": "Population, total",
    # financial inclusion (Global Findex, via WDI)
    "FX.OWN.TOTL.FE.ZS": "Account ownership, female (% 15+)",
    "FX.OWN.TOTL.MA.ZS": "Account ownership, male (% 15+)",
    "FX.OWN.TOTL.ZS": "Account ownership, total (% 15+)",
    # digital
    "IT.NET.USER.ZS": "Individuals using the Internet (% of population)",
    "IT.CEL.SETS.P2": "Mobile cellular subscriptions (per 100 people)",
}

BASE = "https://api.worldbank.org/v2/country/{c}/indicator/{i}"


def main() -> None:
    out = RAW / "wdi"
    out.mkdir(parents=True, exist_ok=True)
    rows, failed = [], []
    cstr = ";".join(COUNTRIES)
    for code, label in INDICATORS.items():
        url = BASE.format(c=cstr, i=code)
        try:
            js = get_json(url, params={"format": "json", "per_page": 20000,
                                       "date": "1960:2025"}, timeout=180)
        except Exception as e:  # noqa: BLE001
            failed.append((code, str(e)[:120]))
            print("FAILED", code, e)
            continue
        (out / f"{code}.json").write_text(json.dumps(js))
        meta, data = js[0], (js[1] or [])
        for d in data:
            if d.get("value") is None:
                continue
            rows.append({
                "country_iso3": d["countryiso3code"] or d["country"]["id"],
                "country": d["country"]["value"],
                "indicator": code,
                "indicator_name": d["indicator"]["value"],
                "year": int(d["date"]),
                "value": d["value"],
                "last_updated": meta.get("lastupdated", ""),
            })
        print(f"{code}: {len(data)} records")
    df = pd.DataFrame(rows).sort_values(["indicator", "country_iso3", "year"])
    df.to_csv(PROC / "wdi_long.csv", index=False)
    # convenience: Bolivia wide table
    bol = df[df.country_iso3 == "BOL"].pivot_table(index="year", columns="indicator", values="value")
    bol.to_csv(PROC / "wdi_bolivia_wide.csv")
    if failed:
        pd.DataFrame(failed, columns=["indicator", "error"]).to_csv(PROC / "wdi_failed.csv", index=False)
    print("wrote", PROC / "wdi_long.csv", len(df), "rows; failed:", failed)


if __name__ == "__main__":
    main()
