"""World Bank Poverty and Inequality Platform (PIP) API: Bolivia survey-year
Gini, mean, deciles and poverty headcounts at the 2021-PPP lines
($3.00, $4.20, $8.30 per person per day).

Raw JSON -> data/raw/pip/, tidy CSV -> data/processed/pip_bolivia.csv
Note: Bolivia's PIP welfare aggregate is per-capita household INCOME
(Encuesta de Hogares), not consumption; `survey_comparability` changes mark
breaks in comparability (compare only within a comparable spell).

Run: python3 data/fetch/fetch_pip.py
API docs: https://pip.worldbank.org/api
"""
from __future__ import annotations

import json

import pandas as pd

from common import PROC, RAW, get_json

URL = "https://api.worldbank.org/pip/v1/pip"
LINES = [3.0, 4.2, 8.3]
COUNTRIES = ["BOL"]
KEEP = ["country_code", "reporting_year", "survey_year", "survey_acronym", "welfare_type",
        "reporting_level", "survey_comparability", "comparable_spell", "ppp_version" ,
        "poverty_line", "headcount", "poverty_gap", "gini", "mean", "median", "mld",
        "decile1", "decile2", "decile3", "decile4", "decile5", "decile6", "decile7",
        "decile8", "decile9", "decile10", "reporting_pop", "distribution_type",
        "estimation_type", "is_interpolated"]


def main() -> None:
    out = RAW / "pip"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in COUNTRIES:
        for pl in LINES:
            js = get_json(URL, params={"country": c, "year": "all", "povline": pl,
                                       "format": "json"}, timeout=180)
            (out / f"pip_{c}_{pl}.json").write_text(json.dumps(js))
            for r in js:
                rows.append({k: r.get(k) for k in KEEP})
    df = pd.DataFrame(rows).sort_values(["country_code", "poverty_line", "reporting_year", "reporting_level"])
    df.to_csv(PROC / "pip_bolivia.csv", index=False)
    print(df[df.poverty_line == 3.0][["reporting_year", "reporting_level", "welfare_type",
                                      "survey_comparability", "gini", "headcount"]].to_string())


if __name__ == "__main__":
    main()
