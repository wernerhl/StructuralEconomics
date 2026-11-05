"""ILOSTAT (ILO) SDMX API: Bolivia employment by sex and economic activity and
by status in employment, from household surveys (Encuesta de Hogares /
Encuesta Continua de Empleo, harmonised by ILO; thousands of persons).

Dataflows:
  DF_EMP_TEMP_SEX_ECO_NB : employment by sex and economic activity (ISIC Rev.4
                           sections ECO_ISIC4_*, and aggregate sectors)
  DF_EMP_TEMP_SEX_STE_NB : employment by sex and status in employment
Endpoint: https://sdmx.ilo.org/rest/data/ILO,<DF>,1.0/BOL.A...?format=csvfile

Outputs: data/raw/ilostat/*.csv; data/processed/ilostat_bolivia_emp_eco_sex.csv
         (year, sex, eco, value_thousands, share_of_total_employment,
         female_share_of_sector, source, notes) and
         data/processed/ilostat_bolivia_emp_status.csv
Run: python3 data/fetch/fetch_ilostat.py
"""
from __future__ import annotations

import pandas as pd

from common import PROC, RAW, download

URL = "https://sdmx.ilo.org/rest/data/ILO,{df},1.0/BOL.A...?format=csvfile"


def main():
    d = RAW / "ilostat"
    eco = pd.read_csv(download(URL.format(df="DF_EMP_TEMP_SEX_ECO_NB"), d / "DF_EMP_TEMP_SEX_ECO_NB_BOL.csv"))
    ste = pd.read_csv(download(URL.format(df="DF_EMP_TEMP_SEX_STE_NB"), d / "DF_EMP_TEMP_SEX_STE_NB_BOL.csv"))
    e = eco[eco.ECO.str.startswith("ECO_ISIC4_") | (eco.ECO == "ECO_SECTOR_TOTAL")].copy()
    e = e.rename(columns={"TIME_PERIOD": "year", "SEX": "sex", "ECO": "eco", "OBS_VALUE": "value_thousands",
                          "SOURCE": "source", "NOTE_SOURCE": "note_source"})
    # one source per year: prefer the source with the most records that year
    e["n"] = e.groupby(["year", "source"]).value_thousands.transform("size")
    best = e.sort_values("n", ascending=False).drop_duplicates("year")[["year", "source"]]
    e = e.merge(best, on=["year", "source"])
    tot = e[(e.eco == "ECO_ISIC4_TOTAL") | (e.eco == "ECO_SECTOR_TOTAL")].drop_duplicates(["year", "sex"]).set_index(
        ["year", "sex"]).value_thousands
    e["share_of_sex_total"] = e.apply(lambda r: r.value_thousands / tot.get((r.year, r.sex), float("nan")), axis=1)
    fem = e[e.sex == "SEX_F"].set_index(["year", "eco"]).value_thousands
    allx = e[e.sex == "SEX_T"].set_index(["year", "eco"]).value_thousands
    e["female_share_of_sector"] = e.apply(lambda r: fem.get((r.year, r.eco), float("nan")) /
                                          allx.get((r.year, r.eco), float("nan")), axis=1)
    e[["year", "sex", "eco", "value_thousands", "share_of_sex_total", "female_share_of_sector", "source",
       "note_source", "NOTE_CLASSIF", "NOTE_INDICATOR"]].sort_values(["year", "sex", "eco"]).to_csv(
        PROC / "ilostat_bolivia_emp_eco_sex.csv", index=False)
    s = ste.rename(columns={"TIME_PERIOD": "year", "SEX": "sex", "STE": "status", "OBS_VALUE": "value_thousands",
                            "SOURCE": "source"})
    s[["year", "sex", "status", "value_thousands", "source"]].sort_values(["year", "sex", "status"]).to_csv(
        PROC / "ilostat_bolivia_emp_status.csv", index=False)
    print(e.year.unique(), s.year.unique())


if __name__ == "__main__":
    main()
