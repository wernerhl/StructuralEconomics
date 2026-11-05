"""Scalar facts for the chapters, derived from the fetched series.

Inputs : processed/wdi_long.csv, processed/pip_bolivia.csv,
         processed/ine_poverty_gini_long.csv, processed/bolivia_io_sector_facts.csv,
         processed/ilostat_bolivia_emp_eco_sex.csv, processed/ilostat_bolivia_emp_status.csv
Output : processed/bolivia_facts.csv  (fact_id, value, unit, year(s), source, derivation)
Run    : python3 data/compute/compute_bolivia_facts.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROC = Path(__file__).resolve().parents[1] / "processed"
BS_PER_USD = 6.96  # official BCB rate, fixed since Nov 2011 (used only for Bs conversions; stated)


def main():
    W = pd.read_csv(PROC / "wdi_long.csv")

    def w(c, i, y):
        s = W[(W.country_iso3 == c) & (W.indicator == i) & (W.year == y)].value
        return float(s.iloc[0]) if len(s) else np.nan

    def wavg_growth(c, y0, y1):
        g = W[(W.country_iso3 == c) & (W.indicator == "NY.GDP.MKTP.KD.ZG") & W.year.between(y0, y1)].value
        return float(np.mean(g)), float((np.prod(1 + g / 100) ** (1 / len(g)) - 1) * 100), len(g)

    F = []

    def add(fid, v, unit, years, src, how):
        F.append(dict(fact_id=fid, value=v, unit=unit, years=years, source=src, derivation=how))

    for (y0, y1) in [(2006, 2019), (2006, 2014), (1990, 2019), (1986, 2019), (2006, 2024)]:
        a, gm, n = wavg_growth("BOL", y0, y1)
        add(f"BOL_gdp_growth_mean_{y0}_{y1}", round(a, 2), "% per year (arithmetic mean)", f"{y0}-{y1}",
            "WDI NY.GDP.MKTP.KD.ZG", f"mean of {n} annual rates")
        add(f"BOL_gdp_growth_geo_{y0}_{y1}", round(gm, 2), "% per year (compound)", f"{y0}-{y1}",
            "WDI NY.GDP.MKTP.KD.ZG", "geometric mean")
    for y in (2019, 2023, 2024):
        add(f"BOL_gdp_usd_{y}", round(w("BOL", "NY.GDP.MKTP.CD", y) / 1e9, 2), "US$ billion (current)", y,
            "WDI NY.GDP.MKTP.CD", "")
    for c in ("BOL", "KOR"):
        add(f"{c}_gdppc_ppp_current_2019", round(w(c, "NY.GDP.PCAP.PP.CD", 2019)), "current intl $ (2021 ICP vintage)", 2019,
            "WDI NY.GDP.PCAP.PP.CD", "")
    for c in ("BOL", "DEU"):
        add(f"{c}_labour_force_2019_mn", round(w(c, "SL.TLF.TOTL.IN", 2019) / 1e6, 2), "million persons", 2019,
            "WDI SL.TLF.TOTL.IN (ILO modelled)", "")
    for y in (2012, 2019, 2024):
        add(f"BOL_population_{y}_mn", round(w("BOL", "SP.POP.TOTL", y) / 1e6, 2), "million", y, "WDI SP.POP.TOTL", "")
    for y in (2007, 2019, 2023, 2024):
        r = w("BOL", "BX.TRF.PWKR.CD.DT", y)
        add(f"BOL_remittances_{y}", round(r / 1e9, 3), "US$ billion", y, "WDI BX.TRF.PWKR.CD.DT", "")
        add(f"BOL_remittances_{y}_Bs", round(r * BS_PER_USD / 1e9, 2), "Bs billion", y, "WDI BX.TRF.PWKR.CD.DT",
            f"x {BS_PER_USD} Bs/US$")
        add(f"BOL_remittances_pct_gdp_{y}", round(w("BOL", "BX.TRF.PWKR.DT.GD.ZS", y), 2), "% of GDP", y,
            "WDI BX.TRF.PWKR.DT.GD.ZS", "")
    for y in (2006, 2013, 2014, 2019, 2023, 2024):
        add(f"BOL_fuel_exports_share_{y}", round(w("BOL", "TX.VAL.FUEL.ZS.UN", y), 1), "% of merchandise exports", y,
            "WDI TX.VAL.FUEL.ZS.UN (UN Comtrade)", "")
        add(f"BOL_ores_metals_exports_share_{y}", round(w("BOL", "TX.VAL.MMTL.ZS.UN", y), 1), "% of merchandise exports",
            y, "WDI TX.VAL.MMTL.ZS.UN", "")
    for y in (2006, 2013, 2019, 2021):
        add(f"BOL_natgas_rents_pct_gdp_{y}", round(w("BOL", "NY.GDP.NGAS.RT.ZS", y), 2), "% of GDP", y,
            "WDI NY.GDP.NGAS.RT.ZS", "")
    for y in (2012, 2019, 2023, 2024):
        add(f"BOL_internet_users_{y}", round(w("BOL", "IT.NET.USER.ZS", y), 1), "% of population", y,
            "WDI IT.NET.USER.ZS (ITU)", "")
    for y in (2017, 2021, 2024):
        for s_, code in (("female", "FX.OWN.TOTL.FE.ZS"), ("male", "FX.OWN.TOTL.MA.ZS")):
            v = w("BOL", code, y)
            if not np.isnan(v):
                add(f"BOL_account_ownership_{s_}_{y}", round(v, 1), "% of age 15+", y, f"WDI {code} (Global Findex)", "")
    for y in (2019, 2024):
        add(f"BOL_self_employed_{y}", round(w("BOL", "SL.EMP.SELF.ZS", y), 1), "% of employment (ILO modelled)", y,
            "WDI SL.EMP.SELF.ZS", "")
        add(f"BOL_female_lfp_{y}", round(w("BOL", "SL.TLF.CACT.FE.ZS", y), 1), "% of women 15+ (ILO modelled)", y,
            "WDI SL.TLF.CACT.FE.ZS", "")

    # INE official poverty and Gini
    P = pd.read_csv(PROC / "ine_poverty_gini_long.csv")
    nat = P[P.area.str.upper().str.startswith("BOLIVIA")]

    def ine(vint, ind, y):
        s = nat[(nat.series_vintage == vint) & nat.indicator.str.startswith(ind) & (nat.year == y)].value
        return float(s.iloc[0]) if len(s) else np.nan
    for y in (2005, 2006, 2018):
        add(f"BOL_INE_extreme_poverty_oldline_{y}", round(ine("old_lines_2005_2018", "Incidencia de pobreza extrema", y), 1),
            "% of population", y, "INE cuadro 3.06.01.02 (EH 2005-2018, old lines)", "")
        add(f"BOL_INE_poverty_oldline_{y}", round(ine("old_lines_2005_2018", "Incidencia de pobreza (FGT0)", y), 1),
            "% of population", y, "INE cuadro 3.06.01.01 (EH 2005-2018, old lines)", "")
    for y in (2016, 2019, 2024, 2025):
        add(f"BOL_INE_extreme_poverty_EPF_{y}", round(ine("EPF2015_lines_2016_2025_rev25", "Incidencia de pobreza extrema", y), 1),
            "% of population", y, "INE cuadro 3.06.01.02-rev25 (EPF 2015-16 lines)", "")
        add(f"BOL_INE_poverty_EPF_{y}", round(ine("EPF2015_lines_2016_2025_rev25", "Incidencia de pobreza (FGT0)", y), 1),
            "% of population", y, "INE cuadro 3.06.01.01-rev25 (EPF 2015-16 lines)", "")
    g = P[(P.series_vintage == "gini_2005_2025_rev25") & (P.area == "Bolivia")].set_index("year").value
    for y in (2005, 2006, 2018, 2019, 2024, 2025):
        add(f"BOL_INE_gini_{y}", round(float(g[y]), 3), "Gini, per-capita household income", y,
            "INE cuadro 3.06.01.03-rev25", "")
    pip = pd.read_csv(PROC / "pip_bolivia.csv")
    pp = pip[(pip.poverty_line == 3.0) & (pip.reporting_level == "national")].set_index("reporting_year")
    for y in (2000, 2005, 2006, 2019, 2024):
        add(f"BOL_PIP_gini_{y}", round(float(pp.gini[y]), 3), "Gini (income), PIP", y, "World Bank PIP API", "")
        add(f"BOL_PIP_pov300_{y}", round(float(pp.headcount[y]) * 100, 1), "% below $3.00/day 2021 PPP", y,
            "World Bank PIP API", "")

    # INE SUT: hydrocarbons and gas-sector structure
    S = pd.read_csv(PROC / "bolivia_io_sector_facts.csv")
    for (b, y) in [(1990, 2000), (1990, 2006), (1990, 2013), (2017, 2017), (2017, 2019), (2017, 2023), (2017, 2024)]:
        r = S[(S.base == b) & (S.year == y) & (S.book12 == "Hydrocarbons")]
        if len(r):
            r = r.iloc[0]
            add(f"BOL_hydrocarbons_va_share_{y}_b{b}", round(r.va_share * 100, 2), "% of total gross value added (basic prices)",
                y, f"INE SUT base {b}", "industry 05 (b2017) / 6 (b1990) VA over total industry VA")
            add(f"BOL_hydrocarbons_output_share_{y}_b{b}", round(r.output_share * 100, 2), "% of gross output", y,
                f"INE SUT base {b}", "")
            add(f"BOL_hydrocarbons_export_share_of_output_{y}_b{b}", round(r.export_share_of_output * 100, 1),
                "% of sector output exported", y, f"INE SUT base {b}", "model-D industry exports / output")
            add(f"BOL_hydrocarbons_dom_inputs_per_output_{y}_b{b}", round(r.dom_int_inputs_per_output, 3),
                "domestic intermediate inputs per unit of output", y, f"INE SUT base {b}", "")
            if not np.isnan(r.comp_share_of_va):
                add(f"BOL_hydrocarbons_comp_share_of_va_{y}_b{b}", round(r.comp_share_of_va * 100, 1),
                    "% of sector VA paid as compensation of employees", y, f"INE SUT base {b}", "D1 / B1b")

    # ILOSTAT employment shares
    E = pd.read_csv(PROC / "ilostat_bolivia_emp_eco_sex.csv")
    e = E[(E.year == 2019) & (E.sex == "SEX_T")].set_index("eco")
    for code, lab in [("ECO_ISIC4_A", "agriculture"), ("ECO_ISIC4_B", "mining_quarrying_incl_hydrocarbons"),
                      ("ECO_ISIC4_C", "manufacturing"), ("ECO_ISIC4_F", "construction"), ("ECO_ISIC4_G", "commerce"),
                      ("ECO_ISIC4_H", "transport"), ("ECO_ISIC4_O", "public_administration"),
                      ("ECO_ISIC4_T", "domestic_services"), ("ECO_ISIC4_K", "finance")]:
        add(f"BOL_emp_share_{lab}_2019", round(float(e.share_of_sex_total[code]) * 100, 1), "% of employment", 2019,
            "ILOSTAT DF_EMP_TEMP_SEX_ECO_NB", f"{code}")
        add(f"BOL_female_share_{lab}_2019", round(float(e.female_share_of_sector[code]) * 100, 1),
            "% of sector employment that is female", 2019, "ILOSTAT DF_EMP_TEMP_SEX_ECO_NB", f"{code}")
    pd.DataFrame(F).to_csv(PROC / "bolivia_facts.csv", index=False)
    print(pd.DataFrame(F)[["fact_id", "value", "unit"]].to_string())


if __name__ == "__main__":
    main()
