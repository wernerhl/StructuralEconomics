"""Structural velocity v_Gamma (Chapter 23) and SCI from real Gamma series.

Gamma = 1 - |lambda_2| of W^P (input shares, diagonal kept, domestic flows, 12-sector
book aggregation).  Velocities are computed WITHIN one data vintage only:
  * Bolivia: INE base 1990 (1988-2014) and INE base 2017 (2017-2024) separately;
  * comparators: OECD ICIO 2021 edition (1995-2014) and 2023 edition (2011-2020) separately.
v_Gamma per decade = OLS slope of Gamma on year x 10 (endpoint differences are also given).
SCI uses Gamma_ref = 0.50 as in ch23:330; it is undefined when Gamma_t0 >= Gamma_ref.

Output: data/processed/gamma_velocity.csv
"""
from pathlib import Path

import numpy as np
import pandas as pd

PROC = Path(__file__).resolve().parents[1] / "processed"


def vel(g):
    g = g.dropna().sort_index()
    y = g.index.values.astype(float)
    slope = np.polyfit(y, g.values, 1)[0] * 10 if len(g) > 2 else np.nan
    return dict(y0=int(y[0]), y1=int(y[-1]), gamma_y0=g.iloc[0], gamma_y1=g.iloc[-1],
                v_endpoint_per_decade=(g.iloc[-1] - g.iloc[0]) / (y[-1] - y[0]) * 10,
                v_ols_per_decade=slope, gamma_min=g.min(), gamma_max=g.max(),
                sci_ref050=((g.iloc[-1] - g.iloc[0]) / (0.5 - g.iloc[0])) if g.iloc[0] < 0.5 else np.nan)


def main():
    rows = []
    B = pd.read_csv(PROC / "bolivia_io_summary.csv")
    B = B[(B.flows == "domestic") & (B.aggregation == "book12") & (B.normalization == "input_shares")]
    for base, lab, rng in [(1990, "INE base 1990", (1988, 2014)), (1990, "INE base 1990 excl. 2014", (1988, 2013)),
                           (2017, "INE base 2017", (2017, 2024)), (1990, "INE base 1990, 1990-2013", (1990, 2013))]:
        g = B[(B.base == base) & B.year.between(*rng)].set_index("year").Gamma
        rows.append(dict(country="BOL", vintage=lab, **vel(g)))
    S = pd.read_csv(PROC / "oecd_io_summary.csv")
    S = S[(S.aggregation == "book12_isic") & (S.normalization == "input_shares")]
    for (ed, c), g in S.groupby(["edition", "country"]):
        rows.append(dict(country=c, vintage=ed, **vel(g.set_index("year").Gamma)))
    out = pd.DataFrame(rows)
    out.to_csv(PROC / "gamma_velocity.csv", index=False)
    print(out[out.country.isin(["BOL", "KOR", "CHN", "CHL", "ARG", "DEU", "USA", "MEX"])].round(3).to_string())


if __name__ == "__main__":
    main()
