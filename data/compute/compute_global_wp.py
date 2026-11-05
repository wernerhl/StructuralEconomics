"""Global production weight matrix from bilateral goods trade (Chapter 16).

W^P_global[i, j] = share of importer i's goods imports that come from exporter j
(row-stochastic; "imports only", as in Exercise ch16:ex:build).  Diagonal = 0
(no domestic absorption: trade data have no own-country flow).
Source: Growth Lab HS92 bilateral totals (fetch/fetch_atlas_bilateral.py).

Output: data/processed/global_wp_summary.csv  (year, N, Gamma, Neff, power_gini,
          top10_share, pi of USA/CHN/DEU/BOL, max/min ratio, ...)
        data/processed/global_wp_pi.csv       (year, country, pi, colsum, import share)
Only the largest strongly connected set of countries is kept (iteratively
dropping countries with zero imports or zero exports inside the set).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import iolib  # noqa: E402

PROC = Path(__file__).resolve().parents[1] / "processed"


def main():
    df = pd.read_csv(PROC / "atlas_bilateral_selected_years.csv.gz")
    summ, pis = [], []
    for year, g in df.groupby("year"):
        K = g.pivot_table(index="importer", columns="exporter", values="trade_value", aggfunc="sum").fillna(0)
        cs = sorted(set(K.index) | set(K.columns))
        Kv = K.reindex(index=cs, columns=cs).fillna(0).values.copy()
        np.fill_diagonal(Kv, 0)
        K = pd.DataFrame(Kv, index=cs, columns=cs)
        keep = np.ones(len(cs), bool)
        while True:
            M = K.values[np.ix_(keep, keep)]
            bad = (M.sum(0) <= 0) | (M.sum(1) <= 0)
            if not bad.any():
                break
            keep[np.where(keep)[0][bad]] = False
        cc = [c for c, k in zip(cs, keep) if k]
        M = K.loc[cc, cc].values
        W = iolib.row_normalize(M)
        d = iolib.diagnostics(W, cc)
        p = d.pop("pi")
        order = np.argsort(-p)
        s = pd.Series(p, index=cc)
        world_exp = M.sum(0) / M.sum()
        summ.append({"year": int(year), **d, "top10_share": float(p[order[:10]].sum()),
                     "top10": ";".join(cc[i] for i in order[:10]),
                     **{f"pi_{c}": float(s.get(c, np.nan)) for c in ["USA", "CHN", "DEU", "JPN", "BOL", "KOR", "CHL", "LAO"]},
                     "pi_max_over_min": float(p.max() / p.min()), "argmin": cc[int(p.argmin())],
                     "pi_USA_over_LAO": float(s.get("USA", np.nan) / s.get("LAO", np.nan)),
                     "pi_DEU_over_BOL": float(s.get("DEU", np.nan) / s.get("BOL", np.nan))})
        pis.append(pd.DataFrame({"year": int(year), "country": cc, "pi": p, "colsum_W": W.sum(0),
                                 "world_export_share": world_exp}))
        print(year, len(cc), {k: summ[-1][k] for k in ("Gamma", "Neff", "power_gini", "top10_share", "pi_USA", "pi_CHN", "pi_DEU", "pi_BOL")})
    pd.DataFrame(summ).to_csv(PROC / "global_wp_summary.csv", index=False)
    pd.concat(pis).to_csv(PROC / "global_wp_pi.csv", index=False)


if __name__ == "__main__":
    main()
