"""Comparator economies: the book's diagnostics on OECD ICIO 2023 domestic blocks.

Input : data/processed/oecd_icio_domestic/<ISO3>.csv.gz, data/processed/oecd_icio_output.csv.gz
        (fetch/fetch_oecd_icio.py)
Output: data/processed/oecd_io_summary.csv   (country, year, aggregation, normalization, Gamma, ...)
        data/processed/oecd_io_linkages.csv  (book12_isic only)
        data/processed/oecd_io_pi.csv        (book12_isic only)
        data/processed/gamma_comparison.csv  Bolivia (INE) vs comparators, like-for-like
                                             (domestic flows, input shares, diagonal kept)

Construction: ICIO domestic block Z_cc = intermediate flows from domestic
industry i to domestic industry j, USD million, basic prices, industry x
industry, 45 ISIC Rev.4 industries (2023 edition).  Output x = 'OUT' row.
Aggregations: 'detail45' and 'book12_isic' (mapping below; the same 12 groups
as the Bolivia 'book12_isic' aggregation in compute_bolivia_io.py).
Caveat: Bolivia's matrix is built from INE SUTs with import proportionality
and the industry-technology assumption; ICIO domestic blocks are built by the
OECD from national SUTs with their own (better) import split.  Differences of a
few hundredths in Gamma are within construction noise.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import iolib  # noqa: E402

DATA = Path(__file__).resolve().parents[1]
PROC = DATA / "processed"
BOOK12 = ["Hydrocarbons", "Mining", "Public administration", "Manufacturing",
          "Construction", "Commerce", "Transport", "Financial services",
          "Agriculture", "Other services", "Personal services", "Domestic services"]


def _num_to_letter(code: str) -> str:
    """2021-edition numeric ISIC division codes ('01T02', '05T06', '45T47') -> ISIC section letter."""
    n = int(code[:2])
    for lo, hi, L in [(1, 3, "A"), (5, 9, "B"), (10, 33, "C"), (35, 35, "D"), (36, 39, "E"), (41, 43, "F"),
                      (45, 47, "G"), (49, 53, "H"), (55, 56, "I"), (58, 63, "J"), (64, 66, "K"), (68, 68, "L"),
                      (69, 75, "M"), (77, 82, "N"), (84, 84, "O"), (85, 85, "P"), (86, 88, "Q"), (90, 93, "R"),
                      (94, 96, "S"), (97, 98, "T")]:
        if lo <= n <= hi:
            return L
    raise KeyError(code)


def icio_map(code: str) -> str:
    if code[:1].isdigit():  # 2021 edition
        L = _num_to_letter(code)
        if code.startswith(("05", "09")):
            return "Hydrocarbons"
        code = L + code
    if code.startswith("A"):
        return "Agriculture"
    if code in ("B05_06", "B09"):          # energy-producing mining + mining support services
        return "Hydrocarbons"
    if code.startswith("B"):
        return "Mining"
    if code.startswith("C"):
        return "Manufacturing"
    if code.startswith(("D", "E", "J", "L", "M", "N", "P", "Q")):
        return "Other services"
    if code.startswith("F"):
        return "Construction"
    if code.startswith("G"):
        return "Commerce"
    if code.startswith("H"):
        return "Transport"
    if code.startswith("I") or code.startswith(("R", "S")):
        return "Personal services"
    if code.startswith("K"):
        return "Financial services"
    if code.startswith("O"):
        return "Public administration"
    if code.startswith("T"):
        return "Domestic services"
    raise KeyError(code)


def main():
    summ, links, pis = [], [], []
    for edition, sfx in (("ICIO2023", ""), ("ICIO2021", "2021")):
        if not (PROC / f"oecd_icio{sfx}_output.csv.gz").exists():
            continue
        run_edition(edition, sfx, summ, links, pis)
    finish(summ, links, pis)


def run_edition(edition, sfx, summ, links, pis):
    out = pd.read_csv(PROC / f"oecd_icio{sfx}_output.csv.gz")
    for f in sorted((PROC / f"oecd_icio{sfx}_domestic").glob("*.csv.gz")):
        c = f.name[:3]
        D = pd.read_csv(f)
        for year, g in D.groupby("year"):
            xo = out[(out.country == c) & (out.year == year)].set_index("industry").output
            inds = list(xo.index)
            Z = g.pivot_table(index="row", columns="col", values="value", aggfunc="sum").reindex(
                index=inds, columns=inds).fillna(0).values
            x = xo.values.astype(float)
            pos = x > 0
            Zp, xp, ip = Z[np.ix_(pos, pos)], x[pos], [i for i, k in zip(inds, pos) if k]
            tag = dict(edition=edition, country=c, year=int(year), aggregation=f"detail{len(ip)}", flows="domestic")
            r, _, _ = iolib.full_report(Zp, xp, ip, tag)
            summ += r
            Zg, xg, groups, _ = iolib.aggregate(Zp, xp, ip, {i: icio_map(i) for i in ip}, order=BOOK12)
            keep = xg > 0
            Zg, xg, groups = Zg[np.ix_(keep, keep)], xg[keep], [gname for gname, k in zip(groups, keep) if k]
            tag = dict(edition=edition, country=c, year=int(year), aggregation="book12_isic", flows="domestic")
            r, lk, p = iolib.full_report(Zg, xg, groups, tag)
            summ += r; links.append(lk); pis.append(p)
        print(edition, c, "done")


def finish(summ, links, pis):
    S = pd.DataFrame(summ)
    S.to_csv(PROC / "oecd_io_summary.csv", index=False)
    pd.concat(links).to_csv(PROC / "oecd_io_linkages.csv", index=False)
    pd.concat(pis).to_csv(PROC / "oecd_io_pi.csv", index=False)

    # like-for-like comparison table
    B = pd.read_csv(PROC / "bolivia_io_summary.csv")
    B = B[(B.flows == "domestic") & (B.base == 2017)].assign(country="BOL (INE base 2017)", edition="INE")
    cols = ["edition", "country", "year", "aggregation", "normalization", "N", "Gamma", "abs_lambda2", "Neff",
            "Neff_over_N", "power_gini", "pi_max", "pi_argmax"]
    C = pd.concat([B[cols], S[cols]])
    C = C[C.aggregation.isin(["book12_isic"]) | C.aggregation.str.startswith("detail")]
    C.to_csv(PROC / "gamma_comparison.csv", index=False)
    k = S[(S.aggregation == "book12_isic") & (S.normalization == "input_shares") & (S.year.isin([1995, 2000, 2005, 2010, 2011, 2014, 2015, 2019, 2020]))]
    print(k.pivot_table(index="country", columns=["edition", "year"], values="Gamma").round(3).to_string())


if __name__ == "__main__":
    main()
