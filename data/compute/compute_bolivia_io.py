"""Bolivia: from INE supply-use tables to the book's diagnostics.

Input : data/processed/ine_sut_long.csv  (fetch/fetch_ine_sut.py)
Output: data/processed/bolivia_io_summary.csv      Gamma, pi concentration, ... per
                                                   (base, year, aggregation, flows, normalization)
        data/processed/bolivia_io_linkages.csv     Rasmussen BL/FL per sector
        data/processed/bolivia_io_pi.csv           Perron vector per sector
        data/processed/bolivia_io_sector_facts.csv output/VA/compensation/export/import
                                                   shares by sector (12-sector aggregate)
        data/processed/bolivia_Z_<base>_<year>_<agg>.csv  domestic flow matrices

CONSTRUCTION (state these with any number taken from here)
 * Tables: INE current-price SUTs; base 2017 (2017-2024; 39 industries x 42
   products; 2023-2024 preliminary) and base 1990 (1988-2009, 2011-2014; 35x35).
   The two bases use different classifications and SNA vintages: do not splice.
 * Valuation: purchasers' -> basic prices by removing, product by product, the
   product's net taxes and its trade/transport margins in proportion to every
   use (uniform-rate assumption); removed margins are reallocated to the margin
   products (commerce, transport: products with negative margin entries in the
   supply table) in proportion to their margin supply.
 * Imports: the use tables do not split domestic/imported.  Import
   proportionality: in every non-export use of product p, the imported share is
   m_p / (total basic-price non-export use of p); exports are assumed domestic
   up to domestic output (excess exports are treated as re-exported imports).
 * Industry x industry: fixed product-sales-structure assumption (Eurostat
   model D), Z = D U_dom with D[i,p] = V[p,i]/q_p.  For base-1990 years without a
   published make matrix (all except 2002, 2003, 2006, 2011-2014) V is taken as
   diagonal (in 2006 the make matrix is 98.6% diagonal).
 * Base 1990: the FISIM dummy column ("imputacion bancaria", purchases of financial
   services by a fictitious industry, 1968-SNA style) is allocated to the 35
   industries in proportion to their intermediate consumption.
 * Coefficients: data convention A = Z / x (columns), L = (I - A)^-1 = the book's
   (I - A_book')^-1.  W^P = D^-1 A_book (input shares).  See iolib.py.
 * Aggregations: 'detail' (all industries) and 'book12' (the 12 sectors of
   Example ch17:ex:bolivia_gamma; mapping below).  Gamma depends on both.
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

MAP17 = {**{f"{i:02d}": "Agriculture" for i in (1, 2, 3, 4)},
         "05": "Hydrocarbons", "06": "Mining",
         **{f"{i:02d}": "Manufacturing" for i in range(7, 23)},
         "23": "Other services",            # electricity, water, waste
         "24": "Construction", "25": "Construction",
         "26": "Commerce", "27": "Transport",
         "28": "Personal services",         # accommodation and food
         "29": "Other services", "30": "Financial services",
         "31": "Other services", "32": "Other services",
         "33": "Public administration", "35": "Public administration",
         "37": "Public administration",     # non-market (government) education/health
         "34": "Other services", "36": "Other services",
         "38": "Personal services", "39": "Domestic services"}
# ISIC-comparable variant (used for comparison with OECD ICIO, whose P85/Q86T88 mix
# market and government education/health): non-market education/health -> Other services
MAP17_ISIC = {**MAP17, "35": "Other services", "37": "Other services"}
MAP90 = {**{str(i): "Agriculture" for i in range(1, 6)},
         "6": "Hydrocarbons", "7": "Mining",
         **{str(i): "Manufacturing" for i in range(8, 24)},
         "24": "Other services", "25": "Construction", "26": "Commerce",
         "27": "Transport", "28": "Other services", "29": "Financial services",
         "30": "Other services", "31": "Other services",
         "32": "Personal services", "33": "Personal services",
         "34": "Domestic services", "35": "Public administration"}


def _col(df, pat, optional=False):
    m = df.columns[df.columns.str.contains(pat, case=False, regex=True)]
    if optional and len(m) == 0:
        return pd.Series(0.0, index=df.index)
    if len(m) != 1:
        raise KeyError(f"{pat!r} matched {list(m)}")
    return df[m[0]]


def load_b17(d, year):
    s = d[(d.base == 2017) & (d.year == year) & (d.table == "supply")]
    u = d[(d.base == 2017) & (d.year == year) & (d.table == "use")]
    prods = sorted(c for c in s.row_code.unique() if isinstance(c, str) and c.isdigit())
    inds = sorted(c for c in u.col_code.dropna().unique() if c.isdigit())
    S = s[s.row_code.isin(prods)].pivot_table(index="row_code", columns="col_label", values="value", aggfunc="sum")
    Vt = s[s.row_code.isin(prods) & s.col_code.isin(inds)].pivot_table(
        index="row_code", columns="col_code", values="value", aggfunc="sum").reindex(index=prods, columns=inds).fillna(0)
    U = u[u.row_code.isin(prods) & u.col_code.isin(inds)].pivot_table(
        index="row_code", columns="col_code", values="value", aggfunc="sum").reindex(index=prods, columns=inds).fillna(0)
    F = u[u.row_code.isin(prods) & (u.col_code.isna() | (u.col_code == ""))].pivot_table(
        index="row_code", columns="col_label", values="value", aggfunc="sum").reindex(index=prods).fillna(0)
    S = S.reindex(prods).fillna(0)
    q = _col(S, r"^VALOR BRUTO")
    m = _col(S, r"^IMPORTACIONES")
    tax = sum(_col(S, p, optional=True) for p in [r"^IMPUESTO AL VALOR", r"^DERECHOS", r"^OTROS IMPUESTOS",
                                                    r"^IMPUESTOS IMPLICITOS", r"^SUBVENCIONES$",
                                                    r"^SUBVENCIONES IMPLICITAS"])
    mg = _col(S, r"^M.RGENES")
    ts = _col(S, r"^OFERTA TOTAL")
    e = _col(F, r"^EXPORTACIONES")
    fd = sum(_col(F, p, optional=True) for p in [r"^Hogares$", r"^ISFLHS$", r"^Gobierno Individual$",
                                                   r"^Gobierno Colectivo$", r"^Formaci.n Bruta de Capital Fijo$",
                                                   r"^Variaciones de Existencias$", r"^Objetos valiosos$"])
    td = _col(F, r"^DEMANDA TOTAL")
    va_rows = u[u.col_code.isin(inds) & u.row_code.isin(["P1", "P2", "B1b", "D1", "D2", "D3", "B2b", "B3b"])]
    VA = va_rows.pivot_table(index="row_code", columns="col_code", values="value", aggfunc="sum").reindex(columns=inds)
    labels = (u[u.col_code.isin(inds)].drop_duplicates("col_code").set_index("col_code").col_label.reindex(inds))
    return dict(prods=prods, inds=inds, V=Vt.values, U=U.values, q=q.values, m=m.values, tax=tax.values,
                mg=mg.values, ts=ts.values, e=e.values, fd=fd.values, td=td.values, VA=VA,
                labels=list(labels.values), has_make=True)


def load_b90(d, year):
    x = d[(d.base == 1990) & (d.year == year)]
    prods = [str(i) for i in range(1, 36)]
    inds = prods
    s = x[x.table == "supply"]
    # aggfunc="first": the 2011-2014 sheets repeat the VBP column under the same label
    S = s[s.row_code.isin(prods)].pivot_table(index="row_code", columns="col_label", values="value",
                                              aggfunc="first").reindex(prods).fillna(0)
    q = _col(S, r"^VALOR BRUTO")
    m = _col(S, r"^IMPORTACI")
    tax = sum(_col(S, p) for p in [r"^DERECHOS", r"^IMPUESTO AL VALOR", r"^IMPUESTO A LAS TRANS"])
    mg = _col(S, r"^M.RGENES")
    ts = _col(S, r"^OFERTA TOTAL")
    mk = x[x.table.isin(["make", "supply"]) & x.row_code.isin(prods) & x.col_code.isin(inds)]
    has_make = len(mk) > 0
    if has_make:
        V = mk.pivot_table(index="row_code", columns="col_code", values="value", aggfunc="sum").reindex(
            index=prods, columns=inds).fillna(0).values
    else:
        V = np.diag(q.values)
    u = x[x.table == "use"]
    U = u[u.row_code.isin(prods) & u.col_code.isin(inds)].pivot_table(
        index="row_code", columns="col_code", values="value", aggfunc="sum").reindex(index=prods, columns=inds).fillna(0)
    F = u[u.row_code.isin(prods) & ~u.col_code.isin(inds)].pivot_table(
        index="row_code", columns="col_label", values="value", aggfunc="first").reindex(prods).fillna(0)
    # FISIM ("imputacion bancaria") dummy industry: allocate its purchases (financial
    # services) to the 35 industries in proportion to their intermediate consumption
    fisim = _col(F, r"^IMPUTACI.N BANCARIA", optional=True).values
    U = U + np.outer(fisim, U.values.sum(0) / U.values.sum()) if fisim.sum() else U
    e = _col(F, r"EXPORTACI")
    fd = (_col(F, r"HOGARES$") + _col(F, r"^ADMINISTRACI.N P.BLICA$") + _col(F, r"^FORMACI.N BRUTA")
          + _col(F, r"^VARIACI.N DE EXIST"))
    td = _col(F, r"^DEMANDA TOTAL")
    va = u[u.col_code.isin(inds) & ~u.row_code.isin(prods)]
    VA = va.pivot_table(index="row_label", columns="col_code", values="value", aggfunc="sum").reindex(columns=inds)
    VA.index = ["P1" if i.upper().startswith("PRODUCCI") else "B1b" if i.upper().startswith("VALOR AGREGADO")
                else "P2" if i.upper().startswith("CONSUMO INTERMEDIO") else i for i in VA.index]
    labels = (u[u.col_code.isin(inds)].drop_duplicates("col_code").set_index("col_code").col_label.reindex(inds))
    U = pd.DataFrame(np.asarray(U), index=prods, columns=inds)
    return dict(prods=prods, inds=inds, V=V, U=U.values, q=q.values, m=m.values, tax=tax.values, mg=mg.values,
                ts=ts.values, e=e.values, fd=fd.values, td=td.values, VA=VA, labels=list(labels.values),
                has_make=has_make)


def sut_to_iot(t):
    U, q, m, tax, mg, ts, e, fd = (t[k].astype(float) for k in ("U", "q", "m", "tax", "mg", "ts", "e", "fd"))
    V = t["V"].astype(float)
    checks = {}
    # total demand vs total supply (purchasers' prices)
    tot_use = U.sum(1) + fd + e
    gaps = np.abs(tot_use - ts) / np.maximum(ts, 1)
    checks["max_rel_supply_use_gap"] = float(gaps.max())
    checks["argmax_gap_product"] = t["prods"][int(gaps.argmax())]
    with np.errstate(divide="ignore", invalid="ignore"):
        tau = np.where(ts != 0, tax / ts, 0.0)
        mu = np.where(ts != 0, np.maximum(mg, 0) / ts, 0.0)
    keep = 1 - tau - mu
    Ub, fdb, eb = U * keep[:, None], fd * keep, e * keep
    neg = mg < 0
    w = np.where(neg, mg, 0.0)
    w = w / w.sum() if w.sum() != 0 else w
    Ub = Ub + np.outer(w, (U * mu[:, None]).sum(0))
    fdb = fdb + w * float((fd * mu).sum())
    eb = eb + w * float((e * mu).sum())
    nonexp = Ub.sum(1) + fdb
    # domestic output q is split between exports and non-export uses; exports are
    # domestic up to q (any excess is re-exported imports, as for machinery and
    # basic metals in some base-1990 years).  When eb <= q this is exactly import
    # proportionality: r = m / nonexp.
    ed = np.minimum(eb, q)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.clip(np.where(nonexp > 0, 1 - (q - ed) / nonexp, 0.0), 0, 1)
    Ud = Ub * (1 - r)[:, None]
    fdd = fdb * (1 - r) + ed
    checks["n_products_reexport"] = int((eb > q + 1e-9).sum())
    with np.errstate(divide="ignore", invalid="ignore"):
        Dm = np.where(q[None, :] > 0, V.T / q[None, :], 0.0)  # industries x products
    Zd = Dm @ Ud
    Zt = Dm @ Ub
    x = V.sum(0)  # industry output (basic prices)
    fdi = Dm @ fdd
    checks["min_Z"] = float(Zd.min())
    checks["row_balance_rel_err"] = float(np.max(np.abs(Zd.sum(1) + fdi - Dm @ (Ud.sum(1) + fdd)) / np.maximum(x, 1)))
    ide = np.abs(Zd.sum(1) + fdi - x) / np.maximum(x, 1)
    checks["identity_rel_err_max"] = float(ide.max())
    checks["identity_argmax"] = t["inds"][int(ide.argmax())]
    checks["identity_rel_err_2nd"] = float(np.sort(ide)[-2])
    checks["import_share_of_intermediate"] = float(1 - Ud.sum() / Ub.sum())
    Zd, Zt = np.clip(Zd, 0, None), np.clip(Zt, 0, None)
    exports_i = Dm @ ed
    imports_used_i = (Ub - Ud).sum(0)  # imported intermediates by using industry
    return dict(Zd=Zd, Zt=Zt, x=x, fd=fdi, exports=exports_i, imp_int=imports_used_i, checks=checks)


def sector_facts(t, io, labels, mapping, codes, base, year):
    VA = t["VA"]
    g = io["x"]
    rows = []
    for c, lab in zip(codes, labels):
        def get(k):
            return float(VA.loc[k, c]) if k in VA.index and c in VA.columns and pd.notna(VA.loc[k, c]) else np.nan
        rows.append(dict(base=base, year=year, code=c, industry=lab, book12=mapping[c], output=g[codes.index(c)],
                         va=get("B1b"), comp=get("D1"), ic_purch=get("P2"),
                         exports=io["exports"][codes.index(c)], dom_int_inputs=io["Zd"][:, codes.index(c)].sum(),
                         imp_int_inputs=io["imp_int"][codes.index(c)]))
    return rows


def main():
    d = pd.read_csv(PROC / "ine_sut_long.csv", dtype={"row_code": str, "col_code": str})
    summ, links, pis, facts, checks = [], [], [], [], []
    years = sorted(set(map(tuple, d[["base", "year"]].drop_duplicates().values)))
    for base, year in years:
        t = load_b17(d, year) if base == 2017 else load_b90(d, year)
        io = sut_to_iot(t)
        mapping = MAP17 if base == 2017 else MAP90
        codes = t["inds"]
        labels = [f"{c} {l}"[:60] for c, l in zip(codes, t["labels"])]
        checks.append(dict(base=base, year=year, has_make=t["has_make"], **io["checks"]))
        facts += sector_facts(t, io, t["labels"], mapping, codes, base, year)
        variants = {"domestic": io["Zd"], "total": io["Zt"]}
        for flows, Z in variants.items():
            # detail
            tag = dict(base=base, year=year, aggregation=f"detail{len(codes)}", flows=flows)
            r, lk, p = iolib.full_report(Z, io["x"], labels, tag)
            summ += r; links.append(lk); pis.append(p)
            aggs = {"book12": mapping}
            if base == 2017:
                aggs["book12_isic"] = MAP17_ISIC
            for aname, amap in aggs.items():
                Zg, xg, groups, _ = iolib.aggregate(Z, io["x"], codes, amap, order=BOOK12)
                tag = dict(base=base, year=year, aggregation=aname, flows=flows)
                r, lk, p = iolib.full_report(Zg, xg, groups, tag)
                summ += r; links.append(lk); pis.append(p)
                if flows == "domestic":
                    pd.DataFrame(Zg, index=groups, columns=groups).to_csv(
                        PROC / "bolivia_Z" / f"Z_{base}_{year}_{aname}.csv")
        print(base, year, "make" if t["has_make"] else "diagV", {k: (round(v, 4) if isinstance(v, float) else v) for k, v in io["checks"].items()})
    S = pd.DataFrame(summ)
    S.to_csv(PROC / "bolivia_io_summary.csv", index=False)
    pd.concat(links).to_csv(PROC / "bolivia_io_linkages.csv", index=False)
    pd.concat(pis).to_csv(PROC / "bolivia_io_pi.csv", index=False)
    F = pd.DataFrame(facts)
    F.to_csv(PROC / "bolivia_io_industry_facts.csv", index=False)
    G = F.groupby(["base", "year", "book12"])[["output", "va", "comp", "exports", "dom_int_inputs",
                                                "imp_int_inputs"]].sum(min_count=1).reset_index()
    tot = G.groupby(["base", "year"])[["output", "va"]].transform("sum")
    G["output_share"] = G.output / tot.output
    G["va_share"] = G.va / tot.va
    G["comp_share_of_va"] = G.comp / G.va
    G["export_share_of_output"] = G.exports / G.output
    G["dom_int_inputs_per_output"] = G.dom_int_inputs / G.output
    G["imported_share_of_int_inputs"] = G.imp_int_inputs / (G.imp_int_inputs + G.dom_int_inputs)
    G.to_csv(PROC / "bolivia_io_sector_facts.csv", index=False)
    pd.DataFrame(checks).to_csv(PROC / "bolivia_io_checks.csv", index=False)
    key = S[(S.flows == "domestic") & (S.normalization == "input_shares")]
    print(key[["base", "year", "aggregation", "Gamma", "abs_lambda2", "Neff", "power_gini", "pi_argmax"]].to_string())


if __name__ == "__main__":
    (PROC / "bolivia_Z").mkdir(exist_ok=True)
    main()
