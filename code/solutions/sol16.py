"""Solutions manual, Chapter 16 (The Global Weight Matrix): computes and prints every
number used in solutions/sol16.tex, with asserts on the key results.
Run from the repo root:  python3 code/solutions/sol16.py
numpy + pandas (pandas reads data/processed and the tracked Growth Lab ranking
file data/raw/atlas/growth_proj_eci_rankings.csv)."""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=120)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROC = os.path.join(ROOT, "data", "processed")
RAW = os.path.join(ROOT, "data", "raw", "atlas")


def gap(W):
    return 1 - np.sort(np.abs(np.linalg.eigvals(W)))[-2]


def perron_eig(W):
    ev, vl = np.linalg.eig(W.T)
    p = np.real(vl[:, np.argmin(np.abs(ev - 1))])
    return p / p.sum()


def power_iteration(W, tol=1e-14, maxit=100000):
    """pi'_{k+1} = pi'_k W, from the uniform vector (Appendix 6A); W aperiodic here."""
    p = np.ones(len(W)) / len(W)
    for k in range(maxit):
        q = p @ W
        if np.abs(q - p).sum() < tol:
            return q / q.sum(), k + 1
        p = q
    raise RuntimeError("no convergence")


def neff(p):
    return 1 / np.sum(p ** 2)


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def build_K(df, year):
    g = df[df.year == year]
    K = g.pivot_table(index="importer", columns="exporter", values="trade_value", aggfunc="sum").fillna(0)
    cs = sorted(set(K.index) | set(K.columns))
    Kv = K.reindex(index=cs, columns=cs).fillna(0).values.copy()
    np.fill_diagonal(Kv, 0)
    keep = np.ones(len(cs), bool)
    while True:                                    # largest set with positive imports and exports
        M = Kv[np.ix_(keep, keep)]
        bad = (M.sum(0) <= 0) | (M.sum(1) <= 0)
        if not bad.any():
            break
        keep[np.where(keep)[0][bad]] = False
    cc = [c for c, k in zip(cs, keep) if k]
    return Kv[np.ix_(keep, keep)], cc


# =====================================================================
print("== Exercise 16.1: the global W^P from Growth Lab bilateral goods trade (HS92), imports only ==")
bil = pd.read_csv(os.path.join(PROC, "atlas_bilateral_selected_years.csv.gz"))
summ = pd.read_csv(os.path.join(PROC, "global_wp_summary.csv")).set_index("year")
RES = {}
for year in (2018, 2019):
    K, cc = build_K(bil, year)
    W = K / K.sum(1, keepdims=True)
    assert np.allclose(W.sum(1), 1) and np.all(np.diag(W) == 0) and W.min() >= 0
    p_eig = perron_eig(W)
    p, its = power_iteration(W)
    assert np.abs(p - p_eig).max() < 1e-10
    col = W.sum(0)
    s = K.sum(0) / K.sum()                         # world export shares (within the set)
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    RES[year] = dict(K=K, cc=cc, W=W, p=p, col=col, s=s, G=1 - ev[1], ev=ev)
    order = np.argsort(-p)
    print(f"{year}: N = {len(cc)}, row sums in [{W.sum(1).min():.12f}, {W.sum(1).max():.12f}], "
          f"power iteration converged in {its} steps")
    print(f"  Gamma = {1-ev[1]:.4f} (|l2| = {ev[1]:.4f}, |l3| = {ev[2]:.4f}); Neff = {neff(p):.2f}; "
          f"power Gini = {gini(p):.3f}; top-10 share = {p[order[:10]].sum():.3f}")
    print("  top 10 pi:", [(cc[i], round(p[i], 4)) for i in order[:10]])
    oc = np.argsort(-col)
    print("  largest column sums:", [(cc[i], round(col[i], 2)) for i in oc[:5]],
          " smallest:", [(cc[i], round(col[i], 6)) for i in oc[-3:]])
    print(f"  mean column sum = {col.mean():.6f}; BOL column sum = {col[cc.index('BOL')]:.4f}, pi_BOL = {p[cc.index('BOL')]:.5f}")
    print(f"  corr(pi, s) = {np.corrcoef(p, s)[0,1]:.4f}; 0.5||pi - s||_1 = {0.5*np.abs(p-s).sum():.4f}; "
          f"Neff(s) = {neff(s):.2f}; 0.5||pi - 1/N||_1 = {0.5*np.abs(p-1/len(p)).sum():.4f}")
    for c in ("CHN", "USA", "DEU", "BOL"):
        print(f"    {c}: pi = {p[cc.index(c)]:.5f}, export share s = {s[cc.index(c)]:.5f}")
    if year == 2019:
        assert abs((1 - ev[1]) - summ.loc[2019, "Gamma"]) < 1e-9 and abs(neff(p) - summ.loc[2019, "Neff"]) < 1e-6
        assert round(1 - ev[1], 3) == 0.452 and round(p[cc.index("BOL")], 5) == 0.00053
    else:
        assert abs((1 - ev[1]) - summ.loc[2018, "Gamma"]) < 1e-9
r18 = RES[2018]
K, cc, p, W = r18["K"], r18["cc"], r18["p"], r18["W"]
# 50-node aggregation: the 49 largest traders (imports + exports) plus a Rest-of-World node
tot = K.sum(0) + K.sum(1)
top = list(np.argsort(-tot)[:49])
rest = [i for i in range(len(cc)) if i not in top]
idx = top + [None]


def agg(K, keep_row_diag=True):
    n = len(top) + 1
    A = np.zeros((n, n))
    groups = [[i] for i in top] + [rest]
    for a, ga in enumerate(groups):
        for b, gb in enumerate(groups):
            A[a, b] = K[np.ix_(ga, gb)].sum()
    if not keep_row_diag:
        A[-1, -1] = 0
    return A


for keep in (True, False):
    Ka = agg(K, keep)
    Wa = Ka / Ka.sum(1, keepdims=True)
    pa = perron_eig(Wa)
    lab = "ROW internal trade kept on the diagonal" if keep else "ROW diagonal set to zero"
    print(f"50-node aggregate ({lab}): Gamma = {gap(Wa):.4f}, Neff = {neff(pa):.2f}, "
          f"pi_ROW = {pa[-1]:.4f}, W_ROW,ROW = {Wa[-1,-1]:.4f}, pi_USA = {pa[top.index(cc.index('USA'))]:.4f}, "
          f"pi_CHN = {pa[top.index(cc.index('CHN'))]:.4f}")
    if keep:
        AG = dict(G=gap(Wa), Neff=neff(pa), piROW=pa[-1], WRR=Wa[-1, -1])
    else:
        AG0 = dict(G=gap(Wa), Neff=neff(pa))
print(f"ROW node = {len(rest)} economies; share of world trade in the 49 = {tot[top].sum()/tot.sum():.3f}")
print(f"aggregation changes: Gamma {r18['G']:.4f} -> {AG['G']:.4f}; Neff {neff(p):.2f} -> {AG['Neff']:.2f} "
      f"(Neff/N {neff(p)/len(cc):.4f} -> {AG['Neff']/50:.4f})")
# pi vs GDP per capita: only the WDI comparators in data/processed
wdi = pd.read_csv(os.path.join(PROC, "wdi_long.csv"))
gdp = wdi[(wdi.indicator == "NY.GDP.PCAP.PP.CD") & (wdi.year == 2018)].set_index("country_iso3").value
gdp = gdp.drop([c for c in ("LCN", "WLD") if c in gdp.index])
cs_ = [c for c in gdp.index if c in cc]
pv = np.array([p[cc.index(c)] for c in cs_]); gv = gdp.loc[cs_].values
r_pi_gdp = np.corrcoef(pv, gv)[0, 1]; r_lpi_lgdp = np.corrcoef(np.log(pv), np.log(gv))[0, 1]
print(f"pi vs GDP per capita PPP 2018, n = {len(cs_)} WDI comparators {cs_}: corr = {r_pi_gdp:.3f}; "
      f"corr(log, log) = {r_lpi_lgdp:.3f}")
# size, not income: corr(log pi, log GDP) where GDP = level
gdpl = wdi[(wdi.indicator == "NY.GDP.MKTP.CD") & (wdi.year == 2018)].set_index("country_iso3").value
gl = gdpl.loc[cs_].values
r_lpi_lsize = np.corrcoef(np.log(pv), np.log(gl))[0, 1]
print(f"corr(log pi, log GDP level, current US$) = {r_lpi_lsize:.3f}")
G16_1 = dict(r_pi_gdp=r_pi_gdp, r_log=r_lpi_lgdp, r_size=r_lpi_lsize, n=len(cs_))

# =====================================================================
print("\n== Exercise 16.2: ECI (Growth Lab HS92, release v12 rankings file) ==")
rk = pd.read_csv(os.path.join(RAW, "growth_proj_eci_rankings.csv"))
r18e = rk[(rk.year == 2018) & (rk.in_rankings)].dropna(subset=["eci_hs92"]).sort_values("eci_hs92", ascending=False)
print(f"2018: {len(r18e)} ranked economies")
print("top 10:", [(a, b) for a, b in zip(r18e.country_iso3_code.head(10), r18e.eci_hs92.head(10))])
print("bottom 10:", [(a, b) for a, b in zip(r18e.country_iso3_code.tail(10), r18e.eci_hs92.tail(10))])
e18 = r18e.set_index("country_iso3_code")
bol, kor = e18.loc["BOL"], e18.loc["KOR"]
print(f"BOL: ECI {bol.eci_hs92} (rank {int(bol.eci_rank_hs92)}); KOR: ECI {kor.eci_hs92} (rank {int(kor.eci_rank_hs92)}); "
      f"gap = {kor.eci_hs92 - bol.eci_hs92:.2f}")
assert round(bol.eci_hs92, 2) == -1.46
cs2 = [c for c in cs_ if c in e18.index]
ev_ = e18.loc[cs2, "eci_hs92"].values; gv2 = gdp.loc[cs2].values
r_eci_gdp = np.corrcoef(ev_, gv2)[0, 1]; r_eci_lgdp = np.corrcoef(ev_, np.log(gv2))[0, 1]
pv2 = np.array([p[cc.index(c)] for c in cs2])
print(f"ECI vs GDP pc 2018, same n = {len(cs2)}: corr = {r_eci_gdp:.3f}; corr(ECI, log GDP pc) = {r_eci_lgdp:.3f}; "
      f"for comparison corr(log pi, log GDP pc) = {np.corrcoef(np.log(pv2), np.log(gv2))[0,1]:.3f}")
# all ranked countries: ECI vs log pi (both available for every country)
common = [c for c in e18.index if c in cc]
r_eci_lpi = np.corrcoef(e18.loc[common, "eci_hs92"].values, np.log([p[cc.index(c)] for c in common]))[0, 1]
print(f"ECI vs log pi over all {len(common)} ranked economies in the trade W: corr = {r_eci_lpi:.3f}")
# ECI vs Gini: latest Gini 2014-2019 for the WDI comparators
gin = wdi[(wdi.indicator == "SI.POV.GINI") & (wdi.year.between(2014, 2019))].dropna(subset=["value"])
gin = gin.sort_values("year").groupby("country_iso3").last()
cs3 = [c for c in gin.index if c in e18.index]
print("Gini (latest 2014-2019):", {c: (round(gin.loc[c, 'value'], 1), int(gin.loc[c, 'year'])) for c in cs3})
r_eci_gini = np.corrcoef(e18.loc[cs3, "eci_hs92"].values, gin.loc[cs3, "value"].values)[0, 1]
print(f"ECI 2018 vs Gini, n = {len(cs3)}: corr = {r_eci_gini:.3f}")
print("  (ECI, Gini):", {c: (e18.loc[c, 'eci_hs92'], round(gin.loc[c, 'value'], 1)) for c in cs3})
G16_2 = dict(gapBK=kor.eci_hs92 - bol.eci_hs92, r_eci_gdp=r_eci_gdp, r_eci_gini=r_eci_gini)

# =====================================================================
print("\n== Exercise 16.4: Bolivia across the boom ==")
ef = pd.read_csv(os.path.join(PROC, "eci_facts.csv"))
eb = ef[ef.country == "BOL"].set_index("year")
print("Bolivia ECI (HS92):", {y: eb.loc[y, "eci_hs92"] for y in (2000, 2003, 2006, 2009, 2012, 2014, 2016, 2019)})
d_boom = eb.loc[2014, "eci_hs92"] - eb.loc[2003, "eci_hs92"]
print(f"ECI change 2003 -> 2014 = {d_boom:+.3f}; 1995 -> 2003 = {eb.loc[2003,'eci_hs92']-eb.loc[1995,'eci_hs92']:+.3f}; "
      f"2014 -> 2019 = {eb.loc[2019,'eci_hs92']-eb.loc[2014,'eci_hs92']:+.3f}")
gp = pd.read_csv(os.path.join(PROC, "global_wp_pi.csv"))
pb = gp[gp.country == "BOL"].set_index("year")
print("pi_BOL by year:", {int(y): round(v, 5) for y, v in pb.pi.items()})
print("BOL world export share by year:", {int(y): round(v, 5) for y, v in pb.world_export_share.items()})
fuel = wdi[(wdi.indicator == "TX.VAL.FUEL.ZS.UN") & (wdi.country_iso3 == "BOL")].set_index("year").value
print("fuel exports % of merchandise exports:", {y: round(fuel.get(y, np.nan), 1) for y in (2000, 2006, 2010, 2014, 2018, 2019)})
# partner decomposition of the BOL column of W between 2006 and 2014 (values only)
dec = {}
for y in (2006, 2014):
    Ky, ccy = build_K(bil, y)
    Wy = Ky / Ky.sum(1, keepdims=True)
    j = ccy.index("BOL")
    dec[y] = {c: Wy[i, j] for i, c in enumerate(ccy) if Wy[i, j] > 0}
both = set(dec[2006]) & set(dec[2014])
new = set(dec[2014]) - set(dec[2006]); lost = set(dec[2006]) - set(dec[2014])
cs06, cs14 = sum(dec[2006].values()), sum(dec[2014].values())
intensive = sum(dec[2014][c] - dec[2006][c] for c in both)
ext = sum(dec[2014][c] for c in new) - sum(dec[2006][c] for c in lost)
print(f"BOL column sum of W: 2006 = {cs06:.4f}, 2014 = {cs14:.4f}; partners {len(dec[2006])} -> {len(dec[2014])} "
      f"(new {len(new)}, lost {len(lost)}); intensive margin {intensive:+.4f}, extensive margin {ext:+.4f}")
top_part = sorted(((dec[2014].get(c, 0) - dec[2006].get(c, 0), c) for c in set(dec[2006]) | set(dec[2014])), reverse=True)[:4]
print("largest increases in BOL's share of a partner's imports:", [(c, round(v, 4)) for v, c in top_part])
assert abs(intensive + ext - (cs14 - cs06)) < 1e-12
G16_4 = dict(d_boom=d_boom, intensive=intensive, ext=ext)

# =====================================================================
print("\n== Exercise 16.6: domestic Gamma and Neff (OECD ICIO 2023, 2019, book12_isic, domestic, input shares) ==")
osum = pd.read_csv(os.path.join(PROC, "oecd_io_summary.csv"))
o19 = osum[(osum.edition == "ICIO2023") & (osum.year == 2019) & (osum.aggregation == "book12_isic")
           & (osum.normalization == "input_shares")].set_index("country")
bsum = pd.read_csv(os.path.join(PROC, "bolivia_io_summary.csv"))
b19 = bsum[(bsum.year == 2019) & (bsum.base == 2017) & (bsum.aggregation == "book12_isic") & (bsum.flows == "domestic")
           & (bsum.normalization == "input_shares")].iloc[0]
tab = o19[["N", "Gamma", "Neff", "power_gini"]].copy()
tab.loc["BOL"] = [b19.N, b19.Gamma, b19.Neff, b19.power_gini]
tab = tab.sort_values("Gamma")
print(tab.round(3).to_string())
print(f"n = {len(tab)}; Gamma range {tab.Gamma.min():.3f}-{tab.Gamma.max():.3f}; Neff range {tab.Neff.min():.2f}-{tab.Neff.max():.2f}; "
      f"corr(Gamma, Neff) = {np.corrcoef(tab.Gamma, tab.Neff)[0,1]:.3f}")
assert round(tab.loc["BOL", "Gamma"], 3) == 0.475 and round(tab.loc["KOR", "Gamma"], 3) == 0.511
# (b) no function of Gamma bounds Var(pi)
N = 11
for gam_ in (0.05, 0.5, 1.0):
    for p1 in (1 / N, 0.5, 0.95):
        pi_ = np.r_[p1, np.full(N - 1, (1 - p1) / (N - 1))]
        Wg = gam_ * np.outer(np.ones(N), pi_) + (1 - gam_) * np.eye(N)
        assert abs(gap(Wg) - gam_) < 1e-10 and np.allclose(perron_eig(Wg), pi_)
        print(f"  gamma = {gam_:.2f}, pi_1 = {p1:.3f}: Gamma = {gap(Wg):.3f}, Var(pi) = {np.mean((pi_-1/N)**2):.5f}")
print(f"  Var(pi) max for N = 11: (N-1)/N^2 = {(N-1)/N**2:.5f}")
# (c) regressions on the small sample with Gini, ECI and Gamma
e19 = rk[(rk.year == 2019) & (rk.in_rankings)].set_index("country_iso3_code")
cs4 = [c for c in gin.index if c in tab.index and c in e19.index]
dat = pd.DataFrame({"gini": gin.loc[cs4, "value"], "eci": e19.loc[cs4, "eci_hs92"],
                    "Gamma": tab.loc[cs4, "Gamma"], "Neff": tab.loc[cs4, "Neff"],
                    "lgdp": np.log(wdi[(wdi.indicator == "NY.GDP.PCAP.PP.CD") & (wdi.year == 2019)]
                                   .set_index("country_iso3").value.loc[cs4])})
print(dat.round(3).to_string())


def ols(y, Xcols):
    X = np.column_stack([np.ones(len(dat))] + [dat[c].values for c in Xcols])
    b, *_ = np.linalg.lstsq(X, dat[y].values, rcond=None)
    e = dat[y].values - X @ b
    n_, k_ = X.shape
    s2 = e @ e / (n_ - k_)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    r2 = 1 - e @ e / np.sum((dat[y] - dat[y].mean()) ** 2)
    return b, se, r2


REG = {}
for spec in (["eci"], ["eci", "lgdp"], ["eci", "lgdp", "Gamma", "Neff"]):
    b, se, r2 = ols("gini", spec)
    REG[tuple(spec)] = (b, se, r2)
    print(f"gini ~ {' + '.join(spec)} (n = {len(dat)}): " +
          ", ".join(f"{nm} {bb:+.2f} ({ss:.2f})" for nm, bb, ss in zip(["const"] + spec, b, se)) + f"; R2 = {r2:.3f}")
b1 = REG[("eci", "lgdp")][0][1]; b2 = REG[("eci", "lgdp", "Gamma", "Neff")][0][1]
print(f"ECI coefficient (controlling log GDP pc) {b1:+.2f} -> {b2:+.2f} after adding Gamma and Neff")
# (d) residuals from gini ~ eci + lgdp
b, se, r2 = REG[("eci", "lgdp")]
X = np.column_stack([np.ones(len(dat)), dat.eci, dat.lgdp])
res = dat.gini - X @ b
print("residuals (Gini minus predicted):", res.round(1).to_dict())
G16_6 = dict(n=len(dat), b1=b1, b2=b2)

# =====================================================================
print("\n== Exercise 16.7: Bolivia's structural position (inputs) ==")
gs = summ
print(f"pi_BOL 2019 = {gs.loc[2019,'pi_BOL']:.5f} (1/N = {1/gs.loc[2019,'N']:.5f}, ratio {gs.loc[2019,'pi_BOL']*gs.loc[2019,'N']:.3f}); "
      f"ECI 2019 = {e19.loc['BOL','eci_hs92']} (rank {int(e19.loc['BOL','eci_rank_hs92'])} of {int(rk[(rk.year==2019)&rk.in_rankings].shape[0])})")
print(f"domestic Gamma 2019 (book12) = {bsum[(bsum.year==2019)&(bsum.base==2017)&(bsum.aggregation=='book12')&(bsum.flows=='domestic')&(bsum.normalization=='input_shares')].Gamma.iloc[0]:.3f}")
fac = pd.read_csv(os.path.join(PROC, "bolivia_facts.csv")).set_index("fact_id").value
print({k: fac[k] for k in fac.index if "gini" in k.lower()} )
# 20-year ECI changes across ranked economies, 1999 -> 2019 (to calibrate what a plan can expect)
e99 = rk[(rk.year == 1999) & (rk.in_rankings)].set_index("country_iso3_code").eci_hs92
d20 = (e19.eci_hs92 - e99).dropna()
q = d20.quantile([0.5, 0.75, 0.9, 0.95])
print(f"20-year ECI change 1999-2019 over {len(d20)} economies: median {q[0.5]:+.2f}, p75 {q[0.75]:+.2f}, "
      f"p90 {q[0.9]:+.2f}, p95 {q[0.95]:+.2f}; max {d20.max():+.2f} ({d20.idxmax()}); BOL {d20['BOL']:+.2f}; "
      f"share of economies gaining >= 0.5: {(d20 >= 0.5).mean():.3f}")
gdp19 = wdi[(wdi.indicator == "NY.GDP.MKTP.CD") & (wdi.year == 2019) & (wdi.country_iso3 == "BOL")].value.iloc[0]
print(f"BOL GDP 2019 = {gdp19/1e9:.1f} bn US$; 1% of GDP for 20 years = {0.01*gdp19*20/1e9:.1f} bn US$")
G16_7 = dict(q=q, d20B=d20["BOL"])
print("\nall Chapter 16 solution checks passed")
