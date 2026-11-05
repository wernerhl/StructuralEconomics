"""Solutions manual, Chapter 23 (Tracking Structural Change Over Time).  Computes and
prints every number used in solutions/sol23.tex, with asserts.

Sourced inputs (data/processed, retrieved 2026-11-05): INE domestic flow matrices
(bolivia_Z/), bolivia_io_summary.csv, bolivia_io_pi.csv, bolivia_io_sector_facts.csv,
oecd_io_summary.csv (OECD ICIO 2021 and 2023 editions), eci_facts.csv (Growth Lab).
Run from the repo root:  python3 code/solutions/sol23.py   (numpy + pandas).
"""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=150)
PROC = os.path.join("data", "processed")


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def neff(p):
    return 1 / np.sum(p ** 2)


def tau1(W):
    n = len(W)
    return max(0.5 * np.abs(W[i] - W[j]).sum() for i in range(n) for j in range(n))


def rownorm(K):
    return K / K.sum(1, keepdims=True)


def read_Z(base, year):
    df = pd.read_csv(os.path.join(PROC, "bolivia_Z", f"Z_{base}_{year}_book12.csv"), index_col=0)
    return df.values.astype(float), list(df.columns)


def wp(Z, labels, keep_labels=None):
    if keep_labels is None:
        keep = np.ones(len(labels), bool)
        while True:
            Zs = Z[np.ix_(keep, keep)]
            bad = (Zs.sum(0) <= 0) | (Zs.sum(1) <= 0)
            if not bad.any():
                break
            keep[np.where(keep)[0][bad]] = False
    else:
        keep = np.array([l in keep_labels for l in labels])
    Zk = Z[np.ix_(keep, keep)]
    return rownorm(Zk.T), [l for l, k in zip(labels, keep) if k]


SUM = pd.read_csv(os.path.join(PROC, "bolivia_io_summary.csv"))
SUM = SUM[(SUM.aggregation == "book12") & (SUM.flows == "domestic")]
NORMS = ["input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag"]
def series(base, norm, col="Gamma"):
    x = SUM[(SUM.base == base) & (SUM.normalization == norm)].set_index("year")[col].sort_index()
    return x

YRS90 = sorted(SUM[SUM.base == 1990].year.unique()); YRS17 = sorted(SUM[SUM.base == 2017].year.unique())

# =====================================================================
print("=" * 70); print("Exercise 23.1  Bolivia's trajectory 2000-2019"); print("=" * 70)
TR = {}
for base, yrs in ((1990, [y for y in YRS90 if y >= 2000]), (2017, [2017, 2018, 2019])):
    for norm in NORMS:
        g = series(base, norm).loc[yrs]; G = series(base, norm, "power_gini").loc[yrs]
        TR[(base, norm)] = (g, G)
print("input shares, diagonal kept (headline):")
for base in (1990, 2017):
    g, G = TR[(base, "input_shares")]
    for y in g.index:
        print(f"  {base} {y}: Gamma = {g[y]:.3f}, G(pi) = {G[y]:.3f}")
CHG = {}
for norm in NORMS:
    for base, y0, y1 in ((1990, 2000, 2013), (1990, 2000, 2014), (2017, 2017, 2019)):
        g, G = TR[(base, norm)]
        CHG[(norm, base, y0, y1)] = (g[y1] - g[y0], G[y1] - G[y0])
for norm in NORMS:
    s = "; ".join(f"{y0}-{y1}: dGamma {CHG[(norm,b,y0,y1)][0]:+.3f}, dG {CHG[(norm,b,y0,y1)][1]:+.3f}"
                  for (n_, b, y0, y1) in CHG if n_ == norm)
    print(f"  {norm:24s} {s}")
signs = {k: (np.sign(v[0]), np.sign(v[1])) for k, v in CHG.items()}
print("sign of dGamma 2000-2013 by normalisation:", {n: int(signs[(n, 1990, 2000, 2013)][0]) for n in NORMS})
print("sign of dG(pi) 2000-2013 by normalisation:", {n: int(signs[(n, 1990, 2000, 2013)][1]) for n in NORMS})
print("sign of dGamma 2017-2019 by normalisation:", {n: int(signs[(n, 2017, 2017, 2019)][0]) for n in NORMS})
print("sign of dG(pi) 2017-2019 by normalisation:", {n: int(signs[(n, 2017, 2017, 2019)][1]) for n in NORMS})
# year-to-year co-movement of Gamma and G(pi), input shares, base 1990, 2000-2014
g, G = TR[(1990, "input_shares")]
print(f"corr(Gamma_t, G(pi)_t), base 1990 2000-2014 = {np.corrcoef(g.values, G.values)[0,1]:+.3f}")
CORR_GG = np.corrcoef(g.values, G.values)[0, 1]
# TikZ coordinates for the (G(pi), Gamma) plot: x = 20 (G - 0.30), y = 20 (Gamma - 0.30)
for base in (1990, 2017):
    g, G = TR[(base, "input_shares")]
    print(f"tikz base {base}:", " ".join(f"({20*(G[y]-0.30):.2f},{20*(g[y]-0.30):.2f})" for y in g.index))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.2  Structural break detection"); print("=" * 70)
g90 = series(1990, "input_shares"); g17 = series(2017, "input_shares")
dG90 = {y: g90[y] - g90[y - 1] for y in YRS90 if y - 1 in g90.index}
dG17 = {y: g17[y] - g17[y - 1] for y in YRS17 if y - 1 in g17.index}
quiet = [dG90[y] for y in dG90 if y <= 1997]
sd = float(np.std(quiet, ddof=1))
print(f"sigma_hat from {len(quiet)} quiet changes 1989-1997 = {sd:.4f}; thresholds 2 sigma = {2*sd:.4f}, 3 sigma = {3*sd:.4f}")
flag2_90 = sorted(y for y in dG90 if abs(dG90[y]) > 2 * sd)
flag2_17 = sorted(y for y in dG17 if abs(dG17[y]) > 2 * sd)
print("base 1990 dGamma:", {y: round(v, 3) for y, v in dG90.items()})
print("base 2017 dGamma:", {y: round(v, 3) for y, v in dG17.items()})
print(f"flagged at k = 2: base 1990 {flag2_90}; base 2017 {flag2_17}")
# bootstrap: multiplicative lognormal noise on every Z cell (CV 10%), independent across years
rng = np.random.default_rng(23)
CV = 0.10; B = 300
sig_ln = np.sqrt(np.log(1 + CV ** 2))
def boot_gamma(base, year):
    Z, lab = read_Z(base, year)
    out = []
    for _ in range(B):
        Zb = Z * rng.lognormal(-sig_ln ** 2 / 2, sig_ln, Z.shape)
        out.append(gap(wp(Zb, lab)[0]))
    return np.array(out)
BOOT = {(1990, y): boot_gamma(1990, y) for y in YRS90}
BOOT.update({(2017, y): boot_gamma(2017, y) for y in YRS17})
band = {}
for base, dG in ((1990, dG90), (2017, dG17)):
    for y in dG:
        d = BOOT[(base, y)] - BOOT[(base, y - 1)]
        band[(base, y)] = np.quantile(d, [0.025, 0.975])
sdb = np.mean([BOOT[k].std() for k in BOOT])
print(f"bootstrap (CV {CV:.0%} per cell, B = {B}): mean s.d. of Gamma_t = {sdb:.4f}")
BREAK = []
for base, dG in ((1990, dG90), (2017, dG17)):
    for y in sorted(dG):
        lo_, hi_ = band[(base, y)]
        excl0 = (lo_ > 0) or (hi_ < 0)
        if abs(dG[y]) > 2 * sd:
            print(f"  {base} {y}: dGamma = {dG[y]:+.3f}; bootstrap 95% band of dGamma [{lo_:+.3f}, {hi_:+.3f}]; excludes 0: {excl0}")
        if abs(dG[y]) > 2 * sd and excl0:
            BREAK.append((base, y))
print("flagged by both criteria:", BREAK)
# tau_1 and realised-product gap on a common 10-sector node set (base 1990)
common = None
for y in YRS90:
    Z, lab = read_Z(1990, y)
    _, s = wp(Z, lab)
    common = set(s) if common is None else common & set(s)
common = sorted(common)
print(f"common node set, base 1990: {len(common)} sectors (drops {sorted(set(read_Z(1990,1988)[1]) - set(common))})")
W90 = {}
for y in YRS90:
    Z, lab = read_Z(1990, y)
    W90[y], lab_common = wp(Z, lab, keep_labels=common)
T1 = {y: tau1(W90[y]) for y in YRS90}
Gc = {y: gap(W90[y]) for y in YRS90}
G5 = {}
for y in YRS90:
    win = list(range(y, y + 5))
    if all(t in W90 for t in win):
        P = np.eye(len(common))
        for t in win:
            P = W90[t] @ P
        G5[y] = 1 - abs(eig_sorted(P)[1]) ** (1 / 5)
print("year: Gamma (common set) | tau_1 | Gamma^(5) starting that year")
for y in YRS90:
    print(f"  {y}: {Gc[y]:.3f} | {T1[y]:.3f} | {G5[y]:.3f}" if y in G5 else f"  {y}: {Gc[y]:.3f} | {T1[y]:.3f} |   -")
print(f"tau_1 range {min(T1.values()):.3f}-{max(T1.values()):.3f}; 1 - tau_1 range {1-max(T1.values()):.3f}-{1-min(T1.values()):.3f}")
assert all(Gc[y] >= 1 - T1[y] - 1e-12 for y in YRS90)       # Prop. 23.dobrushin(ii)
G5v = np.array(list(G5.values()))
print(f"Gamma^(5) range {G5v.min():.3f}-{G5v.max():.3f}; mean per-period Gamma over the same windows "
      f"{np.mean([np.mean([Gc[t] for t in range(y, y+5)]) for y in G5]):.3f}")
print(f"Gamma^(5): 1990 window {G5[1990]:.3f}, 1996 window (incl. 1999-2000) {G5[1996]:.3f}, 2001 window {G5[2001]:.3f}, 2005 window {G5[2005]:.3f}")
# does Gamma^(5) flag the 1999/2000 moves?  windows containing 1999 vs not
w99 = [G5[y] for y in G5 if y <= 1999 <= y + 4]; wno = [G5[y] for y in G5 if not (y <= 1999 <= y + 4)]
print(f"mean Gamma^(5) over windows containing 1999: {np.mean(w99):.3f}; others: {np.mean(wno):.3f}")
print(f"Gamma (common set) 1999 = {Gc[1999]:.3f}, 2000 = {Gc[2000]:.3f}, 2014 = {Gc[2014]:.3f}")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.3  Structural velocity"); print("=" * 70)
OS = pd.read_csv(os.path.join(PROC, "oecd_io_summary.csv"))
OS = OS[(OS.aggregation == "book12_isic") & (OS.normalization == "input_shares") & (OS.flows == "domestic")]
ECI = pd.read_csv(os.path.join(PROC, "eci_facts.csv"))
def eci(c, y):
    v = ECI[(ECI.country == c) & (ECI.year == y)]["eci_hs92"]
    return float(v.iloc[0])
def vel(ser, y0, y1):
    s = ser.loc[[y for y in ser.index if y0 <= y <= y1]]
    ep = (s.loc[y1] - s.loc[y0]) / (y1 - y0) * 10
    ols = np.polyfit(s.index.values.astype(float), s.values, 1)[0] * 10
    return ep, ols, len(s)
VEL = {}
for c in ("CHL", "PER"):
    for ed, y0, y1 in (("ICIO2021", 1995, 2014), ("ICIO2023", 2011, 2020)):
        ser = OS[(OS.country == c) & (OS.edition == ed)].set_index("year").Gamma.sort_index()
        ep, ols, n_ = vel(ser, y0, y1)
        de = eci(c, y1) - eci(c, y0)
        VEL[(c, y0, y1)] = (ep, ols, de, ser.loc[y0], ser.loc[y1])
        print(f"{c} {ed} {y0}-{y1}: Gamma {ser.loc[y0]:.3f} -> {ser.loc[y1]:.3f}; v endpoint {ep:+.3f}, OLS {ols:+.3f} per decade "
              f"(n = {n_}); dECI = {de:+.3f}")
for lab_, base, y0, y1 in (("1995-2014", 1990, 1995, 2014), ("1995-2013 (excl. 2014)", 1990, 1995, 2013),
                           ("2011-2014", 1990, 2011, 2014), ("2017-2020", 2017, 2017, 2020), ("2017-2024", 2017, 2017, 2024)):
    ser = series(base, "input_shares")
    ep, ols, n_ = vel(ser, y0, y1)
    de = eci("BOL", y1) - eci("BOL", y0)
    VEL[("BOL", y0, y1)] = (ep, ols, de, ser.loc[y0], ser.loc[y1])
    print(f"BOL base {base} {lab_}: Gamma {ser.loc[y0]:.3f} -> {ser.loc[y1]:.3f}; v endpoint {ep:+.3f}, OLS {ols:+.3f} (n = {n_}); dECI = {de:+.3f}")
for c in ("CHL", "PER"):
    a, b_ = VEL[(c, 1995, 2014)], VEL[(c, 2011, 2020)]
    print(f"{c}: sign kept across windows? endpoint {np.sign(a[0]) == np.sign(b_[0])}, OLS {np.sign(a[1]) == np.sign(b_[1])}")
gv = pd.read_csv(os.path.join(PROC, "gamma_velocity.csv"))
for c, ed in (("CHL", "ICIO2021"), ("PER", "ICIO2021"), ("CHL", "ICIO2023"), ("PER", "ICIO2023")):
    r = gv[(gv.country == c) & (gv.vintage == ed)].iloc[0]
    key = (c, 1995, 2014) if ed == "ICIO2021" else (c, 2011, 2020)
    assert abs(r.v_endpoint_per_decade - VEL[key][0]) < 1e-9, (c, ed)
print("endpoint velocities match data/processed/gamma_velocity.csv")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.4  Structural-cyclical decomposition of pi_gas (base 1990)"); print("=" * 70)
PI = pd.read_csv(os.path.join(PROC, "bolivia_io_pi.csv"))
ph = PI[(PI.base == 1990) & (PI.aggregation == "book12") & (PI.flows == "domestic") & (PI.normalization == "input_shares")
        & (PI.sector == "Hydrocarbons")].set_index("year").pi.sort_index()
print("pi_gas:", {y: round(v, 3) for y, v in ph.items()})
def hp_missing(years_obs, yobs, lam):
    """HP trend on the full annual grid; the fit term uses observed years only."""
    grid = np.arange(min(years_obs), max(years_obs) + 1)
    n = len(grid)
    S = np.zeros((n, n)); yv = np.zeros(n)
    for y, v in zip(years_obs, yobs):
        k = y - grid[0]; S[k, k] = 1; yv[k] = v
    D2 = np.zeros((n - 2, n))
    for i in range(n - 2):
        D2[i, i:i + 3] = [1, -2, 1]
    tau = np.linalg.solve(S + lam * D2.T @ D2, S @ yv)
    return dict(zip(grid, tau))
HP = {}
for drop14 in (False, True):
    yrs = [y for y in ph.index if not (drop14 and y == 2014)]
    for lam in (6.25, 100, 400):
        tr = hp_missing(yrs, ph.loc[yrs].values, lam)
        share = (tr[2013] - tr[2003]) / (ph[2013] - ph[2003])
        HP[(drop14, lam)] = (tr[2003], tr[2013], share, tr)
        print(f"lambda = {lam:7.2f}, {'without' if drop14 else 'with   '} 2014: trend 2003 = {tr[2003]:.4f}, 2013 = {tr[2013]:.4f}; "
              f"structural share of the 2003-2013 rise = {share:.3f}; trend 2010 (interpolated) = {tr[2010]:.4f}")
print(f"observed rise 2003-2013 = {ph[2013]-ph[2003]:+.4f} ({ph[2003]:.4f} -> {ph[2013]:.4f})")
# HP with missing obs equals ordinary HP when nothing is missing: check on a gap-free window
yy = list(range(1988, 2010)); tr1 = hp_missing(yy, ph.loc[yy].values, 100)
n = len(yy); D2 = np.zeros((n - 2, n))
for i in range(n - 2):
    D2[i, i:i + 3] = [1, -2, 1]
tr2 = np.linalg.solve(np.eye(n) + 100 * D2.T @ D2, ph.loc[yy].values)
assert np.allclose([tr1[y] for y in yy], tr2)
# composition check: number of sectors with a defined row each year
NS = SUM[(SUM.base == 1990) & (SUM.normalization == "input_shares")].set_index("year").N
print("N by year (base 1990):", {y: int(v) for y, v in NS.items()})
# pi_gas on the common 10-sector node set
ph10 = {}
for y in YRS90:
    Wy = W90[y]; ph10[y] = perron(Wy)[lab_common.index("Hydrocarbons")]
print("pi_gas on the common 10-sector set:", {y: round(v, 3) for y, v in ph10.items()})
tr10 = hp_missing(list(ph10), np.array(list(ph10.values())), 100)
sh10 = (tr10[2013] - tr10[2003]) / (ph10[2013] - ph10[2003])
print(f"common-set series, lambda = 100: structural share 2003-2013 = {sh10:.3f}")
p17 = PI[(PI.base == 2017) & (PI.aggregation == "book12") & (PI.flows == "domestic") & (PI.normalization == "input_shares")
         & (PI.sector == "Hydrocarbons")].set_index("year").pi.sort_index()
print(f"base 2017 pi_gas: {dict((y, round(v,3)) for y, v in p17.items())}; 2014 (base 1990) {ph[2014]:.3f} vs 2017 (base 2017) {p17[2017]:.3f}")
# HP cutoff periods (gain 1/2)
for lam in (6.25, 100, 400):
    w = np.arccos(1 - 1 / (2 * np.sqrt(lam)))
    print(f"HP gain-1/2 cutoff period, lambda = {lam}: {2*np.pi/w:.1f} years")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.5  Commodity boom counterfactual (2019 W^P)"); print("=" * 70)
Z19, lab19 = read_Z(2017, 2019)
SF = pd.read_csv(os.path.join(PROC, "bolivia_io_sector_facts.csv"))
SF19 = SF[(SF.base == 2017) & (SF.year == 2019)].set_index("book12")
iH = lab19.index("Hydrocarbons")
xH = SF19.loc["Hydrocarbons", "output"]; vaH = SF19.loc["Hydrocarbons", "va"]
x12 = np.array([SF19.loc[s, "output"] for s in lab19])
def diag(Z):
    W, s = wp(Z, lab19)
    p = perron(W)
    A = Z / x12[None, :]
    L = np.linalg.inv(np.eye(12) - A)
    return dict(G=gap(W), pi=p[s.index("Hydrocarbons")], int=Z[:, iH].sum() / xH, BL=L.sum(0)[iH], neff=neff(p))
base = diag(Z19)
print(f"2019: domestic input intensity {base['int']:.3f}, pi_hydro {base['pi']:.4f}, Gamma {base['G']:.4f}, BL_hydro {base['BL']:.3f}")
print(f"scale: hydrocarbons output {xH:.0f} M Bs, VA {vaH:.0f} M Bs; 30% of VA = {0.3*vaH:.0f} M Bs")
CF = {}
for target in (0.30, 0.381):
    extra = (target - base["int"]) * xH
    Za = Z19.copy(); Za[:, iH] *= target / base["int"]                       # (i) proportional scale-up
    Zb = Z19.copy()                                                          # (ii) new purchases from agr + mfg
    for s, sh in (("Agriculture", 0.5), ("Manufacturing", 0.5)):
        Zb[lab19.index(s), iH] += sh * extra
    CF[("i", target)] = diag(Za); CF[("ii", target)] = diag(Zb)
    print(f"target intensity {target}: extra purchases {extra:.0f} M Bs/yr ({extra/(0.3*vaH):.2f} x 30% of VA)")
    for k in ("i", "ii"):
        d = CF[(k, target)]
        print(f"  ({k}) intensity {d['int']:.3f}, pi_hydro {d['pi']:.4f} ({d['pi']-base['pi']:+.4f}), Gamma {d['G']:.4f} "
              f"({d['G']-base['G']:+.4f}), BL_hydro {d['BL']:.3f}")
    assert abs(CF[("i", target)]["G"] - base["G"]) < 1e-12 and abs(CF[("i", target)]["pi"] - base["pi"]) < 1e-12
# (iii) forward linkage: domestic sectors buy gas inputs worth 10% of their imported intermediates
Zc = Z19.copy()
for k, s in enumerate(lab19):
    if s != "Hydrocarbons" and np.isfinite(SF19.loc[s, "imp_int_inputs"]):
        Zc[iH, k] += 0.10 * SF19.loc[s, "imp_int_inputs"]
CF["iii"] = diag(Zc)
d = CF["iii"]
print(f"  (iii) forward linkage: intensity {d['int']:.3f}, pi_hydro {d['pi']:.4f} ({d['pi']-base['pi']:+.4f}), Gamma {d['G']:.4f} ({d['G']-base['G']:+.4f})")
# first-order sensitivity of lambda_2 to the hydrocarbons row (perturbation lemma, eigenvalue part)
W19, s19 = wp(Z19, lab19)
wr, vr = np.linalg.eig(W19); wl, vl = np.linalg.eig(W19.T)
k2 = np.argsort(-abs(wr))[1]; u2 = vl[:, np.argmin(abs(wl - wr[k2]))].real; v2 = vr[:, k2].real
u2 = u2 / (u2 @ v2)
jH = s19.index("Hydrocarbons")
print(f"|u_2| at hydrocarbons = {abs(u2[jH]):.4f} vs max |u_2| = {abs(u2).max():.4f}; u_2 weight share of hydrocarbons "
      f"= {abs(u2[jH])/abs(u2).sum():.4f}")
U2H = abs(u2[jH]) / abs(u2).sum()

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.6  SCI for Latin American ICIO economies"); print("=" * 70)
LA = ["ARG", "BRA", "CHL", "COL", "CRI", "MEX", "PER"]
SCI = {}
for ed, y0, y1 in (("ICIO2021", 1995, 2014), ("ICIO2023", 2011, 2020)):
    for c in LA:
        ser = OS[(OS.country == c) & (OS.edition == ed)].set_index("year").Gamma.sort_index()
        if len(ser) == 0:
            print(f"{ed} {c}: not available")
            continue
        g0, g1 = ser.loc[y0], ser.loc[y1]
        sci = (g1 - g0) / (0.50 - g0) if g0 < 0.50 else np.nan
        try:
            de = eci(c, y1) - eci(c, y0)
        except IndexError:
            de = np.nan
        SCI[(ed, c)] = (g0, g1, sci, de)
for ed in ("ICIO2021", "ICIO2023"):
    rows = [(c, *SCI[(ed, c)]) for c in LA if (ed, c) in SCI]
    nund = sum(np.isnan(r[3]) for r in rows)
    print(f"{ed}: {len(rows)} economies, SCI undefined for {nund}")
    for c, g0, g1, sci, de in sorted(rows, key=lambda r: -np.nan_to_num(r[3], nan=-1e9)):
        print(f"  {c}: Gamma {g0:.3f} -> {g1:.3f} (dGamma {g1-g0:+.3f}); SCI = {'undefined' if np.isnan(sci) else f'{sci:.2f}'}; "
              f"dECI = {'n/a' if np.isnan(de) else f'{de:+.3f}'}")
    ok = [(g1 - g0, de) for c, g0, g1, sci, de in rows if not np.isnan(de)]
    a = np.array(ok)
    rk = lambda v: np.argsort(np.argsort(v))
    rho = np.corrcoef(rk(a[:, 0]), rk(a[:, 1]))[0, 1]
    print(f"  Spearman rank correlation of dGamma and dECI over {len(ok)} economies: {rho:+.3f}")
    SCI[(ed, "rho")] = (rho, len(ok), nund, len(rows))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 23.7  Trajectory projection of Bolivia's ECI"); print("=" * 70)
b95, b19, b24 = eci("BOL", 1995), eci("BOL", 2019), eci("BOL", 2024)
k95, k19 = eci("KOR", 1995), eci("KOR", 2019)
slopes = {"(a) 1995-2019 trend": (b19 - b95) / 24, "(b) 2019-2024 trend": (b24 - b19) / 5, "(c) Korea's 1995-2019 pace": (k19 - k95) / 24}
print(f"BOL ECI 1995 {b95:+.3f}, 2019 {b19:+.3f}, 2024 {b24:+.3f}; KOR 1995 {k95:+.3f}, 2019 {k19:+.3f}")
print(f"dECI BOL 1995-2019 = {b19-b95:+.3f}; BOL 2019-2024 = {b24-b19:+.3f}; KOR 1995-2019 = {k19-k95:+.3f}")
PROJ = {}
for k, sl in slopes.items():
    path = {y: b24 + sl * (y - 2024) for y in range(2025, 2061)}
    hit = 2024 - b24 / sl if sl > 0 else np.inf
    PROJ[k] = (sl, path, hit)
    print(f"{k}: slope {sl:+.4f}/yr; 2030 {path[2030]:+.3f}, 2040 {path[2040]:+.3f}, 2060 {path[2060]:+.3f}; "
          f"reaches 0 in {'never' if not np.isfinite(hit) else f'{hit:.1f}'}")
print("tikz (x = (year-2020)/5, y = 2 ECI):")
for k, (sl, path, hit) in PROJ.items():
    pts = [(2024, b24)] + [(y, path[y]) for y in (2030, 2040, 2050, 2060)]
    print(f"  {k}: " + " ".join(f"({(y-2020)/5:.2f},{2*v:.2f})" for y, v in pts))
hist = [(y, eci("BOL", y)) for y in (1995, 2000, 2005, 2010, 2015, 2019, 2024)]
print("  history: " + " ".join(f"({(y-2020)/5:.2f},{2*v:.2f})" for y, v in hist))
nC = ECI[ECI.year == 2019].shape[0]
print(f"ECI rank of Bolivia: 2019 {int(ECI[(ECI.country=='BOL')&(ECI.year==2019)].eci_rank_hs92.iloc[0])}, "
      f"2024 {int(ECI[(ECI.country=='BOL')&(ECI.year==2024)].eci_rank_hs92.iloc[0])}")
print("\nall Chapter 23 solution numbers computed; asserts passed")
