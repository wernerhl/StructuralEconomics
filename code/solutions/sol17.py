"""Solutions manual, Chapter 17 (Domestic Structural Diagnostics): computes and prints
every number used in solutions/sol17.tex, with asserts on the key results.
Run from the repo root:  python3 code/solutions/sol17.py
numpy + pandas (pandas reads data/processed, and the tracked Growth Lab ranking file
data/raw/atlas/growth_proj_eci_rankings.csv for one ECI value)."""
import math
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=120)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROC = os.path.join(ROOT, "data", "processed")


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    p = np.real(vl[:, np.argmin(np.abs(ev - 1))])
    return p / p.sum()


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def read_Z(fname):
    Z = pd.read_csv(os.path.join(PROC, "bolivia_Z", fname), index_col=0)
    assert list(Z.index) == list(Z.columns)
    return Z.values.astype(float), list(Z.index)


def wp(Z, labels, sales=False, zero_diag=False):
    """input shares: row s = Z[:, s]/sum (W^P = D^-1 A); sales shares: row s = Z[s, :]/sum.
    Sectors with no intermediate purchases or sales are dropped iteratively."""
    Z = Z.copy()
    if zero_diag:
        np.fill_diagonal(Z, 0)
    keep = np.ones(len(labels), bool)
    while True:
        Zs = Z[np.ix_(keep, keep)]
        bad = (Zs.sum(0) <= 0) | (Zs.sum(1) <= 0)
        if not bad.any():
            break
        keep[np.where(keep)[0][bad]] = False
    Zk = Z[np.ix_(keep, keep)]
    K = Zk if sales else Zk.T
    return K / K.sum(1, keepdims=True), [l for l, k in zip(labels, keep) if k]


SUM = pd.read_csv(os.path.join(PROC, "bolivia_io_summary.csv"))


def summ(year, agg="book12", norm="input_shares", flows="domestic"):
    m = SUM[(SUM.year == year) & (SUM.aggregation == agg) & (SUM.normalization == norm) & (SUM.flows == flows)]
    assert len(m) == 1
    return m.iloc[0]


# =====================================================================
print("== Exercise 17.1: Bolivia's Gamma from INE 2019 (base 2017), domestic, basic prices ==")
Z12, lab12 = read_Z("Z_2017_2019_book12.csv")
W, s11 = wp(Z12, lab12)
ev = eig_sorted(W)
p = perron(W)
print(f"book12 input shares: N = {len(s11)} (dropped: {set(lab12) - set(s11)})")
print("eigenvalues (sorted by modulus):", np.round(ev, 4))
print("pi:", dict(zip(s11, np.round(p, 3))))
print(f"Gamma = {gap(W):.4f}; G(pi) = {gini(p):.4f}; Neff = {1/np.sum(p**2):.3f}; lambda_2 = {ev[1]:.4f}, |lambda_3| = {abs(ev[2]):.4f}")
assert round(gap(W), 3) == 0.467 and round(gini(p), 3) == 0.445
V17 = {}
for agg, fn in (("book12", "Z_2017_2019_book12.csv"), ("book12_isic", "Z_2017_2019_book12_isic.csv")):
    Zx, lx = read_Z(fn)
    for norm, sales, zd in (("input_shares", False, False), ("input_shares_zero_diag", False, True),
                            ("sales_shares", True, False), ("sales_shares_zero_diag", True, True)):
        Wx, sx = wp(Zx, lx, sales, zd)
        px = perron(Wx)
        g_, G_ = gap(Wx), gini(px)
        V17[(agg, norm)] = (g_, G_, 1 / np.sum(px ** 2), len(sx))
        csv = summ(2019, agg, norm)
        assert abs(g_ - csv.Gamma) < 1e-9 and abs(G_ - csv.power_gini) < 1e-9
        print(f"  {agg:12s} {norm:24s}: N = {len(sx)}, Gamma = {g_:.3f}, G(pi) = {G_:.3f}, Neff = {1/np.sum(px**2):.2f}, "
              f"pi_max = {sx[int(np.argmax(px))]} {px.max():.3f}")
for norm in ("input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag"):
    c = summ(2019, "detail39", norm)
    V17[("detail39", norm)] = (c.Gamma, c.power_gini, c.Neff, c.N)
    print(f"  detail39     {norm:24s}: N = {int(c.N)}, Gamma = {c.Gamma:.3f}, G(pi) = {c.power_gini:.3f} (from CSV)")
gs = [v[0] for k, v in V17.items() if k[0] in ("book12", "detail39")]
Gs = [v[1] for k, v in V17.items() if k[0] in ("book12", "detail39")]
print(f"range over book12 + detail39 x 4 normalizations: Gamma {min(gs):.3f}-{max(gs):.3f}; G(pi) {min(Gs):.3f}-{max(Gs):.3f}")
assert round(min(gs), 2) == 0.33 and round(max(gs), 2) == 0.69
# the sales-share matrix: who are the hubs
Ws, ss = wp(Z12, lab12, sales=True)
ps = perron(Ws)
print("sales-share pi:", dict(zip(ss, np.round(ps, 3))))
# an 8-sector aggregation built here from the 12-sector table
MAP8 = {"Hydrocarbons": "Extractive", "Mining": "Extractive", "Manufacturing": "Manufacturing",
        "Construction": "Construction", "Agriculture": "Agriculture", "Commerce": "Trade & transport",
        "Transport": "Trade & transport", "Financial services": "Finance & other services",
        "Other services": "Finance & other services", "Public administration": "Government",
        "Personal services": "Personal & domestic", "Domestic services": "Personal & domestic"}
g8 = sorted(set(MAP8.values()))
S = np.array([[1.0 if MAP8[l] == g else 0.0 for g in g8] for l in lab12])
Z8 = S.T @ Z12 @ S
W8, s8 = wp(Z8, g8)
p8 = perron(W8)
print(f"8-sector aggregation: N = {len(s8)}, Gamma = {gap(W8):.3f}, G(pi) = {gini(p8):.3f}, pi = {dict(zip(s8, np.round(p8, 3)))}")
W8s, _ = wp(Z8, g8, sales=True)
print(f"8-sector, sales shares: Gamma = {gap(W8s):.3f}, G(pi) = {gini(perron(W8s)):.3f}")
A8 = (gap(W8), gini(p8), gap(W8s), gini(perron(W8s)))

# =====================================================================
print("\n== Exercise 17.2: between-block variance (illustrative inputs of Example 17.5) ==")
w3 = np.array([0.08, 0.35, 0.57]); y3 = np.array([18500, 5200, 1800]); VT = 32e6
m = w3 @ y3
VB3 = w3 @ (y3 - m) ** 2
w2 = np.array([w3[0] + w3[1], w3[2]]); y2 = np.array([(w3[:2] @ y3[:2]) / w3[:2].sum(), y3[2]])
VB2 = w2 @ (y2 - m) ** 2
print(f"mean = {m:.0f}; 3 blocks: V_B = {VB3/1e6:.2f}M, sigma_B = {VB3/VT:.3f}")
print(f"2 blocks (merge 1 and 2): formal mean = {y2[0]:.1f}; V_B = {VB2/1e6:.2f}M, sigma_B = {VB2/VT:.3f}")
print(f"  difference = within-coarse-block variance of fine means = "
      f"{w3[:2].sum() * np.sum(w3[:2]/w3[:2].sum() * (y3[:2]-y2[0])**2)/1e6:.2f}M")
assert VB3 >= VB2 and abs(VB3 - VB2 - w3[:2].sum() * np.sum(w3[:2] / w3[:2].sum() * (y3[:2] - y2[0]) ** 2)) < 1e-3
rng = np.random.default_rng(172)
for _ in range(500):
    n = 80; y = rng.lognormal(0, 0.8, n); lab = rng.integers(0, 3, n)
    def VB(l):
        return sum((l == b).mean() * (y[l == b].mean() - y.mean()) ** 2 for b in np.unique(l))
    fine = lab * 3 + rng.integers(0, 3, n)
    assert VB(fine) >= VB(lab) - 1e-12
    assert abs(VB(np.arange(n)) - y.var()) < 1e-9          # singletons: sigma_B = 1
print("refinement never lowers V_B: 500 random checks; singleton partition gives sigma_B = 1")

# =====================================================================
print("\n== Exercise 17.3: power-employment ratio (INE 2019 + ILOSTAT 2019) ==")
EMP_MAP = {"Hydrocarbons and mining": "B", "Public administration": "O", "Manufacturing": "C",
           "Construction": "F", "Commerce": "G", "Transport": "H", "Financial services": "K",
           "Agriculture": "A", "Other services": "DEJLMNPQ", "Personal services": "IRS",
           "Domestic services": "T"}
ilo = pd.read_csv(os.path.join(PROC, "ilostat_bolivia_emp_eco_sex.csv"))
ilo = ilo[ilo.year == 2019]
I = {(r.sex, r.eco): r.value_thousands for r in ilo.itertuples()}
tot = I[("SEX_T", "ECO_ISIC4_TOTAL")]
emp_k = {s: sum(I[("SEX_T", "ECO_ISIC4_" + c)] for c in cs) for s, cs in EMP_MAP.items()}
emp = {s: v / tot for s, v in emp_k.items()}
pd_ = dict(zip(s11, p))
pi_e = {s: (pd_["Hydrocarbons"] + pd_["Mining"]) if s == "Hydrocarbons and mining" else pd_.get(s, np.nan) for s in EMP_MAP}
ratio = {s: pi_e[s] / emp[s] for s in EMP_MAP}
sf = pd.read_csv(os.path.join(PROC, "bolivia_io_sector_facts.csv"))
sf = sf[(sf.year == 2019) & (sf.base == 2017)].set_index("book12")
va = {s: (sf.loc["Hydrocarbons", "va"] + sf.loc["Mining", "va"]) if s == "Hydrocarbons and mining" else sf.loc[s, "va"]
      for s in EMP_MAP}
vapw = {s: va[s] / emp_k[s] for s in EMP_MAP}          # million Bs per thousand workers = thousand Bs per worker
print(f"{'sector':24s} emp share   pi   pi/share  VA/worker (thousand Bs)")
for s in EMP_MAP:
    print(f"  {s:24s} {emp[s]:.3f}  {pi_e[s]:.3f}  {ratio[s]:6.2f}  {vapw[s]:8.1f}")
sec = [s for s in EMP_MAP if not np.isnan(ratio[s])]
lr = np.log([ratio[s] for s in sec]); lv = np.log([vapw[s] for s in sec])
X = np.column_stack([np.ones(len(sec)), lv])
b, *_ = np.linalg.lstsq(X, lr, rcond=None)
res = lr - X @ b
R2 = 1 - res @ res / np.sum((lr - lr.mean()) ** 2)
rk = lambda x: np.argsort(np.argsort(-np.asarray(x)))
spear = np.corrcoef(rk(lr), rk(lv))[0, 1]
print(f"log(pi/share) = {b[0]:+.3f} + {b[1]:.3f} log(VA/worker): R2 = {R2:.3f}; Spearman = {spear:.3f}; n = {len(sec)}")
print("residuals:", {s: round(r, 2) for s, r in zip(sec, res)})
# log(pi/share) = log(pi/VA share) + log(VA per worker relative): decomposition
vash = {s: va[s] / sum(va[t] for t in EMP_MAP) for s in EMP_MAP}
print("pi / VA share:", {s: round(pi_e[s] / vash[s], 2) for s in sec})
assert round(ratio["Financial services"], 1) == 23.4 and round(ratio["Agriculture"], 2) == 0.33
REG3 = dict(b=b, R2=R2, spear=spear)
# (d) per-worker allocation, agriculture -> manufacturing; N = 11 sectors, phi = Phi/Omega
N11 = 11
for ph in (0.0, 0.25, 0.5, 1.0):
    rel = {s: (1 / N11 + (pd_[s] - 1 / N11) * ph) / emp[s] for s in ("Agriculture", "Manufacturing")}
    print(f"  Phi/Omega = {ph:.2f}: per-worker allocation / national mean: agriculture {rel['Agriculture']:.3f}, "
          f"manufacturing {rel['Manufacturing']:.3f}, ratio {rel['Manufacturing']/rel['Agriculture']:.2f}")
rel1 = {s: pd_[s] / emp[s] for s in ("Agriculture", "Manufacturing")}
assert abs(rel1["Manufacturing"] / rel1["Agriculture"] - 4.3) < 0.05

# =====================================================================
print("\n== Exercise 17.4: N_eff and refinement of the node set ==")
s4 = np.array([0.34, 0.37, 0.12, 0.17])
print(f"four ecosystems: Neff = {1/np.sum(s4**2):.3f}")
for k in (2, 5):
    so = np.repeat(s4 / k, k)
    print(f"  each ecosystem split into {k} equal outlets: Neff = {1/np.sum(so**2):.3f} (= {k} x {1/np.sum(s4**2):.3f})")
so_uneq = np.r_[0.34 * np.array([0.8, 0.2]), 0.37 * np.array([0.5, 0.5]), 0.12 * np.array([1.0]), 0.17 * np.array([0.6, 0.4])]
print(f"  an unequal split (0.8/0.2, 0.5/0.5, 1, 0.6/0.4): Neff = {1/np.sum(so_uneq**2):.3f}")
assert 1 / np.sum(so_uneq ** 2) >= 1 / np.sum(s4 ** 2)

# =====================================================================
print("\n== Exercise 17.5: trajectory, every available year (book12, domestic, input shares) ==")
tr = SUM[(SUM.aggregation == "book12") & (SUM.flows == "domestic") & (SUM.normalization == "input_shares")].copy()
PI = pd.read_csv(os.path.join(PROC, "bolivia_io_pi.csv"))
ph_ = PI[(PI.aggregation == "book12") & (PI.flows == "domestic") & (PI.normalization == "input_shares")
         & (PI.sector == "Hydrocarbons")].set_index(["base", "year"]).pi
tr["pi_h"] = [ph_.loc[(b_, y_)] for b_, y_ in zip(tr.base, tr.year)]
for r in tr.itertuples():
    print(f"  {r.base} {r.year}: N = {r.N}, Gamma = {r.Gamma:.3f}, G(pi) = {r.power_gini:.3f}, pi_hydro = {r.pi_h:.3f}")
b90 = tr[tr.base == 1990].set_index("year")
z = (b90.Gamma - b90.Gamma.mean()) / b90.Gamma.std()
print("base-1990 z-scores of Gamma beyond |z| > 1.5:", z[abs(z) > 1.5].round(2).to_dict())
# 2014: which sector carries the slow mode?
Z14, l14 = read_Z("Z_1990_2014_book12.csv")
W14, s14 = wp(Z14, l14)
ev14, V14 = np.linalg.eig(W14)
o = np.argsort(-np.abs(ev14)); v2 = np.real(V14[:, o[1]]); v2 /= np.abs(v2).max()
Z13, l13 = read_Z("Z_1990_2013_book12.csv")
W13, s13 = wp(Z13, l13)
print("2014 slow mode v2:", dict(zip(s14, np.round(v2, 2))))
print("diagonal of W^P 2013 vs 2014:", {s: (round(W13[s13.index(s), s13.index(s)], 3), round(W14[i, i], 3)) for i, s in enumerate(s14)})
dmax = max(s14, key=lambda s: abs(W14[s14.index(s), s14.index(s)] - W13[s13.index(s), s13.index(s)]))
print(f"largest diagonal change 2013 -> 2014: {dmax}")
# structural break tests on the base-1990 annual series (Chow, trend model), 1988-2013
yrs = b90.index[b90.index <= 2013].values; g = b90.loc[yrs, "Gamma"].values


def fcdf(x, d1, d2, n=20000):
    """CDF of the F distribution by numerical integration (Simpson), no scipy."""
    if x <= 0:
        return 0.0
    lb = math.lgamma((d1 + d2) / 2) - math.lgamma(d1 / 2) - math.lgamma(d2 / 2)
    def pdf(t):
        return math.exp(lb + (d1 / 2) * math.log(d1 / d2) + (d1 / 2 - 1) * math.log(t)
                        - ((d1 + d2) / 2) * math.log(1 + d1 * t / d2)) if t > 0 else 0.0
    # integrate in u = t/(1+t) to handle the tail
    us = np.linspace(1e-12, x / (1 + x), n + 1)
    f = np.array([pdf(u / (1 - u)) / (1 - u) ** 2 for u in us])
    h = us[1] - us[0]
    return float(h / 3 * (f[0] + f[-1] + 4 * f[1:-1:2].sum() + 2 * f[2:-1:2].sum()))


assert abs(fcdf(3.4928, 2, 20) - 0.95) < 2e-3            # F(2,20) 5% critical value
BREAK = {}
for excl in (False, True):
    mask = np.ones(len(yrs), bool)
    if excl:
        mask &= yrs != 1999
    yy, gg = yrs[mask], g[mask]
    X0 = np.column_stack([np.ones(len(yy)), yy - 2000])
    rss0 = np.sum((gg - X0 @ np.linalg.lstsq(X0, gg, rcond=None)[0]) ** 2)
    for tb in (2006, 2009):
        D = (yy >= tb).astype(float)
        X1 = np.column_stack([X0, D, D * (yy - 2000)])
        rss1 = np.sum((gg - X1 @ np.linalg.lstsq(X1, gg, rcond=None)[0]) ** 2)
        F = ((rss0 - rss1) / 2) / (rss1 / (len(yy) - 4))
        pval = 1 - fcdf(F, 2, len(yy) - 4)
        BREAK[(excl, tb)] = (F, pval, len(yy))
        print(f"  Chow test (level + trend) at {tb}, n = {len(yy)}{' excl. 1999' if excl else ''}: F = {F:.3f}, p = {pval:.3f}")
    slope = np.linalg.lstsq(X0, gg, rcond=None)[0][1]
    print(f"  OLS trend{' excl. 1999' if excl else ''}: {10*slope:+.4f} per decade")
# (d) sigma_B projection with the illustrative trajectory 0.68 (1990) -> 0.62 (2019)
rate = (0.62 - 0.68) / 29
t40 = 2019 + (0.40 - 0.62) / rate
print(f"illustrative sigma_B trend: {rate*10:+.4f} per decade; reaches 0.40 in {t40:.0f}")

# =====================================================================
print("\n== Exercise 17.6: dashboard comparison Bolivia, Chile, Costa Rica (2019, book12_isic) ==")
osum = pd.read_csv(os.path.join(PROC, "oecd_io_summary.csv"))
o = osum[(osum.edition == "ICIO2023") & (osum.year == 2019) & (osum.aggregation == "book12_isic")
         & (osum.normalization == "input_shares")].set_index("country")
bi = summ(2019, "book12_isic")
D6 = {"BOL": (bi.Gamma, bi.power_gini, bi.Neff), "CHL": (o.loc["CHL", "Gamma"], o.loc["CHL", "power_gini"], o.loc["CHL", "Neff"]),
      "CRI": (o.loc["CRI", "Gamma"], o.loc["CRI", "power_gini"], o.loc["CRI", "Neff"])}
for c, (g_, G_, ne) in D6.items():
    print(f"  {c}: Gamma = {g_:.3f}, G(pi) = {G_:.3f}, Neff(W^P) = {ne:.2f}")
opi = pd.read_csv(os.path.join(PROC, "oecd_io_pi.csv"))
for c in ("CHL", "CRI"):
    q = opi[(opi.edition == "ICIO2023") & (opi.country == c) & (opi.year == 2019) & (opi.aggregation == "book12_isic")
            & (opi.normalization == "input_shares")].sort_values("pi", ascending=False)
    print(f"  {c} top pi:", list(zip(q.sector.head(3), q.pi.head(3).round(3))))
ec = pd.read_csv(os.path.join(PROC, "eci_facts.csv"))
print("  ECI 2019:", {c: ec[(ec.country == c) & (ec.year == 2019)].eci_hs92.iloc[0] for c in ("BOL", "CHL")},
      "(CRI not in the processed ECI file)")
rkr = pd.read_csv(os.path.join(ROOT, "data", "raw", "atlas", "growth_proj_eci_rankings.csv"))
eci_cri = rkr[(rkr.country_iso3_code == "CRI") & (rkr.year == 2019)].eci_hs92.iloc[0]
print(f"  ECI 2019 CRI (Growth Lab rankings file, data/raw/atlas) = {eci_cri}")
wdi = pd.read_csv(os.path.join(PROC, "wdi_long.csv"))
gi = wdi[(wdi.indicator == "SI.POV.GINI") & (wdi.year.between(2014, 2019))].dropna(subset=["value"]).sort_values("year").groupby("country_iso3").last()
print("  WDI Gini latest 2014-2019:", {c: (gi.loc[c, "value"], int(gi.loc[c, "year"])) for c in ("BOL", "CHL") if c in gi.index})
# rescaled dashboard bars (benchmark = 1), measured instruments only
for c, (g_, G_, ne) in D6.items():
    print(f"  bars {c}: Gamma {g_:.3f}, 1 - G(pi) {1-G_:.3f}, Neff/N {ne/11:.3f}")

# =====================================================================
print("\n== Exercise 17.7(e): the Chapter 10 three-block model, 30 periods ==")
n = 3; N = 9
blk = np.repeat([0, 1, 2], n)
same = blk[:, None] == blk[None, :]
WP = np.where(same, 0.94 / n, 0.01)
WI0 = np.where(same, 0.85 / n, 0.025)
alpha = np.array([0.06, 0.04, 0.03])[blk]; beta = np.array([0.10, 0.07, 0.04])[blk]; gam = 0.02
pib = np.array([0.15, 0.11, 1 / 3 - 0.15 - 0.11])[blk]
omega, Phi, phi = 0.5, 4.5, 1.0
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.1
P = np.eye(N) - np.ones((N, N)) / N
C = np.diag(beta) + gam * WP


def alloc(th):
    return omega + pib * Phi + phi * (th - th.mean())


def T(th, WI, rho_=rho, mu_=mu):
    x = alloc(th)
    r = alpha + beta * x + gam * WP @ x
    return rho_ * th + (1 - rho_) * th0 + mu_ * (WI @ r - rstar)


def steady(WI, rho_=rho, mu_=mu):
    J = rho_ * np.eye(N) + mu_ * phi * WI @ C @ P
    c = T(np.zeros(N), WI, rho_, mu_)
    return np.linalg.solve(np.eye(N) - J, c), J


th_a, J0 = steady(WI0)
print(f"Chapter 10 check: theta* by block = {th_a[::3].round(4)}, rho(J) = {max(abs(np.linalg.eigvals(J0))):.3f}")
assert np.allclose(th_a[::3], [1.125, 1.041, 0.987], atol=6e-4) and round(max(abs(np.linalg.eigvals(J0))), 3) == 0.866
# information operator: cross-block weight per pair 0.025 -> 0.05, own-block weights cut equally
WI1 = np.where(same, 0.70 / n, 0.05)
assert np.allclose(WI1.sum(1), 1)
th_b, J1 = steady(WI1)


def VBx(th):
    x = alloc(th)
    return np.mean((x - x.mean()) ** 2)


print(f"new steady state theta* = {th_b[::3].round(4)}; gap theta1 - theta3: {th_a[0]-th_a[-1]:.4f} -> {th_b[0]-th_b[-1]:.4f}")
print(f"allocation variance V(x): {VBx(th_a):.5f} -> {VBx(th_b):.5f} (power part alone: {np.mean((pib*Phi - pib.mean()*Phi)**2):.5f})")
print(f"Gamma(W^I): {gap(WI0):.3f} -> {gap(WI1):.3f}; Gamma(W^P) = {gap(WP):.3f} (unchanged); Neff(pi) = {1/np.sum(pib**2):.2f} (unchanged)")
for lab_, rho_ in (("generational periods (rho = 0.85)", rho), ("annual periods (rho = 0.85^(1/25))", 0.85 ** (1 / 25))):
    mu_ = mu * (1 - rho_) / (1 - rho)                      # same formation multiplier M = mu/(1-rho)
    tha, _ = steady(WI0, rho_, mu_); thb, Jb = steady(WI1, rho_, mu_)
    th = tha.copy(); traj = []
    for t in range(31):
        traj.append((th[0] - th[-1], VBx(th)))
        th = T(th, WI1, rho_, mu_)
    g0, g_inf = tha[0] - tha[-1], thb[0] - thb[-1]
    frac = [(g0 - traj[t][0]) / (g0 - g_inf) for t in (1, 10, 30)]
    print(f"  {lab_}: mu = {mu_:.5f}, rho(J) = {max(abs(np.linalg.eigvals(Jb))):.4f}; gap at t = 0, 10, 30: "
          f"{traj[0][0]:.4f}, {traj[10][0]:.4f}, {traj[30][0]:.4f} (long run {g_inf:.4f}); share of adjustment done "
          f"t=1: {frac[0]:.3f}, t=10: {frac[1]:.3f}, t=30: {frac[2]:.3f}")
print("\nall Chapter 17 solution checks passed")
