"""Solutions manual, Chapter 19 (Scale Invariance).  Computes and prints every
number used in solutions/sol19.tex, with asserts on the key results.

Sourced inputs are read from data/processed (retrieved 2026-11-05; see
data/README.md).  Everything else is illustrative and labelled as such.
Run from the repo root:  python3 code/solutions/sol19.py   (numpy + pandas).
"""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=150)
PROC = os.path.join("data", "processed")


# ---------------------------------------------------------------- helpers
def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def neff(p):
    return 1 / np.sum(p ** 2)


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def rownorm(K):
    return K / K.sum(1, keepdims=True)


def stoch_complement(W, idx):
    rest = [j for j in range(W.shape[0]) if j not in idx]
    return W[np.ix_(idx, idx)] + W[np.ix_(idx, rest)] @ np.linalg.solve(
        np.eye(len(rest)) - W[np.ix_(rest, rest)], W[np.ix_(rest, idx)])


def read_Z(fname):
    df = pd.read_csv(os.path.join(PROC, "bolivia_Z", fname), index_col=0)
    assert list(df.index) == list(df.columns)
    return df.values.astype(float), list(df.columns)


def wp_from_Z(Z, labels, sales=False):
    """W^P = D^-1 A (input shares) after iteratively dropping sectors with no
    intermediate purchases or sales (data/README.md, section 1)."""
    keep = np.ones(len(labels), bool)
    while True:
        Zs = Z[np.ix_(keep, keep)]
        bad = (Zs.sum(0) <= 0) | (Zs.sum(1) <= 0)
        if not bad.any():
            break
        keep[np.where(keep)[0][bad]] = False
    Zk = Z[np.ix_(keep, keep)]
    K = Zk if sales else Zk.T
    return rownorm(K), [l for l, k in zip(labels, keep) if k]


FACTS = pd.read_csv(os.path.join(PROC, "bolivia_facts.csv")).set_index("fact_id")["value"]
WDI = pd.read_csv(os.path.join(PROC, "wdi_long.csv"))
GPI = pd.read_csv(os.path.join(PROC, "global_wp_pi.csv"))
GPI19 = GPI[GPI.year == 2019].set_index("country")
GSUM = pd.read_csv(os.path.join(PROC, "global_wp_summary.csv")).set_index("year").loc[2019]
BSUM = pd.read_csv(os.path.join(PROC, "bolivia_io_summary.csv"))
OSUM = pd.read_csv(os.path.join(PROC, "oecd_io_summary.csv"))
ILO = pd.read_csv(os.path.join(PROC, "ilostat_bolivia_emp_eco_sex.csv"))
ILO19 = ILO[ILO.year == 2019].pivot(index="eco", columns="sex", values="value_thousands")


def wdi(c, ind, y):
    v = WDI[(WDI.country_iso3 == c) & (WDI.indicator == ind) & (WDI.year == y)]["value"]
    assert len(v) == 1, (c, ind, y)
    return float(v.iloc[0])


def bsum(year, agg="book12", norm="input_shares", flows="domestic"):
    m = BSUM[(BSUM.year == year) & (BSUM.aggregation == agg) & (BSUM.normalization == norm) & (BSUM.flows == flows)]
    assert len(m) == 1
    return m.iloc[0]


def osum(c, year, ed="ICIO2023", agg="book12_isic", norm="input_shares"):
    m = OSUM[(OSUM.country == c) & (OSUM.year == year) & (OSUM.edition == ed) & (OSUM.aggregation == agg)
             & (OSUM.normalization == norm) & (OSUM.flows == "domestic")]
    assert len(m) == 1
    return m.iloc[0]


# =====================================================================
print("=" * 70); print("Exercise 19.1  Three-scale comparison"); print("=" * 70)
# --- global column: recompute from the bilateral trade file -------------
bil = pd.read_csv(os.path.join(PROC, "atlas_bilateral_selected_years.csv.gz"))
b = bil[bil.year == 2019]
K = b.pivot_table(index="importer", columns="exporter", values="trade_value", aggfunc="sum").fillna(0)
cs = sorted(set(K.index) | set(K.columns))
Kv = K.reindex(index=cs, columns=cs).fillna(0).values.copy()
np.fill_diagonal(Kv, 0)
keep = np.ones(len(cs), bool)
while True:
    M_ = Kv[np.ix_(keep, keep)]
    bad = (M_.sum(0) <= 0) | (M_.sum(1) <= 0)
    if not bad.any():
        break
    keep[np.where(keep)[0][bad]] = False
cc = [c for c, k in zip(cs, keep) if k]
KG = Kv[np.ix_(keep, keep)]
WG = rownorm(KG)
pG = perron(WG); evG = eig_sorted(WG); GamG = gap(WG)
iG = {c: k for k, c in enumerate(cc)}
print(f"global W (imports only, 2019): N = {len(cc)}, lambda_2 = {evG[1]:.4f}, Gamma_G = {GamG:.4f}")
assert len(cc) == int(GSUM["N"]) == 230 and abs(GamG - float(GSUM["Gamma"])) < 1e-9
assert abs(pG[iG["BOL"]] - GPI19.loc["BOL", "pi"]) < 1e-12
argmax, argmin = cc[int(pG.argmax())], cc[int(pG.argmin())]
print(f"pi_max = {pG.max():.4f} ({argmax}); pi_min = {pG.min():.3e} ({argmin}); max/min = {pG.max()/pG.min():.3e}")
print(f"pi_USA/pi_LAO = {pG[iG['USA']]/pG[iG['LAO']]:.1f}; Neff = {neff(pG):.2f} (Neff/N = {neff(pG)/len(cc):.3f}); G(pi) = {gini(pG):.3f}")
assert argmax == "CHN" and argmin == "HMD"
assert round(pG[iG["USA"]] / pG[iG["LAO"]]) == 193
# max/min over economies with a labour force in WDI is not available for all 230; report
# the 10th/90th percentile ratio instead (robust to uninhabited territories)
q10, q90 = np.quantile(pG, [0.1, 0.9])
print(f"pi P90/P10 = {q90/q10:.1f}")
# --- national column -------------------------------------------------
Z12, lab12 = read_Z("Z_2017_2019_book12.csv")
WN, secN = wp_from_Z(Z12, lab12)
pN = perron(WN); evN = eig_sorted(WN); GamN = gap(WN)
print(f"national W^P (INE 2019, book12, domestic, input shares): N = {len(secN)}, lambda_2 = {evN[1].real:.4f} "
      f"(imag {abs(evN[1].imag):.1e}), |lambda_3| = {abs(evN[2]):.3f}, Gamma_N = {GamN:.4f}")
assert len(secN) == 11 and round(GamN, 3) == 0.467
jmax, jmin = int(pN.argmax()), int(pN.argmin())
print(f"pi_max = {pN[jmax]:.4f} ({secN[jmax]}), pi_min = {pN[jmin]:.4f} ({secN[jmin]}), max/min = {pN[jmax]/pN[jmin]:.1f}")
assert secN[jmax] == "Other services" and secN[jmin] == "Public administration" and round(pN[jmax] / pN[jmin]) == 53
print(f"Neff = {neff(pN):.2f} of 11 (Neff/N = {neff(pN)/11:.3f}); G(pi) = {gini(pN):.3f}")
print(f"Gamma_N/Gamma_G = {GamN/GamG:.3f}")
VAR = {}
for agg in ("book12", "book12_isic", "detail39"):
    for norm in ("input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag"):
        VAR[(agg, norm)] = float(bsum(2019, agg, norm)["Gamma"])
for k, v in VAR.items():
    print(f"  Gamma_N {k[0]:12s} {k[1]:24s} = {v:.3f}")
print(f"range over the 12 variants: {min(VAR.values()):.3f} - {max(VAR.values()):.3f}")
# --- lumpability: exact (Prop. 19.lumpability) and a non-lumpable counterexample ---
print("\n-- Aggregation and Gamma --")
for a in (0.01, 0.25):
    W3 = np.array([[1 - a, a, 0], [0, 0, 1], [1, 0, 0]])
    p3 = perron(W3)
    blocks = [[0], [1, 2]]
    Wt = np.array([[p3[A] @ W3[np.ix_(A, B)].sum(1) / p3[A].sum() for B in blocks] for A in blocks])
    Wu = np.array([[W3[np.ix_(A, B)].sum() / W3[A].sum() for B in blocks] for A in blocks])
    assert np.allclose(perron(Wt), [p3[A].sum() for A in blocks])
    assert np.isclose(gap(W3), 1 - np.sqrt(a)) and np.isclose(gap(Wt), 1 - abs(0.5 - a)) and np.allclose(Wt, Wu)
    lump = np.allclose(W3[1, [1, 2]].sum(), W3[2, [1, 2]].sum())
    print(f"a = {a}: pi = {np.round(p3,4)}, Gamma_fine = 1-sqrt(a) = {gap(W3):.3f}, "
          f"Gamma_coarse = 1-|1/2-a| = {gap(Wt):.3f}; exactly lumpable: {lump}")
print(f"real data, aggregation detail39 -> book12 (INE 2019): {VAR[('detail39','input_shares')]:.3f} -> {VAR[('book12','input_shares')]:.3f}")
for c in ("KOR", "DEU", "USA", "MEX", "BRA"):
    print(f"  ICIO 2019 {c}: detail45 {float(osum(c, 2019, agg='detail45')['Gamma']):.3f} -> book12_isic {float(osum(c, 2019)['Gamma']):.3f}")
# --- individual scale: ARD sample size -------------------------------
z = 1.959964
for pbar, h, deff in ((0.2, 0.02, 2.0),):
    n_srs = z ** 2 * pbar * (1 - pbar) / h ** 2
    print(f"sample size to estimate a cross-group tie share p = {pbar} within +/-{h} (95%): SRS {n_srs:.0f}; "
          f"with design effect {deff}: {n_srs*deff:.0f} respondents per group")
    n_per_group = n_srs * deff
assert round(n_per_group) == 3073

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.2  Power nesting: three Bolivian agents"); print("=" * 70)
p_of = dict(zip(secN, pN))
EMP_MAP = {"Hydrocarbons and mining": "B", "Agriculture": "A", "Personal services": "IRS"}
emp_tot = ILO19.loc["ECO_ISIC4_TOTAL", "SEX_T"]
emp = {s: sum(ILO19.loc["ECO_ISIC4_" + c, "SEX_T"] for c in cs_) / emp_tot for s, cs_ in EMP_MAP.items()}
fem = {s: sum(ILO19.loc["ECO_ISIC4_" + c, "SEX_F"] for c in cs_) / sum(ILO19.loc["ECO_ISIC4_" + c, "SEX_T"] for c in cs_)
       for s, cs_ in EMP_MAP.items()}
pis = {"Hydrocarbons and mining": p_of["Hydrocarbons"] + p_of["Mining"], "Agriculture": p_of["Agriculture"],
       "Personal services": p_of["Personal services"]}
sf = {s: pis[s] / emp[s] for s in pis}
for s in pis:
    print(f"{s:24s}: pi = {pis[s]:.4f}, employment share = {emp[s]:.4f}, sector factor pi/e = {sf[s]:.3f}, female share = {fem[s]:.3f}")
assert round(sf["Hydrocarbons and mining"], 1) == 5.7 and round(sf["Agriculture"], 2) == 0.33 and round(sf["Personal services"], 3) == 0.198
r_eth, r_gen = 0.57, 0.58                      # Chapter 18, illustrative per-capita ratios
w = {"a": 1.0, "b": r_eth, "c": r_eth * r_gen}
print(f"within-sector factors (relative to a non-indigenous man; Ch.18 illustrative): a 1, b {w['b']:.2f}, c {w['c']:.4f}")
P = {"a": sf["Hydrocarbons and mining"] * w["a"], "b": sf["Agriculture"] * w["b"], "c": sf["Personal services"] * w["c"]}
for k in P:
    print(f"  per-capita power index of agent {k} (country factor common) = {P[k]:.4f}")
pairs = [("a", "b"), ("a", "c"), ("b", "c")]
RAT = {}
for i, j in pairs:
    RAT[(i, j)] = P[i] / P[j]
    s_ratio = {"a": sf["Hydrocarbons and mining"], "b": sf["Agriculture"], "c": sf["Personal services"]}
    lsec, lwit = np.log(s_ratio[i] / s_ratio[j]), np.log(w[i] / w[j])
    print(f"  {i}/{j} = {RAT[(i,j)]:.2f} = sector {s_ratio[i]/s_ratio[j]:.2f} x within {w[i]/w[j]:.3f}; "
          f"log share of sector factor = {lsec/(lsec+lwit):.3f}")
assert round(RAT[("a", "c")]) == 88 and round(RAT[("a", "b")], 1) == 30.3 and round(RAT[("b", "c")], 2) == 2.90
# sensitivity: hydrocarbons alone, employment unknown but < mining+hydrocarbons employment
lb = p_of["Hydrocarbons"] / emp["Hydrocarbons and mining"]
print(f"hydrocarbons alone: pi_h / e_B = {p_of['Hydrocarbons']:.4f}/{emp['Hydrocarbons and mining']:.4f} = {lb:.2f} "
      f"(a lower bound on the hydrocarbons factor); a/c lower bound = {lb*w['a']/(sf['Personal services']*w['c']):.1f}")
# mean-one normalisation of the within factor, using gender composition (ILOSTAT) only
msec = {s: (1 - fem[s]) * 1 + fem[s] * r_gen for s in fem}
for s in msec:
    print(f"  mean within factor (gender only) in {s}: {msec[s]:.3f}")
Pn = {"a": P["a"] / msec["Hydrocarbons and mining"], "b": P["b"] / msec["Agriculture"], "c": P["c"] / msec["Personal services"]}
for i, j in pairs:
    print(f"  with mean-one normalisation: {i}/{j} = {Pn[i]/Pn[j]:.2f}")
RATn = {(i, j): Pn[i] / Pn[j] for i, j in pairs}

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.3  The global comparison"); print("=" * 70)
piG = {c: float(GPI19.loc[c, "pi"]) for c in ("DEU", "BRA", "BOL")}
L = {c: wdi(c, "SL.TLF.TOTL.IN", 2019) for c in ("DEU", "BRA", "BOL")}
for c in piG:
    print(f"{c}: pi_G = {piG[c]:.6f}, labour force 2019 = {L[c]/1e6:.2f} M, pi_G/L = {piG[c]/L[c]:.3e} per worker")
cf = {c: (piG[c] / L[c]) / (piG["BOL"] / L["BOL"]) for c in piG}
print(f"country per-worker factor relative to Bolivia: DEU {cf['DEU']:.2f}, BRA {cf['BRA']:.3f}")
assert round(cf["DEU"], 1) == 21.2 and round(cf["BRA"], 2) == 1.45
s_vendor, w_vendor = sf["Personal services"], w["c"]
agents = {"German engineer": (cf["DEU"], 1.0, 1.0), "Brazilian factory worker": (cf["BRA"], 1.0, 1.0),
          "Potosi vendor": (1.0, s_vendor, w_vendor)}
for k, (c_, s_, w_) in agents.items():
    print(f"  {k:26s}: country {c_:.3f} x sector {s_:.3f} x within {w_:.4f} = {c_*s_*w_:.4f}")
Q = {k: np.prod(v) for k, v in agents.items()}
R3 = {"DEU/vendor": Q["German engineer"] / Q["Potosi vendor"], "BRA/vendor": Q["Brazilian factory worker"] / Q["Potosi vendor"],
      "DEU/BRA": Q["German engineer"] / Q["Brazilian factory worker"]}
for k, v in R3.items():
    print(f"  ratio {k} = {v:.1f}")
assert round(R3["DEU/vendor"]) == 324 and round(R3["BRA/vendor"], 1) == 22.2 and round(R3["DEU/BRA"], 1) == 14.6
for k in ("German engineer", "Brazilian factory worker"):
    c_ = agents[k][0]
    logs = np.log([c_, 1 / s_vendor, 1 / w_vendor])
    print(f"  log shares of {k}/vendor gap: country {logs[0]/logs.sum():.3f}, sector {logs[1]/logs.sum():.3f}, within {logs[2]/logs.sum():.3f}")
noL = {c: piG[c] / piG["BOL"] for c in piG}
print(f"forgetting the workforce: DEU/vendor = {noL['DEU']/s_vendor/w_vendor:.0f} (x{L['DEU']/L['BOL']:.2f}), "
      f"BRA/vendor = {noL['BRA']/s_vendor/w_vendor:.0f} (x{L['BRA']/L['BOL']:.2f}); DEU/BRA = {noL['DEU']/noL['BRA']:.2f} instead of {R3['DEU/BRA']:.1f}")
noLr = {"DEU": noL["DEU"] / s_vendor / w_vendor, "BRA": noL["BRA"] / s_vendor / w_vendor}
# sensitivity: apply Chapter 18's illustrative gender factor to the German engineer
print(f"if the German engineer's within factor is 0.58 (Ch.18 gender ratio): DEU/vendor = {R3['DEU/vendor']*0.58:.0f}")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.4  Five processes in one 6x6 matrix (illustrative)"); print("=" * 70)
lab6 = ["C1", "C2", "S1", "S2", "P1", "P2"]
W6 = np.array([
    [0.50, 0.46, 0.03, 0.00, 0.01, 0.00],   # core: dense, buys little from the rest
    [0.46, 0.50, 0.00, 0.03, 0.01, 0.00],
    [0.04, 0.00, 0.60, 0.36, 0.00, 0.00],   # semi-periphery: own block, a little core
    [0.00, 0.04, 0.36, 0.60, 0.00, 0.00],
    [0.30, 0.20, 0.00, 0.00, 0.45, 0.05],   # periphery enclave: buys mainly from the core
    [0.15, 0.10, 0.05, 0.00, 0.20, 0.50]])  # subsistence node: nobody buys from it but P1, a little
assert np.allclose(W6.sum(1), 1) and (W6 >= 0).all()
assert np.linalg.matrix_power(W6 + np.eye(6), 6).min() > 0      # irreducible
ev6 = eig_sorted(W6); p6 = perron(W6)
print("eigenvalues:", np.round(ev6, 4))
print(f"Gamma (whole economy) = {gap(W6):.4f}")
print("pi =", dict(zip(lab6, np.round(p6, 4))), f" Neff = {neff(p6):.2f} of 6")
print("column sums:", dict(zip(lab6, np.round(W6.sum(0), 3))))
blk = {"core": [0, 1], "semi": [2, 3], "periphery": [4, 5]}
Gr = {}
for k, ix in blk.items():
    Wr = rownorm(W6[np.ix_(ix, ix)])
    Gr[k] = gap(Wr)
    print(f"  region {k}: renormalised block = {np.round(Wr,3).tolist()}, Gamma_r = {Gr[k]:.4f}, block power = {p6[ix].sum():.4f}")
varG = np.var(list(Gr.values()))
print(f"Var_r(Gamma_r) = {varG:.4f}")
cross = {(a, b_): W6[np.ix_(blk[a], blk[b_])].sum(1).mean() for a in blk for b_ in blk if a != b_}
print("mean row weight on other blocks:", {f"{a}->{b_}": round(v, 3) for (a, b_), v in cross.items()})
print(f"dependence: periphery->core {cross[('periphery','core')]:.3f} vs core->periphery {cross[('core','periphery')]:.3f}; "
      f"pi_periphery = {p6[4:].sum():.4f} vs 2/N = {2/6:.3f}")
print(f"marginality: pi_P2 = {p6[5]:.4f} (1/N = {1/6:.3f}); column sum of P2 = {W6[:,5].sum():.3f}, excluding its self-loop {W6[:,5].sum()-W6[5,5]:.3f}")
# conductance between the two weakly linked blocks core+periphery vs semi
Sset = [2, 3]
Qflow = sum(p6[i] * W6[i, j] for i in Sset for j in range(6) if j not in Sset)
Phi = Qflow / min(p6[Sset].sum(), 1 - p6[Sset].sum())
print(f"flow out of the semi block = {Qflow:.5f}; conductance Phi = {Phi:.4f}; 1 - lambda_2 = {1-ev6[1].real:.4f} <= 2 Phi = {2*Phi:.4f}")
assert 1 - ev6[1].real <= 2 * Phi + 1e-12
assert gap(W6) < 0.1 and Gr["core"] > 0.9 and Gr["periphery"] < Gr["core"] and p6[4:].sum() < 1 / 3 and p6[5] < 0.03
# development: raise core-semi coupling (a path that raises Gamma)
def couple(t):
    X = W6.copy()
    for i, j in ((0, 2), (1, 3), (2, 0), (3, 1)):
        X[i, j] += t
        X[i, i] -= t
    return X
for t in (0.0, 0.05, 0.10, 0.20):
    print(f"  core-semi coupling +{t:.2f}: Gamma = {gap(couple(t)):.3f}, Neff = {neff(perron(couple(t))):.2f}")
G_couple = {t: gap(couple(t)) for t in (0.0, 0.05, 0.10, 0.20)}

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.5  Cross-scale feedback simulation (illustrative)"); print("=" * 70)
NC, NS, NG = 3, 3, 2
NN = NC * NS * NG


def build(e, lam, att, gx=0.05):
    """Raw intensities: own 1; same sector, other type lam_c; other sector in the same country
    e_c (x lam_c across types); other country gx * att_c' (x lam_c' if the partner is type D)."""
    Kx = np.zeros((NN, NN))
    for i in range(NN):
        ci, si, gi = i // 6, (i // 2) % 3, i % 2
        for j in range(NN):
            cj, sj, gj = j // 6, (j // 2) % 3, j % 2
            if ci == cj:
                Kx[i, j] = (1.0 if gi == gj else lam[ci]) * (1.0 if si == sj else e[ci])
            else:
                Kx[i, j] = gx * att[cj] * (1.0 if gj == 0 else lam[cj])
    return rownorm(Kx)


sig = lambda x: 1 / (1 + np.exp(-x))


def measures(Wx):
    px = perron(Wx)
    piGx = np.array([px[6 * c:6 * c + 6].sum() for c in range(NC)])
    GNx = np.array([gap(stoch_complement(Wx, list(range(6 * c, 6 * c + 6)))) for c in range(NC)])
    return px, piGx, GNx


def step(state, beta):
    e, lam, att = state
    Wx = build(e, lam, att)
    px, piGx, GNx = measures(Wx)
    e_t = 0.02 + 0.98 * sig(beta * ((NC * piGx - 1) + (lam - 0.5)))   # G->N and I->N
    l_t = 0.05 + 0.95 * sig(beta * (GNx - 0.5))                         # N->I
    a_t = 0.05 + 0.95 * sig(beta * (GNx - 0.5))                         # N->G
    eta = 0.2
    return (e + eta * (e_t - e), lam + eta * (l_t - lam), att + eta * (a_t - att)), piGx, GNx


def run(beta, s0, T=400):
    s = tuple(np.array(x, float) for x in s0)
    for _ in range(T):
        s, piGx, GNx = step(s, beta)
    return s, piGx, GNx


S0 = ([0.9, 0.5, 0.05], [0.9, 0.5, 0.1], [0.9, 0.5, 0.1])        # core, semi, periphery
SIM = {}
for beta in (1, 2, 3, 4, 6, 8):
    (e, lam, att), piGx, GNx = run(beta, S0)
    SIM[beta] = (e, lam, att, piGx, GNx)
    print(f"beta = {beta}: e = {np.round(e,3)}, lambda = {np.round(lam,3)}, att = {np.round(att,3)}, "
          f"NC*pi_G = {np.round(NC*piGx,3)}, Gamma_N = {np.round(GNx,3)}")
# critical feedback strength: bisection on whether the periphery ends low
def trapped(beta):
    (e, lam, att), piGx, GNx = run(beta, S0)
    return e[2] < 0.5 * e[0] and NC * piGx[2] < 0.5


lo, hi = 2.0, 4.0
assert not trapped(lo) and trapped(hi)
for _ in range(14):
    mid = (lo + hi) / 2
    if trapped(mid):
        hi = mid
    else:
        lo = mid
beta_star = (lo + hi) / 2
print(f"critical feedback strength beta* = {beta_star:.3f} (unique equilibrium below, periphery trapped above)")
# one-shot pushes at beta = 6: minimal push on a single scale or on all three, for the periphery
BETA = 6
(e6, l6, a6), _, _ = run(BETA, S0)
def escapes(de, dl, da):
    s0 = (e6.copy(), l6.copy(), a6.copy())
    s0[0][2] = min(1, s0[0][2] + de); s0[1][2] = min(1, s0[1][2] + dl); s0[2][2] = min(1, s0[2][2] + da)
    (e, lam, att), piGx, _ = run(BETA, s0, T=400)
    return e[2] > 0.5
PUSH = {}
for name, vec in (("national e only", (1, 0, 0)), ("individual lambda only", (0, 1, 0)),
                  ("global attractiveness only", (0, 0, 1)), ("all three scales", (1, 1, 1))):
    if not escapes(*[v * 1.0 for v in vec]):
        PUSH[name] = None
        print(f"  push on {name}: no push up to 1.0 breaks the trap")
        continue
    lo_, hi_ = 0.0, 1.0
    for _ in range(14):
        mid = (lo_ + hi_) / 2
        if escapes(*[v * mid for v in vec]):
            hi_ = mid
        else:
            lo_ = mid
    PUSH[name] = hi_
    print(f"  push on {name}: minimal size = {hi_:.3f}")
print(f"periphery state before the push: e = {e6[2]:.3f}, lambda = {l6[2]:.3f}, att = {a6[2]:.3f}")
# nesting check at the trapped state: pi_i = pi_G,c * pi_N,s|c * pi_within
Wf = build(e6, l6, a6)
pf = perron(Wf)
countries = [list(range(6 * c, 6 * c + 6)) for c in range(NC)]
recon = np.zeros(NN)
for c, idx in enumerate(countries):
    Sc = stoch_complement(Wf, idx)
    sc_p = perron(Sc)
    secs = [[0, 1], [2, 3], [4, 5]]
    sw = [perron(stoch_complement(Sc, sb)) for sb in secs]
    Cc = np.array([[sw[a] @ Sc[np.ix_(secs[a], secs[b_])].sum(1) for b_ in range(3)] for a in range(3)])
    xiN = perron(Cc)
    for k_, sb in enumerate(secs):
        recon[[idx[i] for i in sb]] = pf[idx].sum() * xiN[k_] * sw[k_]
assert np.allclose(recon, pf)
print("nesting factorisation pi_I = pi_G * pi_N|c * pi_within verified exactly at the trapped state")
dvsA = [pf[6 * c + 1::2][:3].sum() / pf[6 * c::2][:3].sum() for c in range(NC)]
print("per-capita power of type D relative to type A, by country:", np.round(dvsA, 3))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.6  The binding scale: Bolivia, Korea, Uruguay"); print("=" * 70)
LW = wdi("WLD", "SL.TLF.TOTL.IN", 2019)
BIND = {}
for c in ("BOL", "KOR", "URY"):
    pg = float(GPI19.loc[c, "pi"]); Lc = wdi(c, "SL.TLF.TOTL.IN", 2019)
    rel = pg / (Lc / LW)
    BIND[c] = (pg, Lc, rel)
    print(f"{c}: pi_G = {pg:.5f}, labour force = {Lc/1e6:.2f} M, pi_G per worker relative to world = {rel:.3f}")
print(f"world labour force (WDI WLD) = {LW/1e9:.3f} bn; Gamma_G (one number for the whole matrix) = {GamG:.3f}")
GN = {"BOL (book12_isic)": float(bsum(2019, "book12_isic")["Gamma"]), "BOL (book12)": GamN,
      "KOR": float(osum("KOR", 2019)["Gamma"])}
for k, v in GN.items():
    print(f"Gamma_N {k} = {v:.3f}")
assert "URY" not in set(OSUM.country)
print("URY: not in OECD ICIO -> Gamma_N not available")
# Gamma_I proxy for Bolivia: gender exposure across ISIC sections (ILOSTAT 2019)
secs_isic = [s for s in ILO19.index if s.startswith("ECO_ISIC4_") and s != "ECO_ISIC4_TOTAL"]
f_s = ILO19.loc[secs_isic, "SEX_F"].values; m_s = ILO19.loc[secs_isic, "SEX_M"].values
n_s = f_s + m_s; F, Mn = f_s.sum(), m_s.sum()
duncan = 0.5 * np.abs(f_s / F - m_s / Mn).sum()
Q_FF = np.sum(f_s / F * f_s / n_s); Q_MM = np.sum(m_s / Mn * m_s / n_s)
Qg = np.array([[Q_FF, 1 - Q_FF], [1 - Q_MM, Q_MM]])
lam2g = Q_FF + Q_MM - 1
print(f"Bolivia 2019, {len(secs_isic)} ISIC sections: Duncan index D = {duncan:.3f}; isolation Q_FF = {Q_FF:.3f}, Q_MM = {Q_MM:.3f}")
print(f"gender exposure matrix {np.round(Qg,3).tolist()}: lambda_2 = {lam2g:.3f}, Gamma = {1-abs(lam2g):.3f}")
pQ = perron(Qg)
assert np.allclose(pQ, [F / (F + Mn), Mn / (F + Mn)])
print(f"its Perron vector = sex shares of employment ({pQ[0]:.3f}, {pQ[1]:.3f})")
assert round(duncan, 3) == 0.312 and round(1 - abs(lam2g), 3) == 0.802
# equivalent lambda_g in Prop. 18.assortative(a) with a = c: lambda_2 = (1-l)/(1+l)
lam_eq = (1 - lam2g) / (1 + lam2g)
print(f"equal-group equivalent: (1-l)/(1+l) = {lam2g:.3f}  =>  lambda_g = {lam_eq:.3f}")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 19.7  Multi-scale plan for Bolivia"); print("=" * 70)
PCI = pd.read_csv(os.path.join(PROC, "pci_selected_products.csv"))
pci19 = PCI[PCI.year == 2019].set_index("product_hs92_code")
for code, row in pci19.iterrows():
    print(f"  HS {code} {row['description'][:45]:45s} PCI 2019 = {row['pci']:+.3f}")
eci_bol19 = -1.295
# gravity: pi_G tracks world export shares
g19 = GPI19
print(f"corr(pi_G, world export share) across {len(g19)} economies, 2019 = {np.corrcoef(g19.pi, g19.world_export_share)[0,1]:.4f}")
# global: doubling the world's imports from Bolivia, and replacing Bolivia's import row
Kd = KG.copy(); Kd[:, iG["BOL"]] *= 2
pD = perron(rownorm(Kd))
print(f"doubling every economy's imports from Bolivia: pi_BOL {pG[iG['BOL']]:.6f} -> {pD[iG['BOL']]:.6f} "
      f"(x{pD[iG['BOL']]/pG[iG['BOL']]:.3f}); Gamma_G {GamG:.4f} -> {gap(rownorm(Kd)):.4f}")
GG_double = gap(rownorm(Kd)); piBOL_double = pD[iG["BOL"]]
print(f"change in Gamma_G = {GG_double - GamG:+.4f}")
# national operators on the real 2019 W^P (book12, 11 sectors)
SF = pd.read_csv(os.path.join(PROC, "bolivia_io_sector_facts.csv"))
SF19 = SF[(SF.base == 2017) & (SF.year == 2019)].set_index("book12")
iH = lab12.index("Hydrocarbons")
def diag_Z(Zx):
    Wx, sx = wp_from_Z(Zx, lab12)
    px = perron(Wx)
    return gap(Wx), px[sx.index("Hydrocarbons")], neff(px), dict(zip(sx, px))
base = diag_Z(Z12)
print(f"baseline 2019: Gamma = {base[0]:.4f}, pi_hydro = {base[1]:.4f}, Neff = {base[2]:.3f}")
# (i) information operator: W^P unchanged
# (ii) forward linkage / ISI: each sector replaces 10% of its imported intermediate inputs by domestic hydrocarbons
Z2 = Z12.copy()
for k, s in enumerate(lab12):
    if s in SF19.index and s != "Hydrocarbons" and np.isfinite(SF19.loc[s, "imp_int_inputs"]):
        Z2[iH, k] += 0.10 * SF19.loc[s, "imp_int_inputs"]
op2 = diag_Z(Z2)
# (iii) backward linkage: hydrocarbons doubles domestic purchases, the new half from agriculture and manufacturing
Z3 = Z12.copy()
extra = Z12[:, iH].sum()
for s, sh in (("Agriculture", 0.5), ("Manufacturing", 0.5)):
    Z3[lab12.index(s), iH] += sh * extra
op3 = diag_Z(Z3)
Z4 = Z2.copy(); Z4[:, iH] += Z3[:, iH] - Z12[:, iH]
op4 = diag_Z(Z4)
OPS = {"(i) information": base, "(ii) forward linkage (10% import substitution by gas)": op2,
       "(iii) backward linkage (gas doubles domestic purchases, agr+mfg)": op3, "(ii)+(iii)": op4}
for k, v in OPS.items():
    print(f"  {k:62s}: Gamma = {v[0]:.4f} ({v[0]-base[0]:+.4f}), pi_hydro = {v[1]:.4f} ({v[1]-base[1]:+.4f}), Neff = {v[2]:.3f}")
dpi4 = {s: op4[3][s] - base[3][s] for s in base[3]}
print("  change in pi by sector under (ii)+(iii):", {s: round(v, 4) for s, v in dpi4.items()})
print(f"  hydrocarbons domestic inputs per unit of output: {Z12[:, iH].sum()/SF19.loc['Hydrocarbons','output']:.3f} -> "
      f"{Z3[:, iH].sum()/SF19.loc['Hydrocarbons','output']:.3f}")
imp_sum = sum(0.10 * SF19.loc[s, "imp_int_inputs"] for s in lab12 if s in SF19.index and s != "Hydrocarbons"
              and np.isfinite(SF19.loc[s, "imp_int_inputs"]))
print(f"  size of (ii): {imp_sum:.0f} million Bs = {imp_sum/SF19.loc['Hydrocarbons','output']:.3f} of hydrocarbons output")
# individual scale: Chapter 18 formulas
def ratio_gender(lg, aF=0.454, aM=1.0, c=1.0):
    return (aF + lg * c) / (aM + lg * c)
def ratio_minority(p, l):
    return (p + l * (1 - p)) / ((1 - p) + l * p)
IND = {}
for lg, lr in ((0.30, 0.35), (0.50, 0.50), (0.70, 0.70)):
    IND[(lg, lr)] = (ratio_gender(lg), ratio_minority(0.41, lr), (1 - lg) / (1 + lg), 2 * lg / (1 + lg))
    print(f"  lambda_g = {lg:.2f}, lambda_r = {lr:.2f}: women/men pi ratio = {IND[(lg,lr)][0]:.3f}; "
          f"indigenous/non-indig. homophily ratio = {IND[(lg,lr)][1]:.3f}; aspiration-gap factor (1-l)/(1+l) = {IND[(lg,lr)][2]:.3f}; "
          f"Gamma (two equal groups, a=c) = {IND[(lg,lr)][3]:.3f}")
assert round(IND[(0.30, 0.35)][0], 2) == 0.58 and round(IND[(0.30, 0.35)][1], 2) == 0.84
print(f"individual-scale two-group Gamma change, lambda_g 0.3 -> 0.5: {IND[(0.50,0.50)][3]-IND[(0.30,0.35)][3]:+.3f}")
rho = 0.85
print(f"formation half-life ln2/ln(1/rho) at rho = {rho}: {np.log(2)/np.log(1/rho):.2f} generations")
vel = pd.read_csv(os.path.join(PROC, "gamma_velocity.csv"))
vB = vel[(vel.country == "BOL") & (vel.vintage == "INE base 1990, 1990-2013")].iloc[0]
print(f"Bolivia structural velocity 1990-2013 = {vB.v_endpoint_per_decade:+.3f} per decade")
gain = sorted([s for s, v in dpi4.items() if v > 0]); lose = sorted([s for s, v in dpi4.items() if v < 0])
pi_gain = sum(base[3][s] for s in gain); pi_lose = sum(base[3][s] for s in lose)
print(f"(ii)+(iii) gainers {gain}: current pi sum = {pi_gain:.3f}; losers {lose}: current pi sum = {pi_lose:.3f}")
g3 = ["Agriculture", "Hydrocarbons", "Manufacturing"]
print(f"pi of the three main gainers (agr, hydro, mfg) = {sum(base[3][s] for s in g3):.3f}")
assert round(pi_gain, 3) == 0.327 and round(pi_lose, 3) == 0.673
print("\nall Chapter 19 solution numbers computed; asserts passed")
