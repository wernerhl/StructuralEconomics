"""Solutions manual, Chapter 22 (Measuring the Utility Network).  Computes and prints
every number used in solutions/sol22.tex, with asserts.

Illustrative inputs: the chapter's transfer matrix (Example ch22:ex:transfers), the
per-capita alignment positions (Example ch22:ex:alignment), Chapter 21's audience
table, Chapter 5's layer weights and Chapter 2's 4-agent calibration.  Sourced: WDI
remittances and Gini (data/processed, retrieved 2026-11-05).
Run from the repo root:  python3 code/solutions/sol22.py   (numpy + pandas).
"""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=150)
PROC = os.path.join("data", "processed")


def perron(W):
    w, v = np.linalg.eig(W.T)
    q = np.real(v[:, np.argmin(abs(w - 1))])
    return q / q.sum()


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def neff(q):
    return 1 / np.sum(q ** 2)


def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])


G4 = ['Urban NI', 'Urban Indig', 'Rural NI', 'Rural Indig']
raw = np.array([[.45, .08, .05, .02], [.10, .40, .05, .25], [.08, .03, .50, .10], [.03, .15, .07, .55]])
TU = raw / raw.sum(1, keepdims=True)                   # Example ch22:ex:transfers (rows: senders)
FACT = pd.read_csv(os.path.join(PROC, "bolivia_facts.csv")).set_index("fact_id")["value"]
WDI = pd.read_csv(os.path.join(PROC, "wdi_long.csv"))

# =====================================================================
print("=" * 70); print("Exercise 22.1  Transfer-revealed W^U"); print("=" * 70)
print("transfer-share matrix:\n", np.round(TU, 4))
evU = eig_sorted(TU); GU = gap(TU); piU = perron(TU)
print("eigenvalues:", np.round(evU, 4))
print(f"Gamma = {GU:.4f}; pi = {np.round(piU,4)}; Neff = {neff(piU):.3f} of 4 (Neff/n = {neff(piU)/4:.3f})")
assert round(GU, 3) == 0.300 and round(neff(piU), 2) == 3.85
NIc, INc = [0, 2], [1, 3]; URc, RUc = [0, 1], [2, 3]
cross_e = [TU[i, j] for i in range(4) for j in range(4) if (i in NIc) != (j in NIc)]
within_e = [TU[i, j] for i in range(4) for j in range(4) if (i in NIc) == (j in NIc)]
lam_e = np.mean(cross_e) / np.mean(within_e)
cross_g = [TU[i, j] for i in range(4) for j in range(4) if (i in URc) != (j in URc)]
within_g = [TU[i, j] for i in range(4) for j in range(4) if (i in URc) == (j in URc)]
lam_g = np.mean(cross_g) / np.mean(within_g)
print(f"lambda^U_ethnicity = mean of 8 cross cells {np.mean(cross_e):.4f} / mean of 8 within cells {np.mean(within_e):.4f} = {lam_e:.3f}")
print(f"lambda^U_geography = {np.mean(cross_g):.4f} / {np.mean(within_g):.4f} = {lam_g:.3f}")
assert round(lam_e, 2) == 0.20 and round(lam_g, 2) == 0.28
print(f"column sums (incoming weight): {np.round(TU.sum(0),4)}")
rem, rem_bs = FACT["BOL_remittances_2019"], FACT["BOL_remittances_2019_Bs"]
print(f"remittances received 2019: US$ {rem:.3f} bn = {rem_bs:.2f} bn Bs; {FACT['BOL_remittances_pct_gdp_2019']:.2f}% of GDP "
      f"(peak {FACT['BOL_remittances_pct_gdp_2007']:.2f}% in 2007)")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.2  Directed dictator game: treatments and sample size"); print("=" * 70)
z, zb = 1.959964, 0.841621
n_treat = 2 * 2 * 2
print(f"recipient profiles: 2 (ethnicity) x 2 (community) x 2 (gender) = {n_treat}, plus anonymous stranger = {n_treat+1} decisions")
sd_d, delta = 0.15, 0.05
n_main = ((z + zb) * sd_d / delta) ** 2
sd_int = 0.25
n_int = ((z + zb) * sd_int / delta) ** 2
m_c, icc = 10, 0.05
deff = 1 + (m_c - 1) * icc
n_bs = 2 * (z + zb) ** 2 * 0.20 ** 2 / delta ** 2
print(f"within-subject main effect (delta = {delta}, SD of paired difference {sd_d}): n = {n_main:.1f}")
print(f"two-way interaction contrast (SD {sd_int}): n = {n_int:.1f}; design effect (m = {m_c}, ICC = {icc}) = {deff:.2f} -> {np.ceil(n_int)*deff:.0f}")
print(f"between-subject alternative (SD of shares 0.20): {n_bs:.0f} per arm x {n_treat+1} arms = {np.ceil(n_bs)*(n_treat+1):.0f}")
N22 = dict(main=n_main, inter=n_int, deff=deff, total=np.ceil(n_int) * deff, bs=n_bs, bs_tot=np.ceil(n_bs) * (n_treat + 1))
assert np.ceil(n_main) == 71 and np.ceil(n_int) == 197
# delta method for W = s/(1-s)
for s_ in (0.12, 0.25, 0.38):
    se_s = 0.20 / np.sqrt(np.ceil(n_int))
    print(f"  s = {s_:.2f}: W = {s_/(1-s_):.3f}, SE(s) = {se_s:.4f} -> SE(W) = {se_s/(1-s_)**2:.4f}")
SEW = {s_: (s_ / (1 - s_), (0.20 / np.sqrt(np.ceil(n_int))) / (1 - s_) ** 2) for s_ in (0.12, 0.25, 0.38)}

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.3  Fiscal W^U (illustrative incidence of the chapter)"); print("=" * 70)
inc = {"poor households": 0.12, "non-poor households": 0.08, "reinvested in the sector": 0.25}
print(f"illustrative IDH incidence {inc}; unallocated = {1-sum(inc.values()):.2f}")
UNALLOC = 1 - sum(inc.values())

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.4  The alignment test on four groups"); print("=" * 70)
S6 = np.array([[.45, .15, .02, .30, .08], [.40, .18, .03, .28, .11], [.20, .30, .15, .20, .15],
               [.25, .20, .10, .15, .30], [.08, .22, .30, .10, .30], [.05, .15, .25, .08, .47]])
p6 = np.array([0.17, 0.18, 0.25, 0.10, 0.15, 0.15])
m6 = p6 @ S6
W6 = S6 @ np.diag(1 / m6) @ S6.T @ np.diag(p6)
map6 = [[0, 1], [2], [3], [4, 5]]                        # Urban NI, Urban Indig, Rural NI, Rural Indig
pi6 = perron(W6)
WI4 = np.array([[pi6[A] @ W6[np.ix_(A, B)].sum(1) / pi6[A].sum() for B in map6] for A in map6])
piI = perron(WI4)
p21 = np.array([p6[A].sum() for A in map6])
assert np.allclose(piI, p21)
print("W^I aggregated to four groups (pi-weighted):\n", np.round(WI4, 4))
print(f"pi^I = {np.round(piI,4)} = Chapter 21 population shares aggregated; Gamma(W^I_4) = {gap(WI4):.4f}")
p22 = np.array([0.35, 0.25, 0.15, 0.25])
print(f"Chapter 22 population shares: {p22}; Chapter 21 aggregated: {p21}")
# illustrative production positions (Example ch22:ex:alignment)
rP = np.array([1.50, 0.80, 1.00, 0.46]); rP[2] = (1 - (p22 * rP)[[0, 1, 3]].sum()) / p22[2]
rI_ex = np.array([1.40, 0.90, 1.20, 0.35]); rI_ex[2] = (1 - (p22 * rI_ex)[[0, 1, 3]].sum()) / p22[2]
rU_ex = np.array([1.35, 1.00, 1.15, 0.30]); rU_ex[2] = (1 - (p22 * rU_ex)[[0, 1, 3]].sum()) / p22[2]
piP = rP * p22
print(f"pi^P (illustrative, from r^P = {np.round(rP,3)}): {np.round(piP,4)}")
print(f"pi^U (transfer matrix): {np.round(piU,4)}; per capita r^U = pi^U/p22 = {np.round(piU/p22,3)}")
print(f"per capita r^I = pi^I/p21 = {np.round(piI/p21,3)} (identically 1); with p22: {np.round(piI/p22,3)}")
C_lev = {"P,U": corr(piP, piU), "P,I": corr(piP, piI), "I,U": corr(piI, piU)}
C_pc = {"P,U": corr(rP, piU / p22)}
print("correlations of pi levels:", {k: round(v, 3) for k, v in C_lev.items()})
print("correlation of per-capita positions r^P, r^U:", {k: round(v, 3) for k, v in C_pc.items()})
C_ex = {"P,U": corr(rP * p22, rU_ex * p22), "P,I": corr(rP * p22, rI_ex * p22), "I,U": corr(rI_ex * p22, rU_ex * p22)}
print("Example 22.alignment's own correlations:", {k: round(v, 2) for k, v in C_ex.items()})
assert C_pc["P,U"] < 0

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.5  Misalignment as policy"); print("=" * 70)
alpha = np.ones(3) / 3
def comp_r(rp, ri, ru):
    return alpha @ np.vstack([rp, ri, ru])
r0 = comp_r(rP, rI_ex, rU_ex)
pi0 = r0 * p22
q = np.array([0, p22[1], 0, p22[3]]); q = q / q.sum()        # indigenous population, normalised
print(f"rank-one layers, equal weights. before: composite r = {np.round(r0,3)}, Var(pi) = {pi0.var():.5f}, Neff = {neff(pi0):.3f}")
MIS = {}
for s_ in (0.25, 0.50):
    piU2 = (1 - s_) * rU_ex * p22 + s_ * q
    rU2 = piU2 / p22
    r1 = comp_r(rP, rI_ex, rU2); pi1 = r1 * p22
    MIS[s_] = (rU2, r1, pi1)
    print(f"s = {s_}: r^U = {np.round(rU2,3)}; composite r = {np.round(r1,3)}; Var(pi) = {pi1.var():.5f}; Neff = {neff(pi1):.3f}; "
          f"corr(pi^P, pi^U) = {corr(rP*p22, piU2):+.3f}")
    assert np.isclose(r1[3] - r0[3], alpha[2] * (rU2[3] - rU_ex[3]))
print(f"composite change = alpha_U x change in r^U (rank-one): verified")
# full matrices: W^U transfers (4x4), W^I aggregated media (4x4), W^P rank-one
WP = np.outer(np.ones(4), piP)
def composite(Wu):
    Wc = (WP + WI4 + Wu) / 3
    return Wc, perron(Wc)
Wc0, pc0 = composite(TU)
TU2 = TU.copy()
for i in (1, 3):                     # indigenous rows move 0.10 from NI recipients to indigenous recipients
    moved = 0.10
    ni_mass = TU2[i, [0, 2]].sum()
    TU2[i, [0, 2]] -= moved * TU2[i, [0, 2]] / ni_mass
    TU2[i, [1, 3]] += moved * TU2[i, [1, 3]] / TU2[i, [1, 3]].sum()
assert np.allclose(TU2.sum(1), 1) and (TU2 >= 0).all()
Wc1, pc1 = composite(TU2)
print(f"full layers: pi^U {np.round(piU,4)} -> {np.round(perron(TU2),4)}; Gamma(W^U) {gap(TU):.4f} -> {gap(TU2):.4f}")
print(f"composite pi: {np.round(pc0,4)} -> {np.round(pc1,4)}; per capita (p22): {np.round(pc0/p22,3)} -> {np.round(pc1/p22,3)}")
print(f"composite Gamma {gap(Wc0):.4f} -> {gap(Wc1):.4f}; Neff {neff(pc0):.3f} -> {neff(pc1):.3f}; Var {pc0.var():.5f} -> {pc1.var():.5f}")
FULL = (pc0, pc1, gap(Wc0), gap(Wc1), perron(TU2), gap(TU2))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.6  Inequality inputs (WDI Gini)"); print("=" * 70)
GINI = {}
for c in ("BOL", "CHL", "URY"):
    g = WDI[(WDI.country_iso3 == c) & (WDI.indicator == "SI.POV.GINI") & (WDI.year <= 2019)].sort_values("year").iloc[-1]
    GINI[c] = (int(g.year), float(g.value))
    print(f"{c}: Gini {g.value:.1f} ({int(g.year)})")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 22.7  The complete measurement (illustrative)"); print("=" * 70)
aU, aP, aI = 0.2, 0.5, 0.3                               # Chapter 5, Example ch05:ex:LaPaz
Wc = aU * TU + aP * WP + aI * WI4
assert np.allclose(Wc.sum(1), 1)
pc = perron(Wc); Gc = gap(Wc)
print("composite W:\n", np.round(Wc, 4))
print(f"composite Gamma = {Gc:.4f}; pi = {np.round(pc,4)}; Neff = {neff(pc):.3f}; eigenvalues {np.round(eig_sorted(Wc),4)}")
N, Omega, kappa, dbar, gamU = 4, 20.0, 0.1, 3.0, 0.3
Phi = Omega - N * dbar
U = lambda x: x - kappa / 2 * x ** 2
xB = Omega / N + (pc - 1 / N) * Phi
assert (kappa * xB < 1).all()
WU0 = TU.copy(); np.fill_diagonal(WU0, 0); WU0 = WU0 / WU0.sum(1, keepdims=True)
c = WU0.sum(0)
print(f"x^B = {np.round(xB,4)} (kappa x < 1: max {kappa*xB.max():.3f}); c_j (zero-diagonal W^U column sums) = {np.round(c,4)}, sum {c.sum():.1f}")
Sfun = lambda x: np.sum((1 + gamU * c) * U(x))
# first best: (1 + g c_j)(1 - kappa x_j) = mu, sum x = Omega
wts = 1 + gamU * c
mu = (N - kappa * Omega) / np.sum(1 / wts)
xFB = (1 - mu / wts) / kappa
assert np.isclose(xFB.sum(), Omega) and (kappa * xFB < 1).all()
xE = np.full(N, Omega / N)
L_tot = Sfun(xFB) - Sfun(xB)
L_pow = kappa * N / 2 * Phi ** 2 * pc.var()
L_ext = Sfun(xFB) - Sfun(xE)
L_int = gamU * np.sum(c * (U(Omega / N) - U(xB)))
L_int2 = -gamU * N * (1 - kappa * Omega / N) * Phi * np.mean((c - c.mean()) * (pc - 1 / N)) + gamU * kappa / 2 * Phi ** 2 * np.sum(c * (pc - 1 / N) ** 2)
assert np.isclose(L_tot, L_pow + L_ext + L_int) and np.isclose(L_int, L_int2)
assert np.isclose(L_pow, np.sum(U(Omega / N) - U(xB)))
print(f"x^FB = {np.round(xFB,4)}")
print(f"L = {L_tot:.5f} = L_power {L_pow:.5f} + L_externality {L_ext:.5f} + L_interaction {L_int:.5f}")
fr = {k: v / L_tot for k, v in (("power", L_pow), ("externality", L_ext), ("interaction", L_int))}
print("fractions:", {k: f"{v:.3f}" for k, v in fr.items()})
print(f"Cov(c, pi) = {np.mean((c - c.mean()) * (pc - 1/N)):.6f}")
print(f"S(x^FB) = {Sfun(xFB):.4f}, S(x^E) = {Sfun(xE):.4f}, S(x^B) = {Sfun(xB):.4f}")
WD = dict(pc=pc, Gc=Gc, xB=xB, c=c, xFB=xFB, L=L_tot, Lp=L_pow, Le=L_ext, Li=L_int, fr=fr)
print("\nall Chapter 22 solution numbers computed; asserts passed")
