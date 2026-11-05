"""Solutions manual, Chapter 18 (The Individual Scale): computes and prints every number
used in solutions/sol18.tex, with asserts on the key results.
Run from the repo root:  python3 code/solutions/sol18.py
numpy + pandas (pandas reads data/processed only, Exercise 18.2)."""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=120)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROC = os.path.join(ROOT, "data", "processed")


def gap(W):
    return 1 - np.sort(np.abs(np.linalg.eigvals(W)))[-2]


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    p = np.real(vl[:, np.argmin(np.abs(ev - 1))])
    return p / p.sum()


def rownorm(K):
    return K / K.sum(1, keepdims=True)


# =====================================================================
print("== Exercise 18.1: assortative mixing and Gamma (3 men, 3 women) ==")
grp = np.array([0, 0, 0, 1, 1, 1])
same = grp[:, None] == grp[None, :]
G1 = {}
for lam in (1.0, 0.7, 0.4, 0.1, 0.0):
    K = np.where(same, 1.0, lam)                     # K^max = 1 for all pairs, self-loops included
    W = rownorm(K)
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    g = gap(W)
    formula = 2 * lam / (1 + lam)
    Kn = K.copy(); np.fill_diagonal(Kn, 0)           # variant without self-loops
    gn = gap(rownorm(Kn))
    evn = np.sort(np.real(np.linalg.eigvals(rownorm(Kn))))[::-1]
    G1[lam] = (g, gn)
    print(f"lambda = {lam:.1f}: |eigs| = {np.round(ev, 4)}, Gamma = {g:.4f} (2l/(1+l) = {formula:.4f}); "
          f"no self-loops: Gamma = {gn:.4f}, eigs = {np.round(evn, 4)}")
    assert abs(g - formula) < 1e-12
lams = np.linspace(0, 1, 11)
gl = 2 * lams / (1 + lams)
print("Gamma on a grid:", np.round(gl, 4))
print(f"slope at 0 = 2, at 1 = 0.5; chord value at 0.5 = 0.5 vs Gamma(0.5) = {2*0.5/1.5:.4f} (concave)")
# without self-loops: closed form lambda_2 = (2 - 3l)/(2 + 3l), other eigenvalue -1/(2+3l)
for lam in (1.0, 0.7, 0.4, 0.1):
    l2 = (2 - 3 * lam) / (2 + 3 * lam); l3 = -1 / (2 + 3 * lam)
    assert abs(G1[lam][1] - (1 - max(abs(l2), abs(l3)))) < 1e-12

# =====================================================================
print("\n== Exercise 18.2: Bolivia's gender power gap (ILOSTAT 2019 tabulation of the EH; INE 2019) ==")
EMP_MAP = {"Hydrocarbons and mining": "B", "Public administration": "O", "Manufacturing": "C",
           "Construction": "F", "Commerce": "G", "Transport": "H", "Financial services": "K",
           "Agriculture": "A", "Other services": "DEJLMNPQ", "Personal services": "IRS",
           "Domestic services": "T"}
ilo = pd.read_csv(os.path.join(PROC, "ilostat_bolivia_emp_eco_sex.csv"))
ilo = ilo[ilo.year == 2019]
I = {(r.sex, r.eco): r.value_thousands for r in ilo.itertuples()}
E = {s: sum(I[("SEX_T", "ECO_ISIC4_" + c)] for c in cs) for s, cs in EMP_MAP.items()}
F = {s: sum(I[("SEX_F", "ECO_ISIC4_" + c)] for c in cs) for s, cs in EMP_MAP.items()}
fsh = {s: F[s] / E[s] for s in EMP_MAP}
Ftot, Mtot = I[("SEX_F", "ECO_ISIC4_TOTAL")], I[("SEX_M", "ECO_ISIC4_TOTAL")]
pis = pd.read_csv(os.path.join(PROC, "bolivia_io_pi.csv"))
pis = pis[(pis.base == 2017) & (pis.year == 2019) & (pis.aggregation == "book12") & (pis.flows == "domestic")
          & (pis.normalization == "input_shares")].set_index("sector").pi
lk = pd.read_csv(os.path.join(PROC, "bolivia_io_linkages.csv"))
lk = lk[(lk.year == 2019) & (lk.base == 2017) & (lk.aggregation == "book12") & (lk.flows == "domestic")].set_index("sector")
pi_e = {s: (pis["Hydrocarbons"] + pis["Mining"]) if s == "Hydrocarbons and mining" else pis.get(s, 0.0) for s in EMP_MAP}
bl = {}
for s in EMP_MAP:
    if s == "Hydrocarbons and mining":
        w = lk.loc[["Hydrocarbons", "Mining"], "output"]
        bl[s] = float((lk.loc[["Hydrocarbons", "Mining"], "BL"] * w).sum() / w.sum())
    else:
        bl[s] = float(lk.loc[s, "BL"])
cls = {s: ("male" if fsh[s] < 0.4 else "female" if fsh[s] > 0.6 else "mixed") for s in EMP_MAP}
print(f"{'sector':24s} female share  class   BL     pi")
for s in EMP_MAP:
    print(f"  {s:24s} {fsh[s]:.3f}   {cls[s]:6s} {bl[s]:.3f}  {pi_e[s]:.3f}")
GR = {}
for c in ("male", "mixed", "female"):
    ss = [s for s in EMP_MAP if cls[s] == c]
    simple = np.mean([bl[s] for s in ss])
    empw = sum(bl[s] * E[s] for s in ss) / sum(E[s] for s in ss)
    GR[c] = (ss, simple, empw)
    print(f"{c:6s}-dominated: {ss}; mean BL = {simple:.3f} (employment-weighted {empw:.3f})")


def ratio_FM(fs):
    piF = sum(pi_e[s] * fs[s] for s in EMP_MAP) / sum(fs[s] * E[s] for s in EMP_MAP)
    piM = sum(pi_e[s] * (1 - fs[s]) for s in EMP_MAP) / sum((1 - fs[s]) * E[s] for s in EMP_MAP)
    return piF / piM


r0 = ratio_FM(fsh)
print(f"pi_F / pi_M (sector pi shared equally among workers; domestic service gets 0) = {r0:.4f}")
assert round(r0, 2) == 0.86
male = [s for s in EMP_MAP if cls[s] == "male"]
fs_male = sum(F[s] for s in male) / sum(E[s] for s in male)
print(f"female share in male-dominated sectors {male}: {fs_male:.3f}")
# counterfactual: female share 40% in every male-dominated sector, sector employment fixed;
# the extra women come from (and the displaced men go to) the other sectors in proportion to employment
extra = sum(0.4 * E[s] - F[s] for s in male if fsh[s] < 0.4)
others = [s for s in EMP_MAP if s not in male]
fs1 = dict(fsh)
for s in male:
    fs1[s] = max(0.4, fsh[s])
for s in others:
    fs1[s] = (F[s] - extra * E[s] / sum(E[t] for t in others)) / E[s]
Fmap = sum(F.values())
assert abs(sum(fs1[s] * E[s] for s in EMP_MAP) - Fmap) < 1e-6
r1 = ratio_FM(fs1)
# variant: the women come only from female-dominated sectors
fem_s = [s for s in EMP_MAP if cls[s] == "female"]
fs2 = dict(fsh)
for s in male:
    fs2[s] = 0.4
for s in fem_s:
    fs2[s] = (F[s] - extra * E[s] / sum(E[t] for t in fem_s)) / E[s]
assert abs(sum(fs2[s] * E[s] for s in EMP_MAP) - Fmap) < 1e-6
r2 = ratio_FM(fs2)
print(f"women moved into male-dominated sectors: {extra:.1f} thousand ({extra/Fmap:.2%} of female employment)")
print(f"counterfactual ratio (from all other sectors) = {r1:.4f}; gap 1 - r closes from {1-r0:.4f} to {1-r1:.4f} "
      f"({(r1-r0)/(1-r0):.1%} of the gap)")
print(f"counterfactual ratio (from female-dominated sectors only) = {r2:.4f} ({(r2-r0)/(1-r0):.1%} of the gap)")
# which sector drives it: contribution of each male-dominated sector
for s in male:
    fsx = dict(fs1); fsx[s] = fsh[s]
    print(f"  without the {s} change: ratio {ratio_FM(fsx):.4f}")
G2 = dict(r0=r0, r1=r1, r2=r2, extra=extra)

# =====================================================================
print("\n== Exercise 18.3: intersectionality, four group-agents ==")
lg, lr = 0.3, 0.35
attrs = [(0, 0), (1, 0), (0, 1), (1, 1)]           # (female?, indigenous?)
phi = np.array([[(1 if a[0] == b[0] else lg) * (1 if a[1] == b[1] else lr) for b in attrs] for a in attrs])
print("phi =\n", phi)
Wa = rownorm(phi * 0.25)
pa = perron(Wa)
print(f"(b) W^max = 11'/4: row sums of phi = {phi.sum(1)}; pi = {pa}")
assert np.allclose(pa, 0.25) and np.allclose(phi.sum(1), (1 + lg) * (1 + lr))
nvec = np.array([0.295, 0.295, 0.205, 0.205])
Kc = phi * nvec[None, :]
Wc = rownorm(Kc)
pc = perron(Wc)
d = Kc.sum(1)
print(f"(c) pi = {pc}; per-capita pi/n = {pc/nvec}; relative to agent 1: {(pc/nvec)/(pc[0]/nvec[0])}")
assert np.allclose(pc, nvec * d / (nvec @ d))             # pi_i proportional to n_i d_i (reversible)
assert abs(pc[0] / nvec[0] - pc[1] / nvec[1]) < 1e-12 and abs(pc[2] / nvec[2] - pc[3] / nvec[3]) < 1e-12
r41 = (pc[3] / nvec[3]) / (pc[0] / nvec[0])
p_ = 0.41
fml = (p_ + lr * (1 - p_)) / ((1 - p_) + lr * p_)
print(f"(d) (pi4/n4)/(pi1/n1) = {r41:.4f}; lambda_g lambda_r = {lg*lr:.3f}; formula (p = 0.41, lambda_r) = {fml:.4f}")
assert abs(r41 - fml) < 1e-12 and round(r41, 2) == 0.84
print(f"Gamma of W (c) = {gap(Wc):.4f}; eigenvalues = {np.round(np.sort(np.real(np.linalg.eigvals(Wc)))[::-1], 4)}")
G3 = dict(r41=r41, pc=pc)

# =====================================================================
print("\n== Exercise 18.4: the formation trap with endogenous returns ==")
rho, mu, thb, rs = 0.85, 0.3, 1.0, 0.1


def steady4(lam, cD=0.03):
    w = 1 / (1 + lam)
    # theta = rho theta + (1-rho) thb + mu (w r_own + (1-w) r_other - rs), r_A = 0.05 + .1 th_A, r_D = cD + .1 th_D
    M = np.array([[1 - rho - mu * w * 0.1, -mu * (1 - w) * 0.1],
                  [-mu * (1 - w) * 0.1, 1 - rho - mu * w * 0.1]])
    c = np.array([(1 - rho) * thb + mu * (w * 0.05 + (1 - w) * cD - rs),
                  (1 - rho) * thb + mu * (w * cD + (1 - w) * 0.05 - rs)])
    th = np.linalg.solve(M, c)
    J = np.array([[rho + mu * w * 0.1, mu * (1 - w) * 0.1], [mu * (1 - w) * 0.1, rho + mu * w * 0.1]])
    return th, J


def dtheta_closed(lam):
    q = (1 - lam) / (1 + lam)
    return 0.02 * mu * q / ((1 - rho) - 0.1 * mu * q)


for lam in (0.25, 0.5):
    th, J = steady4(lam)
    rA, rD = 0.05 + 0.1 * th[0], 0.03 + 0.1 * th[1]
    print(f"lambda = {lam}: own weight {1/(1+lam):.4f}; theta*_A = {th[0]:.5f}, theta*_D = {th[1]:.5f}, "
          f"d theta* = {th[0]-th[1]:.5f} (closed form {dtheta_closed(lam):.5f}); r_A = {rA:.5f}, r_D = {rD:.5f}; "
          f"rho(J) = {max(abs(np.linalg.eigvals(J))):.4f}")
    assert abs(th[0] - th[1] - dtheta_closed(lam)) < 1e-12
d25, d50 = dtheta_closed(0.25), dtheta_closed(0.5)
print(f"d theta* shrinks by {d25-d50:.5f} ({1-d50/d25:.1%})")
# fixed-returns benchmark (Prop. ch18:thm:attribute_trap): M q (r_A - r_D) with r_A - r_D = 0.02
Mf = mu / (1 - rho)
print(f"formation multiplier M = {Mf:.3f}; fixed-returns gap at lambda = 0.25: {Mf*0.6*0.02:.5f}; "
      f"amplification from endogenous returns 1/(1 - 0.1 M q) = {1/(1-0.1*Mf*0.6):.4f}")
# (d) threshold
qstar = 0.01 * (1 - rho) / (0.02 * mu + 0.01 * 0.1 * mu)
lstar = (1 - qstar) / (1 + qstar)
print(f"d theta* < 0.01 iff q < {qstar:.6f} iff lambda > {lstar:.4f}")
assert abs(dtheta_closed(lstar) - 0.01) < 1e-12 and dtheta_closed(lstar + 1e-6) < 0.01
# (e) identical opportunity structure
for lam in (0.0, 0.25, 0.5, 0.9):
    th, _ = steady4(lam, cD=0.05)
    assert abs(th[0] - th[1]) < 1e-12
print("with r^D = 0.05 + 0.1 theta^D: d theta* = 0 for lambda in {0, 0.25, 0.5, 0.9}")
th_e, _ = steady4(0.25, cD=0.05)
print(f"common steady state theta* = {th_e[0]:.5f}")
G4 = dict(d25=d25, d50=d50, lstar=lstar)

# =====================================================================
print("\n== Exercise 18.5: the cost of discrimination as a share of GDP ==")
kappa = 0.5


def U(x):
    return x - kappa / 2 * x ** 2


def per_capita_loss(r):
    xM, xF = 2 / (1 + r), 2 * r / (1 + r)               # Omega/N = 1, Phi/N = 1: x = 1 + (N pi - 1)
    N2V = ((1 - r) / (1 + r)) ** 2
    L = kappa / 2 * N2V
    assert abs(L - (U(1) - (U(xM) + U(xF)) / 2)) < 1e-12 and kappa * max(xM, xF) < 1
    return L, xM, xF


r_03 = (0.454 + 0.3) / (1 + 0.3)                  # Prop. homophily_power (iii): a_F = 0.454 a_M, c = a_M
r_1 = (0.454 + 1.0) / (1 + 1.0)
print(f"calibrated ratios: r(0.3) = {r_03:.4f}, r(1) = {r_1:.4f}")
L58, xM, xF = per_capita_loss(r_03)
L73, _, _ = per_capita_loss(r_1)
dL = L58 - L73
Psi = 2.3
print(f"L_power/N: r = {r_03:.3f}: {L58:.5f} (x_M = {xM:.4f}, x_F = {xF:.4f}); r = {r_1:.3f}: {L73:.5f}; "
      f"discrimination (power channel) = {dL:.5f}; with (1 + Psi): {(1+Psi)*dL:.5f}")
assert round(L58, 4) == 0.0177 and round(L73, 4) == 0.0062 and round(dL, 4) == 0.0114
# (i) linear money metric: divide by marginal utility at the mean, U'(1) = 1 - kappa
mu_bar = 1 - kappa
lin = (1 + Psi) * dL / mu_bar
lin_pow = dL / mu_bar
# (ii) exact equally-distributed-equivalent income: U(x_ede) = U(1) - L/N
def ede(Lpc):
    # solve x - k/2 x^2 = U(1) - Lpc on x < 1/kappa
    a, b, c = -kappa / 2, 1.0, -(U(1) - Lpc)
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)
for nm, l1, l0 in (("power channel only", L58, L73), ("with (1+Psi)", (1 + Psi) * L58, (1 + Psi) * L73)):
    e1, e0 = ede(l1), ede(l0)
    print(f"  {nm}: EDE income {e1:.5f} at lambda = 0.3 vs {e0:.5f} at lambda = 1 -> loss = {e0-e1:.5f} of per-capita GDP")
EDE = {nm: ede(l0) - ede(l1) for nm, l1, l0 in (("pow", L58, L73), ("tot", (1 + Psi) * L58, (1 + Psi) * L73))}
print(f"  linear approximation: power channel {lin_pow:.4f}, total {lin:.4f} of per-capita GDP")
G5 = dict(dL=dL, lin=lin, lin_pow=lin_pow, EDE=EDE)

# =====================================================================
print("\n== Exercise 18.6: three interventions (illustrative design) ==")
aF, aM, c = 0.454, 1.0, 1.0


def pratio(lam):
    return (aF + lam * c) / (aM + lam * c)


def gam_layer(lam):                                    # Prop 18.? (a) with a = c: Gamma = 2 lam / (1 + lam)
    return 2 * lam / (1 + lam)


def asp_factor(lam):                                   # (1 - lam)/(1 + lam): aspiration-gap factor
    return (1 - lam) / (1 + lam)


base = 0.3
print(f"base lambda = 0.3: pi_F/pi_M = {pratio(0.3):.4f}; Gamma(layer) = {gam_layer(0.3):.4f}; aspiration factor = {asp_factor(0.3):.4f}")
assert round(pratio(0.3), 2) == 0.58 and round(pratio(1.0), 2) == 0.73
INT = {"quota (W^P)": ("P", 0.45), "mentoring (W^I)": ("I", 0.40), "media campaign (W^I)": ("I", 0.35)}
rows6 = {}
for nm, (layer, lam1) in INT.items():
    dG = gam_layer(lam1) - gam_layer(base)
    if layer == "P":
        E_ = pratio(lam1) - pratio(base)                # enforcement: change in power ratio at fixed preferences
        tot = (1 + Psi) * E_
        dasp = 0.0
        rows6[nm] = (lam1, dG, E_, tot, dasp)
        print(f"  {nm}: lambda_g^P 0.30 -> {lam1:.2f}; dGamma(W^P) = {dG:+.4f}; enforcement d(pi_F/pi_M) = {E_:+.4f}; "
              f"total (1 + Psi) E = {tot:+.4f}")
    else:
        dasp = asp_factor(lam1) - asp_factor(base)     # change in the aspiration-gap factor (relative)
        rel = dasp / asp_factor(base)
        rows6[nm] = (lam1, dG, 0.0, None, rel)
        print(f"  {nm}: lambda_g^I 0.30 -> {lam1:.2f}; dGamma(W^I) = {dG:+.4f}; enforcement = 0 (Psi undefined); "
              f"aspiration gap changes by {rel:+.1%}")
# power calculation for a cluster-randomized evaluation of the mentoring programme
z = 1.959964 + 0.841621                                # alpha = 0.05 two-sided, power 0.80
p0, mde, m, icc = 0.12, 0.04, 20, 0.05
n_ind = z ** 2 * 2 * p0 * (1 - p0) / mde ** 2           # per arm, individual randomization (approx., equal p)
deff = 1 + (m - 1) * icc
n_cl = n_ind * deff
J_arm = int(np.ceil(n_cl / m))
print(f"power calc: z-sum = {z:.4f}; n per arm (individual) = {n_ind:.1f}; design effect = {deff:.2f}; "
      f"n per arm (clustered) = {n_cl:.1f}; clusters per arm (m = {m}) = {J_arm}")
G6 = dict(n_ind=n_ind, deff=deff, n_cl=n_cl, J=J_arm, rows=rows6)
print("\nall Chapter 18 solution checks passed")
