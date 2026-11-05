"""Solutions manual, Chapter 2: every number used in solutions/sol02.tex.
Run from the repository root:  python3 code/solutions/sol02.py
numpy throughout; pandas only to read data/processed (Exercise 2.5)."""
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)


def hdr(s):
    print("\n== " + s)


def U(x, k):
    x = np.asarray(x, float)
    return x - k / 2 * x ** 2


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eigs(W):
    e = np.linalg.eigvals(W)
    return e[np.argsort(-abs(e))]


def gap(W):
    return 1 - abs(eigs(W)[1])


def neff(p):
    return 1 / np.sum(p ** 2)


def var(p):
    return np.mean((p - 1 / len(p)) ** 2)


# ---------------------------------------------------------------- 2.1
hdr("2.1 Doubly stochastic check")
mats = {
    "a": np.ones((3, 3)) / 3,
    "b": np.array([[.6, .3, .1], [.1, .6, .3], [.3, .1, .6]]),
    "c": np.array([[.8, .2], [.5, .5]]),
    "d": np.array([[.4, .3, .3], [.3, .4, .3], [.3, .3, .4]]),
}
for k_, W in mats.items():
    p = perron(W)
    print("(%s) row sums %s col sums %s  pi %s  Var(pi) %.5f  N_eff %.4f  eig %s  Gamma %.4f"
          % (k_, W.sum(1), W.sum(0), p.round(4), var(p), neff(p), np.round(eigs(W), 4), gap(W)))
    assert np.allclose(W.sum(1), 1)
p = perron(mats["c"]); assert np.allclose(p, [5 / 7, 2 / 7])
print("(c) pi = (5/7, 2/7) =", [5 / 7, 2 / 7], " Var = (5/7-1/2)^2 =", (5 / 7 - .5) ** 2)
print("(b) |lambda_2| = |0.6 + 0.3 w + 0.1 w^2|, w = exp(2 pi i/3):",
      abs(.6 + .3 * np.exp(2j * np.pi / 3) + .1 * np.exp(4j * np.pi / 3)).round(4))

# ---------------------------------------------------------------- 2.2
hdr("2.2 Is there a clawback? (Example 2.15: N=4, Omega=20, dbar=3, Phi=8)")
N, Om, dbar, kap = 4, 20., 3., 0.1
Phi = Om - N * dbar


def Wb(b, scale=1.0):
    g = (1 - 2 * b * scale) / 2
    t = (1 - 2 * .02 * scale) / 2
    return np.array([[g, g, b * scale, b * scale], [g, g, b * scale, b * scale],
                     [.02 * scale, .02 * scale, t, t], [.02 * scale, .02 * scale, t, t]])


W0 = Wb(.01); pi0 = perron(W0)
print("pi0 =", pi0.round(4), " Phi =", Phi)
for tau in [0.5, 1, 2, 3]:
    d = np.array([dbar, dbar - tau, dbar + tau, dbar])
    x = d + pi0 * (Om - d.sum())
    print("tau=%.1f  d=%s  x=%s  x3 - x3(0) = %.4f  max kappa x=%.3f" % (tau, d, x.round(4), x[2] - (dbar + pi0[2] * Phi), (kap * x).max()))
    assert np.isclose(x[2] - (dbar + pi0[2] * Phi), tau) and np.isclose(x.sum(), Om)
for scale in [1, 5]:
    print("-- cross-block entries x%d" % scale)
    for b in [.01, .0125, .015, .02]:
        W = Wb(b, scale); p = perron(W); x = dbar + p * Phi
        dpi = 0.5 * np.abs(p - 1 / N).sum()
        e1, e2 = 2 * b * scale, 2 * .02 * scale
        print("b=%.4f  eps1=%.3f eps2=%.3f  pi=%s  block mass %.4f (formula %.4f)  x=%s  Phi*d_pi=%.4f  formula 4(0.02-b)/(b+0.02)=%.4f  Gamma=%.4f  eps1+eps2=%.3f  N_eff=%.3f  L_power=%.4f"
              % (b, e1, e2, p.round(4), p[:2].sum(), e2 / (e1 + e2), x.round(4), Phi * dpi, 4 * (.02 - b) / (b + .02),
                 gap(W), e1 + e2, neff(p), kap * N / 2 * Phi ** 2 * var(p)))
        assert np.isclose(p[:2].sum(), e2 / (e1 + e2)) and np.isclose(Phi * dpi, 4 * (.02 - b) / (b + .02))
        assert np.isclose(gap(W), e1 + e2)
        print("     spectrum:", np.round(eigs(W), 4))

# ---------------------------------------------------------------- 2.3
hdr("2.3 Welfare decomposition: N=3, gamma=0.5, kappa=0.5, Omega=3, Phi=1")
N, gam, kap, Om, Phi = 3, .5, .5, 3., 1.
m = Om / N


def decomposition(W, label):
    p = perron(W)
    c = W.sum(0) - np.diag(W)
    alpha = 1 + gam * c
    H = N / np.sum(1 / alpha)
    Lext = (1 - kap * m) ** 2 / (2 * kap) * np.sum((alpha - H) ** 2 / alpha)
    Lext2 = N * (1 - kap * m) ** 2 / (2 * kap) * gam ** 2 * np.var(c)
    Lpow = kap * N / 2 * Phi ** 2 * var(p)
    cov = np.mean((p - p.mean()) * (c - c.mean()))
    Lint = kap * gam * Phi ** 2 / 2 * np.sum(c * (p - 1 / N) ** 2) - gam * Phi * N * (1 - kap * m) * cov
    xs = m + (1 - kap * m) * (alpha - H) / (kap * alpha)
    xb = m + (p - 1 / N) * Phi
    Ltot = np.sum(alpha * U(xs, kap)) - np.sum(alpha * U(xb, kap))
    print("[%s] pi=%s  Var(pi)=%.6f  c=%s  Var(c)=%.6f  Cov(pi,c)=%.6f  alpha=%s  H=%.6f"
          % (label, p.round(4), var(p), c.round(4), np.var(c), cov, alpha.round(4), H))
    print("     x_social=%s  x_bargain=%s  max kappa x=%.3f" % (xs.round(4), xb.round(4), kap * max(xs.max(), xb.max())))
    print("     L_ext=%.6f (2nd-order %.6f)  L_power=%.6f  L_int=%.6f [first term %.6f, cov term %.6f]  sum=%.6f  direct L=%.6f"
          % (Lext, Lext2, Lpow, Lint, kap * gam * Phi ** 2 / 2 * np.sum(c * (p - 1 / N) ** 2),
             -gam * Phi * N * (1 - kap * m) * cov, Lext + Lpow + Lint, Ltot))
    assert np.isclose(Lext + Lpow + Lint, Ltot) and kap * max(xs.max(), xb.max()) < 1
    return p, c


decomposition(np.ones((3, 3)) / 3, "a: W^CE")
Wd = np.array([[.5, .4, .1], [.2, .5, .3], [.1, .3, .6]])
p, c = decomposition(Wd, "b")
print("column sums incl. diagonal:", Wd.sum(0), " Gamma(W)=%.4f  N_eff=%.4f" % (gap(Wd), neff(p)))

# ---------------------------------------------------------------- 2.4
hdr("2.4 Spectral gap across topologies, N = 8")
N = 8
T = {}
T["CE"] = np.ones((N, N)) / N
R = np.zeros((N, N))
for i in range(N):
    R[i, (i + 1) % N] = R[i, (i - 1) % N] = .5
T["ring"] = R
B = np.zeros((N, N)); B[:4, :4] = .25; B[4:, 4:] = .25
T["blocks"] = B
S = np.zeros((N, N)); S[0, 1:] = 1 / 7; S[1:, 0] = 1
T["star"] = S
for lazy in [False, True]:
    for name, W in T.items():
        M = (np.eye(N) + W) / 2 if lazy else W
        ev = eigs(M)
        n1 = int(np.sum(np.isclose(ev, 1)))
        msg = "%-7s %-6s |lambda| sorted: %s  #eig=1: %d  Gamma=%.4f" % (name, "lazy" if lazy else "", np.round(abs(ev), 4), n1, gap(M))
        if n1 == 1:
            p = perron(M)
            msg += "  pi=%s  N_eff=%.4f  Var(pi)=%.5f" % (p.round(4), neff(p), var(p))
        print(msg)
print("two blocks: N_eff of t-mixture = 4/(t^2+(1-t)^2): t=0 ->", 4 / (0 + 1), " t=1/2 ->", 4 / (.25 + .25))
print("ring eigenvalues cos(2 pi k/8):", np.round(np.cos(2 * np.pi * np.arange(8) / 8), 4))
print("lazy ring lambda_2 = (1+cos(pi/4))/2 =", round((1 + np.cos(np.pi / 4)) / 2, 4))
print("star: pi_h = 1/2, spokes 1/14 =", round(1 / 14, 4), " N_eff = 1/(1/4+7/196) =", round(1 / (.25 + 7 / 196), 4))
assert np.isclose(gap(T["ring"]), 0) and np.isclose(gap(T["star"]), 0) and np.isclose(gap(T["blocks"]), 0)

# ---------------------------------------------------------------- 2.5
hdr("2.5 Block disconnection with data (INE SUT 2019 VA; ILOSTAT 2019 employment)")
f = pd.read_csv("data/processed/bolivia_io_sector_facts.csv")
f = f[(f.base == 2017) & (f.year == 2019)].set_index("book12")
emp = pd.read_csv("data/processed/ilostat_bolivia_emp_eco_sex.csv")
emp = emp[(emp.year == 2019) & (emp.sex == "SEX_T")].set_index("eco").value_thousands
# book12 -> ISIC sections (approximate; public administration and other services merged
# because INE puts non-market education/health in public administration)
groups = {
    "Hydrocarbons+Mining": (["Hydrocarbons", "Mining"], ["B"]),
    "Agriculture": (["Agriculture"], ["A"]),
    "Manufacturing": (["Manufacturing"], ["C"]),
    "Construction": (["Construction"], ["F"]),
    "Commerce": (["Commerce"], ["G"]),
    "Transport": (["Transport"], ["H"]),
    "Financial services": (["Financial services"], ["K"]),
    "Personal services": (["Personal services"], ["I", "R", "S"]),
    "Domestic services": (["Domestic services"], ["T"]),
    "Public adm + other services": (["Public administration", "Other services"], ["D", "E", "J", "L", "M", "N", "O", "P", "Q", "U"]),
}
rows = []
for g, (secs, isic) in groups.items():
    va = f.loc[secs, "va"].sum()                      # million Bs
    n = sum(emp["ECO_ISIC4_" + s] for s in isic)      # thousand workers
    rows.append((g, va, n, va / n * 1000))            # Bs per worker per year
D = pd.DataFrame(rows, columns=["group", "va_mBs", "emp_k", "va_per_worker_Bs"]).set_index("group")
print(D.round(1).to_string())
tot_emp = emp["ECO_ISIC4_TOTAL"]
print("employment total (ILOSTAT, thousands): %.3f  sum of groups: %.3f" % (tot_emp, D.emp_k.sum()))
assert np.isclose(D.emp_k.sum(), tot_emp, atol=0.01)
y = D.va_per_worker_Bs.values; n = D.emp_k.values
ybar = np.sum(n * y) / n.sum()
blk = np.array([g == "Hydrocarbons+Mining" for g in D.index])
y1 = np.sum(n[blk] * y[blk]) / n[blk].sum(); y2 = np.sum(n[~blk] * y[~blk]) / n[~blk].sum()
w1 = n[blk].sum() / n.sum()
between = w1 * (y1 - ybar) ** 2 + (1 - w1) * (y2 - ybar) ** 2
within2 = np.sum(n[~blk] * (y[~blk] - y2) ** 2) / n.sum()
total = np.sum(n * (y - ybar) ** 2) / n.sum()
print("block-2 employment (thousands): %.1f" % n[~blk].sum())
print("mean VA/worker: economy %.0f Bs, block 1 %.0f, block 2 %.0f, ratio %.2f, block-1 employment share %.4f"
      % (ybar, y1, y2, y1 / y2, w1))
print("between var %.4e  within (sector-level) %.4e  total %.4e  between share %.4f  between/within %.4f"
      % (between, within2, total, between / total, between / within2))
assert np.isclose(between + within2, total)
# log version (scale-free)
ly = np.log(y); lbar = np.sum(n * ly) / n.sum()
l1 = ly[blk][0]; l2 = np.sum(n[~blk] * ly[~blk]) / n[~blk].sum()
lb = w1 * (l1 - lbar) ** 2 + (1 - w1) * (l2 - lbar) ** 2
lw = np.sum(n[~blk] * (ly[~blk] - l2) ** 2) / n.sum()
print("log VA/worker: between %.4f within %.4f between share %.4f" % (lb, lw, lb / (lb + lw)))
kap = 1e-6
print("kappa=1e-6 per Bs: max kappa*y = %.3f;  L_block per worker = %.0f ; within per worker = %.0f (utility units)"
      % (kap * y.max(), kap / 2 * between, kap / 2 * within2))
assert kap * y.max() < 1
print("hydrocarbons comp/VA = %.4f, mining comp/VA = %.4f" % (f.loc["Hydrocarbons", "comp_share_of_va"], f.loc["Mining", "comp_share_of_va"]))
comp1 = f.loc[["Hydrocarbons", "Mining"], "comp"].sum() / n[blk].sum() * 1000
print("compensation per worker, block 1: %.0f Bs/yr" % comp1)

# ---------------------------------------------------------------- 2.6
hdr("2.6 El Alto vs Zona Sur")
nA, nS = 5, 3; wA, wS = 800., 4500.; kap = 1e-4
wbar = (nA * wA + nS * wS) / 8
Lb = kap / 2 * (nA * (wA - wbar) ** 2 + nS * (wS - wbar) ** 2)
print("wbar=%.2f  terms %.2f %.2f  L_block=%.4f  kappa*omega: %.2f %.2f" % (wbar, nA * (wA - wbar) ** 2, nS * (wS - wbar) ** 2, Lb, kap * wA, kap * wS))
for eps in [.02, .15]:
    W = np.zeros((8, 8))
    W[:5, :5] = (1 - 3 * eps) / 5; W[:5, 5:] = eps
    W[5:, 5:] = (1 - 5 * eps) / 3; W[5:, :5] = eps
    p = perron(W)
    print("eps=%.2f  spectrum %s  Gamma=%.4f (8 eps=%.2f, 1-|1-8eps|=%.2f)  col sums %s  pi %s  N_eff %.3f"
          % (eps, np.round(eigs(W).real, 4), gap(W), 8 * eps, 1 - abs(1 - 8 * eps), W.sum(0).round(4), p.round(4), neff(p)))
    hl = np.log(2) / -np.log(abs(1 - 8 * eps))
    print("    half-life of cross-block deviations: %.3f periods" % hl)
    assert np.isclose(gap(W), 1 - abs(1 - 8 * eps)) and np.allclose(p, 1 / 8)

# ---------------------------------------------------------------- 2.7
hdr("2.7 Furtado: a fixed point of underdevelopment")
wM, wT, cbar, delta, sig, e0, Imax = 8., 3.2, 3.5, .1, .05, .02, .5


def F(e):
    xT = wT + e * (wM - wT)
    inv = sig * min(max(xT - cbar, 0.), Imax)
    return (1 - delta) * e + delta * e0 + inv


eth = (cbar - wT) / (wM - wT)
# fixed points: (i) low: e = e0 if e0 < eth ; (ii) middle: delta (e - e0) = sig (wT + e dw - cbar)
dw = wM - wT
emid = (sig * (cbar - wT) - delta * e0) / (sig * dw - delta)
ehigh = e0 + sig * Imax / delta
esat = (cbar + Imax - wT) / dw          # coupling at which investment saturates
print("threshold eps_th=%.4f  low fixed point %.4f  middle %.4f  high %.4f  saturation at %.4f" % (eth, e0, emid, ehigh, esat))
for e in [e0, emid, ehigh]:
    print("  F(%.4f) = %.6f" % (e, F(e)))
    assert np.isclose(F(e), e)
print("slopes: below threshold %.2f, between threshold and saturation %.2f, above saturation %.2f"
      % (1 - delta, 1 - delta + sig * dw, 1 - delta))
for e_start in [.05, .09, .095, .15]:
    e = e_start
    for t in range(200):
        e = F(e)
    print("  start %.3f -> after 200 periods %.4f" % (e_start, e))
print("L_block of Example 2.12 (8 agents, means 8 and 3.2, kappa=0.1): %.2f" % (0.1 / 2 * (3 * (8 - 5) ** 2 + 5 * (3.2 - 5) ** 2)))
xT0 = wT + e0 * dw
print("traditional income at the low fixed point: %.4f (subsistence %.2f)" % (xT0, cbar))
print("minimal one-off jump in coupling from e0 to escape: %.4f ; largest jump one period of investment can give: sig*Imax = %.4f" % (emid - e0, sig * Imax))
print("modern endowment above which the low equilibrium disappears (eps0>0): omega_M > %.2f" % (wT + (cbar - wT) / e0))
print("permanent rise in omega_T that removes the low equilibrium: %.4f" % (cbar - wT - e0 * dw))
# one-off transfer of income tau to T in one period: raises coupling by sig*tau (below saturation)
for tau in [0.5, 1.0, 1.5, 2.0]:
    e = (1 - delta) * e0 + delta * e0 + sig * min(max(wT + tau + e0 * dw - cbar, 0), Imax)
    for t in range(300):
        e = F(e)
    print("  one-off transfer tau=%.1f -> long run %.4f" % (tau, e))
print("\nAll Chapter 2 solution checks passed.")
