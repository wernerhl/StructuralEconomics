"""Solutions manual, Chapter 7 (The Bargaining Equilibrium).
Computes and prints every number used in solutions/sol07.tex.
Run from the repository root:  python3 code/solutions/sol07.py
numpy only, except Exercise 7.7, which reads data/processed with pandas."""
import pathlib
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)
ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"


def perron(W):
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmin(np.abs(vals - 1))
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def gap(W):
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


def neff(p):
    return 1 / np.sum(p ** 2)


def varpi(p):
    return np.mean((p - 1 / len(p)) ** 2)


def rule(pi, d, Omega):
    """TU asymmetric Nash (Prop ch07:prop:allocation-rule): x_i = d_i + pi_i (Omega - sum d)."""
    return d + pi * (Omega - d.sum())


def Uq(x, kappa):
    return x - kappa / 2 * x ** 2


def gini(x):
    x = np.sort(x); n = len(x)
    return np.sum((2 * np.arange(1, n + 1) - n - 1) * x) / (n * x.sum())


# ===================================================================== 7.1
print("=" * 72); print("Exercise 7.1: landlord-farmer with Y = A K^0.3 L^0.7")
K, Lab, Y = 10.0, 2.0, 20.0
dL, dF = 10.0, 2.0
S = Y - dL - dF
for a in (0.5, 0.3):
    A = Y / (K ** a * Lab ** (1 - a))
    MPK, MPL = a * Y / K, (1 - a) * Y / Lab
    payL, payF = MPK * K, MPL * Lab
    piF_needed = (payF - dF) / S
    print(f"capital exponent {a}: A = {A:.4f}, MPK = {MPK:.2f}, MPL = {MPL:.2f}, "
          f"competitive (L, F) = ({payL:.1f}, {payF:.1f}); pi_F needed = {piF_needed:.3f}")
    assert np.isclose(A * K ** a * Lab ** (1 - a), Y)
for pL in (0.7,):
    print(f"bargaining pi=(0.7,0.3): x = ({dL + pL * S:.1f}, {dF + (1 - pL) * S:.1f}) (independent of the technology)")
# range of competitive Farmer payments reproducible by some pi_F in [0,1]
print("reproducible Farmer payments: [d_F, Y - d_L] = [", dF, ",", Y - dL, "];"
      " capital share a must satisfy a*Y >= d_L, i.e. a >=", dL / Y)

# ===================================================================== 7.2
print("=" * 72); print("Exercise 7.2: Omega = 3, d = (0.6, 0.4, 0.2), pi = (0.5, 0.3, 0.2), kappa = 0.25")
Om = 3.0; d = np.array([.6, .4, .2]); pi = np.array([.5, .3, .2]); kap = 0.25
print("(a) endowment economy: omega = (1.5, 1.0, 0.5), sum =", 1.5 + 1.0 + 0.5)


def solve_foc(pi, d, Om, Sfun, dSfun, x0, iters=100):
    """Newton on F(x, mu) = [pi_i dS_i/dx_i / S_i - mu ; sum x - Om] = 0."""
    n = len(pi); z = np.r_[x0, 1.0]
    for _ in range(iters):
        x, mu = z[:n], z[n]
        Sx, dS = Sfun(x), dSfun(x)
        F = np.r_[pi * dS / Sx - mu, x.sum() - Om]
        # derivative of pi dS/S wrt x_i: pi (d2S S - dS^2)/S^2
        d2S = dS2fun(x)
        J = np.zeros((n + 1, n + 1))
        J[np.arange(n), np.arange(n)] = pi * (d2S * Sx - dS ** 2) / Sx ** 2
        J[:n, n] = -1; J[n, :n] = 1
        z = z - np.linalg.solve(J, F)
        if np.max(np.abs(F)) < 1e-14:
            break
    return z[:n], z[n]


# TU: S_i = x_i - d_i
dS2fun = lambda x: np.zeros_like(x)
xTU, muTU = solve_foc(pi, d, Om, lambda x: x - d, lambda x: np.ones_like(x), np.full(3, 1.0))
xcl = rule(pi, d, Om)
print("(c) Newton on FOC (TU): x =", np.round(xTU, 6), " mu =", round(muTU, 6), " 1/mu =", round(1 / muTU, 6))
print("    closed form       : x =", xcl, " surplus Omega - sum d =", Om - d.sum())
assert np.allclose(xTU, xcl)
print("    kappa x =", kap * xcl, " max =", kap * xcl.max())
assert np.all(kap * xcl < 1)
# second-order / global check: random feasible points never beat x
rng = np.random.default_rng(7)
f = lambda z: np.sum(pi * np.log(z - d))
for _ in range(5000):
    z = d + rng.dirichlet(np.ones(3)) * (Om - d.sum())
    assert f(z) <= f(xcl) + 1e-12
print("    ok: no random feasible allocation beats x (5000 draws); Nash product value", round(f(xcl), 6))
# (d) equal power
xeq = rule(np.ones(3) / 3, d, Om)
print("(d) equal power x =", np.round(xeq, 4), " difference x(pi) - x(eq) =", np.round(xcl - xeq, 4),
      " = (pi - 1/3)*1.8 =", np.round((pi - 1 / 3) * 1.8, 4))
for nm, x in (("pi", xcl), ("equal", xeq), ("equal division 1,1,1", np.ones(3))):
    print(f"    utilities U(x) at {nm}: {np.round(Uq(x, kap), 6)}  sum = {Uq(x, kap).sum():.6f}")
# NTU variant: Nash over utilities, S_i = U(x_i) - U(d_i)
dS2fun = lambda x: np.full_like(x, -kap)
xN, muN = solve_foc(pi, d, Om, lambda x: Uq(x, kap) - Uq(d, kap), lambda x: 1 - kap * x, xcl.copy())
print("    NTU variant S_i = U(x_i) - U(d_i): x =", np.round(xN, 4), " FOC residuals",
      np.round(pi * (1 - kap * xN) / (Uq(xN, kap) - Uq(d, kap)) - muN, 12))
dS2fun = lambda x: np.zeros_like(x)

# ===================================================================== 7.3
print("=" * 72); print("Exercise 7.3: welfare decomposition, Economies A and B")
kappa, Omega, N, dbar = 0.1, 20.0, 4, 3.0
Phi = Omega - N * dbar
WA = np.array([[.6, .2, .1, .1], [.2, .5, .2, .1], [.1, .3, .4, .2], [.1, .1, .2, .6]])
WB = np.array([[.48, .50, .01, .01], [.40, .58, .01, .01], [.03, .03, .64, .30], [.03, .03, .30, .64]])
U = lambda z: Uq(z, kappa)


def fb(w, kappa, Omega):
    """argmax sum_j w_j U(x_j) s.t. sum x = Omega: w_j (1 - kappa x_j) = mu."""
    mu = (len(w) - kappa * Omega) / np.sum(1 / w)
    return (1 - mu / w) / kappa


res = {}
for name, W in (("A", WA), ("B", WB)):
    p = perron(W)
    x = Omega / N + (p - 1 / N) * Phi
    assert np.all(kappa * x < 1)
    Lp = np.sum(U(np.full(N, Omega / N)) - U(x))
    assert np.isclose(Lp, kappa * N / 2 * Phi ** 2 * varpi(p))
    res[name] = (p, x, Lp)
    print(f"Economy {name}: pi = {np.round(p, 4)}, Var(pi) = {varpi(p):.6f}, x = {np.round(x, 4)}, L_power = {Lp:.4f}")
print("ratio L_B / L_A =", round(res["B"][2] / res["A"][2], 2))
g = 0.3
dec = {}
for name, W in (("A", WA), ("B", WB)):
    p, x, Lp = res[name]
    WU = W - np.diag(np.diag(W))
    c = WU.sum(0)
    w = 1 + g * c
    SW = lambda z: np.sum(w * U(z))
    xfb = fb(w, kappa, Omega)
    assert np.all(kappa * xfb < 1) and np.isclose(xfb.sum(), Omega)
    Lext = SW(xfb) - SW(np.full(N, Omega / N))
    Lint = g * np.sum(c * (U(Omega / N) - U(x)))
    Ltot = SW(xfb) - SW(x)
    assert np.isclose(Ltot, Lp + Lext + Lint)
    up = 1 - kappa * Omega / N
    cov = np.mean((c - c.mean()) * (p - 1 / N))
    lin = -g * N * up * Phi * cov
    quad = g * kappa / 2 * Phi ** 2 * np.sum(c * (p - 1 / N) ** 2)
    assert np.isclose(Lint, lin + quad)
    corr = np.corrcoef(c, p)[0, 1]
    dec[name] = (Lp, Lext, Lint, Ltot, cov)
    print(f"Economy {name}: W^U row sums {np.round(WU.sum(1), 2)}; c = {np.round(c, 2)}; weights 1+g c = {np.round(w, 3)}")
    print(f"   x_FB = {np.round(xfb, 4)}  S(x_FB) = {SW(xfb):.4f}  S(x_E) = {SW(np.full(N, 5.0)):.4f}  S(x_B) = {SW(x):.4f}")
    print(f"   L_power = {Lp:.6f}  L_ext = {Lext:.6f}  L_int = {Lint:.6f} (linear {lin:.6f}, quadratic {quad:.6f})  "
          f"total = {Ltot:.6f}")
    print(f"   Cov(c, pi) = {cov:.6f}  corr(c, pi) = {corr:.3f}  U'(Omega/N) = {up}")

# ===================================================================== 7.4
print("=" * 72); print("Exercise 7.4: transfer vs structural policy in Economy B")
pB, xB, LB = res["B"]
SWu = lambda z: np.sum(U(z))
print("baseline: x =", np.round(xB, 4), " sum U =", round(SWu(xB), 4), " loss vs equal =", round(LB, 4),
      " Gini(x) =", round(gini(xB), 4), " max/min =", round(xB.max() / xB.min(), 4))
dT = np.array([2., 3, 3, 4])
xT = rule(pB, dT, Omega)
print("(a) transfer d' = (2,3,3,4): Omega - sum d' =", Omega - dT.sum(), " x =", np.round(xT, 4),
      " sum =", round(xT.sum(), 6))
print("    sum U =", round(SWu(xT), 4), " gain =", round(SWu(xT) - SWu(xB), 4),
      " loss vs equal =", round(np.sum(U(np.full(N, 5.)) - U(xT)), 4), " Gini(x) =", round(gini(xT), 4),
      " max/min =", round(xT.max() / xT.min(), 4))
assert np.allclose(xT - xB, dT - dbar)          # pass-through: dx = dd when sum d fixed
WB2 = WB.copy()
for i in range(4):
    for j in range(4):
        if (i < 2) != (j < 2):
            WB2[i, i] -= 0.05 - WB2[i, j]; WB2[i, j] = 0.05
assert np.allclose(WB2.sum(1), 1)
p2 = perron(WB2); x2 = Omega / N + (p2 - 1 / N) * Phi
print("(b) structural W^P =\n", WB2)
print("    pi =", np.round(p2, 4), " block-1 mass =", round(p2[:2].sum(), 4), " (eps1=eps2=0.10 -> 0.5)",
      " Gamma =", round(gap(WB), 4), "->", round(gap(WB2), 4), " N_eff =", round(neff(pB), 4), "->", round(neff(p2), 4))
print("    x =", np.round(x2, 4), " sum U =", round(SWu(x2), 4), " gain =", round(SWu(x2) - SWu(xB), 4),
      " loss vs equal =", round(np.sum(U(np.full(N, 5.)) - U(x2)), 4), " Gini(x) =", round(gini(x2), 4),
      " max/min =", round(x2.max() / x2.min(), 4))
assert np.isclose(p2[:2].sum(), 0.5)
# eigenvalues of WB2
print("    eigenvalues new:", np.round(np.sort(np.linalg.eigvals(WB2).real)[::-1], 4))
# combined policy
xboth = rule(p2, dT, Omega)
print("    both policies: x =", np.round(xboth, 4), " sum U =", round(SWu(xboth), 4))

# ===================================================================== 7.5
print("=" * 72); print("Exercise 7.5: escape threshold, star N=5, a=0.9, s=0.5")


def star(N, a, s):
    W = np.zeros((N, N))
    W[0, 0] = s; W[0, 1:] = (1 - s) / (N - 1)
    for i in range(1, N):
        W[i, 0] = a; W[i, i] = 1 - a
    return W


Ws = star(5, 0.9, 0.5); ps = perron(Ws)
print("pi =", np.round(ps, 4), " pi_h = a/(a+1-s) =", round(0.9 / 1.4, 4), " Gamma =", round(gap(Ws), 4),
      " spectrum", np.round(np.sort(np.linalg.eigvals(Ws).real), 4), " N_eff =", round(neff(ps), 4))
target = 1.5 * ps[1]
pih_new = 1 - 4 * target
a_new = pih_new * 0.5 / (1 - pih_new)            # pi_h = a/(a + 1 - s)  =>  a = pi_h (1-s)/(1-pi_h)
print("spoke target =", round(target, 4), " => pi_h =", round(pih_new, 4), " => a =", round(a_new, 4),
      " decrease =", round(0.9 - a_new, 4), " (", round(100 * (1 - a_new / 0.9), 1), "%)")
assert np.isclose(perron(star(5, a_new, 0.5))[1], target)
Wn = star(5, a_new, 0.5)
print("new Gamma =", round(gap(Wn), 4), " new spectrum", np.round(np.sort(np.linalg.eigvals(Wn).real), 4))
# general closed form: spoke power = (1-pi_h)/(N-1) = (1-s)/((N-1)(a+1-s))
# one spoke alone lowers her dependence (others keep a = 0.9)
for ai in (0.6, 0.43, 0.3, 0.1):
    W = Ws.copy(); W[1, 0] = ai; W[1, 1] = 1 - ai
    q = perron(W)
    print(f"  only spoke 1 lowers a_1 to {ai}: pi_1 = {q[1]:.4f} ({100*(q[1]/ps[1]-1):+.1f}%), other spokes {q[2]:.4f}, hub {q[0]:.4f}")
# exact single-spoke threshold for +50%: bisection
lo, hi = 0.0, 0.9
for _ in range(100):
    mid = (lo + hi) / 2
    W = Ws.copy(); W[1, 0] = mid; W[1, 1] = 1 - mid
    if perron(W)[1] >= target:
        lo = mid
    else:
        hi = mid
print("  single spoke reaches +50% at a_1 =", round(lo, 4))
# first-order sensitivity d pi_spoke / d a (all spokes) at a = 0.9
da = 1e-6
dps = (perron(star(5, 0.9 + da, 0.5))[1] - perron(star(5, 0.9 - da, 0.5))[1]) / (2 * da)
print("  linear prediction of the spoke gain:", round(-dps * (0.9 - a_new), 4), " required gain:", round(target - ps[1], 4))
print("  d pi_spoke / d a at a=0.9 =", round(dps, 4), " closed form -(1-s)/((N-1)(a+1-s)^2) =", round(-0.5 / (4 * 1.4 ** 2), 4))

# ===================================================================== 7.6
print("=" * 72); print("Exercise 7.6: Bolivia 4x4 (Gas, Government, Urban informal, Rural traditional)")
lab = ["Gas", "Gov", "Urban", "Rural"]
WP = np.array([[.50, .30, .15, .05],
               [.30, .40, .20, .10],
               [.15, .10, .50, .25],
               [.10, .05, .25, .60]])
WU = np.array([[.60, .30, .05, .05],
               [.30, .50, .10, .10],
               [.10, .10, .60, .20],
               [.05, .05, .20, .70]])
WI = np.array([[.60, .30, .05, .05],
               [.30, .60, .05, .05],
               [.20, .30, .40, .10],
               [.10, .40, .10, .40]])
alpha = np.array([.4, .3, .3])
for M in (WP, WU, WI):
    assert np.allclose(M.sum(1), 1)


def composite(WP, WU, WI):
    return alpha[0] * WP + alpha[1] * WU + alpha[2] * WI


W = composite(WP, WU, WI)
pc = perron(W)
for nm, M in (("W^P", WP), ("W^U", WU), ("W^I", WI), ("composite", W)):
    q = perron(M)
    print(f"  pi({nm}) = {np.round(q, 4)}  N_eff = {neff(q):.3f}  Gamma = {gap(M):.4f}")
print("composite W =\n", np.round(W, 4))
d6 = np.array([4., 3., 2., 1.]); Om6 = 20.0; k6 = 0.05
x6 = rule(pc, d6, Om6)
print("(b) d =", d6, " Omega =", Om6, " surplus =", Om6 - d6.sum(), " x =", np.round(x6, 4),
      " kappa x max =", round(k6 * x6.max(), 3), " sum U =", round(Uq(x6, k6).sum(), 4))
assert np.all(k6 * x6 < 1)
# (i) transfer of 1 unit of outside option from Gas to Rural
d6i = d6 + np.array([-1, 0, 0, 1.])
xi = rule(pc, d6i, Om6)
# (ii) plant: Rural buys fertilizer from Gas (0.10 of Rural's W^P row from self to Gas); output +2
WPii = WP.copy(); WPii[3, 0] += 0.10; WPii[3, 3] -= 0.10
pii = perron(composite(WPii, WU, WI))
xii_net = rule(pii, d6, Om6)                 # network effect only
xii = rule(pii, d6, Om6 + 2.0)               # with productivity gain
# (ii') plant that buys Rural output (Gas's W^P row moves 0.10 from self to Rural), output +2
WPiib = WP.copy(); WPiib[0, 3] += 0.10; WPiib[0, 0] -= 0.10
piib = perron(composite(WPiib, WU, WI))
xiib = rule(piib, d6, Om6 + 2.0)
# (iii) rural internet in W^I: Rural listens more to Urban; Urban and Government listen more to Rural
WIiii = WI.copy()
WIiii[3, 3] -= 0.10; WIiii[3, 2] += 0.10           # Rural: own 0.4 -> 0.3, Urban 0.1 -> 0.2
WIiii[2, 2] -= 0.10; WIiii[2, 3] += 0.10           # Urban: own 0.4 -> 0.3, Rural 0.1 -> 0.2
WIiii[1, 1] -= 0.05; WIiii[1, 3] += 0.05           # Gov: own 0.6 -> 0.55, Rural 0.05 -> 0.10
assert np.allclose(WIiii.sum(1), 1)
piii = perron(composite(WP, WU, WIiii))
xiii = rule(piii, d6, Om6)
# (iii') one-way internet: Rural only listens more
WIone = WI.copy(); WIone[3, 3] -= 0.10; WIone[3, 2] += 0.10
pone = perron(composite(WP, WU, WIone)); xone = rule(pone, d6, Om6)
for nm, q, x in (("baseline", pc, x6), ("(i) transfer", pc, xi), ("(ii) plant, network only", pii, xii_net),
                 ("(ii) plant, +2 output", pii, xii), ("(ii') plant buying rural output, +2", piib, xiib),
                 ("(iii) two-way internet", piii, xiii), ("(iii') one-way internet", pone, xone)):
    print(f"  {nm:38s} pi = {np.round(q, 4)}  x = {np.round(x, 4)}  dx_R = {x[3]-x6[3]:+.4f}  sum U = {Uq(x, k6).sum():.4f}")
print("  structural policy beats a 1-unit transfer iff d pi_R * surplus > 1, i.e. d pi_R >", 1 / (Om6 - d6.sum()))

# ===================================================================== 7.7
print("=" * 72); print("Exercise 7.7: labor share (INE SUT, base 2017) and implied pi_labor")
d = pd.read_csv(DATA / "ine_sut_long.csv", dtype={"row_code": str, "col_code": str})
u = d[(d.base == 2017) & (d.table == "use") & (d.row_code.isin(["B1b", "D1", "D2", "D3", "B2b", "B3b"]))]
u = u[u.col_code.str.fullmatch(r"\d\d", na=False)]          # the 39 industries
t = u.groupby(["year", "row_code"]).value.sum().unstack()
ident = (t.D1 + t.D2 - t.D3 + t.B2b + t.B3b) / t.B1b
assert np.allclose(ident, 1, atol=1e-6), ident
t["sL"] = t.D1 / t.B1b
# Gollin-type adjustment: attribute to labor the employee share of mixed income
t["sL_adj"] = (t.D1 + t.B3b * t.D1 / (t.D1 + t.B2b)) / t.B1b
t["mixed"] = t.B3b / t.B1b
print(t[["B1b", "D1", "B2b", "B3b", "sL", "sL_adj", "mixed"]].round(4).to_string())
s19, s19a = float(t.loc[2019, "sL"]), float(t.loc[2019, "sL_adj"])
wdi = pd.read_csv(DATA / "wdi_long.csv")
self_emp = wdi[(wdi.indicator == "SL.EMP.SELF.ZS") & (wdi.year == 2019) & wdi.country_iso3.isin(["BOL", "CHL", "MEX"])]
print("self-employment % of employment, 2019 (WDI SL.EMP.SELF.ZS):",
      {r.country_iso3: round(r.value, 1) for r in self_emp.itertuples()})


def pi_labor(sL, dLab, dK):
    return (sL - dLab) / (1 - dLab - dK)


print("implied pi_labor = (s_L - d_L)/(1 - d_L - d_K), Omega = GVA = 1:")
for sL, nm in ((s19, "s_L = D1/B1b"), (s19a, "s_L adjusted")):
    for dLab, dK in ((0, 0), (0.2, 0.2), (0.3, 0.1), (0.1, 0.3), (0.35, 0.2)):
        print(f"   {nm:14s} {sL:.4f}  d_L={dLab:.2f} d_K={dK:.2f}:  pi_L = {pi_labor(sL, dLab, dK):.4f}")
print("derivatives at d=(0.2,0.2), s_L=", round(s19, 4), ": dpi/d d_L =", round((s19 - 1 + 0.2) / (1 - 0.4) ** 2, 4),
      " dpi/d d_K =", round((s19 - 0.2) / (1 - 0.4) ** 2, 4))
# identified set when outside options are only known to be nonnegative and individually rational
print("pi_L ranges over (0,1) as d_L ranges over [0, s_L) and d_K over [0, 1 - s_L): not identified")
print("all assertions passed")
