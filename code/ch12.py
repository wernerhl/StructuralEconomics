"""Chapter 12: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch12.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def perron(W):
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmax(vals.real)
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def gap(W):
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


def v(z, kappa=0.5):
    """Reference-dependent utility v(x - a) = z - (kappa/2) z^2, increasing for kappa*z < 1."""
    return z - 0.5 * kappa * z ** 2


# ---------------------------------------------------------------------------
# Example ch12:ex:satisfaction: Chapter 10's Bolivia allocations, aspirations a = M W^I x
N = 9
blk = np.repeat([0, 1, 2], 3)
xb = np.array([1.2486, 0.9852, 0.7662])          # block means from code/ch10.py
x = np.repeat(xb, 3)
WI = np.where(blk[:, None] == blk[None, :], 0.85 / 3, 0.025)
Wce = np.ones((N, N)) / N
Ma = 0.8                                          # aspiration multiplier mu/(1-rho), e.g. rho=.75, mu=.2
kappa = 0.5
print("== Example satisfaction")
rows = {}
for name, a in [("endogenous (block W^I)", Ma * WI @ x),
                ("integrated (W^CE)", Ma * Wce @ x),
                ("no aspirations (a=0)", np.zeros(N))]:
    z = x - a
    assert np.all(kappa * z < 1)
    s = v(z, kappa)
    rows[name] = s
    print(f"{name:24s} a={np.round(a[::3],3)}  x-a={np.round(z[::3],3)}  s={np.round(s[::3],3)}  "
          f"gas-rural gap={s[0]-s[6]:.3f}")
g_endo = rows["endogenous (block W^I)"][0] - rows["endogenous (block W^I)"][6]
g_int = rows["integrated (W^CE)"][0] - rows["integrated (W^CE)"][6]
print("compression ratio endogenous/integrated =", g_endo / g_int, " (1 - M) =", 1 - Ma)
print("allocation gap gas-rural =", xb[0] - xb[2])

# Proposition ch12:thm:paradox: two equal blocks, pure block vs integrated -> ratio exactly 1-M
for xM, xT, M_, k_ in [(1.2, 0.8, 0.8, 1.0), (1.0, 0.4, 0.5, 0.8)]:
    sb = v((1 - M_) * xM, k_) - v((1 - M_) * xT, k_)
    xbar = 0.5 * (xM + xT)
    si = v(xM - M_ * xbar, k_) - v(xT - M_ * xbar, k_)
    assert abs(sb / si - (1 - M_)) < 1e-12
print("two-block exact ratio check passed")

# Proposition ch12:prop:endo_linear: sum_i (x_i - M (W x)_i) = Omega - M sum_j c_j x_j
cj = WI.sum(0)
lhs = np.sum(x - Ma * WI @ x)
print("linear-v endogenous total:", lhs, " = Omega - M c'x:", x.sum() - Ma * cj @ x)
# watching the rich lowers measured welfare: everyone puts extra weight on gas block
Wrich = 0.7 * WI + 0.3 * np.tile(np.where(blk == 0, 1 / 3, 0.0), (N, 1))
print("   if 30% of every row moves to the gas block: total =", np.sum(x - Ma * Wrich @ x))

# ---------------------------------------------------------------------------
# Proposition ch12:prop:fixed_eval and Remark (Gamma is not power)
print("\n== Fixed-evaluation ranking and the Gamma counterexample")
piA = np.array([0.5, 0.3, 0.2])
WA = np.outer(np.ones(3), piA)
WB = 0.5 * np.ones((3, 3)) / 3 + 0.5 * np.eye(3)
Om, Phi, kap = 3.0, 1.0, 0.5
for nm, W in [("A", WA), ("B", WB)]:
    p = perron(W)
    xa = Om / 3 + (p - 1 / 3) * Phi
    assert np.all(kap * xa < 1)
    U = np.sum(xa - 0.5 * kap * xa ** 2)
    Vp = np.mean((p - 1 / 3) ** 2)
    print(f"W_{nm}: Gamma={gap(W):.3f}, pi={np.round(p,3)}, N_eff={1/np.sum(p**2):.3f}, Var(pi)={Vp:.5f}, "
          f"x={np.round(xa,3)}, sum U={U:.5f}, loss (kN/2)Phi^2 Var = {0.5*kap*3*Phi**2*Vp:.5f}")
U0 = 3 * (1 - 0.5 * kap * 1)
pA = perron(WA)
assert abs((U0 - np.sum((1 + (pA - 1 / 3)) - 0.5 * kap * (1 + (pA - 1 / 3)) ** 2))
           - 0.5 * kap * 3 * Phi ** 2 * np.mean((pA - 1 / 3) ** 2)) < 1e-12

# ---------------------------------------------------------------------------
# Capability: e_i = omega_i + pi_i Phi, O_i = [0, c_i e_i]
print("\n== Capability")
pi2 = np.array([0.7, 0.3])
e = 5 + pi2 * 2
print("Ex 12.3: entitlements e =", e, " with c_2 = 0.8: capability bounds", np.array([1.0, 0.8]) * e)
rng = np.random.default_rng(3)
for _ in range(1000):
    p = rng.dirichlet(np.ones(4))
    assert np.min(1 + p * 2) <= 1 + 0.25 * 2 + 1e-12     # min entitlement maximized at uniform pi
    assert abs(np.sum(1 + p * 2) - (4 + 2)) < 1e-12        # sum is invariant to pi

# ---------------------------------------------------------------------------
# Proposition ch12:prop:consent (Pareto nesting under endogenous aspirations)
print("\n== Consent / Pareto nesting")
W2 = np.array([[0.5, 0.5], [0.5, 0.5]])
M2 = 0.8
x0 = np.array([1.0, 1.0])
x1 = np.array([1.0, 1.1])
s0 = x0 - M2 * W2 @ x0
s1 = x1 - M2 * W2 @ x1
print("raise x_2 by 0.1: change in endogenous satisfaction (linear v):", s1 - s0)
# (I - M W)^{-1} >= 0 for M < 1
for _ in range(500):
    W = rng.dirichlet(np.ones(4), size=4)
    assert np.all(np.linalg.inv(np.eye(4) - 0.8 * W) >= -1e-12)

# ---------------------------------------------------------------------------
# Example ch12:ex:education_policy: rural aspirations under Ch10 Policy C at fixed x
WIc = WI.copy()
for i in np.where(blk == 2)[0]:
    WIc[i] = np.where(blk == 0, 0.075, np.where(blk == 2, 0.70 / 3, 0.025))
a_old = Ma * WI @ x
a_new = Ma * WIc @ x
print("\n== Education policy: rural aspiration", a_old[6], "->", a_new[6],
      " satisfaction", v(x[6] - a_old[6]), "->", v(x[6] - a_new[6]))

# ---------------------------------------------------------------------------
print("\n== Exercises")
# 12.1: satisfaction gap vs aspiration multiplier
for M_ in [0.0, 0.4, 0.8]:
    s = v(x - M_ * WI @ x)
    print(f"12.1 M={M_}: block satisfaction {np.round(s[::3],3)}, gas-rural gap {s[0]-s[6]:.3f}")
# 12.2: laundered with W^CE at the current allocation
ahat = Ma * x.mean()
print("12.2 laundered aspiration", ahat, " satisfaction by block", np.round(v(xb - ahat), 3))
# 12.7 Sen vs Bentham
k7 = 0.1
for pi_ in [np.array([0.8, 0.2]), np.array([0.5, 0.5])]:
    x7 = 5 + (pi_ - 0.5) * 4
    a7 = np.array([4.0, 1.5])
    endo = v(x7 - a7, k7)
    laund = v(x7 - 3.0, k7)
    print(f"12.7 pi={pi_}: x={x7}, Bentham(endog) {endo.sum():.4f} {np.round(endo,4)}, "
          f"laundered {laund.sum():.4f} {np.round(laund,4)}, Rawls(endog) {endo.min():.4f}, "
          f"capability min e {x7.min():.2f}")
# Appendix 12B identity: U = theta x - (k/2) x^2 at x = (theta + d)/k equals theta^2/(2k) - d^2/(2k)
for th_, d_, k_ in [(1.1, 0.07, 0.5), (0.9, -0.2, 1.3)]:
    xx = (th_ + d_) / k_
    assert abs(th_ * xx - 0.5 * k_ * xx ** 2 - (th_ ** 2 - d_ ** 2) / (2 * k_)) < 1e-12
print("Appendix 12B identity verified")
