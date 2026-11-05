"""Chapter 18: reproduces every number in the worked examples, viz panel and
exercise answers of chapters/chapter18.tex.  numpy only.
Run from the repo root:  python3 code/ch18.py
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def perron(W):
    """Left Perron vector of an irreducible row-stochastic W (sums to 1)."""
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def gap(W):
    """Gamma = 1 - |lambda_2| (docs/conventions.md, item 3)."""
    m = np.sort(abs(np.linalg.eigvals(W)))[::-1]
    return 1 - m[1]


def rownorm(K):
    return K / K.sum(1, keepdims=True)


print("=== Prop. two-block (equitable partition) ===")
# within-group intensity a, unconstrained cross intensity c <= a, constraint lambda
n, a, c = 3, 1.0, 1.0
for lam in [0.1, 0.3, 0.5, 1.0]:
    K = np.full((2 * n, 2 * n), lam * c)
    K[:n, :n] = a
    K[n:, n:] = a
    W = rownorm(K)
    G = gap(W)
    formula = 1 - abs(a - lam * c) / (a + lam * c)
    assert np.isclose(G, formula)
    assert np.isclose(G, 2 * lam * c / (a + lam * c))
    assert np.allclose(perron(W), 1 / (2 * n))  # pi uniform for every lambda > 0
    print(f"lambda={lam:.1f}: Gamma={G:.3f}  (2*lam/(1+lam)={2*lam/(1+lam):.3f}), pi uniform")

print("--- counterexample when c > a (cross stronger than within) ---")
a, c = 0.05, 1.0
for lam in [1.0, 0.5, 0.2, 0.05]:
    K = np.full((2 * n, 2 * n), lam * c)
    K[:n, :n] = a
    K[n:, n:] = a
    G = gap(rownorm(K))
    assert np.isclose(G, 1 - abs(a - lam * c) / (a + lam * c))
    print(f"lambda={lam:.2f}: Gamma={G:.3f}")

print("\n=== Prop. homophily and group power (pi proportional to degree) ===")
def degree_ratio(p, lam):
    return (p + lam * (1 - p)) / ((1 - p) + lam * p)

r_bol = degree_ratio(0.41, 0.35)
print(f"minority p=0.41, lambda_r=0.35: per-capita pi ratio = {r_bol:.3f}")
assert abs(r_bol - 0.84) < 0.005
assert np.isclose(degree_ratio(0.5, 0.3), 1.0)
# individual-level check: 41 indigenous, 59 non-indigenous agents, uniform W^max
nI, nM, lam = 41, 59, 0.35
g = np.array([1] * nI + [0] * nM)
K = np.where(g[:, None] == g[None, :], 1.0, lam)
p = perron(rownorm(K))
assert np.allclose(p, K.sum(1) / K.sum())   # pi_i = d_i / sum d
print(f"  agent-level check: {p[:nI].mean()/p[nI:].mean():.4f}")

# asymmetric densities, equal sizes: ratio = (aF + lam c)/(aM + lam c)
lam_g, aM, c = 0.3, 1.0, 1.0
aF = 0.58 * (aM + lam_g * c) - lam_g * c
print(f"calibration: aF/aM = {aF:.3f} gives pi_F/pi_M = {(aF+lam_g*c)/(aM+lam_g*c):.3f} at lambda_g=0.3")
r_open = (aF + c) / (aM + c)
print(f"  same densities at lambda_g = 1: pi_F/pi_M = {r_open:.3f}")

print("\n=== Example: Bolivia ethnic structure (illustrative per-capita index) ===")
pI, pM, uI, uM = 0.41, 0.59, 0.008, 0.014
share = pI * uI / (pI * uI + pM * uM)
print(f"indigenous power share = {share:.3f}; power/population = {share/pI:.3f}; "
      f"per-capita ratio = {uI/uM:.3f}; population ratio = {pI/pM:.3f}")
print(f"indigenous woman vs non-indigenous man: {uI/uM:.2f} x 0.58 = {uI/uM*0.58:.2f}")

print("\n=== Intersectional separable baseline ===")
print(f"lambda_r*lambda_g = 0.35*0.30 = {0.35*0.30:.3f}; x lambda_c 0.4 = {0.3*0.35*0.4:.3f}")

print("\n=== Viz panel: 8 agents, within M=55, within F=35, cross=10 ===")
K = np.full((8, 8), 10.0)
K[:4, :4] = 55
K[4:, 4:] = 35
W = rownorm(K)
pv = perron(W)
print("pi =", pv, f" ratio F/M = {pv[4:].mean()/pv[:4].mean():.3f}  Gamma = {gap(W):.3f}")
print("bar lengths (pi*18):", (pv * 18).round(3), " 1/N*18 =", 18 / 8)

print("\n=== Exercise 18.3 (population-weighted W^max) ===")
lg, lr = 0.3, 0.35
attrs = [(0, 0), (1, 0), (0, 1), (1, 1)]  # (female?, indigenous?)
nsz = np.array([0.295, 0.295, 0.205, 0.205])
phi = np.array([[(1 if a_[0] == b_[0] else lg) * (1 if a_[1] == b_[1] else lr)
                 for b_ in attrs] for a_ in attrs])
print("phi =\n", phi)
Wu = rownorm(phi * 0.25)
print("uniform W^max: pi =", perron(Wu), "(uniform: every row of phi sums to", phi.sum(1)[0], ")")
Kp = phi * nsz[None, :]
Wp = rownorm(Kp)
pp = perron(Wp)
pc = pp / nsz
print("population-weighted: pi =", pp, " per-capita pi/n =", pc / pc[0])
print(f"  pi4/pi1 per capita = {pc[3]/pc[0]:.4f}  vs lambda_g*lambda_r = {lg*lr:.3f}")
assert np.isclose(pc[3] / pc[0], degree_ratio(0.41, 0.35))

print("\n=== Prop. formation gap and Exercise 18.4 ===")
rho, mu, tbar, rstar = 0.85, 0.3, 1.0, 0.1
def steady(lam):
    om = lam / (1 + lam)          # cross-group weight after row-normalisation
    # theta_g = rho theta_g + (1-rho) tbar + mu((1-om) r_g + om r_other - r*)
    # r_A = 0.05 + 0.1 theta_A ; r_D = 0.03 + 0.1 theta_D
    b0 = np.array([0.05, 0.03])
    Wi = np.array([[1 - om, om], [om, 1 - om]])
    M = (1 - rho) * np.eye(2) - mu * 0.1 * Wi
    rhs = (1 - rho) * tbar + mu * (Wi @ b0 - rstar)
    return np.linalg.solve(M, rhs)
for lam in [0.25, 0.5]:
    th = steady(lam)
    print(f"lambda={lam}: own weight {1/(1+lam):.3f}; theta*_A={th[0]:.4f} theta*_D={th[1]:.4f} "
          f"Delta={th[0]-th[1]:.4f}")
x_star = 0.0015 / 0.0063
lam_star = (1 - x_star) / (1 + x_star)
print(f"Delta < 0.01 iff (1-lam)/(1+lam) < {x_star:.4f}, i.e. lambda > {lam_star:.3f}")
assert steady(lam_star + 1e-6)[0] - steady(lam_star + 1e-6)[1] < 0.01 < steady(lam_star - 1e-3)[0] - steady(lam_star - 1e-3)[1]
# symmetric version of the original exercise gives Delta = 0
def steady_sym(lam):
    om = lam / (1 + lam); Wi = np.array([[1 - om, om], [om, 1 - om]])
    M = (1 - rho) * np.eye(2) - mu * 0.1 * Wi
    return np.linalg.solve(M, (1 - rho) * tbar + mu * (Wi @ np.array([.05, .05]) - rstar))
assert np.isclose(*steady_sym(0.25))

print("\n=== Welfare example (per-capita, two equal groups) ===")
kappa, phi_pc, omega_pc = 0.5, 1.0, 1.0
def loss_pc(r):
    dev = (1 - r) / (1 + r)            # N*pi_M - 1 = 1 - N*pi_F
    N2Var = dev ** 2                   # N^2 Var(pi)
    xM, xF = omega_pc + dev * phi_pc, omega_pc - dev * phi_pc
    assert kappa * xM < 1
    return kappa / 2 * phi_pc ** 2 * N2Var, xM, xF, N2Var
L03, xM, xF, nv = loss_pc(0.58)
L1, *_ = loss_pc(r_open)
U = lambda x: x - kappa / 2 * x ** 2
assert np.isclose(L03, U(1.0) - (U(xM) + U(xF)) / 2)
print(f"r=0.58: N^2 Var(pi)={nv:.4f}, x_M={xM:.3f}, x_F={xF:.3f}, L/N={L03:.4f} "
      f"({L03/U(1.0)*100:.2f}% of U(Omega/N)={U(1.0):.2f})")
print(f"lambda_g=1 (r={r_open:.3f}): L/N={L1:.4f};  Delta L/N = {L03-L1:.4f} "
      f"({(L03-L1)/L03*100:.0f}% of the power loss)")
N = 11e6
print(f"Var(pi) maximum at N=11e6: (N-1)/N^2 = {(N-1)/N**2:.2e}")
print(f"(1+Psi) with Psi=2.3: {1+2.3:.1f};  enforcement share 1/(1+Psi) = {1/3.3:.2f}")

print("\n=== Half-lives (case study) ===")
for r in [0.85, 0.90]:
    h = np.log(2) / np.log(1 / r)
    print(f"rho={r}: {h:.2f} generations = {25*h:.0f}-{30*h:.0f} years (25-30y generations); "
          f"500 years = {500/(30*h):.1f}-{500/(25*h):.1f} half-lives")
print(f"residual after 2.5-4.7 half-lives: {0.5**4.7:.3f}-{0.5**2.5:.3f}")

print("\n=== Prop. assortative (b): Cheeger bound on random symmetric K^max ===")
rng = np.random.default_rng(18)
for trial in range(200):
    n1, n2 = rng.integers(2, 6, size=2)
    Km = rng.random((n1 + n2, n1 + n2)); Km = (Km + Km.T) / 2
    lam = rng.random()
    grp = np.array([0] * n1 + [1] * n2)
    same = grp[:, None] == grp[None, :]
    K = Km * np.where(same, 1.0, lam)
    C = Km[np.ix_(grp == 0, grp == 1)].sum()
    A0 = Km[np.ix_(grp == 0, grp == 0)].sum(); A1 = Km[np.ix_(grp == 1, grp == 1)].sum()
    bound = 2 * lam * C / (min(A0, A1) + lam * C)
    assert gap(rownorm(K)) <= bound + 1e-12
print("Cheeger bound holds in 200 random trials; tight in the equal-group example:",
      np.isclose(2 * 0.3 / 1.3, 1 - abs(1 - 0.3) / 1.3))
