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

print("\n=== Example (Bolivia's gender structure): ILOSTAT 2019 and Chapter 17 sector power per worker ===")
# Imported from code/ch17.py (Example ch17:ex:bolivia_pi), which reads INE and ILOSTAT from
# data/processed, so Ch17 and Ch18 use the same sector mapping and cannot drift.
import contextlib, io, os, runpy
with contextlib.redirect_stdout(io.StringIO()):
    _c17 = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch17.py"))
_r, _fem, _emp, _pie = _c17["ratio"], _c17["fem_share"], _c17["emp"], _c17["pi_emp"]
_chk = _c17["check"]
# female share of employment, table of Example ch18:ex:bolivia_gender (percent, one decimal)
FEM_TEXT = {"Hydrocarbons and mining": 6.5, "Manufacturing": 44.1, "Construction": 3.1, "Commerce": 61.4,
            "Agriculture": 46.0, "Domestic services": 93.9, "Public administration": 34.3}
_facts = _c17["FACTS"]
_fid = {"Hydrocarbons and mining": "mining_quarrying_incl_hydrocarbons", "Manufacturing": "manufacturing",
        "Construction": "construction", "Commerce": "commerce", "Agriculture": "agriculture",
        "Domestic services": "domestic_services", "Public administration": "public_administration"}
for _s, _v in FEM_TEXT.items():
    _chk("female share " + _s, 100 * _fem[_s], _v, 1)
    _chk("male share " + _s, 100 * (1 - _fem[_s]), round(100 - _v, 1), 1)
    assert _facts[f"BOL_female_share_{_fid[_s]}_2019"] == _v          # data/processed/bolivia_facts.csv
    print(f"  {_s:24s} female {100*_fem[_s]:5.1f}%  pi/share = {_r[_s]:.2f}")
for _s, _v in {"Personal services": 69, "Financial services": 47}.items():
    _chk("female share " + _s, 100 * _fem[_s], _v, 0)
for _s, _v in {"Personal services": 0.20, "Hydrocarbons and mining": 5.7, "Financial services": 23.4,
               "Construction": 0.37, "Public administration": 0.19}.items():
    _chk("pi/share " + _s, _r[_s], _v, 1 if _v >= 5 else 2)
# per-capita power by sex: each sector's pi shared equally among its workers; domestic service
# (no W^P row) gets none.  pi_F / pi_M = [sum_s pi_s f_s / F] / [sum_s pi_s (1 - f_s) / M]
_ILO = _c17["ILO"]
_F, _M = _ILO[("SEX_F", "ECO_ISIC4_TOTAL")], _ILO[("SEX_M", "ECO_ISIC4_TOTAL")]
_piF = sum(0 if np.isnan(_pie[s]) else _pie[s] * _fem[s] for s in _pie) / _F
_piM = sum(0 if np.isnan(_pie[s]) else _pie[s] * (1 - _fem[s]) for s in _pie) / _M
print(f"per-capita pi_F / pi_M from sector composition = {_piF/_piM:.3f}")
_chk("pi_F/pi_M sectoral", _piF / _piM, 0.86, 2)
# Exercise 18.2(d): female share in the table's male-dominated (>60% male) sectors, "about 10%"
_md = [s for s in FEM_TEXT if _fem[s] < 0.4]
_fs = sum(_fem[s] * _c17["emp_k"][s] for s in _md) / sum(_c17["emp_k"][s] for s in _md)
print(f"male-dominated sectors {_md}: female share {100*_fs:.1f}%")
_chk("female share in male-dominated sectors", 100 * _fs, 10, 0)
# Findex account ownership (WDI), Example ch18:ex:gender_aspiration
for _k, _v in {"female_2021": 63.3, "male_2021": 74.0, "female_2024": 53.5, "male_2024": 60.3}.items():
    assert _facts[f"BOL_account_ownership_{_k}"] == _v, _k
print("Findex account ownership 2021 F/M = 63.3/74.0; 2024 F/M = 53.5/60.3 (bolivia_facts.csv)")

print("\n=== Proof pass: Prop. assortative (b) with three classes (constrained A_S) ===")
for trial in range(300):
    sizes = rng.integers(1, 5, size=3)
    grp = np.repeat(np.arange(3), sizes)
    n_ = len(grp)
    Km = rng.random((n_, n_)); Km = (Km + Km.T) / 2
    lam = rng.random()
    same = grp[:, None] == grp[None, :]
    K = Km * np.where(same, 1.0, lam)
    for g0 in range(3):
        S = grp == g0
        C = Km[np.ix_(S, ~S)].sum()
        AS, ASc = K[np.ix_(S, S)].sum(), K[np.ix_(~S, ~S)].sum()
        assert gap(rownorm(K)) <= 2 * lam * C / (min(AS, ASc) + lam * C) + 1e-12
print("bound holds for every class S in 300 random three-class economies")

print("\n=== Proof pass: Prop. gender formation / attribute trap (closed form vs iteration) ===")
for trial in range(200):
    rho_, mu_, tb, rs = rng.uniform(0, 0.99), rng.uniform(0, 1), rng.normal(), rng.normal()
    lam = rng.uniform(0, 1); rF, rM = rng.normal(size=2)
    Wg = np.array([[1, lam], [lam, 1]]) / (1 + lam)
    sbar = Wg @ np.array([rF, rM])
    th = rng.normal(size=2) * 5
    for t in range(4000):
        th = rho_ * th + (1 - rho_) * tb + mu_ * (sbar - rs)
    closed = tb + mu_ / (1 - rho_) * (sbar - rs)
    assert np.allclose(th, closed, atol=1e-8)
    assert np.isclose(th[0] - th[1], mu_ / (1 - rho_) * (1 - lam) / (1 + lam) * (rF - rM), atol=1e-8)
    assert np.isclose(np.sort(np.linalg.eigvals(Wg).real)[0], (1 - lam) / (1 + lam))
print("steady-state gap M (1-lam)/(1+lam) (r_F - r_M) verified in 200 random draws")

print("\n=== Proof pass: Prop. homophily (ii)-(iii) at agent level ===")
for trial in range(100):
    nA, nB = rng.integers(1, 30, size=2); lam = rng.uniform(0.01, 1)
    g_ = np.array([1] * nA + [0] * nB)
    K = np.where(g_[:, None] == g_[None, :], 1.0, lam)
    pv = perron(rownorm(K)); pA = nA / (nA + nB)
    assert np.isclose(pv[:nA].mean() / pv[nA:].mean(), degree_ratio(pA, lam))
    aF, aM_, c_ = rng.uniform(0.1, 2, size=3); n_ = rng.integers(1, 10)
    K = np.full((2 * n_, 2 * n_), lam * c_); K[:n_, :n_] = aF; K[n_:, n_:] = aM_
    pv = perron(rownorm(K))
    assert np.isclose(pv[:n_].mean() / pv[n_:].mean(), (aF + lam * c_) / (aM_ + lam * c_))
print("degree-ratio formulas verified in 100 random draws")
