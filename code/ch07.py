"""Chapter 7: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch07.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)

# ---- Sourced data (data/processed, retrieved 2026-11-05; see data/README.md) ----
import csv, pathlib
DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "processed"
def rows(name):
    with open(DATA / name, newline="", encoding="utf-8") as f: return list(csv.DictReader(f))
def fact(fid):
    return float(next(r["value"] for r in rows("bolivia_facts.csv") if r["fact_id"] == fid))
def bol_summary(year, agg, norm="input_shares", flows="domestic", base=None):
    base = base or ("2017" if int(year) >= 2017 else "1990")
    return next(r for r in rows("bolivia_io_summary.csv") if (r["base"], r["year"], r["aggregation"], r["flows"], r["normalization"]) == (base, str(year), agg, flows, norm))
def oecd_summary(iso, year=2019, agg="book12_isic", edition="ICIO2023"):
    return next(r for r in rows("oecd_io_summary.csv") if (r["edition"], r["country"], r["year"], r["aggregation"], r["flows"], r["normalization"]) == (edition, iso, str(year), agg, "domestic", "input_shares"))
def bol_pi(sector, year=2019, agg="book12"):
    return float(next(r["pi"] for r in rows("bolivia_io_pi.csv") if (r["base"], r["year"], r["aggregation"], r["flows"], r["normalization"], r["sector"]) == ("2017", str(year), agg, "domestic", "input_shares", sector)))
def sector_fact(sector, col, year=2019):
    return float(next(r[col] for r in rows("bolivia_io_sector_facts.csv") if (r["base"], r["year"], r["book12"]) == ("2017", str(year), sector)))
def check(label, value, text, nd):
    """Assert that the sourced value rounds to the number printed in the text (nd decimals)."""
    print("  %-58s %10.4f  text: %s" % (label, value, text)); assert round(value, nd) == round(text, nd), (label, value, text)


def perron(W):
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmax(vals.real)
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
    """TU asymmetric Nash: x_i = d_i + pi_i (Omega - sum d)."""
    return d + pi * (Omega - d.sum())


print("=" * 70)
print("Example ch07:ex:landlord")
A = np.sqrt(20)                        # TFP so that A*10^.5*2^.5 = 20
Y = A * 10 ** 0.5 * 2 ** 0.5
MPK, MPL = 0.5 * Y / 10, 0.5 * Y / 2
print("A =", round(A, 3), " Y =", round(Y, 3), " MPK =", MPK, " MPL =", MPL,
      " payments L,F =", MPK * 10, MPL * 2)
dL, dF, S = 10, 2, 20 - 12
for pL in (0.7, 0.5, 0.3):
    print(f"pi_L={pL}: x_L={dL + pL * S:.1f} x_F={dF + (1 - pL) * S:.1f}")
# Exercise 7.1: Y = A K^.3 L^.7 with Y = 20
print("Ex 7.1 competitive with K^.3 L^.7 (Y=20): L gets", 0.3 * 20, ", F gets", 0.7 * 20,
      "; needed pi_F =", (14 - 2) / 8, "> 1 -> impossible; L's payment 6 < d_L = 10")
xF = dF + 0.3 * S
print("Ex 7.1(b) pi=(0.7,0.3): x =", (dL + 0.7 * S, xF))

print("=" * 70)
print("Prop ch07:prop:allocation-rule check (random) and exact power loss")
rng = np.random.default_rng(0)
for _ in range(200):
    N = rng.integers(2, 7)
    pi = rng.dirichlet(np.ones(N)); d = np.full(N, rng.random()); Om = d.sum() + 1 + rng.random()
    x = rule(pi, d, Om)
    Phi = Om - N * d[0]
    assert np.allclose(x, Om / N + (pi - 1 / N) * Phi)
    # x maximizes sum pi ln(x - d) s.t. sum x = Omega: check against random feasible points
    f = lambda z: np.sum(pi * np.log(z - d))
    for _ in range(20):
        z = d + rng.dirichlet(np.ones(N)) * (Om - d.sum())
        assert f(z) <= f(x) + 1e-12
    k = 0.9 / x.max()
    U = lambda z: z - k / 2 * z ** 2
    assert np.isclose(np.sum(U(np.full(N, Om / N)) - U(x)), k * N / 2 * Phi ** 2 * varpi(pi))
print("ok: rule maximizes the weighted Nash product; L_power = (kN/2) Phi^2 Var(pi) exactly")

print("=" * 70)
print("Prop ch07:prop:nesting: symmetric Nash != Walras (Cobb-Douglas counterexample)")
eA, eB = np.array([.9, .3]), np.array([.1, .7])
dA = .5 * np.log(eA).sum(); dB = .5 * np.log(eB).sum()
t = np.linspace(1e-4, 1 - 1e-4, 2_000_001)            # contract curve x_A = y_A = t
SA, SB = np.log(t) - dA, np.log(1 - t) - dB
ok = (SA > 0) & (SB > 0)
f = np.where(ok, np.log(np.where(ok, SA, 1)) + np.log(np.where(ok, SB, 1)), -np.inf)
tN = t[f.argmax()]
# Walras: p_y = 1, p_x = p; demand x_A = (p*.9+.3)/(2p); clearing gives p = 1
print("symmetric Nash x_A = y_A =", round(tN, 3), "; Walras x_A = y_A = 0.600")
assert abs(tN - 0.631) < 1e-3

print("=" * 70)
print("Section ch07:sec:worked: Economies A and B")
kappa, Omega, N, dbar = 0.1, 20.0, 4, 3.0
Phi = Omega - N * dbar
WA = np.array([[.6, .2, .1, .1], [.2, .5, .2, .1], [.1, .3, .4, .2], [.1, .1, .2, .6]])
WB = np.array([[.48, .50, .01, .01], [.40, .58, .01, .01], [.03, .03, .64, .30], [.03, .03, .30, .64]])
U = lambda z: z - kappa / 2 * z ** 2
res = {}
for name, W in (("A", WA), ("B", WB)):
    assert np.allclose(W.sum(1), 1)
    p = perron(W)
    x = Omega / N + (p - 1 / N) * Phi
    assert np.all(kappa * x < 1)
    L = np.sum(U(np.full(N, Omega / N)) - U(x))
    assert np.isclose(L, kappa * N / 2 * Phi ** 2 * varpi(p))
    res[name] = (p, x, L)
    print(f"Economy {name}: pi={np.round(p, 3)} N_eff={neff(p):.2f} Var(pi)={varpi(p):.5f} Gamma={gap(W):.3f}")
    print(f"   x = {np.round(x, 2)}  L_power = {L:.4f}  col sums = {W.sum(0)}")
    print("   heatmap cells (100*W):", np.round(100 * W).astype(int).tolist())
    print("   allocation bars pct (10*x):", np.round(10 * x).astype(int).tolist())
print("Phi =", Phi, " ratio L_B/L_A =", round(res['B'][2] / res['A'][2], 1))
pB, xB, _ = res['B']
print("enclave power share =", round(pB[:2].sum(), 3), " predicted e2/(e1+e2) =", 0.06 / (0.02 + 0.06))
print("enclave consumption =", round(xB[:2].sum(), 2), "vs 10 (",
      round(100 * (xB[:2].sum() / 10 - 1), 1), "%)")
WBsym = np.array([[.48, .50, .01, .01], [.40, .58, .01, .01], [.01, .01, .68, .30], [.01, .01, .30, .68]])
print("B with symmetric leakage (0.02 both ways): pi =", np.round(perron(WBsym), 3), " Gamma =", round(gap(WBsym), 3))

print("=" * 70)
print("Exercise 7.2 (TU): Omega=3, d=(0.6,0.4,0.2), kappa=0.25")
d = np.array([.6, .4, .2]); pi = np.array([.5, .3, .2])
x = rule(pi, d, 3.0); xe = rule(np.ones(3) / 3, d, 3.0)
print("x(pi) =", x, " x(equal) =", np.round(xe, 3), " max kappa x =", 0.25 * x.max())

print("=" * 70)
print("Theorem ch07:thm:decomposition; Exercise 7.3 (gamma=0.3, W^U = W^P off-diagonal)")


def fb(w, kappa, Omega):
    """argmax sum w_j U(x_j) s.t. sum x = Omega, U = x - k/2 x^2:
    w_j (1 - k x_j) = mu."""
    # x_j = (1 - mu/w_j)/k ; sum = Omega
    mu = (len(w) - kappa * Omega) / np.sum(1 / w)
    return (1 - mu / w) / kappa


for name, W in (("A", WA), ("B", WB)):
    p, x, Lp = res[name]
    WU = W - np.diag(np.diag(W))
    g = 0.3
    c = WU.sum(0)                                  # c_j = sum_i W^U_ij
    w = 1 + g * c
    SW = lambda z: np.sum(w * U(z))
    xfb = fb(w, kappa, Omega)
    assert np.all(kappa * xfb < 1)
    Lext = SW(xfb) - SW(np.full(N, Omega / N))
    Lint = g * np.sum(c * (U(Omega / N) - U(x)))
    Ltot = SW(xfb) - SW(x)
    assert np.isclose(Ltot, Lp + Lext + Lint)
    up = 1 - kappa * Omega / N
    lin = -g * N * up * Phi * np.mean((c - c.mean()) * (p - 1 / N))
    quad = g * kappa / 2 * Phi ** 2 * np.sum(c * (p - 1 / N) ** 2)
    assert np.isclose(Lint, lin + quad)
    print(f"Economy {name}: c={np.round(c, 2)} L_power={Lp:.4f} L_ext={Lext:.4f} L_int={Lint:.4f} total={Ltot:.4f} Cov(c,pi)={np.mean((c-c.mean())*(p-1/N)):.5f}")

print("=" * 70)
print("Exercise 7.4: transfer vs structural change in Economy B")
pB = res['B'][0]
dT = np.array([2., 3, 3, 4])
xT = rule(pB, dT, Omega)
print("transfer d'=(2,3,3,4): x =", np.round(xT, 2), " L_power-type loss =",
      round(np.sum(U(np.full(N, 5.)) - U(xT)), 4))
WB2 = WB.copy()
for i in range(4):
    for j in range(4):
        if (i < 2) != (j < 2):
            WB2[i, i] -= 0.05 - WB2[i, j]; WB2[i, j] = 0.05
assert np.allclose(WB2.sum(1), 1)
p2 = perron(WB2); x2 = Omega / N + (p2 - 1 / N) * Phi
print("structural (cross-block entries -> 0.05): pi =", np.round(p2, 3), " x =", np.round(x2, 2),
      " loss =", round(np.sum(U(np.full(N, 5.)) - U(x2)), 4), " Gamma:", round(gap(WB), 3), "->", round(gap(WB2), 3))

print("=" * 70)
print("Exercise 7.5 / Prop ch07:prop:escape_hard: star N=5, a=0.9, s=0.5")


def star(N, a, s):
    W = np.zeros((N, N))
    W[0, 0] = s; W[0, 1:] = (1 - s) / (N - 1)
    for i in range(1, N):
        W[i, 0] = a; W[i, i] = 1 - a
    return W


Ws = star(5, 0.9, 0.5); ps = perron(Ws)
print("pi =", np.round(ps, 4), " Gamma =", round(gap(Ws), 3))
target = 1.5 * ps[1]
a_new = (1 - 4 * target) * 0.5 / (1 - (1 - 4 * target))
print("spoke power +50% requires pi_h =", round(1 - 4 * target, 4), " i.e. a =", round(a_new, 3))
assert np.isclose(perron(star(5, a_new, 0.5))[1], target)

print("first-order escape bound: |d pi_P| <= ||E||_inf ||Z||_inf")
for eps in (0.2, 0.05, 0.01):
    h = 2
    Wb = np.full((4, 4), eps / h); Wb[:h, :h] = Wb[h:, h:] = (1 - eps) / h
    p = perron(Wb); Z = np.linalg.inv(np.eye(4) - Wb + np.outer(np.ones(4), p))
    print(f"  weak blocks eps={eps}: Gamma={gap(Wb):.3f} ||Z||_inf={np.abs(Z).sum(1).max():.1f}"
          f"  d(block-2 mass)/d(eps_1) = eps2/(eps1+eps2)^2 = {eps / (2 * eps) ** 2:.1f}")
# verify the bound on random perturbations
for _ in range(500):
    W = rng.random((5, 5)); W /= W.sum(1, keepdims=True)
    E = rng.standard_normal((5, 5)); E -= E.mean(1, keepdims=True)
    p = perron(W); Z = np.linalg.inv(np.eye(5) - W + np.outer(np.ones(5), p))
    dp = p @ E @ Z
    assert np.abs(dp).sum() <= np.abs(E).sum(1).max() * np.abs(Z).sum(1).max() + 1e-12
print("ok: bound holds on 500 random draws")

print("=" * 70)
print("Proof checks: Thm ch07:thm:decomposition (random c, pi) and Prop ch07:prop:comp_statics")
for _ in range(300):
    N = rng.integers(2, 7)
    pi = rng.dirichlet(np.ones(N))
    c = rng.random(N) if rng.random() < 0.7 else np.full(N, rng.random())   # also test equal c_j
    g = rng.random(); dbar = rng.random(); Om = N * dbar + 0.5 + rng.random()
    Phi = Om - N * dbar
    x = rule(pi, np.full(N, dbar), Om)
    w = 1 + g * c
    k = 0.5 / x.max()
    xfb = fb(w, k, Om)
    if np.any(k * xfb >= 1):
        continue
    U = lambda z: z - k / 2 * z ** 2
    SW = lambda z: np.sum(w * U(z))
    xE = np.full(N, Om / N)
    for _ in range(10):                            # x^FB maximizes S on the budget line
        z = xfb + rng.standard_normal(N) * 0.1; z += (Om - z.sum()) / N
        assert SW(z) <= SW(xfb) + 1e-12
    Lp = k * N / 2 * Phi ** 2 * varpi(pi)
    Lext = SW(xfb) - SW(xE)
    cov = np.mean((c - c.mean()) * (pi - 1 / N))
    Lint = -g * N * (1 - k * Om / N) * Phi * cov + g * k / 2 * Phi ** 2 * np.sum(c * (pi - 1 / N) ** 2)
    assert np.isclose(SW(xfb) - SW(x), Lp + Lext + Lint)
    assert Lext >= -1e-12
    if np.allclose(g * c, g * c[0]):               # zero iff all gamma c_j equal
        assert abs(Lext) < 1e-10 and np.allclose(xfb, xE)
    else:
        assert Lext > 0 and not np.allclose(xfb, xE)
print("ok: L = L_power + L_ext + L_int exactly, L_ext = 0 iff gamma*c_j all equal (300 draws)")
for _ in range(200):
    N = rng.integers(2, 6)
    pi = rng.dirichlet(np.ones(N)); d = rng.random(N); Om = d.sum() + 1 + rng.random()
    h = 1e-6
    J = np.array([(rule(pi, d + h * np.eye(N)[j], Om) - rule(pi, d - h * np.eye(N)[j], Om)) / (2 * h)
                  for j in range(N)]).T          # J[i, j] = dx_i / dd_j
    assert np.allclose(np.diag(J), 1 - pi) and np.allclose(J - np.diag(np.diag(J)),
                                                            -np.outer(pi, np.ones(N)) + np.diag(pi))
    assert np.allclose((rule(pi, d, Om + h) - rule(pi, d, Om - h)) / (2 * h), pi)
    W = rng.random((N, N)); W /= W.sum(1, keepdims=True)
    E = rng.standard_normal((N, N)); E -= E.mean(1, keepdims=True)
    p = perron(W); Z = np.linalg.inv(np.eye(N) - W + np.outer(np.ones(N), p))
    dx = (rule(perron(W + h * E), d, Om) - rule(perron(W - h * E), d, Om)) / (2 * h)
    assert np.allclose(dx, (Om - d.sum()) * (p @ E @ Z), atol=1e-6)
print("ok: dx_i/dd_i = 1 - pi_i, dx_i/dd_j = -pi_i, dx/dOmega = pi, dx/dt = (Omega - sum d) pi'EZ")
print("\n== Sourced values quoted after Economy B (Sec. worked economies)")
s19 = bol_summary(2019, "book12"); assert int(s19["N"]) == 11
check("pi hydrocarbons 2019, book12 W^P", bol_pi("Hydrocarbons"), 0.039, 3)
check("uniform share 1/11", 1 / 11, 0.091, 3); assert bol_pi("Hydrocarbons") < 1 / 11
check("hydrocarbons % of output exported 2019", 100 * sector_fact("Hydrocarbons", "export_share_of_output"), 70, 0)
check("power Gini 2019 Bolivia, book12_isic (INE)", float(bol_summary(2019, "book12_isic")["power_gini"]), 0.435, 3)
check("power Gini 2019 KOR, ICIO 2023 book12_isic", float(oecd_summary("KOR")["power_gini"]), 0.658, 3)
check("power Gini 2019 DEU, ICIO 2023 book12_isic", float(oecd_summary("DEU")["power_gini"]), 0.611, 3)
print("all assertions passed")
