"""Chapter 6: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch06.py
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
    """Left Perron vector of an irreducible nonnegative W, normalized to sum 1."""
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmax(vals.real)
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def gap(W):
    """Gamma = 1 - |lambda_2| (Definition ch06:def:gap)."""
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


def neff(p):
    return 1 / np.sum(p ** 2)


def varpi(p):
    N = len(p)
    return np.mean((p - 1 / N) ** 2)


def gini(p):
    p = np.sort(p)
    N = len(p)
    return np.sum((2 * np.arange(1, N + 1) - N - 1) * p) / (N * p.sum())


def lazy(W):
    return 0.5 * (np.eye(len(W)) + W)


def rownorm(A):
    return A / A.sum(1, keepdims=True)


print("=" * 70)
print("Section 6.1: five-agent economy (eq. ch06:eq:W5), row-stochastic")
W5 = np.array([[.4, .4, .1, .1, 0],
               [.3, .3, .2, .2, 0],
               [.1, .2, .3, .2, .2],
               [0, .1, .4, .3, .2],
               [0, .1, .1, .3, .5]])
assert np.allclose(W5.sum(1), 1)
I5 = np.eye(5)
indeg = W5.sum(0)
katz = (np.linalg.inv(I5 - 0.5 * W5.T) - I5) @ np.ones(5)
bon = np.linalg.inv(I5 - 0.3 * W5.T) @ W5.T @ np.ones(5)
katz03 = (np.linalg.inv(I5 - 0.3 * W5.T) - I5) @ np.ones(5)
assert np.allclose(bon, katz03 / 0.3)          # Bonacich(beta) = Katz(beta)/beta
pi5 = perron(W5)
print("in-degree       ", np.round(indeg, 2))
print("Katz delta=0.5  ", np.round(katz, 2))
print("Bonacich beta=.3", np.round(bon, 2), " (= Katz(0.3)/0.3)")
print("pi              ", np.round(pi5, 3))
print("eigenvalues     ", np.round(np.linalg.eigvals(W5), 3))
print("Gamma =", round(gap(W5), 3), " N_eff =", round(neff(pi5), 2))
# normalized Katz -> pi as delta -> 1
for d in (0.9, 0.99, 0.999):
    kz = (1 - d) / 5 * np.linalg.inv(I5 - d * W5.T) @ np.ones(5)
    print(f"  normalized Katz delta={d}:", np.round(kz / kz.sum(), 4))
print("heatmap cells (200*W):", np.round(200 * W5).astype(int).tolist())
print("power bars (pct = 200*pi):", np.round(200 * pi5).astype(int).tolist())
# Exercise 6.1(c): W31 0.1 -> 0.2, taking 0.1 from W33
W5c = W5.copy(); W5c[2, 0] += 0.1; W5c[2, 2] -= 0.1
pic = perron(W5c)
print("Ex 6.1(c) new pi:", np.round(pic, 3), " d pi_1 =", round(pic[0] - pi5[0], 3))
assert pic[0] >= pi5[0]                          # Prop ch06:prop:monotone

print("=" * 70)
print("Prop ch06:prop:equal_power: counterexample for non-stochastic W")
Wn = np.array([[.2, .6], [.4, 0.]])
print("col sums", Wn.sum(0), "row sums", Wn.sum(1), "pi", perron(Wn),
      "rho", max(abs(np.linalg.eigvals(Wn))))
assert np.allclose(perron(Wn), [.5, .5])

print("=" * 70)
print("Prop ch06:prop:monotone: random check")
rng = np.random.default_rng(0)
for _ in range(2000):
    W = rownorm(rng.random((6, 6)))
    j, i, k = rng.choice(6, 3, replace=False)
    if rng.random() < 0.3:
        j = i                                      # also test j = i
    t = rng.random() * W[j, k]
    W2 = W.copy(); W2[j, i] += t; W2[j, k] -= t
    assert perron(W2)[i] >= perron(W)[i] - 1e-12
print("ok: moving weight in row j toward i never lowers pi_i (2000 draws)")

print("=" * 70)
print("Prop ch06:prop:star_power: pi_h = a/(a+1-s); spectrum {1, s-a, 1-a}")


def star(N, a, s):
    W = np.zeros((N, N))
    W[0, 0] = s; W[0, 1:] = (1 - s) / (N - 1)
    for i in range(1, N):
        W[i, 0] = a; W[i, i] = 1 - a
    return W


for N in (5, 10, 50):
    for a, s in ((1, 0), (0.5, 0.5), (0.9, 0.5), (0.5, 0.9)):
        W = star(N, a, s)
        p = perron(W)
        G = 1 - max(abs(s - a), 1 - a) if N > 2 else None
        assert np.isclose(p[0], a / (a + 1 - s))
        assert np.isclose(gap(W), G)
        print(f"N={N:2d} a={a} s={s}: pi_h={p[0]:.3f} Gamma={gap(W):.3f} N_eff={neff(p):.2f}")

print("=" * 70)
print("Prop ch06:prop:block_power: exact block masses under uniform leakage")
A1 = rownorm(rng.random((3, 3))); A2 = rownorm(rng.random((4, 4)))
R12 = rownorm(rng.random((3, 4))); R21 = rownorm(rng.random((4, 3)))
e1, e2 = 0.05, 0.02
WB = np.block([[(1 - e1) * A1, e1 * R12], [e2 * R21, (1 - e2) * A2]])
pB = perron(WB)
print("block-1 mass", round(pB[:3].sum(), 4), " predicted e2/(e1+e2) =", round(e2 / (e1 + e2), 4))
assert np.isclose(pB[:3].sum(), e2 / (e1 + e2))
assert np.any(np.isclose(np.linalg.eigvals(WB), 1 - e1 - e2))
print("lambda = 1-e1-e2 =", 1 - e1 - e2, "is an eigenvalue, so Gamma <=", e1 + e2,
      "; Gamma =", round(gap(WB), 4))

print("=" * 70)
print("Prop ch06:prop:welfare_revisited: W_mkt = 1 s'")
s = np.array([.5, .3, .2])
Wm = np.outer(np.ones(3), s)
print("pi(1s') =", perron(Wm), " Gamma =", round(gap(Wm), 3))

print("=" * 70)
print("Prop ch06:prop:gap-power-independent: W = g 1 pi' + (1-g) I")
for pi_, g in ((np.r_[0.91, np.full(9, 0.01)], 1.0),
               (np.r_[0.91, np.full(9, 0.01)], 0.05),
               (np.full(10, 0.1), 0.05)):
    W = g * np.outer(np.ones(10), pi_) + (1 - g) * np.eye(10)
    assert np.allclose(perron(W), pi_) and np.isclose(gap(W), g)
    print(f"  pi_1={pi_[0]:.2f} g={g}: Gamma={gap(W):.3f} N_eff={neff(pi_):.2f}")
eps = 0.01
Wb = np.full((10, 10), eps / 5); Wb[:5, :5] = Wb[5:, 5:] = (1 - eps) / 5
print("  doubly stochastic weak blocks eps=0.01: Gamma =", round(gap(Wb), 3),
      " pi uniform:", np.allclose(perron(Wb), 0.1))
print("  Var(pi) max (N-1)/N^2 at N=10:", 9 / 100, " Var(e_1) =", varpi(np.eye(10)[0]))

print("=" * 70)
print("Example ch06:ex:gap: topology table (N=10) and viz panel (N=8)")


def topos(N, eps=0.01):
    T = {}
    T['Complete'] = np.ones((N, N)) / N
    S = np.roll(np.eye(N), 1, 1)
    T['Ring (lazy)'] = lazy(0.5 * (S + S.T))
    A = S + S.T
    A[0, N // 2] = A[N // 2, 0] = 1
    A[2, 2 + N // 2] = A[2 + N // 2, 2] = 1
    T['Small-world (lazy)'] = lazy(rownorm(A))
    St = np.zeros((N, N)); St[1:, 0] = 1; St[0, 1:] = 1 / (N - 1)
    T['Star (lazy)'] = lazy(St)
    h = N // 2
    B = np.full((N, N), eps / h); B[:h, :h] = B[h:, h:] = (1 - eps) / h
    T['Weak blocks'] = B
    c = 2
    C = np.zeros((N, N)); C[:c, :c] = 0.9 / c; C[:c, c:] = 0.1 / (N - c)
    for i in range(c, N):
        C[i, i] = 0.5; C[i, :c] = 0.5 / c
    T['Core-periphery'] = C
    p = np.r_[0.91, np.full(N - 1, 0.09 / (N - 1))]
    T['Broadcast hub'] = np.outer(np.ones(N), p)
    return T


for N in (10, 8):
    print(f"N = {N}")
    for k, W in topos(N).items():
        assert np.allclose(W.sum(1), 1)
        p = perron(W)
        print(f"  {k:20s} Gamma={gap(W):.3f}  max pi={p.max():.3f}  N_eff={neff(p):.2f}  Gini={gini(p):.2f}")
        if N == 8:
            print("      pi =", np.round(p, 3))
# non-lazy ring / star are periodic: Gamma = 0
S = np.roll(np.eye(10), 1, 1)
St = np.zeros((10, 10)); St[1:, 0] = 1; St[0, 1:] = 1 / 9
print("non-lazy ring Gamma =", round(gap(0.5 * (S + S.T)), 3), "; non-lazy star Gamma =", round(gap(St), 3))

print("=" * 70)
print("Example ch06:ex:rich_poor")
Wrp = np.array([[.8, .2], [.9, .1]])
print("pi =", perron(Wrp), " = (9/11, 2/11) =", (9 / 11, 2 / 11))

print("=" * 70)
print("Exercise answers")
print("6.2(b) pi =", perron(np.array([[.5, .5], [.8, .2]])))
W62c = np.array([[.5, .5, 0], [0, .5, .5], [.5, 0, .5]])
print("6.2(c) col sums", W62c.sum(0), "pi", perron(W62c))
print("6.2(d) pi", perron(np.array([[.4, .3, .3], [.3, .4, .3], [.3, .3, .4]])))
Wbs = np.array([[0, 0, .5, .5], [0, 0, .5, .5], [.5, .5, 0, 0], [.5, .5, 0, 0]])
I4 = np.eye(4)
for b in (0.3, -0.3):
    print(f"6.4 Bonacich beta={b}:", np.linalg.inv(I4 - b * Wbs.T) @ Wbs.T @ np.ones(4))
print("6.4 pi =", perron(Wbs), " Gamma =", round(gap(Wbs), 3), "(periodic)")
W65 = np.array([[.6, .4], [.2, .8]])
v = np.ones(2) / 2
for t in range(1, 6):
    v = W65.T @ v; v /= v.sum(); print(f"6.5 t={t}", np.round(v, 4))
print("6.5 pi =", perron(W65), " lambda_2 =", np.sort(np.linalg.eigvals(W65))[0], " Gamma =", round(gap(W65), 3))

print("=" * 70)
print("Lemma ch06:lem:perturbation")
Wr = rownorm(rng.random((6, 6)))
E = rng.standard_normal((6, 6)); E -= E.mean(1, keepdims=True)
p0 = perron(Wr)
Z = np.linalg.inv(np.eye(6) - Wr + np.outer(np.ones(6), p0))
t = 1e-6
fd = (perron(Wr + t * E) - perron(Wr - t * E)) / (2 * t)
lem = p0 @ E @ Z
lam, V = np.linalg.eig(Wr); U = np.linalg.inv(V); k1 = np.argmax(lam.real)
eigform = sum((p0 @ E @ V[:, k]) / (1 - lam[k]) * U[k] for k in range(6) if k != k1).real
oldform = sum(U[k] * (U[k] @ E @ V[:, k1]) / (1 - lam[k]) for k in range(6) if k != k1).real
print("finite diff   ", fd)
print("pi'EZ         ", lem)
print("eigen form    ", eigform)
print("old formula   ", oldform, "(identically zero since E1 = 0)")
assert np.allclose(fd, lem, atol=1e-6) and np.allclose(lem, eigform) and np.allclose(oldform, 0)
# complex-safe d|lambda|: a matrix whose lambda_2 is a complex pair
Wc = np.array([[.5, .3, .2], [.1, .6, .3], [.2, .2, .6]])
Ec = np.array([[-1., 1., 0.], [0., -1., 1.], [1., 0., -1.]])     # rows sum to 0
lam, V = np.linalg.eig(Wc); U = np.linalg.inv(V)
o = np.argsort(-np.abs(lam)); l2 = lam[o[1]]; v2 = V[:, o[1]]; u2 = U[o[1]]
dl = u2 @ Ec @ v2
dabs = (np.conj(l2) * dl).real / abs(l2)
wrong = (np.sign(l2.real) * dl).real
fdabs = (np.sort(np.abs(np.linalg.eigvals(Wc + t * Ec)))[-2] - np.sort(np.abs(np.linalg.eigvals(Wc - t * Ec)))[-2]) / (2 * t)
print("old sign(lambda)*dlambda formula:", round(wrong, 6))
print("lambda_2 =", np.round(l2, 4), " d|lambda_2|: formula", round(dabs, 6), " finite diff", round(fdabs, 6))
assert np.isclose(dabs, fdabs, atol=1e-5)
# gradient intuition: two path-like blocks joined by a weak bridge (nodes 4-5)
N = 10
A = np.zeros((N, N))
for b in (0, 5):
    for i in range(5):
        for j in range(5):
            if i != j and abs(i - j) <= 2:
                A[b + i, b + j] = 1
A[4, 5] = A[5, 4] = 0.3
A += np.eye(N)
Wg = rownorm(A)
lam, V = np.linalg.eig(Wg); o = np.argsort(-np.abs(lam)); U = np.linalg.inv(V)
print("bridge example: Gamma =", round(gap(Wg), 4))
print("  |v2| (scaled)", np.round(np.abs(V[:, o[1]].real) / np.abs(V[:, o[1]]).max(), 2))
print("  |u2| (scaled)", np.round(np.abs(U[o[1]].real) / np.abs(U[o[1]]).max(), 2))


def dgap(W, j, i, k, t=1e-7):
    E = np.zeros_like(W); E[j, i] += 1; E[j, k] -= 1
    return (gap(W + t * E) - gap(W)) / t


print("  dGamma, core 0 -> core 9 (from self-loop):", round(dgap(Wg, 0, 9, 0), 3))
print("  dGamma, bridge 4 -> bridge 5 (from self-loop):", round(dgap(Wg, 4, 5, 4), 3))
assert dgap(Wg, 0, 9, 0) > dgap(Wg, 4, 5, 4)

print("=" * 70)
print("Proof checks: Kac formula (Prop ch06:prop:monotone), Birkhoff contraction and")
print("periodic eigenvalues (Appendix 6B), lumpable star / block spectra")
for _ in range(500):
    n = rng.integers(3, 8)
    W = rownorm(rng.random((n, n)) * (rng.random((n, n)) < 0.6) + 0.3 * np.roll(np.eye(n), 1, 1))
    p = perron(W)
    i = rng.integers(n)
    keep = [l for l in range(n) if l != i]
    Q = W[np.ix_(keep, keep)]
    assert max(abs(np.linalg.eigvals(Q))) < 1
    h = np.linalg.solve(np.eye(n - 1) - Q, np.ones(n - 1))
    m_i = 1 + W[i, keep] @ h
    assert np.isclose(p[i], 1 / m_i)                       # Kac: pi_i = 1 / m_i
print("ok: pi_i = 1/m_i on 500 random irreducible chains")


def hilbert(x, y):
    r = x / y
    return np.log(r.max() / r.min())


for _ in range(500):
    n = rng.integers(2, 6)
    Wp = rng.random((n, n)) + 0.01
    phi = min(Wp[i, k] * Wp[j, l] / (Wp[j, k] * Wp[i, l])
              for i in range(n) for j in range(n) for k in range(n) for l in range(n))
    tau = (1 - np.sqrt(phi)) / (1 + np.sqrt(phi))
    assert phi >= (Wp.min() / Wp.max()) ** 2 - 1e-12
    x, y = rng.random(n) + 1e-3, rng.random(n) + 1e-3
    assert hilbert(Wp @ x, Wp @ y) <= tau * hilbert(x, y) + 1e-10
print("ok: Birkhoff d_H(Wx, Wy) <= tau(W) d_H(x, y) on 500 draws")
for d in (2, 3, 4):
    # irreducible with period d: cyclic classes of size 2, edges only C_r -> C_{r+1}
    n = 2 * d
    W = np.zeros((n, n))
    for r in range(d):
        W[2 * r:2 * r + 2, 2 * ((r + 1) % d):2 * ((r + 1) % d) + 2] = rng.random((2, 2)) + 0.1
    W = rownorm(W)
    lam = np.linalg.eigvals(W)
    for k in range(d):
        assert np.min(np.abs(lam - np.exp(2j * np.pi * k / d))) < 1e-9
    assert np.isclose(gap(W), 0)
print("ok: period-d chains have the d-th roots of unity in their spectrum (Gamma = 0)")
for _ in range(200):
    n = rng.integers(3, 9)
    a, s = rng.random() * 0.98 + 0.01, rng.random() * 0.98
    W = star(n, a, s)
    ev = np.sort_complex(np.linalg.eigvals(W))
    target = np.sort_complex(np.r_[1, s - a, np.full(n - 2, 1 - a)].astype(complex))
    assert np.allclose(ev, target, atol=1e-7)
print("ok: star spectrum = {1, s-a, 1-a (N-2 times)} on 200 draws")
print("\n== Sourced Bolivia and comparator values quoted in the Prebisch box (Sec. gap diagnostic) and Exercise 6.x")
check("Gamma 2019 Bolivia, book12_isic (INE)", float(bol_summary(2019, "book12_isic")["Gamma"]), 0.475, 3)
check("power Gini 2019 Bolivia, book12_isic (INE)", float(bol_summary(2019, "book12_isic")["power_gini"]), 0.435, 3)
for iso, g in [("USA", .463), ("MEX", .649), ("KOR", .511)]:
    check("Gamma 2019 %s, ICIO 2023 book12_isic" % iso, float(oecd_summary(iso)["Gamma"]), g, 3)
cmp_g = [float(r["Gamma"]) for r in rows("oecd_io_summary.csv") if (r["edition"], r["year"], r["aggregation"], r["flows"], r["normalization"]) == ("ICIO2023", "2019", "book12_isic", "domestic", "input_shares") and r["country"] in ("KOR", "DEU", "USA", "CHL", "CHN", "MEX")]
assert min(cmp_g) < float(bol_summary(2019, "book12_isic")["Gamma"]) < max(cmp_g)
for iso, g in [("KOR", .658), ("DEU", .611)]:
    check("power Gini 2019 %s, ICIO 2023 book12_isic" % iso, float(oecd_summary(iso)["power_gini"]), g, 3)
    assert float(oecd_summary(iso)["power_gini"]) > float(bol_summary(2019, "book12_isic")["power_gini"])
gs = [float(bol_summary(2019, agg, nm)["Gamma"]) for agg in ("book12", "detail39")
      for nm in ("input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag")]
check("Gamma 2019 range over 8 variants: min", min(gs), 0.33, 2); check("Gamma 2019 range over 8 variants: max", max(gs), 0.69, 2)
check("pi hydrocarbons 2019, book12 W^P", bol_pi("Hydrocarbons"), 0.039, 3); assert int(bol_summary(2019, "book12")["N"]) == 11
check("hydrocarbons % of output exported 2019", 100 * sector_fact("Hydrocarbons", "export_share_of_output"), 70, 0)
check("Gamma 2019 Bolivia, book12 (Exercise)", float(bol_summary(2019, "book12")["Gamma"]), 0.467, 3)
check("N_eff 2019 Bolivia, book12 (Exercise)", float(bol_summary(2019, "book12")["Neff"]), 6.48, 2)
print("all assertions passed")
