"""Chapter 2: reproduces every number in the worked examples, tables, viz panels
and exercise answers.  Run: python3 code/ch02.py  (numpy only)."""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
def hdr(s): print("\n== " + s)
def U(x, k): x = np.asarray(x, float); return x - k / 2 * x ** 2
def perron(W):
    w, v = np.linalg.eig(W.T); p = np.real(v[:, np.argmin(abs(w - 1))]); return p / p.sum()
def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True); return 1 - e[1]
def neff(p): return 1 / np.sum(p ** 2)
def var(p): return np.mean((p - 1 / len(p)) ** 2)
def Lpower(k, N, Phi, p): return k * N / 2 * Phi ** 2 * var(p)
def Lext(alpha, k, m):
    """exact gain of the V-utilitarian optimum over equal division (quadratic U)"""
    Up = 1 - k * m; H = len(alpha) / np.sum(1 / alpha)
    z = Up * (alpha - H) / (k * alpha)             # optimal deviations from m
    return Up ** 2 / (2 * k) * np.sum((alpha - H) ** 2 / alpha), m + z

hdr("Prop doubly stochastic: counterexample without W_ii = 0")
W = np.array([[.8, .1, .1], [.1, .5, .4], [.4, .4, .2]])
print("rows", W.sum(1), "cols", W.sum(0), "c_j (off-diag col sums)", W.sum(0) - np.diag(W))

hdr("Example power_4: N=4, Omega=20, kappa=0.1, dbar=4.5 -> Phi=2")
k, N, Om, Phi = 0.1, 4, 20.0, 2.0
for name, p in [("competitive", np.full(4, .25)), ("mild", np.array([.30, .28, .22, .20])),
                ("star", np.array([.55, .15, .15, .15]))]:
    x = Om / N + (p - 1 / N) * Phi
    L = N * U(Om / N, k) - U(x, k).sum()
    print(name, "Var*1e3 =", round(var(p) * 1e3, 3), "N_eff =", round(neff(p), 3), "x =", x, "L =", round(L, 5),
          "formula =", round(Lpower(k, N, Phi, p), 5), "max kx", (k * x).max())
    assert np.isclose(L, Lpower(k, N, Phi, p)) and np.all(k * x < 1)
print("ratio star/mild =", round(var(np.array([.55, .15, .15, .15])) / var(np.array([.30, .28, .22, .20])), 2))

hdr("Example bond: N=3, gamma=0.5, kappa=0.5, Omega=3")
k, m, g = 0.5, 1.0, 0.5
WA = np.array([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]]); WB = np.array([[0, .6, 0], [.5, 0, 0], [0, .5, 0]])
for name, Wu in [("A", WA), ("B", WB)]:
    c = Wu.sum(0) - np.diag(Wu); alpha = 1 + g * c; Le, xV = Lext(alpha, k, m)
    approx = len(c) * (1 - k * m) ** 2 / (2 * k) * g ** 2 * np.var(c)
    print(name, "c =", c, "Var(c) =", round(np.var(c), 4), "alpha =", alpha, "L_externality exact =", round(Le, 5),
          "2nd-order approx =", round(approx, 5), "x^V =", xV.round(4), "max kx", (k * xV).max().round(3))
    assert np.isclose(xV.sum(), 3) and np.all(k * xV < 1)
    # brute-force check of the exact formula
    best = -1e9
    for a in np.linspace(0, 3, 601):
        for b in np.linspace(0, 3 - a, 301):
            xx = np.array([a, b, 3 - a - b]); best = max(best, np.sum(alpha * U(xx, k)))
    print("   grid check: max sum alpha U - sum alpha U(m) =", round(best - np.sum(alpha * U(np.full(3, m), k)), 4))

hdr("Example two blocks: kappa=0.1")
om_b = np.array([8, 3.2]); Nb = np.array([3, 5]); mbar = (Nb * om_b).sum() / Nb.sum()
bv = (Nb * (om_b - mbar) ** 2).sum() / Nb.sum(); print("mean", mbar, "between var", bv, "L_block", 0.1 / 2 * Nb.sum() * bv)

hdr("Viz: three failures (heatmap weights, rows normalized)")
def mk(f, n=8):
    K = np.array([[f(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)], float); return K / K.sum(1, keepdims=True)
Ws = mk(lambda i, j: 65 if j == 1 else 15 if i == 1 else 20 if i == j else 8)
Wc = mk(lambda i, j: 60 if (i <= 3 and j <= 3 and i != j) else 30 if i == j else 12)
Wb = mk(lambda i, j: 60 if (i <= 4 and j <= 4) or (i > 4 and j > 4) else 3)
p = perron(Ws); print("star: pi1 =", p[0].round(3), "N_eff =", round(neff(p), 2), "Gamma =", round(gap(Ws), 3))
c = Wc.sum(0) - np.diag(Wc); print("cluster: c =", c.round(3), "Var(c) =", round(np.var(c), 4))
p = perron(Wb); print("blocks: Gamma =", round(gap(Wb), 4), "N_eff =", round(neff(p), 2))

hdr("Transfer vs structural example: N=4, Omega=20, kappa=0.1, dbar=3 -> Phi=8")
def enclave(a, b):
    W = np.zeros((4, 4)); W[:2, :2] = (1 - 2 * b) / 2; W[2:, 2:] = (1 - 2 * a) / 2; W[:2, 2:] = b; W[2:, :2] = a; return W
k, N, Om, Phi = 0.1, 4, 20.0, 8.0
rows = {}
for name, (a, b) in [("no policy", (.02, .01)), ("integration x5", (.10, .05)), ("rebalancing", (.02, .015))]:
    W = enclave(a, b); p = perron(W); x = Om / N + (p - 1 / N) * Phi
    dpi = 0.5 * np.abs(p - 1 / N).sum()
    rows[name] = (p, x)
    print(name, "W rows", W[0].round(3), W[2].round(3), "pi =", p.round(4), "Gamma =", round(gap(W), 3),
          "N_eff =", round(neff(p), 3), "x =", x.round(3), "gap =", round(x.max() - x.min(), 3),
          "L_power =", round(Lpower(k, N, Phi, p), 4), "transfer volume Phi*d_pi =", round(Phi * dpi, 3))
    assert np.all(k * x < 1)
p0 = rows["no policy"][0]; tau = (1 / N - p0) * Phi; print("neutralizing transfer tau =", tau.round(3))
d = np.array([3, 1, 5, 3.]); x = d + p0 * (Om - d.sum()); print("asset transfer 2 from agent 2 to 3: x =", x.round(3))

hdr("Decomposition example: N=4, gamma=0.3, kappa=0.1, Omega=20, Phi=8")
k, g, m, Phi, N = 0.1, 0.3, 5.0, 8.0, 4
for name, p, c in [("connected", np.array([.26, .25, .25, .24]), np.array([1.02, 1.01, .99, .98])),
                   ("enclave", np.array([1, 1, .5, .5]) / 3, np.array([1.24, 1.30, .72, .74]))]:
    alpha = 1 + g * c; Le, xV = Lext(alpha, k, m); z = (p - 1 / N) * Phi; x = m + z
    Lp = Lpower(k, N, Phi, p)
    Li = k * g * Phi ** 2 / 2 * np.sum(c * (p - 1 / N) ** 2) - g * Phi * N * (1 - k * m) * np.mean((c - c.mean()) * (p - 1 / N))
    total = (np.sum(alpha * U(xV, k))) - np.sum(alpha * U(x, k))
    print(name, "Cov(pi,c) =", round(np.mean((c - c.mean()) * (p - p.mean())), 5), "L_externality =", round(Le, 5),
          "L_power =", round(Lp, 5), "L_interaction =", round(Li, 5), "sum =", round(Le + Lp + Li, 5), "direct =", round(total, 5))
    assert np.isclose(Le + Lp + Li, total) and np.all(k * xV < 1) and np.all(k * x < 1)

hdr("Gamma and pi independent (W = g 1pi' + (1-g)I)")
pi = np.array([.91] + [.01] * 9)
for gg in [1.0, .5, .1]:
    W = gg * np.outer(np.ones(10), pi) + (1 - gg) * np.eye(10)
    print("gamma", gg, "Gamma", round(gap(W), 4), "pi max", perron(W).max().round(3), "N_eff", round(neff(pi), 3))

hdr("Viz: four structures (heatmap weights, rows normalized)")
pan = {"competitive": lambda i, j: 80,
       "weakly dual": lambda i, j: 65 if (i <= 4 and j <= 4) else 55 if (i > 4 and j > 4) else 30,
       "enclave": lambda i, j: 65 if (i <= 4 and j <= 4) else 55 if (i > 4 and j > 4) else 4,
       "monopoly": lambda i, j: 85 if (j == 1 and i != 1) else 12 if i == 1 else 8 if i == j else 3}
for kk, f in pan.items():
    W = mk(f); p = perron(W); print(kk, "Gamma =", round(gap(W), 3), "N_eff =", round(neff(p), 2), "pi1 =", p[0].round(3),
                                   "L_power (N=8, kappa=0.1, Phi=8) =", round(Lpower(0.1, 8, 8.0, p), 4))

hdr("Exercises")
for M in [np.ones((3, 3)) / 3, np.array([[.6, .3, .1], [.1, .6, .3], [.3, .1, .6]]), np.array([[.8, .2], [.5, .5]]),
          np.array([[.4, .3, .3], [.3, .4, .3], [.3, .3, .4]])]:
    p = perron(M); print("2.1 cols", M.sum(0), "pi", p.round(4), "Var", round(var(p), 5))
W3 = np.array([[.5, .4, .1], [.2, .5, .3], [.1, .3, .6]]); p = perron(W3); c = W3.sum(0) - np.diag(W3)
k, g, m, Phi = 0.5, 0.5, 1.0, 1.0; alpha = 1 + g * c; Le, xV = Lext(alpha, k, m)
Li = k * g * Phi ** 2 / 2 * np.sum(c * (p - 1 / 3) ** 2) - g * Phi * 3 * (1 - k * m) * np.mean((c - c.mean()) * (p - 1 / 3))
print("2.3b pi", p.round(4), "Var", round(var(p), 5), "c", c, "Var(c)", round(np.var(c), 5),
      "Cov", round(np.mean((c - c.mean()) * (p - p.mean())), 5), "L_externality", round(Le, 5), "L_power", round(Lpower(k, 3, Phi, p), 5), "L_interaction", round(Li, 5))
n = 8; R = np.roll(np.eye(n), 1, axis=1); S = np.zeros((n, n)); S[0, 1:] = 1 / 7; S[1:, 0] = 1
B = np.zeros((n, n)); B[:4, :4] = .25; B[4:, 4:] = .25
for name, M in [("W^CE", np.ones((n, n)) / n), ("ring", R), ("blocks", B), ("star", S)]:
    lazy = .5 * (np.eye(n) + M); pl = perron(lazy) if name != "blocks" else None
    print("2.4", name, "Gamma =", round(gap(M), 4), "lazy Gamma =", round(gap(lazy), 4),
          "N_eff(lazy) =", None if pl is None else round(neff(pl), 3))
om = np.array([800.] * 5 + [4500.] * 3); mb = om.mean()
print("2.6 mean", mb, "L_block (kappa=1e-4) =", 1e-4 / 2 * np.sum((om - mb) ** 2), "max kappa*x", 1e-4 * 4500)
def elalto(eps):
    W = np.zeros((8, 8)); W[:5, 5:] = eps; W[5:, :5] = eps
    W[:5, :5] = (1 - 3 * eps) / 5; W[5:, 5:] = (1 - 5 * eps) / 3; return W
for eps in [.02, .15]:
    print("2.6 eps", eps, "Gamma =", round(gap(elalto(eps)), 3))

hdr("Furtado fixed-point example: eps' = F(eps) = ebar eps^2/(eps^2+h^2), ebar=0.3, h=0.1")
eb, h = .3, .1; F = lambda e: eb * e ** 2 / (e ** 2 + h ** 2); dF = lambda e: eb * 2 * e * h ** 2 / (e ** 2 + h ** 2) ** 2
roots = [0.0] + sorted(np.roots([1, -eb, h ** 2]).real)
for r in roots: print("fixed point", round(r, 4), "F'(e*) =", round(dF(r), 4), "stable" if abs(dF(r)) < 1 else "unstable")
for e0 in [.03, .05]:
    e = e0
    for t in range(60): e = F(e)
    print("start", e0, "-> after 60 periods", round(e, 4))

hdr("Proof checks (randomized, seed 2): Props FWT_matrix, doubly_stochastic, power_loss, ext_loss, block_loss, transfer_trap; Thm welfare_loss")
rng = np.random.default_rng(2)
for _ in range(400):
    n = rng.integers(2, 8)
    Wz = rng.random((n, n)); np.fill_diagonal(Wz, 0); Wz /= Wz.sum(1, keepdims=True)      # doubly_stochastic
    c = Wz.sum(0) - np.diag(Wz); assert np.isclose(c.sum(), n) and np.allclose(c, Wz.sum(0))
    P = np.roll(np.eye(n), 1, axis=1); Wd = 0.5 * P + 0.5 * P.T if n > 2 else P            # zero-diagonal doubly stochastic
    cd = Wd.sum(0) - np.diag(Wd); assert np.allclose(cd, 1)
    k = rng.random() * 0.3 + 0.01; m = rng.random() * 0.8 / k; Up = 1 - k * m; g = rng.random()  # ext_loss
    alpha = 1 + g * c; H = n / np.sum(1 / alpha)
    Le, xV = Lext(alpha, k, m); val = lambda x: np.sum(alpha * U(x, k))
    assert np.isclose(xV.sum(), n * m) and np.isclose(val(xV) - val(np.full(n, m)), Le)
    for _ in range(10):
        e = rng.normal(size=n); e -= e.mean(); e *= 0.1
        assert val(xV + e) <= val(xV) + 1e-12
    # FWT_matrix: at the blind optimum (equal division) the social gradient alpha_j U'(m) is equal iff c_j are equal
    grad = alpha * Up; assert np.allclose(grad, grad[0]) == np.allclose(c, c[0])
    pi = rng.random(n); pi /= pi.sum(); Phi = rng.random() * 2; z = (pi - 1 / n) * Phi; x = m + z   # power_loss + decomposition
    assert np.isclose(n * U(m, k) - U(x, k).sum(), Lpower(k, n, Phi, pi))
    Cov = np.mean((c - c.mean()) * (pi - 1 / n))
    Li = k * g * Phi ** 2 / 2 * np.sum(c * (pi - 1 / n) ** 2) - g * Phi * n * Up * Cov
    assert np.isclose(val(xV) - val(x), Le + Lpower(k, n, Phi, pi) + Li)
    nb = rng.integers(1, 4); sizes = rng.integers(1, 5, nb); mb = rng.random(nb) * 5             # block_loss
    om_all = np.repeat(mb, sizes); Wbar = om_all.mean()
    Lb = sizes.sum() * U(Wbar, k) - np.sum(sizes * U(mb, k))
    assert np.isclose(Lb, k / 2 * np.sum(sizes * (mb - Wbar) ** 2))
    d = rng.random(n); Om = d.sum() + rng.random() * 3 + .01; Ph = Om - d.sum(); xb = d + pi * Ph   # transfer_trap
    tau = rng.normal(size=n); tau -= tau.mean(); xt = (d + tau) + pi * (Om - (d + tau).sum())
    assert np.allclose(xt, xb + tau)
    tau_eq = Om / n - d - pi * Ph; assert np.isclose(tau_eq.sum(), 0) and np.allclose(d + tau_eq + pi * Ph, Om / n)
    tau_eq2 = (1 / n - pi) * Ph; assert np.isclose(0.5 * np.abs(tau_eq2).sum(), Ph * 0.5 * np.abs(pi - 1 / n).sum())
# second-order expansion of L_externality: L / gamma^2 -> N (1-km)^2/(2k) Var(c)
c = np.array([.5, 1.1, 0.]); k, m = .5, 1.
for g in [1e-2, 1e-3]:
    Le, _ = Lext(1 + g * c, k, m); print("gamma", g, "L_ext/gamma^2 =", Le / g ** 2, " limit =", 3 * (1 - k * m) ** 2 / (2 * k) * np.var(c))
assert abs(Lext(1 + 1e-4 * c, k, m)[0] / 1e-8 - 3 * (1 - k * m) ** 2 / (2 * k) * np.var(c)) < 1e-3
print("all randomized proof checks passed")
print("\nall assertions passed")
