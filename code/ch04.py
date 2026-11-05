"""Chapter 4: reproduces the numbers used in the propositions, viz panels and
exercise answers.  Run: python3 code/ch04.py  (numpy only)."""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
def hdr(s): print("\n== " + s)
def perron(W):
    w, v = np.linalg.eig(W.T); p = np.real(v[:, np.argmin(abs(w - 1))]); return p / p.sum()
def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True); return 1 - e[1]
def neff(p): return 1 / np.sum(p ** 2)

hdr("Prop A2_growth: comparative static of the steady state, check by finite differences")
alpha, s, dlt, g = np.array([.4, .3]), np.array([.2, .2]), np.array([.05, .05]), 1.0
A = np.array([[.10, .25], [.20, .10]]); Lab = np.ones(2)
f = lambda K: K ** alpha * Lab ** (1 - alpha)
def steady(A):
    K = np.ones(2) * 5
    for _ in range(20000):
        K = np.diag(s) @ np.linalg.inv(np.eye(2) - g * A.T) @ f(K) + (1 - dlt) * K
    return K
K0 = steady(A); L = np.linalg.inv(np.eye(2) - g * A.T); Fp = np.diag(alpha * f(K0) / K0)
dA = np.array([[0, .01], [0, 0]])
pred = np.linalg.solve(np.diag(dlt) - np.diag(s) @ L @ Fp, np.diag(s) @ L @ (g * dA.T) @ L @ f(K0))
fd = steady(A + dA) - K0
print("K* =", K0.round(4), " predicted dK* =", pred.round(4), " finite difference =", fd.round(4))
assert np.allclose(pred, fd, rtol=0.05)

hdr("Prop A3: DeGroot consensus = pi' b0")
W = np.array([[.5, .3, .2], [.1, .7, .2], [.2, .1, .7]]); b = np.array([1., 0., 0.]); p = perron(W)
print("W^100 b0 =", (np.linalg.matrix_power(W, 100) @ b).round(6), " pi'b0 =", round(p @ b, 6))

hdr("Ex 4.2")
print("col sums", W.sum(0), " ||1'W - 1'||_2 =", round(np.linalg.norm(W.sum(0) - 1), 4), " pi =", p.round(4),
      " Var(pi) =", round(np.mean((p - 1 / 3) ** 2), 5), " N_eff =", round(neff(p), 3))

hdr("Market interaction with unequal wealth: pi = s even under A1-A4")
s_w = np.array([.425, .3625, .2125]); Wm = np.outer(np.ones(3), s_w); print("pi(1s') =", perron(Wm), " Gamma =", round(gap(Wm), 4))

hdr("Var(pi) maximum (N-1)/N^2")
for n in [3, 4, 9, 20]: print("N =", n, "max Var =", round((n - 1) / n ** 2, 4))

hdr("Viz: five departures (heatmap weights, rows normalized)")
def mk(fn, n=8):
    K = np.array([[fn(i, j) for j in range(1, n + 1)] for i in range(1, n + 1)], float); return K / K.sum(1, keepdims=True)
pan = {"W^CE": lambda i, j: 80, "A1 cluster": lambda i, j: 95 if (i <= 3 and j <= 3 and i != j) else 50,
       "A2 blocks": lambda i, j: 65 if (i <= 4 and j <= 4) else 55 if (i > 4 and j > 4) else 5,
       "A3-A5 hub": lambda i, j: 70 if j == 1 else 15 if i == 1 else 25 if i == j else 6}
for k, fn in pan.items():
    Wk = mk(fn); pk = perron(Wk); print(k, "Gamma =", round(gap(Wk), 4), "N_eff =", round(neff(pk), 3), "pi1 =", pk[0].round(3))

hdr("Appendix 4A: pi of a composite is continuous in alpha (irreducible layers)")
WPl = mk(pan["A2 blocks"]); WIl = mk(pan["A3-A5 hub"]); prev = None; jump = 0
for a in np.linspace(0, 1, 1001):
    pa = perron(a * WPl + (1 - a) * WIl)
    if prev is not None: jump = max(jump, np.abs(pa - prev).max())
    prev = pa
print("max change of pi between alpha steps of 0.001:", round(jump, 5))

hdr("Proof checks (randomized, seed 4): Props A1_cost, A2_growth (dL formula), A3_beliefs, A4_loss, A5_decomp")
rng = np.random.default_rng(4)
for _ in range(300):
    n = rng.integers(2, 7)
    A_r = rng.random((n, n)); A_r *= rng.random() * 0.9 / max(abs(np.linalg.eigvals(A_r))); g_r = rng.random() * 0.9 + 0.1
    Lr = lambda A: np.linalg.inv(np.eye(n) - g_r * A.T); E = rng.normal(size=(n, n)) * 1e-6        # dL = L (g dA') L
    assert np.allclose(Lr(A_r + E) - Lr(A_r), Lr(A_r) @ (g_r * E.T) @ Lr(A_r), atol=1e-10)
    Wr = rng.random((n, n)) + .01; Wr /= Wr.sum(1, keepdims=True); pr = perron(Wr); P = np.outer(np.ones(n), pr)  # A3
    assert np.allclose(np.linalg.matrix_power(Wr - P, 7), np.linalg.matrix_power(Wr, 7) - P)
    ev = sorted(abs(np.linalg.eigvals(Wr)), reverse=True)
    assert np.isclose(max(abs(np.linalg.eigvals(Wr - P))), ev[1]) and ev[1] < 1
    b0 = rng.normal(size=n); assert np.allclose(np.linalg.matrix_power(Wr, 400) @ b0, (pr @ b0) * np.ones(n))
    lay = [rng.random((n, n)) for _ in range(3)]; lay = [M / M.sum(1, keepdims=True) for M in lay]       # A5
    al = rng.random(3); al /= al.sum(); Wc = sum(a * M for a, M in zip(al, lay)); pc = perron(Wc)
    assert np.allclose(Wc.sum(1), 1) and np.allclose(pc, sum(a * (pc @ M) for a, M in zip(al, lay)))
    k = rng.random() * 0.3 + .01; Phi = rng.random() * 2; m = 1.0                                         # A4
    x = m + (pr - 1 / n) * Phi; Lp = n * (m - k / 2 * m ** 2) - np.sum(x - k / 2 * x ** 2)
    assert np.isclose(Lp, k * n / 2 * Phi ** 2 * np.mean((pr - 1 / n) ** 2))
# A1: L_ext / gamma^2 -> Var(c) Psi  (exact formula from Ch2)
c = np.array([.5, 1.1, 0.]); k, m = .5, 1.; Psi = 3 * (1 - k * m) ** 2 / (2 * k)
def Lext(al):
    H = len(al) / np.sum(1 / al); return (1 - k * m) ** 2 / (2 * k) * np.sum((al - H) ** 2 / al)
for gU in [1e-1, 1e-2, 1e-3]: print("gamma_U", gU, "L_ext/gamma^2 =", round(Lext(1 + gU * c) / gU ** 2, 6), " Var(c) Psi =", round(np.var(c) * Psi, 6))
assert abs(Lext(1 + 1e-4 * c) / 1e-8 - np.var(c) * Psi) < 1e-3
print("all randomized proof checks passed")
print("\nall assertions passed")
