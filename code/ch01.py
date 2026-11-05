"""Chapter 1: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run: python3 code/ch01.py  (numpy only)."""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
def hdr(s): print("\n== " + s)
def U(x, k): return x - k / 2 * np.asarray(x, float) ** 2
def perron(W):
    w, v = np.linalg.eig(W.T); p = np.real(v[:, np.argmin(abs(w - 1))]); return p / p.sum()
def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True); return 1 - e[1]

hdr("Allocation / endowments")
X = np.array([[4, 2], [3, 5], [1, 3.]]); Om = np.array([[6, 1], [1, 6], [1, 3.]])
print("1'X =", X.sum(0), " 1'Omega =", Om.sum(0)); assert np.allclose(X.sum(0), Om.sum(0))

hdr("Cobb-Douglas at p=(1,1) and equilibrium")
w = Om @ np.array([1, 1.]); print("w(p=1,1) =", w, " demand col sums =", (np.outer(w, [.5, .5])).sum(0))
p2 = (16 - 8) / 10; p = np.array([1, p2]); w = Om @ p
Xs = np.column_stack([w / 2, w / (2 * p2)])
print("p2* =", p2, " w =", w, "\nX* =\n", Xs, "\ncol sums", Xs.sum(0))
assert np.allclose(Xs.sum(0), Om.sum(0))
s = w / w.sum(); print("wealth shares s =", s, " N_eff(s) =", 1 / np.sum(s ** 2))
Wmkt = np.outer(np.ones(3), s); print("W^mkt = 1s': pi =", perron(Wmkt), " Gamma =", round(gap(Wmkt), 6))

hdr("One-good example (kappa=0.5, Omega=3, omega=(1.8,.6,.6))")
k = 0.5; om = np.array([1.8, .6, .6]); N = 3
WFB = N * U(1.0, k); Wno = U(om, k).sum()
print("U(omega) =", U(om, k), " W_FB =", WFB, " W_noTrade =", Wno, " L =", WFB - Wno,
      " share lost =", (WFB - Wno) / WFB)
print("marginal utilities at omega:", 1 - k * om)
V = np.mean(om ** 2) - np.mean(om) ** 2; print("Var(omega) =", V, " (kN/2)Var =", k * N / 2 * V)
assert np.isclose(WFB - Wno, k * N / 2 * V) and np.all(k * om < 1)

hdr("TU Nash bargaining: harvest 3, d=(0.6,0.3,0.3)")
Omega = 3.0; d = np.array([.6, .3, .3]); Phi = Omega - d.sum()
def nash(lam, d, Omega): return d + lam * (Omega - d.sum())
# brute-force check that the closed form maximizes sum lam*ln(x-d)
lam = np.array([.6, .2, .2]); g = np.linspace(0, Omega, 601); best = (-1e9, None)
for a in g:
    for b in g:
        x = np.array([a, b, Omega - a - b])
        if np.all(x > d):
            v = np.sum(lam * np.log(x - d))
            if v > best[0]: best = (v, x)
print("surplus Phi =", Phi, " closed form x =", nash(lam, d, Omega), " grid argmax =", best[1])
assert np.allclose(best[1], nash(lam, d, Omega), atol=0.011)
for lam in [np.array([.6, .2, .2]), np.ones(3) / 3]:
    x = nash(lam, d, Omega); print("lam", lam.round(3), "x =", x, "U =", U(x, k).round(4),
                                  "W =", round(U(x, k).sum(), 4), "L vs FB =", round(WFB - U(x, k).sum(), 4))
    assert np.all(k * x < 1)
print("W at disagreement =", U(d, k).sum(), "   lambda* giving FB =", (1 - d) / Phi)

hdr("Allocation rule with equal disagreement payoffs and exact loss")
for Npi in [np.array([.6, .2, .2])]:
    dbar = 0.3; Phi2 = Omega - 3 * dbar; x = Omega / 3 + (Npi - 1 / 3) * Phi2
    L = 3 * U(1, k) - U(x, k).sum(); Vp = np.mean((Npi - 1 / 3) ** 2)
    print("Phi =", Phi2, "x =", x, "Var(pi) =", round(Vp, 5), "L =", L, "(kN/2)Phi^2 Var(pi) =", k * 3 / 2 * Phi2 ** 2 * Vp)
    assert np.isclose(L, k * 3 / 2 * Phi2 ** 2 * Vp)

hdr("Leontief example")
A = np.array([[.10, .05, .10], [.20, .05, .05], [.05, .20, .10]])
L = np.linalg.inv(np.eye(3) - A.T)
print("row sums of A (input per unit output):", A.sum(1), " rho(A) =", max(abs(np.linalg.eigvals(A))).round(4))
print("L = (I-A')^-1 =\n", L.round(3), "\nmultipliers (col sums) =", L.sum(0).round(3))
dd = np.array([0, 0, 100.]); r1 = A.T @ dd; r2 = A.T @ r1
print("cascade round1 =", r1, " round2 =", r2, " total y = L d =", (L @ dd).round(2))
assert np.allclose(np.linalg.inv(np.eye(3) - A.T), sum(np.linalg.matrix_power(A.T, j) for j in range(200)))
m = L.max()   # one colour scale for both panels: 0 .. max(L)
for name, M in [("A", A), ("L", L)]:
    cells = ["%d/%d/%d/%.2f/%s" % (i + 1, j + 1, round(100 * M[i, j] / m), M[i, j],
             "white" if 100 * M[i, j] / m >= 50 else "black") for i in range(3) for j in range(3)]
    print("heatmap %s cells (i/j/colour/value/text):" % name, ",".join(cells))

hdr("W^CE properties")
for n in [3, 8]:
    Wce = np.ones((n, n)) / n
    print("N", n, "eig", np.sort(abs(np.linalg.eigvals(Wce)))[::-1].round(6), "Gamma", round(gap(Wce), 6),
          "pi", perron(Wce).round(4))
    assert np.allclose(Wce @ Wce, Wce)

hdr("Stylized dual-economy heatmap (70 / 55 / 5, blocks 3+5), rows normalized")
Nb = 8; K = np.array([[70 if (i < 3 and j < 3) else 55 if (i >= 3 and j >= 3) else 5 for j in range(Nb)] for i in range(Nb)], float)
Wd = K / K.sum(1, keepdims=True); pd = perron(Wd)
print("pi =", pd.round(4), " block-1 mass", pd[:3].sum().round(3), " Gamma =", round(gap(Wd), 3),
      " N_eff =", round(1 / np.sum(pd ** 2), 2))

hdr("Conventions section: Gamma and pi independent")
for gam in [1.0, 0.3]:
    pi = np.array([.91] + [.01] * 9); Wg = gam * np.outer(np.ones(10), pi) + (1 - gam) * np.eye(10)
    print("gamma", gam, "pi(W) max", perron(Wg).max().round(3), "Gamma", round(gap(Wg), 4),
          "N_eff", round(1 / np.sum(pi ** 2), 3))
for n in [4, 20]:
    print("max Var(pi) N=%d: (N-1)/N^2 =" % n, (n - 1) / n ** 2,
          " check e1:", np.mean((np.eye(n)[0] - 1 / n) ** 2))

hdr("Exercise answers")
# 1.2 CD 1/3,2/3
Om2 = np.array([[4, 1], [1, 3], [3, 2.]]); p2 = (8 * 3 - 8) / 6; w2 = Om2 @ [1, p2]
X2 = np.column_stack([w2 / 3, 2 * w2 / (3 * p2)]); print("1.2 p2* =", p2, "X* =\n", X2.round(3), X2.sum(0))
assert np.allclose(X2.sum(0), Om2.sum(0))
# 1.3 kappa = 0.1
for om5 in [np.ones(5), np.array([3, 1, .5, .25, .25]), np.array([5, 0, 0, 0, 0.])]:
    print("1.3 omega", om5, "L =", round(0.1 / 2 * (np.sum(om5 ** 2) - 25 / 5), 5), "max kappa*x", 0.1 * om5.max())
# 1.4 TU
d2 = np.array([.6, .2]); 
for lam in [np.array([.5, .5]), np.array([.8, .2])]:
    x = nash(lam, d2, 2.0); print("1.4 lam", lam, "x =", x, "W =", U(x, .5).sum().round(4))
print("1.4 FB W =", 2 * U(1, .5), " lambda_1 giving FB =", (1 - .6) / 1.2)
# 1.5
A5 = np.array([[.2, .1], [.3, .1]]); S = np.eye(2); P = np.eye(2)
for j in range(1, 5):
    P = P @ A5.T; S = S + P; print("1.5 / App 1C partial sum k=%d" % j, S.flatten().round(4))
L5 = np.linalg.inv(np.eye(2) - A5.T); print("exact", L5.flatten().round(4), "rel err (1,2) at k=4",
      round(abs(S[0, 1] - L5[0, 1]) / L5[0, 1], 4))
print("eig A5", np.linalg.eigvals(A5).round(4), "multipliers", L5.sum(0).round(4))
# 1.6(d)
W6 = np.array([[1, 1, 1], [2, .5, .5], [0, 1.5, 1.5]]) / 3; print("1.6d rows", W6.sum(1), "cols", W6.sum(0),
      "Wx", W6 @ [3, 6, 9], "Wce x", np.full(3, 6.))
# 1.9
A9 = np.array([[.05, .3, .3], [.01, .05, .05], [.01, .05, .05]]); L9 = np.linalg.inv(np.eye(3) - A9.T)
print("1.9 L =\n", L9.round(4), "\nmultipliers", L9.sum(0).round(4), " response to d=e2", (L9 @ [0, 1, 0]).round(4))

hdr("Proof checks (randomized, seed 1): Props welfare_loss, nash_division, allocation_rule, pi_uniform, W_CE; Thm leontief")
rng = np.random.default_rng(1)
for _ in range(500):
    n = rng.integers(2, 8); kk = rng.random() * 0.5 + 0.01
    om_r = rng.random(n) * 3; Om_r = om_r.sum()
    Lr = n * U(Om_r / n, kk) - U(om_r, kk).sum()                      # welfare_loss
    assert np.isclose(Lr, kk * n / 2 * np.var(om_r)) and Lr >= -1e-12
    lam = rng.random(n) + .01; lam /= lam.sum(); d_r = rng.random(n); Om_r = d_r.sum() + rng.random() * 3 + .01
    xr = nash(lam, d_r, Om_r); obj = lambda x: np.sum(lam * np.log(x - d_r))
    assert np.isclose(xr.sum(), Om_r) and np.all(xr > d_r)              # nash_division: feasible ...
    for _ in range(20):                                                 # ... and beats feasible perturbations
        e = rng.normal(size=n); e -= e.mean(); e *= 0.5 * (xr - d_r).min() / np.abs(e).max()
        assert obj(xr + e) < obj(xr) + 1e-12
    W_r = rng.random((n, n)); W_r /= W_r.sum(1, keepdims=True); pr = perron(W_r)   # allocation_rule
    dbar = rng.random(); Phi_r = rng.random() * 3 + .01; Om_r = n * dbar + Phi_r
    x_rule = Om_r / n + (pr - 1 / n) * Phi_r
    assert np.allclose(x_rule, nash(pr, np.full(n, dbar), Om_r))
    assert np.isclose(n * U(Om_r / n, kk) - U(x_rule, kk).sum(), kk * n / 2 * Phi_r ** 2 * np.mean((pr - 1 / n) ** 2))
    D_r = rng.random((n, n)) + .01                                      # pi_uniform: Sinkhorn -> doubly stochastic
    for _ in range(500): D_r /= D_r.sum(1, keepdims=True); D_r /= D_r.sum(0, keepdims=True)
    D_r /= D_r.sum(1, keepdims=True)
    assert np.allclose(perron(D_r), 1 / n, atol=1e-6)
    assert np.allclose(perron(W_r), 1 / n) == np.allclose(W_r.sum(0), 1)
    A_r = rng.random((n, n)); A_r *= rng.random() * 0.95 / max(abs(np.linalg.eigvals(A_r)))   # leontief
    L_r = np.linalg.inv(np.eye(n) - A_r.T); S_r = sum(np.linalg.matrix_power(A_r.T, k) for k in range(2000))
    assert np.allclose(L_r, S_r, atol=1e-8) and np.all(L_r >= -1e-12) and np.all(L_r - np.eye(n) >= -1e-12)
    Wce_r = np.ones((n, n)) / n; ev = np.sort(abs(np.linalg.eigvals(Wce_r)))[::-1]     # W_CE
    assert np.isclose(ev[0], 1) and np.allclose(ev[1:], 0) and np.allclose(perron(Wce_r), 1 / n)
    s_r = rng.random(n); s_r /= s_r.sum(); assert np.allclose(perron(np.outer(np.ones(n), s_r)), s_r)
print("all randomized proof checks passed")
print("\nall assertions passed")
