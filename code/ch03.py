"""Chapter 3: reproduces every number in the worked examples, tables, viz panel
and exercise answers.  Convention (Ch1): a_ij = input j per unit of i, L = (I-A')^-1,
W^P = D^-1 A (row-stochastic input shares).  Run: python3 code/ch03.py (numpy only)."""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
def hdr(s): print("\n== " + s)
def Linv(A, g=1.0): return np.linalg.inv(np.eye(len(A)) - g * A.T)
def WP(A): return A / A.sum(1, keepdims=True)
def perron(W):
    w, v = np.linalg.eig(W.T); p = np.real(v[:, np.argmin(abs(w - 1))]); return p / p.sum()
def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True); return 1 - e[1]
def cross_mass(A, blocks):
    tot = A.sum(); c = sum(A[i, j] for i in range(len(A)) for j in range(len(A)) if blocks[i] != blocks[j]); return c / tot

hdr("Example two_WP / propagation (rows: Modern, Traditional)")
AA = np.array([[.10, .25], [.20, .10]]); AB = np.array([[.10, .02], [.02, .10]])
for name, A in [("A (connected)", AA), ("B (disconnected)", AB)]:
    L = Linv(A); print(name, "L = (I-A')^-1 =\n", L.round(3), "\n  modern final-demand impulse -> (M, T) =", L[:, 0].round(3),
                       " rho(A) =", round(max(abs(np.linalg.eigvals(A))), 3))
    assert np.allclose(L, sum(np.linalg.matrix_power(A.T, k) for k in range(300)))

hdr("Bolivia stylized example")
Ab = np.array([[.15, .03], [.02, .20]]); Ac = np.array([[.15, .15], [.15, .20]])
for name, A in [("enclave", Ab), ("connected", Ac)]:
    L = Linv(A); col = L[:, 0]; W = WP(A); p = perron(W)
    ev = np.sort(np.linalg.eigvals(A).real)[::-1]
    print(name, "L =\n", L.round(4), "\n  gas impulse ->", col.round(4), "multiplier", col.sum().round(4),
          "share to traditional", round(col[1] / col.sum(), 4), " eig(A) =", ev.round(4), "rho-|l2| (old) =", round(ev[0] - abs(ev[1]), 4))
    print("  W^P = D^-1 A =", W.round(4).tolist(), " Gamma(W^P) =", round(gap(W), 4), " = w12+w21 =", round(W[0, 1] + W[1, 0], 4),
          " pi(W^P) =", p.round(4), " cross-block mass of A =", round(cross_mass(A, [0, 1]), 4))
    assert np.isclose(gap(W), W[0, 1] + W[1, 0])
Lb, Lc = Linv(Ab), Linv(Ac); print("ratio of spillovers connected/enclave =", round(Lc[1, 0] / Lb[1, 0], 2))
m = max(Lb.max(), Lc.max())
for name, L in [("enclave", Lb), ("connected", Lc)]:
    print("heatmap", name, ",".join("%d/%d/%d/%.2f" % (i + 1, j + 1, round(100 * L[i, j] / m), L[i, j]) for i in range(2) for j in range(2)))
A06 = np.diag([0.6, 0.05]); print("diag(0.6,0.05): rho-|l2| =", 0.55, " Gamma of W^P (identity) =", round(gap(WP(A06)), 3))

hdr("Growth-propagation bound")
W = np.array([[.5, .9], [0, .5]]); print("entry 0.9 > rho =", max(abs(np.linalg.eigvals(W))))
C = np.array([[.85, .88, .05], [0, .86, .72], [0, 0, .88]]); rho = max(abs(np.linalg.eigvals(C)))
L = np.linalg.inv(np.eye(3) - C.T); off = (L.sum() - np.trace(L)) / 6
print("3x3 nonnormal: rho =", rho.round(3), "avg off-diag L =", round(off, 2), "old 'bound' ||C||_2/(1-rho) =",
      round(np.linalg.norm(C, 2) / (1 - rho), 2), " ||C||_2 =", round(np.linalg.norm(C, 2), 3), "(new bound needs ||C||<1)")
for name, A in [("Bolivia", Ab), ("connected", Ac), ("A", AA)]:
    n2 = np.linalg.norm(A, 2); L = Linv(A); N = len(A); off = (L.sum() - np.trace(L)) / (N * (N - 1))
    print(name, "avg off-diag =", round(off, 4), "bound g||A||/(1-g||A||) =", round(n2 / (1 - n2), 4)); assert off <= n2 / (1 - n2)
rng = np.random.default_rng(0)
for _ in range(20000):   # random check of the corrected bound
    N = rng.integers(2, 6); A = rng.random((N, N)) * rng.random() * 0.6; n2 = np.linalg.norm(A, 2)
    if n2 >= 0.99: continue
    L = Linv(A); assert (L - np.diag(np.diag(L))).max() <= n2 / (1 - n2) + 1e-12
print("corrected bound holds entrywise on 20000 random draws")

hdr("Exercises")
s = lambda th: th / (1 + th); a = .3
print("3.1 sigma_M, sigma_F =", s(2), s(.5), " sbar =", (s(2) + s(.5)) / 2,
      " K*(sbar) = (sbar(1-a))^(1/(1-a)) ratio vs both theta=2:", round(((s(2) + s(.5)) / 2 / s(2)) ** (1 / (1 - a)), 4))
A2 = np.array([[.10, .20, .05], [.15, .10, .15], [.05, .10, .10]]); L2 = Linv(A2)
print("3.2 L =\n", L2.round(4), "\n  impulse to sector 2 ->", L2[:, 1].round(4), " multipliers", L2.sum(0).round(4),
      " rho(A) =", round(max(abs(np.linalg.eigvals(A2))), 4), " Gamma(W^P) =", round(gap(WP(A2)), 4))
A3 = A2.copy(); A3[2, 0] = A3[2, 1] = A3[0, 2] = A3[1, 2] = .01; L3 = Linv(A3)
print("3.3 L =\n", L3.round(4), "\n  impulse to sector 3 ->", L3[:, 2].round(4), " Gamma(W^P) =", round(gap(WP(A3)), 4),
      " decline factor (spill to 1,2):", (L2[:2, 2] / L3[:2, 2]).round(2))
A4 = np.array([[.15, .10], [.10, .20]]); L4 = Linv(A4); print("3.4 links 0.10: L =", L4.round(4).tolist(), " Gamma(W^P) =", round(gap(WP(A4)), 4))
th = .5
for t in range(1, 41):
    th = .8 * th + .2 * 2 + .1 * (.7 * .1 - .1)
    if t in (1, 5, 10, 20, 40): print("3.7 theta_F(t=%d) = %.4f gap = %.4f" % (t, th, 2 - th))
print("3.7 fixed point with lambda=0.7:", round((.4 + .1 * (.07 - .1)) / .2, 4))
print("\nall assertions passed")
