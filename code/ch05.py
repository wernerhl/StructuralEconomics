"""Chapter 5: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch05.py
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


def neff(p):
    return 1 / np.sum(p ** 2)


def rownorm(K):
    """W = D^{-1} K; an all-zero row becomes a self-loop (the agent attends only to herself)."""
    K = K.astype(float).copy()
    for i in range(len(K)):
        if K[i].sum() == 0:
            K[i, i] = 1.0
    return K / K.sum(1, keepdims=True)


print("=" * 70)
print("Example ch05:ex:vendors: U = x - x^2/2, gamma_U = 0.5")
U = lambda x: x - 0.5 * x ** 2
x = np.array([.6, .4, .8])
KU = np.array([[0, .30, .02], [.25, 0, .05], [.02, .05, 0]])
V = U(x) + 0.5 * KU @ U(x)
print("U =", U(x), " V =", np.round(V, 4), " externality share for Maria =", round((V[0] - U(x)[0]) / U(x)[0], 4))

print("=" * 70)
print("Example ch05:ex:LaPaz: raw layers, row-normalized, composite (0.2, 0.5, 0.3)")
KUl = np.array([[0, .3, 0], [.2, 0, 0], [0, 0, 0]])
KPl = np.array([[.1, .4, .1], [.3, .2, 0], [0, 0, .3]])
WI = np.array([[.3, .5, .2], [.4, .3, .3], [.2, .1, .7]])
WU, WP = rownorm(KUl), rownorm(KPl)
W = 0.2 * WU + 0.5 * WP + 0.3 * WI
assert np.allclose(W.sum(1), 1)
p = perron(W)
print("W^U =\n", WU, "\nW^P =\n", np.round(WP, 4), "\ncomposite W =\n", np.round(W, 4))
print("pi =", np.round(p, 3), " N_eff =", round(neff(p), 2), " Gamma =", round(gap(W), 3))
assert np.allclose(p, 0.2 * p @ WU + 0.5 * p @ WP + 0.3 * p @ WI)   # Prop ch05:prop:composite_power
print("pi(W^I) =", np.round(perron(WI), 3), " col sums W^I =", WI.sum(0))
print("heatmaps (100*W): U", np.round(100 * WU).astype(int).tolist(), " P", np.round(100 * WP).astype(int).tolist(),
      " I", np.round(100 * WI).astype(int).tolist(), " W", np.round(100 * W).astype(int).tolist())
# raw (unnormalized) composite, for comparison with the earlier draft
Wraw = 0.2 * KUl + 0.5 * KPl + 0.3 * WI
print("raw composite W21 =", round(Wraw[1, 0], 2), " pi(raw) =", np.round(perron(Wraw), 3),
      " rho =", round(max(abs(np.linalg.eigvals(Wraw))), 3))
# rescaling a raw layer changes pi: non-invariance of unnormalized composites
print("raw composite with K^U doubled: pi =", np.round(perron(0.2 * 2 * KUl + 0.5 * KPl + 0.3 * WI), 3))
# Exercise 5.7(a): vary alpha_P with alpha_U/alpha_I = 2/3
for aP in (0, .25, .5, .75, .95):
    aU, aI = 0.4 * (1 - aP), 0.6 * (1 - aP)
    print(f"  Ex 5.7 alpha_P={aP}: pi_3 = {perron(aU * WU + aP * WP + aI * WI)[2]:.3f}")

print("=" * 70)
print("Exercise 5.3 / Appendix 5B")
A = np.array([[.5, .3, .2], [.1, .6, .3], [.2, .2, .6]])
v = np.ones(3) / 3; prev = v
for t in range(1, 6):
    v = v @ A
    if t in (1, 2, 5):
        print(f"t={t}: pi={np.round(v, 3)} max change={np.abs(v - prev).max():.3f}")
    prev = v
pA = perron(A)
lam = np.linalg.eigvals(A)
print("pi =", np.round(pA, 3), " eigenvalues", np.round(lam, 3), " |lambda2| =", round(np.sort(np.abs(lam))[-2], 3),
      " Gamma =", round(gap(A), 3))
b0 = np.array([0, .5, 1])
print("Ex 5.3(d): b* = pi'b0 =", round(pA @ b0, 4), "; b(50) =", np.round(np.linalg.matrix_power(A, 50) @ b0, 4))
print("col sums", A.sum(0), " W13, W31 =", A[0, 2], A[2, 0], " W23, W32 =", A[1, 2], A[2, 1])
# per-step changes for the table
v = np.ones(3) / 3; out = []
for t in range(1, 6):
    nv = v @ A; out.append(np.abs(nv - v).max()); v = nv
print("max change by step:", np.round(out, 3))

print("=" * 70)
print("Case study: colonial star W = 1 e_1' (+ small self-weights)")
N = 6
Wc = 0.95 * np.outer(np.ones(N), np.eye(N)[0]) + 0.05 * np.eye(N)
pc = perron(Wc)
print("pi =", np.round(pc, 3), " N_eff =", round(neff(pc), 2), " Gamma =", round(gap(Wc), 3))
print("all assertions passed")
