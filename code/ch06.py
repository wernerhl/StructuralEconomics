"""Chapter 6: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch06.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


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
print("all assertions passed")
