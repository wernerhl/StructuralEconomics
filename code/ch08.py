"""Chapter 8: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch08.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def stationary(W):
    """Left eigenvector for eigenvalue 1, normalized (unique if W has one closed class)."""
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmin(np.abs(vals - 1))
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def gap(W):
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


def neff(p):
    return 1 / np.sum(p ** 2)


def halflife(W):
    l2 = np.sort(np.abs(np.linalg.eigvals(W)))[::-1][1]
    if l2 < 1e-12:
        return 0.0
    return np.log(2) / (-np.log(l2)) if l2 < 1 else np.inf


print("=" * 70)
print("Example ch08:ex:three_WI: b0 = (0.9, 0.8, 0.2, 0.1)")
b0 = np.array([.9, .8, .2, .1])
Wc = np.ones((4, 4)) / 4
Ws = np.array([[.7, .1, .1, .1], [.6, .2, .1, .1], [.6, .1, .2, .1], [.6, .1, .1, .2]])
Wb = np.array([[.6, .4, 0, 0], [.4, .6, 0, 0], [0, 0, .6, .4], [0, 0, .4, .6]])
for name, W in (("Complete", Wc), ("Star", Ws)):
    p = stationary(W)
    print(f"{name}: pi={np.round(p, 3)} b*={p @ b0:.3f} Gamma={gap(W):.2f} "
          f"half-life={halflife(W):.2f} N_eff={neff(p):.2f} eigs={np.round(np.linalg.eigvals(W).real, 3)}")
assert np.allclose(stationary(Ws), [2 / 3, 1 / 9, 1 / 9, 1 / 9])
bt = np.linalg.matrix_power(Wb, 200) @ b0
print("Block: limit beliefs", np.round(bt, 3), " Gamma =", round(gap(Wb), 3), "(reducible)")
print("heatmap cells (100*W): complete", np.round(100 * Wc).astype(int)[0].tolist(),
      " star", np.round(100 * Ws).astype(int).tolist(), " block", np.round(100 * Wb).astype(int).tolist())

print("=" * 70)
print("Prop ch08:thm:wisdom: Var(e) >= sigma_min^2 / N_eff")
rng = np.random.default_rng(0)
for _ in range(1000):
    p = rng.dirichlet(np.ones(6)); s2 = rng.random(6) + 0.1
    assert p @ (p * s2) >= s2.min() / neff(p) - 1e-12
# with heterogeneous precision, equal weights are not variance-minimizing
s2 = np.array([1., 1, 4, 4])
eq = np.full(4, .25); opt = (1 / s2) / np.sum(1 / s2)
print("sigma^2 =", s2, " Var(equal) =", eq @ (eq * s2), " Var(1/sigma^2 weights) =", round(opt @ (opt * s2), 4))

print("=" * 70)
print("N_eff illustrations")
N = 1000
p = np.r_[np.full(50, .9 / 50), np.full(950, .1 / 950)]
print("5% hold 90%: N_eff/N =", round(neff(p) / N, 4), " (", round(N / neff(p), 1), "times fewer voices)")
p = np.r_[np.full(3, .2), np.full(200, .4 / 200)]
print("Washington Consensus illustration (3 x 0.2, 200 share 0.4): N_eff =", round(neff(p), 2))
for p1 in (0.5, 0.9):
    p = np.r_[p1, np.full(99, (1 - p1) / 99)]
    print(f"Ex 8.2: pi_1={p1}: N_eff = {neff(p):.2f}")
# Ex 8.2(d): N_eff < 10 when pi_1 above
grid = np.linspace(0, 1, 100001)
ne = 1 / (grid ** 2 + (1 - grid) ** 2 / 99)
print("Ex 8.2(d): N_eff drops below 10 at pi_1 =", round(grid[np.argmax((ne < 10) & (grid > 0.01))], 4))

print("=" * 70)
print("Prop ch08:thm:manipulation: first-order bound |db*| <= t ||E|| ||Z|| max|b0 - b*|")
for _ in range(500):
    W = rng.random((5, 5)); W /= W.sum(1, keepdims=True)
    E = rng.standard_normal((5, 5)); E -= E.mean(1, keepdims=True)
    b = rng.random(5)
    p = stationary(W); Z = np.linalg.inv(np.eye(5) - W + np.outer(np.ones(5), p))
    db = p @ E @ Z @ b
    t = 1e-6
    fd = (stationary(W + t * E) @ b - stationary(W - t * E) @ b) / (2 * t)
    assert np.isclose(db, fd, atol=1e-6)
    assert abs(db) <= np.abs(E).sum(1).max() * np.abs(Z).sum(1).max() * np.abs(b - p @ b).max() + 1e-12
print("ok (500 draws); Z 1 = 1 holds:", np.allclose(Z @ np.ones(5), 1))
# two-block example: cost of moving the consensus falls as coupling weakens
b = np.array([.9, .9, .1, .1])
for eps in (0.2, 0.05, 0.01):
    W = np.full((4, 4), eps / 2); W[:2, :2] = W[2:, 2:] = (1 - eps) / 2
    # raise block-2 rows' weight on block 1 from eps to eps+de: block-1 mass = (eps+de)/(2eps+de)
    de = 0.01
    W2 = W.copy(); W2[2:, :2] += de / 2; W2[2:, 2:] -= de / 2
    print(f"  eps={eps}: Gamma={gap(W):.3f} b*={stationary(W) @ b:.3f} -> {stationary(W2) @ b:.3f} "
          f"after moving 0.01 of weight per row (||E||_inf = 0.02)")

print("=" * 70)
print("Example ch08:ex:platform")
bC = np.array([.9, .3, .5]); bU = np.array([.8, .5, .2])
b = np.r_[bC, bU]
Wo = np.zeros((6, 6))
Wo[:3, :3] = 0.1; np.fill_diagonal(Wo[:3, :3], 0.8)
Wo[3:, :3] = 1 / 3
p = stationary(Wo)
lim = np.linalg.matrix_power(Wo, 500) @ b
print("organic: pi =", np.round(p, 3), " consensus =", np.round(lim, 3), " Gamma =", round(gap(Wo), 3))
Wa = np.zeros((6, 6)); Wa[:3, :3] = np.eye(3)
Wa[3, [0, 2]] = [.875, .125]; Wa[4, 2] = 1.0; Wa[5, [1, 2]] = [.875, .125]
lim = np.linalg.matrix_power(Wa, 500) @ b
print("algorithmic: limit beliefs =", np.round(lim, 3), " Gamma =", round(gap(Wa), 3))

print("=" * 70)
print("Exercise 8.1")
W = np.array([[.5, .3, .2], [.1, .6, .3], [.3, .2, .5]]); b = np.array([1., 0, .5])
for t in range(1, 4):
    b = W @ b; print(f"b({t}) =", np.round(b, 4))
p = stationary(W); l = np.linalg.eigvals(W)
print("pi =", np.round(p, 4), " b* =", round(p @ np.array([1., 0, .5]), 4), " eigs =", np.round(l, 4),
      " |l2| =", round(np.sort(np.abs(l))[-2], 4), " half-life =", round(halflife(W), 3))

print("=" * 70)
print("Exercise 8.4: W_ij ~ exp(-beta |b_i - b_j|), 4 at 0.8, 4 at 0.2")
b0 = np.r_[np.full(4, .8), np.full(4, .2)]
for beta in (0, 1, 5, 20):
    W = np.exp(-beta * np.abs(b0[:, None] - b0[None, :])); W /= W.sum(1, keepdims=True)
    p = stationary(W)
    assert np.isclose(gap(W), 1 - np.tanh(0.3 * beta))
    print(f"beta={beta:2d}: b*={p @ b0:.3f} Gamma={gap(W):.6f} half-life={halflife(W):.3g} rounds")

print("=" * 70)
print("Aligned layers (case study): common stationary pi is preserved by any mixture")
p = rng.dirichlet(np.ones(4))
Ls = []
for _ in range(3):
    g = rng.random() * 0.8 + 0.1
    Ls.append(g * np.outer(np.ones(4), p) + (1 - g) * np.eye(4))
Wmix = 0.2 * Ls[0] + 0.5 * Ls[1] + 0.3 * Ls[2]
assert np.allclose(stationary(Wmix), p)
print("ok: pi(composite) = common pi =", np.round(p, 3))

print("=" * 70)
print("Proof checks: Thm ch08:thm:degroot, Prop ch08:prop:rate, Prop ch08:prop:misinfo")
for _ in range(300):
    n = rng.integers(3, 8)
    W = rng.random((n, n)) * (rng.random((n, n)) < 0.5) + 0.3 * np.roll(np.eye(n), 1, 1)
    W[0, 0] += 0.2
    W /= W.sum(1, keepdims=True)                 # strongly connected and aperiodic
    p = stationary(W); P = np.outer(np.ones(n), p)
    M = W - P
    l2 = np.sort(np.abs(np.linalg.eigvals(W)))[-2]
    assert np.isclose(max(abs(np.linalg.eigvals(M))), l2, atol=1e-8)
    for t in (1, 2, 5, 17):
        assert np.allclose(np.linalg.matrix_power(M, t), np.linalg.matrix_power(W, t) - P)
    # misinformation: E[b*] - omega = pi_h delta, by linearity of the consensus
    h, omega, delta = rng.integers(n), rng.random(), rng.standard_normal()
    eps = rng.standard_normal((2000, n)) * 0.1
    b0 = omega + eps; b0[:, h] += delta
    bstar = b0 @ p
    assert np.allclose(bstar, omega + p[h] * delta + eps @ p)
    assert np.allclose(np.linalg.matrix_power(W, 4000) @ b0[0], bstar[0], atol=1e-8)
print("ok: (W - 1pi')^t = W^t - 1pi', rho(W - 1pi') = |lambda_2|, b* = omega + pi_h delta + pi'eps")
print("all assertions passed")
