"""Chapter 25: dense subgraph vs bridge, Cheeger bound, truncation numbers.
Run from the repo root: python3 code/ch25.py  (numpy only)."""
import itertools
import numpy as np


def rownorm(K):
    return K / K.sum(1, keepdims=True)


def gap(W):
    ev = np.linalg.eigvals(W)
    ev = ev[np.argsort(-np.abs(ev))]
    return 1 - abs(ev[1])


def goldberg_obj(phi, S):
    S = list(S)
    sub = phi[np.ix_(S, S)]
    return np.triu(sub, 1).sum() / len(S)


def charikar_peel(phi):
    """Greedy peeling (Charikar 2000): 1/2-approximation to max e(D)/|D|."""
    alive = list(range(len(phi)))
    best, bestS = -1, None
    while alive:
        val = goldberg_obj(phi, alive)
        if val > best:
            best, bestS = val, list(alive)
        deg = phi[np.ix_(alive, alive)].sum(1)
        alive.pop(int(np.argmin(deg)))
    return bestS, best


def conductance(phi, S):
    S = np.array(S)
    T = np.setdiff1d(np.arange(len(phi)), S)
    d = phi.sum(1)
    return phi[np.ix_(S, T)].sum() / min(d[S].sum(), d[T].sum())


print("== Case B: D is a dense MODULE (6 periphery products well mixed, 4 core) ==")
nP, nD = 6, 4
n = nP + nD
phiB = np.full((n, n), 0.05)
phiB[:nP, :nP] = 0.3
phiB[nP:, nP:] = 0.9
np.fill_diagonal(phiB, 0)
dens = lambda M, S: M[np.ix_(S, S)].sum() / (len(S) * (len(S) - 1))
D = list(range(nP, n))
print(f"density of D = {dens(phiB, D):.3f}, of all = {dens(phiB, list(range(n))):.3f}, ratio = {dens(phiB, D)/dens(phiB, list(range(n))):.2f}")
S, val = charikar_peel(phiB)
# exact Goldberg optimum by brute force (n = 10)
best = max((goldberg_obj(phiB, c), c) for k in range(1, n + 1) for c in itertools.combinations(range(n), k))
print(f"peeling finds {S} (obj {val:.3f}); brute-force optimum {list(best[1])} (obj {best[0]:.3f})")
gB, gB_noD = gap(rownorm(phiB)), gap(rownorm(phiB[:nP, :nP]))
print(f"Gamma with D = {gB:.3f}; without D = {gB_noD:.3f}  (removing the dense module RAISES Gamma)")
assert gB_noD > gB

print("\n== Case A: D is a BRIDGE between three periphery blocks ==")
blocks, b = 3, 3
m = blocks * b + nD
phiA = np.full((m, m), 0.001)
for k in range(blocks):
    phiA[k * b:(k + 1) * b, k * b:(k + 1) * b] = 0.5
phiA[blocks * b:, :] = 0.3
phiA[:, blocks * b:] = 0.3
phiA[blocks * b:, blocks * b:] = 0.9
np.fill_diagonal(phiA, 0)
gA, gA_noD = gap(rownorm(phiA)), gap(rownorm(phiA[:blocks * b, :blocks * b]))
print(f"Gamma with D = {gA:.3f}; without D = {gA_noD:.3f}  (removing the bridge LOWERS Gamma)")
assert gA_noD < gA
# Cheeger upper bound after removal: Gamma <= 2 * conductance of block-1 cut
rest = phiA[:blocks * b, :blocks * b]
phi_cut = conductance(rest, list(range(b)))
print(f"conductance of block-1 cut without D = {phi_cut:.4f}; bound 2*Phi = {2*phi_cut:.4f} >= Gamma {gA_noD:.4f}")
assert gA_noD <= 2 * phi_cut + 1e-12
# with D present the same block-1 cut has conductance
print(f"conductance of block-1 cut with D present = {conductance(phiA, list(range(b))):.3f}")

# random check of Gamma <= 2 Phi(S) for reversible W = D^-1 phi
rng = np.random.default_rng(25)
for _ in range(300):
    k = 6
    P = rng.uniform(0, 1, (k, k)); P = (P + P.T) / 2; np.fill_diagonal(P, rng.uniform(0, 1, k))
    W = rownorm(P)
    ev = np.sort(np.real(np.linalg.eigvals(W)))[::-1]
    for r_ in range(1, k):
        for Sset in itertools.combinations(range(k), r_):
            assert 1 - ev[1] <= 2 * conductance(P, list(Sset)) + 1e-9
            assert gap(W) <= 1 - ev[1] + 1e-9
print("checked Gamma <= 1 - lambda_2 <= 2*Phi(S) for 300 random reversible W, all cuts")

# Proof check for Proposition ch25:prop:bridge_bound: Rayleigh quotient of the test vector
for _ in range(300):
    k = int(rng.integers(3, 9))
    P = rng.uniform(0, 1, (k, k)); P = (P + P.T) / 2; np.fill_diagonal(P, 0.0)
    d = P.sum(1)
    Sset = rng.choice(k, size=int(rng.integers(1, k)), replace=False)
    T = np.setdiff1d(np.arange(k), Sset)
    vS, vT = d[Sset].sum(), d[T].sum()
    f = np.zeros(k); f[Sset] = 1 / vS; f[T] = -1 / vT
    assert abs(f @ d) < 1e-12
    num = 0.5 * sum(P[p_, q_] * (f[p_] - f[q_]) ** 2 for p_ in range(k) for q_ in range(k))
    den = f @ (d * f)
    kap = 1 / vS + 1 / vT
    cut = P[np.ix_(Sset, T)].sum()
    assert abs(num - cut * kap ** 2) < 1e-10 and abs(den - kap) < 1e-10
    lam2 = np.sort(np.real(np.linalg.eigvals(P / d[:, None])))[::-1][1]
    assert 1 - lam2 <= num / den + 1e-10 <= 2 * conductance(P, list(Sset)) + 2e-10
print("checked: test-vector Rayleigh quotient = cut*(1/vol S + 1/vol S^c) >= 1 - lambda_2, <= 2 Phi(S)")

print("\n== Degenerate definition: average including diagonal phi_pp = 1 ==")
phid = phiB.copy(); np.fill_diagonal(phid, 1.0)
single = phid[np.ix_([0], [0])].mean()
print(f"|D| = 1 gives density {single:.2f} = maximum: the old definition picks a single product")

print("\n== Truncation index (illustrative) ==")
for c, k in (("Bolivia", 2), ("South Korea", 82)):
    print(f"{c}: tau = 1 - {k}/100 = {1 - k/100:.2f}")

print("\n== Relaxation-time penalty (1-G)/G, Bolivia G = 0.13 ==")
f = lambda g: (1 - g) / g
print(f"(1-0.13)/0.13 = {f(0.13):.3f}; (1-0.35)/0.35 = {f(0.35):.3f}; difference = {f(0.13)-f(0.35):.3f}")
for g0 in (0.25, 0.30, 0.35, 0.45, 0.52):
    print(f"counterfactual Gamma(0) = {g0:.2f}: share = {(f(0.13)-f(g0))/f(0.13):.3f}")

print("\n== Bolivia's incremental path timing ==")
lo, hi = 4 * 20, 4 * 30
print(f"4 steps x 20-30 years = {lo}-{hi} years; at 7.5x velocity: {lo/7.5:.1f}-{hi/7.5:.1f} years")
print("\n== Korea path (Ch23 linear): Gamma(1973) =", round(0.15 + 0.0075 * 13, 3), " Gamma(1993) =", round(0.15 + 0.0075 * 33, 3))
