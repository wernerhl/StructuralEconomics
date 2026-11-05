"""Chapter 22: utility-network measurement.  Reproduces every number in
chapters/chapter22.tex.  numpy only.  Run from the repo root:  python3 code/ch22.py
"""
import numpy as np

np.set_printoptions(precision=3, suppress=True, linewidth=150)
rng = np.random.default_rng(22)


def perron(W):
    w, v = np.linalg.eig(W.T)
    q = np.real(v[:, np.argmin(abs(w - 1))])
    return q / q.sum()


def gap(W):
    m = np.sort(abs(np.linalg.eigvals(W)))[::-1]
    return 1 - m[1]


def neff(q):
    return 1 / (q ** 2).sum()


print("=== Example: transfer shares among four domestic groups (illustrative) ===")
g = ['Urban NI', 'Urban Indig', 'Rural NI', 'Rural Indig']
raw = np.array([[.45, .08, .05, .02],      # domestic part of the earlier draft's rows
                [.10, .40, .05, .25],
                [.08, .03, .50, .10],
                [.03, .15, .07, .55]])
T = raw / raw.sum(1, keepdims=True)
print(T.round(3))
G, piT = gap(T), perron(T)
print(f"Gamma = {G:.3f}; pi = {piT.round(3)}; Neff = {neff(piT):.2f} of 4 (Neff/n = {neff(piT)/4:.2f})")
diag = np.diag(T)
print("own-group shares:", diag.round(2))
eth_cross = np.mean([T[0, 1], T[0, 3], T[2, 1], T[2, 3], T[1, 0], T[1, 2], T[3, 0], T[3, 2]])
eth_within = np.mean([T[0, 0], T[0, 2], T[2, 0], T[2, 2], T[1, 1], T[1, 3], T[3, 1], T[3, 3]])
geo_cross = np.mean([T[0, 2], T[0, 3], T[1, 2], T[1, 3], T[2, 0], T[2, 1], T[3, 0], T[3, 1]])
geo_within = np.mean([T[0, 0], T[0, 1], T[1, 0], T[1, 1], T[2, 2], T[2, 3], T[3, 2], T[3, 3]])
print(f"lambda_ethnicity^U = {eth_cross/eth_within:.2f}; lambda_geography^U = {geo_cross/geo_within:.2f}")
print("heat cells (round 100*share):", np.rint(100 * T).astype(int).tolist())

print("\n=== Dictator shares -> altruism weights, W = s/(1-s) (log utility) ===")
s = {'same comm, same eth': .38, 'same comm, diff eth': .25, 'diff comm, same eth': .22,
     'diff comm, diff eth': .12, 'stranger': .18}
Wd = {k: v / (1 - v) for k, v in s.items()}
for k, v in Wd.items():
    print(f"  {k:22s} s={s[k]:.2f}  W={v:.3f}")
lc = (Wd['diff comm, same eth'] / Wd['same comm, same eth'], Wd['diff comm, diff eth'] / Wd['same comm, diff eth'])
le = (Wd['same comm, diff eth'] / Wd['same comm, same eth'], Wd['diff comm, diff eth'] / Wd['diff comm, same eth'])
print(f"lambda_community in [{min(lc):.2f}, {max(lc):.2f}]; lambda_ethnicity in [{min(le):.2f}, {max(le):.2f}]")
# check the log-altruism FOC numerically
E = 1.0
for Wv in [0.613, 0.333]:
    xs = np.linspace(1e-4, E - 1e-4, 200001)
    x = xs[np.argmax(np.log(E - xs) + Wv * np.log(xs))]
    assert abs(x / E - Wv / (1 + Wv)) < 1e-3

print("\n=== Appendix 22A: altruism FOC, CRRA ===")
sigma, ci, cj = 2.0, 3.0, 1.0
W_foc = (cj / ci) ** sigma
print(f"c_i={ci}, c_j={cj}, sigma={sigma}: W_ij = (c_j/c_i)^sigma = {W_foc:.3f} "
      f"(inverted formula would give {(ci/cj)**sigma:.1f})")
up = lambda c: c ** (-sigma)
assert np.isclose(up(ci), W_foc * up(cj))

print("\n=== Prop. alignment (rank-one layers): pi = sum alpha_k pi^k, Var <= max Var_k ===")
for trial in range(2000):
    N = rng.integers(3, 8); K = rng.integers(2, 4)
    P = rng.dirichlet(np.ones(N) * rng.uniform(0.3, 3), size=K)
    a = rng.dirichlet(np.ones(K))
    Wc = sum(ak * np.outer(np.ones(N), pk) for ak, pk in zip(a, P))
    pc = perron(Wc)
    assert np.allclose(pc, a @ P)
    assert pc.var() <= max(pk.var() for pk in P) + 1e-15
print("verified on 2000 random draws")
# non-rank-one counterexample to the old 'Var(composite) >= max Var' claim
pP = np.array([.4, .3, .2, .1]); pI = np.array([.3, .27, .23, .2]); pU = np.array([.28, .26, .24, .22])
def withpi(p):
    n = len(p); W = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                W[i, j] = 0.3 * min(1, p[j] / p[i])
        W[i, i] = 1 - W[i].sum()
    return W
Ws = [withpi(x) for x in (pP, pI, pU)]
for Wl, x in zip(Ws, (pP, pI, pU)):
    assert np.allclose(perron(Wl), x)
pcomp = perron(sum(Ws) / 3)
print(f"aligned non-rank-one layers: Var(composite) = {pcomp.var():.4f} < max Var = {max(x.var() for x in (pP,pI,pU)):.4f}")

print("\n=== Example: alignment on four groups (illustrative per-capita positions) ===")
p = np.array([0.35, 0.25, 0.15, 0.25])      # Urban NI, Urban Indig, Rural NI, Rural Indig
# per-capita positions r = pi/p; last group as in the text; others chosen so pi sums to 1
rP = np.array([1.50, 0.80, 1.00, 0.46]); rI = np.array([1.40, 0.90, 1.20, 0.35]); rU = np.array([1.35, 1.00, 1.15, 0.30])
for name, r in [('P', rP), ('I', rI), ('U', rU)]:
    rr = r.copy()
    rr[2] = (1 - (p * r)[[0, 1, 3]].sum()) / p[2]   # adjust rural NI to normalise
    if name == 'P': rP = rr
    if name == 'I': rI = rr
    if name == 'U': rU = rr
for name, r in [('P', rP), ('I', rI), ('U', rU)]:
    assert np.isclose((p * r).sum(), 1)
    print(f"r^{name} = {r.round(3)}")
alpha = np.ones(3) / 3
rc = alpha @ np.vstack([rP, rI, rU])
print(f"composite per-capita position of rural indigenous (rank-one layers, equal weights): {rc[3]:.3f}")
corr = lambda a, b: np.corrcoef(a, b)[0, 1]
print(f"alignment: corr(P,U) = {corr(rP*p, rU*p):.2f}, corr(P,I) = {corr(rP*p, rI*p):.2f}, corr(I,U) = {corr(rI*p, rU*p):.2f}")
# anti-aligned counterfactual: reverse the utility layer's ordering for the two indigenous groups
rU2 = rU.copy(); rU2[3], rU2[1] = rU[1] * p[1] / p[3], rU[3] * p[3] / p[1]
assert np.isclose((p * rU2).sum(), 1)
rc2 = alpha @ np.vstack([rP, rI, rU2])
print(f"if W^U favoured rural indigenous (r^U = {rU2[3]:.2f}): composite position {rc2[3]:.3f}")

print("\n=== Proof pass: Prop. alignment, (a) shared Perron vector; bound fails off rank one ===")
for trial in range(300):
    N = rng.integers(3, 7); q = rng.dirichlet(np.ones(N))
    def with_pi(qv):                      # Metropolis chain with stationary vector qv
        P = rng.random((len(qv), len(qv))); P = (P + P.T) / 2; P /= P.sum(1).max() * 1.01
        Wm = P * np.minimum(1, qv[None, :] / qv[:, None]); np.fill_diagonal(Wm, 0)
        np.fill_diagonal(Wm, 1 - Wm.sum(1)); return Wm
    Ls = [with_pi(q) for _ in range(3)]; a = rng.dirichlet(np.ones(3))
    assert np.allclose(perron(sum(ak * L for ak, L in zip(a, Ls))), q)
print("(a) verified on 300 random draws (layers sharing a stationary vector)")
W1 = np.array([[.8, .1, .1], [.9, .1, 0], [0, .1, .9]])
W2 = np.array([[.7, .3, 0], [.1, .5, .4], [.7, .2, .1]])
for Wl in (W1, W2):
    assert np.linalg.matrix_power(Wl + np.eye(3), 3).min() > 0          # irreducible
p1, p2, pc = perron(W1), perron(W2), perron((W1 + W2) / 2)
ratio = pc.var() / max(p1.var(), p2.var())
print(f"non-rank-one counterexample: pi1 = {p1.round(3)}, pi2 = {p2.round(3)}, composite = {pc.round(3)}; "
      f"Var(composite)/max Var = {ratio:.2f}")
assert ratio > 1.5
# two rank-one layers at correlation -1: Var = (a1 s1 - a2 s2)^2
d = rng.normal(size=5); d -= d.mean(); q1 = 0.2 + 0.02 * d; q2 = 0.2 - 0.03 * d
a1 = 0.4; qc = a1 * q1 + (1 - a1) * q2
assert np.isclose(qc.var(), (a1 * q1.std() - (1 - a1) * q2.std()) ** 2)
