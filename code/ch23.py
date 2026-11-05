"""Chapter 23: reproduces every number used in the chapter's examples,
tables, viz panel and the inhomogeneous-products section.
Run from the repo root:  python3 code/ch23.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def tau1(W):
    n = len(W)
    return max(0.5 * np.abs(W[i] - W[j]).sum() for i in range(n) for j in range(n))


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmax(np.real(ev))
    p = np.real(vl[:, k])
    return p / p.sum()


print("== Section: products of time-varying W ==")
X = np.array([[1, 0, 0], [0, 0, 1], [0.85, 0, 0.15]])
Y = np.array([[0, 0.95, 0.05], [0, 1, 0], [1, 0, 0]])
for name, M in [("X", X), ("Y", Y), ("XY", X @ Y), ("YX", Y @ X)]:
    assert np.allclose(M.sum(1), 1)
    print(f"{name}: Gamma = {gap(M):.3f}, tau1 = {tau1(M):.3f}, |lambda2| = {abs(eig_sorted(M)[1]):.3f}")
assert abs(gap(X) - 0.850) < 1e-3 and abs(gap(Y) - 0.776) < 1e-3
assert abs(gap(X @ Y) - 0.042) < 1e-3
# realized-product gap per period
print(f"realized two-period gap 1-|l2(XY)|^(1/2) = {1-abs(eig_sorted(X@Y)[1])**0.5:.3f}")

# tau1 properties on random matrices
rng = np.random.default_rng(23)
for _ in range(2000):
    A = rng.dirichlet(np.ones(4) * 0.5, size=4)
    B = rng.dirichlet(np.ones(4) * 0.5, size=4)
    assert tau1(A @ B) <= tau1(A) * tau1(B) + 1e-12
    assert abs(eig_sorted(A)[1]) <= tau1(A) + 1e-12
print("checked: tau1(AB) <= tau1(A) tau1(B) and |lambda2| <= tau1 on 2000 random pairs")

# Proof checks for Proposition ch23:prop:dobrushin
# Step 1: ||delta' A||_1 <= tau1(A) ||delta||_1 on zero-sum delta, equality at e_i - e_j
n_complex = 0
for _ in range(1000):
    n = int(rng.integers(2, 7))
    A = rng.dirichlet(np.ones(n) * 0.5, size=n)
    d = rng.normal(size=n); d -= d.mean()
    assert np.abs(d @ A).sum() <= tau1(A) * np.abs(d).sum() + 1e-12
    i, j = max(((i, j) for i in range(n) for j in range(n)),
               key=lambda ij: np.abs(A[ij[0]] - A[ij[1]]).sum())
    e = np.zeros(n); e[i], e[j] = 1, -1
    assert abs(np.abs(e @ A).sum() - 2 * tau1(A)) < 1e-12
    # (ii) every eigenvalue other than one copy of 1 (including complex ones) has modulus <= tau1
    ev = eig_sorted(A)
    n_complex += int(np.any(np.abs(np.imag(ev)) > 1e-9))
    assert np.all(np.abs(ev[1:]) <= tau1(A) + 1e-10)
# repeated unit eigenvalue (reducible): tau1 = 1, Gamma = 0
R = np.array([[0.6, 0.4, 0, 0], [0.3, 0.7, 0, 0], [0, 0, 0.5, 0.5], [0, 0, 0.2, 0.8]])
assert abs(tau1(R) - 1) < 1e-12 and abs(gap(R)) < 1e-12
# Hajnal bound: prod tau1 <= exp(-sum(1 - tau1))
taus = rng.uniform(0, 1, 200)
assert np.prod(taus) <= np.exp(-np.sum(1 - taus)) + 1e-300
print(f"checked: tau1 operator-norm identity, eigenvalue bound ({n_complex} cases with complex spectrum), reducible case, Hajnal bound")

# weak ergodicity: product of scrambling matrices collapses rows; absolute probability sequence
T = 60
Ws = [rng.dirichlet(np.ones(4), size=4) for _ in range(T)]
P = np.eye(4)
bound = 1.0
for t in range(T):
    P = Ws[t] @ P           # backward product W_{t} ... W_0
    bound *= tau1(Ws[t])
print(f"tau1 of 60-period backward product = {tau1(P):.2e} <= prod tau1 = {bound:.2e}")
assert tau1(P) <= bound + 1e-15
# absolute probability sequence: pi_s' = row of lim P (any row), pi_{t+1}' W_t = pi_t'
pis = [None] * (T + 1)
pis[T] = np.ones(4) / 4     # any terminal vector; backward recursion
for t in range(T - 1, -1, -1):
    pis[t] = pis[t + 1] @ Ws[t]
b0 = rng.uniform(0, 1, 4)
b = b0.copy()
for t in range(T):
    b = Ws[t] @ b
print(f"consensus after 60 periods = {b.mean():.4f}; pi_0'b_0 = {pis[0] @ b0:.4f}")
assert np.allclose(b, pis[0] @ b0, atol=1e-8)

print("\n== Perron-vector sensitivity (perturbation lemma) ==")


def two_block(c, m=3):
    W = np.full((2 * m, 2 * m), c / m)
    W[:m, :m] = (1 - c) / m
    W[m:, m:] = (1 - c) / m
    return W


sens = {}
for G in (0.13, 0.45):
    W = two_block(G / 2)
    p = perron(W)
    Z = np.linalg.inv(np.eye(6) - W + np.outer(np.ones(6), p))
    E = np.zeros((6, 6))
    E[0, 1] = -1.0
    E[0, 4] = 1.0               # move weight from within block to cross block in one row
    dpi = p @ E @ Z             # perturbation lemma: dpi' = pi' E Z
    h = 1e-6
    fd = (perron(W + h * E) - p) / h
    assert np.allclose(dpi, fd, atol=1e-4)
    sens[G] = dpi[:3].sum()
    print(f"Gamma = {gap(W):.2f}: d(pi_block1)/d(delta) = {sens[G]:+.3f}")
print(f"sensitivity ratio Gamma 0.13 vs 0.45 = {sens[0.13]/sens[0.45]:.2f}")

print("\n== HP filter: gain-1/2 cutoff period for annual data ==")
for lam in (6.25, 100, 400, 1600):
    w = np.arccos(1 - 1 / (2 * np.sqrt(lam)))
    print(f"lambda = {lam:7.2f}: cutoff period = {2*np.pi/w:5.1f} years")

print("\n== Bolivia epochs (illustrative levels) ==")
G = {1985: 0.11, 1990: 0.10, 2005: 0.12, 2009: 0.13, 2019: 0.13}
S = {1985: 0.65, 1990: 0.68, 2005: 0.63, 2009: 0.61, 2019: 0.62}
ECI = {2010: -0.85, 2019: -0.80}
print("Epoch 1 dGamma", round(G[1990] - G[1985], 2), "dsigma", round(S[1990] - S[1985], 2))
print("Epoch 2 dGamma", round(G[2005] - G[1990], 2), "dsigma", round(S[2005] - S[1990], 2))
print("Epoch 3 dGamma", round(G[2009] - G[2005], 2), "dsigma", round(S[2009] - S[2005], 2))
print("Epoch 4 dGamma", round(G[2019] - G[2009], 2), "dsigma", round(S[2019] - S[2009], 2), "dECI", round(ECI[2019] - ECI[2010], 2))
print("1985-2019 dGamma =", round(G[2019] - G[1985], 2))

print("\n== Viz panel coordinates: x = 10*Gamma, y = 6*sigma_B ==")
pts = {"BOL1985": (0.11, 0.65), "BOL1990": (0.10, 0.68), "BOL2000": (0.11, 0.65),
       "BOL2006": (0.12, 0.63), "BOL2010": (0.13, 0.61), "BOL2019": (0.13, 0.62),
       "KOR1960": (0.15, 0.50), "KOR1980": (0.30, 1 / 3), "KOR2000": (0.45, 0.20)}
for k, (g, s) in pts.items():
    print(f"{k}: ({10*g:.1f}, {6*s:.2f})")

print("\n== Velocities and SCI (illustrative) ==")
traj = {"South Korea": (0.15, 0.45, 1960, 2000), "China": (0.12, 0.35, 1980, 2020),
        "Chile": (0.18, 0.28, 1990, 2020), "Bolivia": (0.10, 0.13, 1990, 2019)}
target = 0.50
for c, (g0, g1, t0, t1) in traj.items():
    v = (g1 - g0) / ((t1 - t0) / 10)
    sci = (g1 - g0) / (target - g0)
    print(f"{c:12s} v = {v:+.4f}/decade  SCI = {sci:.3f}")
print("SCI Argentina-type (g0=0.20, g1=0.19):", round((0.19 - 0.20) / (target - 0.20), 3), "(negative)")
print("years to Gamma 0.45 at +0.01/decade:", round((0.45 - 0.13) / 0.01 * 10))
print("years to full closure at Bolivia's average SCI rate:", round(29 / 0.075))

print("\n== Gas decomposition (illustrative) ==")
pg = {2003: 0.07, 2014: 0.11, 2019: 0.09}
dS, dC = 0.02, 0.02
assert abs(pg[2014] - pg[2003] - (dS + dC)) < 1e-12 and abs(pg[2019] - pg[2003] - dS) < 1e-12
print("pi_gas 2003->2014 change", round(pg[2014] - pg[2003], 2), "= S", dS, "+ C", dC, "; 2019 keeps S only")
