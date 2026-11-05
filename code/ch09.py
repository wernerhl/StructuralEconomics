"""Chapter 9: reproduces every number used in the worked examples, viz panels
and exercise answers.  Run from the repo root:  python3 code/ch09.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# Example ch09:ex:two_steady  (heterogeneous steady states)
# theta_i = thbar + M (sum_j W_ij r_j - r*),  r_M = .05 + .1 th,  r_T = .03 + .02 th
rho, mu, thbar, rstar = 0.8, 0.3, 1.0, 0.1
M = mu / (1 - rho)
a = np.array([0.05, 0.03])
b = np.array([0.10, 0.02])


def steady(W):
    A = np.eye(2) - M * W @ np.diag(b)
    return np.linalg.solve(A, thbar + M * (W @ a - rstar))


print("== Example two_steady")
print("M = mu/(1-rho) =", M)
for w in [0.0, 0.3, 0.4, 0.5]:
    W = np.array([[1 - w, w], [w, 1 - w]])
    th = steady(W)
    print(f"cross weight w={w}: theta=({th[0]:.4f}, {th[1]:.4f}), gap={th[0]-th[1]:.4f}")
th_block = steady(np.eye(2))
assert abs(th_block[0] - 1.0882) < 1e-4 and abs(th_block[1] - 0.9227) < 1e-4
slopes = rho + mu * b
print("map slopes rho+mu*g' (M,T):", slopes)
print("contraction constants M*g' (M,T):", M * b)
assert np.all(M * b < 1)
# Phase diagram (scaled by 2.5): theta' = slope*theta + intercept
icpt = (1 - rho) * thbar + mu * (a - rstar)
print("intercepts:", icpt, "scaled x2.5:", 2.5 * icpt, "line end at X=5:", 2.5 * icpt + 5 * slopes)
print("fixed points scaled:", 2.5 * icpt / (1 - slopes))
assert np.allclose(icpt / (1 - slopes), th_block)

# ---------------------------------------------------------------------------
# Proposition ch09:thm:hysteresis and Example ch09:ex:hysteresis (persistence)
print("\n== Persistence / half-lives (lambda = rho when Dg = 0)")
for r in [0.95, 0.90, 0.80, 0.70, 0.50]:
    h = np.log(2) / np.log(1 / r)
    print(f"rho={r}: t_half={h:.2f} generations = {25*h:.0f} years")
print("share of gap closed after 1,2,4 generations (lambda=0.9):", [round(1 - 0.9 ** t, 3) for t in (1, 2, 4)])
print("speed ratio rho=0.9 vs 0.7:", (np.log(2) / np.log(1 / 0.9)) / (np.log(2) / np.log(1 / 0.7)))
lam_M = slopes[0]
print("Modern agent lambda =", lam_M, "half-life", np.log(2) / np.log(1 / lam_M), "vs rho-only", np.log(2) / np.log(1 / rho))
# lambda >= rho with nonnegative feedback (Perron-Frobenius)
rng = np.random.default_rng(0)
for _ in range(200):
    W = rng.dirichlet(np.ones(4), size=4)
    D = np.diag(rng.uniform(0, 0.3, 4))
    J = 0.8 * np.eye(4) + 0.3 * W @ D
    assert max(abs(np.linalg.eigvals(J))) >= 0.8 - 1e-12

# ---------------------------------------------------------------------------
# Exercises
print("\n== Exercises")
# 9.1
th = 0.5
traj = []
for t in range(10):
    th = 0.8 * th + 0.2 * 1 + 0.2 * (0.15 - 0.1)
    traj.append(th)
ss = 1 + (0.2 / 0.2) * 0.05
print("9.1 theta_1..10:", np.round(traj, 4), "steady state", ss)
n90 = next(t for t in range(1, 200) if (ss - (ss - (ss - 0.5) * 0.8 ** t)) >= 0 and abs(ss - (ss - (ss - 0.5) * 0.8 ** t)) <= 0.1 * (ss - 0.5))
print("   periods to close 90% of gap (rho=0.8):", n90, "; rho=0.95:", int(np.ceil(np.log(0.1) / np.log(0.95))))
# 9.2(b),(c)
W = np.array([[0.7, 0.3], [0.3, 0.7]])
print("9.2(b) theta:", steady(W), " (c) coincide at w = 0.5")
# 9.3
s3, c3 = 0.85 + 0.1 * 0.3, 0.15 * 0.5 - 0.1 * 0.1
th = 0.2
path = []
for t in range(50):
    if t == 10:
        th = 0.6
    path.append(th)
    th = s3 * th + c3
print("9.3 slope", s3, "unique steady state", c3 / (1 - s3), "theta_9", round(path[9], 4), "theta_49", round(path[49], 4))
# 9.4 smartphone: T observes M with weight 0.3, M unchanged
W = np.array([[1.0, 0.0], [0.3, 0.7]])
new = steady(W)
print("9.4 old T steady state", th_block[1], "new", new, "half-gap generations (lambda_T=0.806 if M fixed):",
      np.log(2) / np.log(1 / (rho + mu * 0.7 * 0.02)))
# 9.6
E = 1 / 0.9 - 1
F = 1.05 / 0.9 - 1 / 0.9
print(f"9.6 enforcement {E:.4f}, formation {F:.4f}, Psi = {F/E:.3f}")
# 9.7
th7 = 0.70 + 4 * (0.4 * 0.85 + 0.6 * 0.70 - 0.70)
w_star = (1 - 0.70) / (4 * (0.85 - 0.70))
print("9.7 theta* =", th7, " exceeds 1 when w >", w_star)
assert abs(th7 - 0.94) < 1e-12 and abs(w_star - 0.5) < 1e-12
# Used in Chapter 12, Appendix 12B: connecting the agents lowers mean theta*
print("\nmean theta*: block", steady(np.eye(2)).mean(), " w=0.3", steady(np.array([[0.7, 0.3], [0.3, 0.7]])).mean(),
      " w=0.5", steady(np.array([[0.5, 0.5], [0.5, 0.5]])).mean())
# Viz panel heatmap colours (normalized to each matrix's max): block (1,0) -> 100/0; connected (0.6,0.4) -> 100/67
print("heatmap colours connected:", round(100 * 0.6 / 0.6), round(100 * 0.4 / 0.6))

# ---------------------------------------------------------------------------
# Proof checks (proof pass)
print("\n== Proof checks")
rng = np.random.default_rng(9)
# Prop ch09:thm:hysteresis: with Dg >= 0, lambda = rho + rho(mu W Dg) exactly
for _ in range(300):
    n = rng.integers(2, 7)
    W = rng.dirichlet(np.ones(n), size=n)
    Dg = rng.uniform(0, 0.5, (n, n)) * (rng.random((n, n)) < 0.5)
    r0, m0 = rng.uniform(0.1, 0.9), rng.uniform(0.01, 0.5)
    A = m0 * W @ Dg
    lam = max(abs(np.linalg.eigvals(r0 * np.eye(n) + A)))
    assert abs(lam - (r0 + max(abs(np.linalg.eigvals(A))))) < 1e-9
    assert max(abs(np.linalg.eigvals(A))) >= max(np.diag(A)) - 1e-12
# lambda = rho with Dg != 0: nobody observes agent 2 (column 2 of W is zero), only g_2 responds
Wn = np.array([[1.0, 0.0], [1.0, 0.0]])
Dgn = np.array([[0.0, 0.0], [0.0, 0.4]])
assert abs(max(abs(np.linalg.eigvals(0.8 * np.eye(2) + 0.3 * Wn @ Dgn))) - 0.8) < 1e-12
print("lambda = rho + rho(mu W Dg) verified; lambda = rho with Dg != 0 example verified")
# Local contraction bound (affine case): ||J^t|| <= C (lambda+eps)^t
J = 0.8 * np.eye(3) + 0.3 * rng.dirichlet(np.ones(3), size=3) @ np.diag([0.1, 0.2, 0.05])
lamJ = max(abs(np.linalg.eigvals(J)))
ratios = [np.linalg.norm(np.linalg.matrix_power(J, t), 2) / (lamJ + 0.01) ** t for t in range(1, 400)]
assert max(ratios[200:]) <= max(ratios[:200]) + 1e-9
# Prop ch09:thm:self_confirming (affine case): unique steady state is block-constant
blocks = [[0, 1, 2], [3, 4]]
Wb = np.zeros((5, 5))
Gb = np.zeros((5, 5))
ab = np.zeros(5)
for B, (alpha, beta) in zip(blocks, [(0.05, 0.1), (0.03, 0.02)]):
    Wb[np.ix_(B, B)] = rng.dirichlet(np.ones(len(B)), size=len(B))   # arbitrary within-block weights
    row = rng.dirichlet(np.ones(len(B))) * beta                         # identical rows within block
    for j in B:
        Gb[j, B] = row
        ab[j] = alpha
thb = np.linalg.solve(np.eye(5) - M * Wb @ Gb, thbar + M * (Wb @ ab - rstar))
assert np.allclose(thb[:3], thb[0]) and np.allclose(thb[3:], thb[3])
assert np.allclose(thb[[0, 3]], th_block)
print("block-wise steady state is block-constant:", np.round(thb, 4))
