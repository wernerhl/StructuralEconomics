"""Solutions manual, Chapter 15 (Operator Design): computes and prints every
number used in solutions/sol15.tex, with asserts on the key results.
Run from the repo root:  python3 code/solutions/sol15.py
numpy throughout; pandas only to read data/processed (Exercise 15.6)."""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=120)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROC = os.path.join(ROOT, "data", "processed")


def eigs(W):
    lam, V = np.linalg.eig(W)
    o = np.argsort(-np.abs(lam))
    lam, V = lam[o], V[:, o]
    return lam, V, np.linalg.inv(V)          # rows of U = V^-1: biorthogonal left vectors


def gap(W):
    return 1 - np.sort(np.abs(np.linalg.eigvals(W)))[-2]


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    p = np.real(vl[:, np.argmin(np.abs(ev - 1))])
    return p / p.sum()


def var_pi(p):
    return np.mean((p - 1 / len(p)) ** 2)


def grad_gamma(W):
    """G_ij = -Re(conj(l2) u_{2i} v_{2j}) / |l2|  (eq. ch15:eq:gamma_sensitivity)."""
    lam, V, U = eigs(W)
    return -np.real(np.conj(lam[1]) * np.outer(U[1], V[:, 1])) / abs(lam[1]), lam


def rownorm(K):
    return K / K.sum(1, keepdims=True)


# =====================================================================
print("== Exercise 15.1: optimal edge selection, 4x4 nearly block-diagonal W ==")
W = np.array([[0.49, 0.49, 0.01, 0.01], [0.49, 0.49, 0.01, 0.01],
              [0.01, 0.01, 0.49, 0.49], [0.01, 0.01, 0.49, 0.49]])
assert np.allclose(W.sum(1), 1)
lam, V, U = eigs(W)
print("eigenvalues:", np.round(np.real(lam), 6))
G, _ = grad_gamma(W)
print(f"lambda_2 = {np.real(lam[1]):.4f}, Gamma = {gap(W):.4f}")
print("G =\n", G)
assert abs(np.real(lam[1]) - 0.96) < 1e-12 and abs(gap(W) - 0.04) < 1e-12
assert np.allclose(np.abs(G), 0.25)
# all pairwise moves (i, j <- l), rate G_ij - G_il
moves = [(G[i, j] - G[i, l], i, j, l) for i in range(4) for j in range(4) for l in range(4) if j != l]
best = max(moves)
print(f"best rate G_ij - G_il = {best[0]:.4f} (e.g. row {best[1]+1}: from col {best[3]+1} to col {best[2]+1})")
assert abs(best[0] - 0.5) < 1e-12
n_best = sum(abs(m[0] - 0.5) < 1e-12 for m in moves)
print(f"number of tied optimal moves: {n_best} (every within -> cross move)")
assert n_best == 16
print("rate per unit of l1 budget (moving B/2): ", best[0] / 2)
# finite-difference check of the gradient for one move, and the exact second-order effect
h = 1e-7
E = np.zeros((4, 4)); E[0, 2] = 1; E[0, 0] = -1
fd = (gap(W + h * E) - gap(W)) / h
print(f"finite difference for move (1,1)->(1,3): {fd:.4f}")
assert abs(fd - 0.5) < 1e-5
for eps in (0.01, 0.05):
    print(f"exact dGamma for eps = {eps}: {gap(W + eps * E) - gap(W):.5f} (first order {0.5 * eps:.5f})")
dG001 = gap(W + 0.01 * E) - gap(W)
# within-block move has zero first-order effect
E2 = np.zeros((4, 4)); E2[0, 1] = 1; E2[0, 0] = -1
fd2 = (gap(W + h * E2) - gap(W)) / h
print(f"within-block move (1,1)->(1,2): finite difference {fd2:.2e}")
assert abs(fd2) < 1e-5
# symmetric spreading: same l1 budget split over all four rows
Es = np.zeros((4, 4))
for i, (jw, jc) in enumerate([(0, 2), (1, 3), (2, 0), (3, 1)]):
    Es[i, jc] += 0.25; Es[i, jw] -= 0.25          # ||Es||_1 = 2, same as E
print(f"spread over 4 rows, same l1 norm: first order {np.sum(G*Es):.4f}; exact at eps=0.01: "
      f"{gap(W + 0.01*Es) - gap(W):.5f}")
dGs001 = gap(W + 0.01 * Es) - gap(W)
assert abs(np.sum(G * Es) - 0.5) < 1e-12 and dGs001 > dG001
# perturbed asymmetric W: argmax still cross-block
rng = np.random.default_rng(151)
cross_count = 0
for _ in range(200):
    Wr = W + rng.uniform(-0.005, 0.005, (4, 4)) * (W > 0.1)
    Wr = rownorm(Wr)
    Gr, lr = grad_gamma(Wr)
    if abs(np.imag(lr[1])) > 1e-12:
        continue
    mv = max((Gr[i, j] - Gr[i, l], i, j, l) for i in range(4) for j in range(4) for l in range(4) if j != l)
    cross_count += ((mv[1] < 2) != (mv[2] < 2))
print(f"random asymmetric perturbations: optimal j is cross-block in {cross_count} of 200")
assert cross_count == 200

# =====================================================================
print("\n== Exercise 15.2: ISI-A vs ISI-B ==")
A = np.array([[0.33, 0.25, 0.01, 0.01], [0.25, 0.33, 0.01, 0.01],
              [0.012, 0.012, 0.15, 0.126], [0.012, 0.012, 0.126, 0.15]])
WP = rownorm(A)
print("row sums of A:", A.sum(1))


def out(A):
    return np.ones(len(A)) @ np.linalg.inv(np.eye(len(A)) - A.T) @ np.ones(len(A))


L = np.linalg.inv(np.eye(4) - A.T)
p0 = perron(WP)
print("W^P =\n", WP)
print(f"base: Gamma = {gap(WP):.4f}, lambda = {np.round(np.linalg.eigvals(WP), 4)}, pi = {p0}, "
      f"Var(pi) = {var_pi(p0):.5f}, total output = {out(A):.4f}")
print("L = (I-A')^-1 =\n", L)
BL, FL = L.sum(0), L.sum(1)
print("BL (column sums):", BL, " FL (row sums):", FL)
res = {}
for name, j in (("A", 1), ("B", 2)):
    A2 = A.copy(); A2[0, j] += 0.05
    W2 = rownorm(A2); p2 = perron(W2)
    res[name] = dict(dG=gap(W2) - gap(WP), dV=var_pi(p2) - var_pi(p0), dY=out(A2) - out(A), pi=p2,
                     G=gap(W2), V=var_pi(p2), Y=out(A2))
    # Sherman-Morrison: dY = eps * BL_j * FL_1 / (1 - eps * L_{1j})  (L indexed (row M1, col j))
    sm = 0.05 * BL[j] * FL[0] / (1 - 0.05 * L[0, j])
    print(f"ISI-{name}: Gamma = {res[name]['G']:.4f} (d {res[name]['dG']:+.4f}), Var(pi) = {res[name]['V']:.5f} "
          f"(d {res[name]['dV']:+.5f}), output = {res[name]['Y']:.4f} (d {res[name]['dY']:+.4f}); "
          f"Sherman-Morrison {sm:.4f}; first order 0.05*BL_j*FL_M1 = {0.05*BL[j]*FL[0]:.4f}; pi = {p2}")
    assert abs(sm - res[name]['dY']) < 1e-10
print(f"L_(M1,M2) = {L[0,1]:.4f}, L_(M1,T1) = {L[0,2]:.4f}")
rA, rB = res["A"], res["B"]
assert round(rA["dG"], 3) == -0.001 and round(rB["dG"], 3) == 0.037
assert round(rA["dV"], 4) == 0.0010 and round(rB["dV"], 4) == -0.0103
assert round(rA["dY"], 3) == 0.312 and round(rB["dY"], 3) == 0.181
assert abs(out(A) - 7.826) < 5e-4


# scaling the traditional recipes (rows T1, T2 of A) by k leaves W^P unchanged
def gains(k):
    Ak = A.copy(); Ak[2:] *= k
    out0 = out(Ak)
    d = []
    for j in (1, 2):
        A2 = Ak.copy(); A2[0, j] += 0.05
        d.append(out(A2) - out0)
    return d, Ak


lo, hi = 1.0, 3.0
assert gains(lo)[0][1] < gains(lo)[0][0] and gains(hi)[0][1] > gains(hi)[0][0]
for _ in range(80):
    mid = (lo + hi) / 2
    if gains(mid)[0][1] > gains(mid)[0][0]:
        hi = mid
    else:
        lo = mid
kstar = hi
dstar, Ak = gains(kstar)
print(f"k* = {kstar:.4f}: traditional row sum = {Ak[2].sum():.4f} (from {A[2].sum():.2f}); "
      f"rho(A) = {max(abs(np.linalg.eigvals(Ak))):.4f}; ISI-A dY = {dstar[0]:.4f}, ISI-B dY = {dstar[1]:.4f}")
assert np.allclose(rownorm(Ak), WP)                       # W^P unchanged by row scaling
Lk = np.linalg.inv(np.eye(4) - Ak.T)
print(f"at k*: BL_M2 = {Lk.sum(0)[1]:.4f}, BL_T1 = {Lk.sum(0)[2]:.4f}; L_(M1,M2) = {Lk[0,1]:.4f}, L_(M1,T1) = {Lk[0,2]:.4f}")
for k in (1.5, 2.5):
    d, Akk = gains(k)
    print(f"k = {k}: traditional row sum {Akk[2].sum():.2f}, ISI-A dY = {d[0]:.4f}, ISI-B dY = {d[1]:.4f}")
d20 = gains(2.5)[0]
assert d20[1] > d20[0]
# alternative: raise only the traditional sectors' purchases from the modern block (a_T,M)
def gains_c(c):
    Ac = A.copy(); Ac[2:, :2] = c
    o = out(Ac); r = []
    for j in (1, 2):
        A2 = Ac.copy(); A2[0, j] += 0.05; r.append(out(A2) - o)
    return r, Ac
rc, Ac = gains_c(0.10)
print(f"variant a_T,M = 0.10: ISI-A dY = {rc[0]:.4f}, ISI-B dY = {rc[1]:.4f}; "
      f"Gamma(W^P) = {gap(rownorm(Ac)):.4f}")

# =====================================================================
print("\n== Exercise 15.3: the smartphone operator in the Chapter 10 formation model (two blocks) ==")
n = 3; N = 2 * n
blk = np.array([0] * n + [1] * n)                      # 0 = modern (M), 1 = traditional (T)
same = blk[:, None] == blk[None, :]
WPf = np.where(same, 0.94 / n, 0.06 / n)
WI0 = np.where(same, 0.85 / n, 0.15 / n)               # cross weight 0.05 per pair
beta = np.where(blk == 0, 0.10, 0.04); gam = 0.02
pi_b = np.where(blk == 0, 0.20, 1 / 3 - 0.20)           # per agent, sums to 1
assert abs(pi_b.sum() - 1) < 1e-12
omega, Phi, phi = 0.5, 3.0, 1.0                         # Omega = 6, Omega/N = 1
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.1
P = np.eye(N) - np.ones((N, N)) / N
C = np.diag(beta) + gam * np.diag(np.ones(N)) @ WPf


def alloc(th):
    return omega + pi_b * Phi + phi * (th - th.mean())


def returns(th, alpha):
    x = alloc(th)
    return alpha + beta * x + gam * WPf @ x


def T(th, alpha, WI):
    return rho * th + (1 - rho) * th0 + mu * (WI @ returns(th, alpha) - rstar)


def steady(alpha, WI):
    J = rho * np.eye(N) + mu * phi * WI @ C @ P
    c = T(np.zeros(N), alpha, WI)
    return np.linalg.solve(np.eye(N) - J, c), J


# calibrate block baselines alpha_M, alpha_T so that theta*_M = 1.09, theta*_T = 0.92 (Example 15.11)
def th_of(aM, aT):
    return steady(np.where(blk == 0, aM, aT), WI0)[0]


base = th_of(0, 0); dM = th_of(1, 0) - base; dT = th_of(0, 1) - base
Mx = np.array([[dM[0], dT[0]], [dM[-1], dT[-1]]])
aM, aT = np.linalg.solve(Mx, np.array([1.09, 0.92]) - np.array([base[0], base[-1]]))
alpha = np.where(blk == 0, aM, aT)
th_ss0, J0 = steady(alpha, WI0)
print(f"calibrated alpha_M = {aM:.4f}, alpha_T = {aT:.4f}; theta* = {th_ss0[0]:.4f} (M), {th_ss0[-1]:.4f} (T); "
      f"rho(J) = {max(abs(np.linalg.eigvals(J0))):.4f}")
assert abs(th_ss0[0] - 1.09) < 1e-10 and abs(th_ss0[-1] - 0.92) < 1e-10
r0 = returns(th_ss0, alpha); x0 = alloc(th_ss0)
print(f"base returns r_M = {r0[0]:.4f}, r_T = {r0[-1]:.4f}; allocations x_M = {x0[0]:.4f}, x_T = {x0[-1]:.4f}")
# the smartphone: rural rows gain +0.05 on each urban/modern column, own-block weights cut equally
WI1 = WI0.copy()
WI1[blk == 1] = np.where(blk[None, :] == 0, 0.15 / n + 0.05, (0.85 - n * 0.05) / n)
assert np.allclose(WI1.sum(1), 1) and np.allclose(WI1[blk == 1][:, blk == 0], 0.10)
th_ss1, J1 = steady(alpha, WI1)
x1 = alloc(th_ss1); r1 = returns(th_ss1, alpha)
print(f"new steady state: theta*_M = {th_ss1[0]:.4f}, theta*_T = {th_ss1[-1]:.4f}; x_M = {x1[0]:.4f}, x_T = {x1[-1]:.4f}")
path = [th_ss0]
for t in range(5):
    path.append(T(path[-1], alpha, WI1))
path = np.array(path)
print("theta path (M, T), generations 0-5:")
for t, th in enumerate(path):
    print(f"  t={t}: theta_M = {th[0]:.4f}, theta_T = {th[-1]:.4f}, gap = {th[0]-th[-1]:.4f}, x_T = {alloc(th)[-1]:.4f}")
gap0 = th_ss0[0] - th_ss0[-1]; gap1 = th_ss1[0] - th_ss1[-1]
dthT = th_ss1[-1] - th_ss0[-1]; dthM = th_ss1[0] - th_ss0[0]
share5 = (path[5, -1] - th_ss0[-1]) / dthT
print(f"steady-state d theta_T = {dthT:+.4f}, d theta_M = {dthM:+.4f}; gap {gap0:.4f} -> {gap1:.4f} "
      f"(closes {1-gap1/gap0:.1%}); after 1 generation the gap is {path[1,0]-path[1,-1]:.4f}; "
      f"share of long-run d theta_T reached after 5 generations = {share5:.3f}")
dxT = x1[-1] - x0[-1]; dxM = x1[0] - x0[0]
print(f"steady-state d x_T = {dxT:+.4f}, d x_M = {dxM:+.4f}; sum of d x = {np.sum(x1-x0):+.2e}")
assert abs(np.sum(x1 - x0)) < 1e-12 and dxT > 0
# enforcement effect at theta^0: x and r do not depend on W^I
E_x = alloc(th_ss0) - alloc(th_ss0)
E_r = returns(th_ss0, alpha) - returns(th_ss0, alpha)
print(f"enforcement effect on x and r at fixed theta: max |E| = {max(abs(E_x).max(), abs(E_r).max()):.1e}")
print(f"aggregate return sum: before {r0.sum():.5f}, after {r1.sum():.5f} (change {r1.sum()-r0.sum():+.5f})")
print(f"rho(J) after = {max(abs(np.linalg.eigvals(J1))):.4f}; half-life of deviations ln2/-ln(rho(J)) = "
      f"{np.log(2)/-np.log(max(abs(np.linalg.eigvals(J1)))):.2f} generations")
SP = dict(aM=aM, aT=aT, dthT=dthT, gap0=gap0, gap1=gap1, dxT=dxT)

# =====================================================================
print("\n== Exercise 15.4: sequencing in a stylized 6-sector model (gas G, urban U, rural R) ==")
lab = ["G1", "G2", "U1", "U2", "R1", "R2"]
bl6 = np.array([0, 0, 1, 1, 2, 2])
K0 = np.where(bl6[:, None] == bl6[None, :], 1.0, 0.0)
cross = {(0, 1): 0.04, (1, 2): 0.04, (0, 2): 0.01}
for (a, b), v in cross.items():
    K0[np.ix_(bl6 == a, bl6 == b)] = v; K0[np.ix_(bl6 == b, bl6 == a)] = v
print(f"base Gamma(W^P) = {gap(rownorm(K0)):.4f}")
G_TARGET = 0.30


def simulate(order, eta, T_=10, starts=(1, 4, 7)):
    """order: list of operators; eta: complementarity exponent of ISI uptake.
    State: raw intensities K (W^P = D^-1 K), rural aspiration ratio a = theta_R/theta_U,
    rural capability c.  Information raises the rural cross-block W^I weight from 0.10 to 0.30,
    and theta_R then evolves by the formation equation with fixed returns.  Education adds
    U-R and G-R raw intensity 0.05 and raises c from 0.5 to 1.  ISI-B adds G-U and G-R
    intensity 0.30 * (a*c)^eta (the enclave buys domestic inputs).  Each operator acts once, at the start of its phase."""
    K = K0.copy(); thR, thU = 0.92, 1.09; c = 0.5; wcross = 0.10
    rU, rR = 0.20, 0.08
    gam_path = [gap(rownorm(K))]
    when = dict(zip(order, starts))
    for t in range(1, T_ + 1):
        for op, s in when.items():
            if s != t:
                continue
            if op == "info":
                wcross = 0.30
            elif op == "educ":
                for (a_, b_) in ((1, 2), (0, 2)):
                    K[np.ix_(bl6 == a_, bl6 == b_)] += 0.05; K[np.ix_(bl6 == b_, bl6 == a_)] += 0.05
                c = 1.0
            elif op == "isi":
                up = (min(1.0, thR / thU) * c) ** eta
                for (a_, b_) in ((0, 1), (0, 2)):
                    K[np.ix_(bl6 == a_, bl6 == b_)] += 0.30 * up; K[np.ix_(bl6 == b_, bl6 == a_)] += 0.30 * up
        sR = (1 - wcross) * rR + wcross * rU
        thR = rho * thR + (1 - rho) * th0 + mu * (sR - rstar)
        gam_path.append(gap(rownorm(K)))
    return np.array(gam_path)


def first_hit(g):
    idx = np.where(g >= G_TARGET)[0]
    return int(idx[0]) if len(idx) else None


SEQ = {"A": ["isi", "educ", "info"], "B": ["info", "educ", "isi"]}
results = {}
for eta in (0.0, 1.0, 2.0):
    for nm, order in SEQ.items():
        g = simulate(order, eta)
        results[(nm, eta)] = g
        print(f"eta = {eta}: sequence {nm} {order}: Gamma_t = {np.round(g, 3)}; first t with Gamma >= 0.30: {first_hit(g)}")
gA0, gB0 = results[("A", 0.0)], results[("B", 0.0)]
assert first_hit(gA0) < first_hit(gB0)                 # no complementarity: ISI first wins
assert abs(gA0[-1] - gB0[-1]) < 1e-12                   # same final W when operators are additive
# rural aspiration path under information
thR = 0.92; ths = []
for t in range(10):
    thR = rho * thR + (1 - rho) * th0 + mu * (0.7 * 0.08 + 0.3 * 0.20 - rstar); ths.append(thR)
thR_noinfo = 0.92; ths0 = []
for t in range(10):
    thR_noinfo = rho * thR_noinfo + (1 - rho) * th0 + mu * (0.9 * 0.08 + 0.1 * 0.20 - rstar); ths0.append(thR_noinfo)
print("theta_R with information from t=1:", np.round(ths, 4))
print("theta_R without information      :", np.round(ths0, 4))
print(f"uptake base a*c for ISI at t=1 in A: {(0.92/1.09*0.5):.4f}")
# exact uptake used in B at t = 7: theta_R after 6 periods of information (applied at t=1)
thR = 0.92
for t in range(1, 7):
    sR = 0.7 * 0.08 + 0.3 * 0.20
    thR = rho * thR + (1 - rho) * th0 + mu * (sR - rstar)
upB = min(1, thR / 1.09) * 1.0
print(f"theta_R entering t=7 under B = {thR:.4f}; uptake base a*c = {upB:.4f}")
# crossover complementarity: smallest eta at which B reaches the target no later than A
etas = np.round(np.arange(0, 6.01, 0.25), 2)
cross_eta = None
for e in etas:
    hA, hB = first_hit(simulate(SEQ["A"], e)), first_hit(simulate(SEQ["B"], e))
    hA = 99 if hA is None else hA; hB = 99 if hB is None else hB
    if hB <= hA:
        cross_eta = e; break
print(f"smallest eta on the 0.25 grid at which B reaches 0.30 no later than A: {cross_eta}")
SEQ_RES = dict(cross_eta=cross_eta)

# =====================================================================
print("\n== Exercise 15.6: Hirschman vs spectral gap on INE 2019 (base 2017, book12, domestic) ==")
Z = pd.read_csv(os.path.join(PROC, "bolivia_Z", "Z_2017_2019_book12.csv"), index_col=0)
lk = pd.read_csv(os.path.join(PROC, "bolivia_io_linkages.csv"))
lk = lk[(lk.year == 2019) & (lk.base == 2017) & (lk.aggregation == "book12") & (lk.flows == "domestic")].set_index("sector")
sec = list(Z.index)
X = lk.loc[sec, "output"].values
Zv = Z.values
Abk = (Zv / X[None, :]).T                           # a_ss' = Z[s', s] / X_s (row s = recipe of s)
print("row sums of A (domestic intermediate input per unit output):", dict(zip(sec, np.round(Abk.sum(1), 3))))
Lb = np.linalg.inv(np.eye(len(sec)) - Abk.T)
BLb, FLb = Lb.sum(0), Lb.sum(1)
for s, b, f in zip(sec, BLb, FLb):
    print(f"  {s:22s} BL = {b:.3f}  FL = {f:.3f}  (CSV BL {lk.loc[s,'BL']:.3f}, FL {lk.loc[s,'FL']:.3f})")
assert np.allclose(BLb, lk.loc[sec, "BL"].values, atol=1e-6) and np.allclose(FLb, lk.loc[sec, "FL"].values, atol=1e-6)
print(f"mean BL = {BLb.mean():.4f} = mean FL = {FLb.mean():.4f}")
assert abs(BLb.mean() - FLb.mean()) < 1e-12 and BLb.min() >= 1 - 1e-12
# W^P on the 11 sectors with intermediate flows (domestic services dropped)
keep = [i for i, s in enumerate(sec) if s != "Domestic services"]
s11 = [sec[i] for i in keep]
A11 = Abk[np.ix_(keep, keep)]
W11 = rownorm(A11)
for gmm in (0.5, 0.9):
    Lg = np.linalg.inv(np.eye(11) - gmm * W11)
    Lgt = np.linalg.inv(np.eye(11) - gmm * W11.T)
    print(f"gamma = {gmm}: row sums of (I - g W^P)^-1 all = {Lg.sum(1).min():.4f}..{Lg.sum(1).max():.4f} "
          f"(1/(1-g) = {1/(1-gmm):.4f}); column sums of (I - g W^P')^-1: {Lgt.sum(0).min():.4f}..{Lgt.sum(0).max():.4f}")
    assert np.allclose(Lg.sum(1), 1 / (1 - gmm)) and np.allclose(Lgt.sum(0), 1 / (1 - gmm))
infl = (1 - 0.99) / 11 * np.ones(11) @ np.linalg.inv(np.eye(11) - 0.99 * W11)
p11 = perron(W11)
print(f"influence vector at gamma = 0.99 vs pi: max abs diff = {np.abs(infl - p11).max():.4f}")
G11, lam11 = grad_gamma(W11)
print(f"lambda_2 = {lam11[1]:.4f}, |lambda_3| = {abs(lam11[2]):.4f}, Gamma = {gap(W11):.4f}")
assert abs(gap(W11) - 0.467) < 5e-4 and abs(np.imag(lam11[1])) < 1e-12
lam_, V_, U_ = eigs(W11)
u2, v2 = np.real(U_[1]), np.real(V_[:, 1])
sgn = np.sign(v2[np.argmax(np.abs(v2))])
print("v2 (right, scaled max|.|=1):", dict(zip(s11, np.round(sgn * v2 / np.abs(v2).max(), 3))))
print("u2 (left, scaled max|.|=1): ", dict(zip(s11, np.round(sgn * u2 / np.abs(u2).max(), 3))))
print("G (rows = buyers i, cols = suppliers j):")
print(pd.DataFrame(G11, index=s11, columns=[s[:6] for s in s11]).round(3).to_string())
# finite-difference check
i0, j0, l0 = 0, 1, 6
E = np.zeros((11, 11)); E[i0, j0] = 1; E[i0, l0] = -1
fd = (gap(W11 + 1e-7 * E) - gap(W11)) / 1e-7
assert abs(fd - (G11[i0, j0] - G11[i0, l0])) < 1e-4
# best feasible pairwise moves: take weight from an entry with W_il > 0
mv = sorted(((G11[i, j] - G11[i, l], i, j, l) for i in range(11) for j in range(11) for l in range(11)
             if j != l and W11[i, l] > 0.01), reverse=True)
print("top 8 moves (rate, buyer, to supplier, from supplier, W_il):")
for r_, i, j, l in mv[:8]:
    print(f"  {r_:.3f}  {s11[i]:22s} -> {s11[j]:22s} from {s11[l]:22s} (W_il = {W11[i,l]:.3f})")
best15 = mv[0]
mv_off = [m for m in mv if m[3] != m[1] and m[2] != m[1]]
print("best move that neither uses nor creates a self-purchase:")
for r_, i, j, l in mv_off[:3]:
    print(f"  {r_:.3f}  {s11[i]:22s} -> {s11[j]:22s} from {s11[l]:22s} (W_il = {W11[i,l]:.3f})")
print(f"diagonal of W^P: {dict(zip(s11, np.round(np.diag(W11), 3)))}")
# exact check of the best move with l1 budget 0.1 (move 0.05)
Eb = np.zeros((11, 11)); Eb[best15[1], best15[2]] = 0.05; Eb[best15[1], best15[3]] = -0.05
print(f"exact dGamma of best move with 0.05 of weight: {gap(W11 + Eb) - gap(W11):.4f} (first order {0.05*best15[0]:.4f})")
best_exact = gap(W11 + Eb) - gap(W11)
# supplier score: best rate over buyers i and donors l of moving weight to supplier j
sup_score = np.array([max(G11[i, j] - G11[i, l] for i in range(11) for l in range(11) if l != j and W11[i, l] > 0.01)
                      for j in range(11)])
buy_score = np.array([max(G11[i, j] - G11[i, l] for j in range(11) for l in range(11) if l != j and W11[i, l] > 0.01)
                      for i in range(11)])
BL11 = BLb[keep]
rk = lambda x: np.argsort(np.argsort(-x))
spear_sup = np.corrcoef(rk(BL11), rk(sup_score))[0, 1]
spear_buy = np.corrcoef(rk(BL11), rk(buy_score))[0, 1]
print("sector: BL, rank BL, supplier score (best dGamma/dt into j), buyer score (best in row i)")
for k, s in enumerate(s11):
    print(f"  {s:22s} {BL11[k]:.3f} {rk(BL11)[k]+1:2d}  {sup_score[k]:+.3f} ({rk(sup_score)[k]+1:2d})  {buy_score[k]:+.3f} ({rk(buy_score)[k]+1:2d})")
print(f"Spearman(BL, supplier score) = {spear_sup:.3f}; Spearman(BL, buyer score) = {spear_buy:.3f}")
H15 = dict(best=best15, best_exact=best_exact, spear_sup=spear_sup, spear_buy=spear_buy)

# =====================================================================
print("\n== Exercise 15.7: political cost, stylized N = 4 (Gas, Government, Urban, Rural) ==")
pi4 = np.array([0.35, 0.30, 0.20, 0.15])
Q = np.full((4, 4), 0.25)
W4 = np.minimum(1, pi4[None, :] / pi4[:, None]) * Q        # Metropolis chain, reversible w.r.t. pi4
np.fill_diagonal(W4, 0)
np.fill_diagonal(W4, 1 - W4.sum(1))
print("W =\n", W4)
assert np.allclose(perron(W4), pi4) and W4.min() > 0.1
lam4 = np.sort(np.abs(np.linalg.eigvals(W4)))[::-1]
print(f"eigenvalue moduli = {lam4}; Gamma = {gap(W4):.4f}; Var(pi) = {var_pi(pi4):.6f}")
Ed = np.outer(np.ones(4), pi4) - W4
nF = np.linalg.norm(Ed); n1 = np.abs(Ed).sum()
print(f"||1 pi' - W||_F = {nF:.4f}, ||.||_1 = {n1:.4f}")
for B_ in (0.1, nF):
    tb = min(1, B_ / nF)
    Wt = W4 + tb * Ed
    print(f"budget {B_:.4f}: t = {tb:.4f}, Gamma = {gap(Wt):.4f} (= 1-(1-t)|l2| = {1-(1-tb)*lam4[1]:.4f}), "
          f"pi = {perron(Wt)}")
    assert np.allclose(perron(Wt), pi4)
tb01 = 0.1 / nF
# (b) Var(pi) under a Frobenius budget B = 0.1.
# Reachability: a target pi~ is the Perron vector of W + E for some E with E1 = 0 and
# ||E||_F <= B iff ||(I - W')pi~|| <= B ||pi~||; the least-norm E is pi~ b'/||pi~||^2, b = (I - W')pi~.
BUD = 0.1
Zf = np.linalg.inv(np.eye(4) - W4 + np.outer(np.ones(4), pi4))


def least_norm_E(pt):
    b = pt - W4.T @ pt
    return np.outer(pt, b) / (pt @ pt)


# check the lemma on random targets: E* works and no random feasible E with the same pi~ is smaller
for _ in range(50):
    pt = pi4 + rng.normal(scale=0.02, size=4); pt -= (pt.sum() - 1) / 4
    Es = least_norm_E(pt)
    assert np.allclose(Es.sum(1), 0) and np.allclose(pt @ (W4 + Es), pt)
    assert abs(np.linalg.norm(Es) - np.linalg.norm(pt - W4.T @ pt) / np.linalg.norm(pt)) < 1e-12
    # any other solution: E* + F with F1 = 0 and F'pt = 0 has larger norm (orthogonality)
    F = rng.normal(size=(4, 4)); F -= F.mean(1, keepdims=True)
    F -= np.outer(pt, pt @ F) / (pt @ pt); F -= F.mean(1, keepdims=True)
    F -= np.outer(pt, pt @ F) / (pt @ pt)
    assert np.linalg.norm(Es + F) >= np.linalg.norm(Es) - 1e-12


def best_target(cons, step=0.001, half=0.10):
    """min Var(pi~) over the simplex, subject to reachability and pi~_k >= pi_k for k in cons."""
    g = np.arange(-half, half + 1e-12, step)
    best = (np.inf, None)
    for d1 in g:
        a, b_ = np.meshgrid(g, g, indexing="ij")
        d = np.stack([np.full_like(a, d1), a, b_, -(d1 + a + b_)], -1).reshape(-1, 4)
        pt = pi4 + d
        ok = (pt > 0).all(1)
        for k in cons:
            ok &= pt[:, k] >= pi4[k] - 1e-12
        bb = pt - pt @ W4                                  # rows: (I - W')pi~ for each target
        ok &= np.linalg.norm(bb, axis=1) <= BUD * np.linalg.norm(pt, axis=1)
        if ok.any():
            v = np.mean((pt[ok] - 0.25) ** 2, 1)
            k = np.argmin(v)
            if v[k] < best[0]:
                best = (v[k], pt[ok][k])
    return best


rows = {}
for nm, S in {"unconstrained": [], "{Gas, Gov}": [0, 1], "{Gas}": [0], "{Gas, Gov, Urban}": [0, 1, 2],
              "{Gov, Urban}": [1, 2], "{Urban, Rural}": [2, 3]}.items():
    v, pt = best_target(S)
    Es = least_norm_E(pt)
    W2 = W4 + Es
    assert W2.min() >= 0 and np.allclose(perron(W2), pt, atol=1e-9)
    rows[nm] = (v, pt, np.linalg.norm(Es))
    print(f"{nm:18s}: min Var = {v:.6f} (Omega_2 = {var_pi(pi4)-v:.6f}); pi~ = {np.round(pt, 4)}; "
          f"||E*||_F = {np.linalg.norm(Es):.4f}; min entry of W+E* = {W2.min():.4f}")
Om_u = var_pi(pi4) - rows["unconstrained"][0]
Om_c = var_pi(pi4) - rows["{Gas, Gov}"][0]
Cpol = Om_u - Om_c
print(f"Omega_2 unconstrained = {Om_u:.6f}, constrained = {Om_c:.6f}, political cost = {Cpol:.6f} "
      f"({Cpol/Om_u:.1%} of the unconstrained gain)")
for nm in rows:
    print(f"  political cost for coalition {nm:18s}: {Om_u - (var_pi(pi4) - rows[nm][0]):.6f}")
pt_c = rows["{Gas, Gov}"][1]
assert abs(pt_c[2] - pt_c[3]) < 2e-3 and abs(pt_c[0] - 0.35) < 1e-9 and abs(pt_c[1] - 0.30) < 1e-9
assert rows["{Urban, Rural}"][0] == rows["unconstrained"][0]
# the unconstrained optimum, and the best uniform target (is pi = 1/N reachable?)
pu = np.ones(4) / 4
print(f"reaching uniform pi needs ||E||_F >= {np.linalg.norm(pu - W4.T @ pu)/np.linalg.norm(pu):.4f}")
Bneed_u = np.linalg.norm(pu - W4.T @ pu) / np.linalg.norm(pu)
# first-order (perturbation-lemma) cross-check of the direction: dVar/dt at t=0 along E* (unconstrained)
Eu = least_norm_E(rows["unconstrained"][1])
dvar_fo = (2 / 4) * (pi4 @ Eu @ Zf) @ (pi4 - 0.25)
dvar_fd = (var_pi(perron(W4 + 1e-6 * Eu)) - var_pi(pi4)) / 1e-6
print(f"perturbation lemma along E*: dVar/dt = {dvar_fo:.6f}, finite difference = {dvar_fd:.6f}")
assert abs(dvar_fo - dvar_fd) < 1e-5
print("\nall Chapter 15 solution checks passed")
