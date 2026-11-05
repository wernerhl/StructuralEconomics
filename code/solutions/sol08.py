"""Solutions manual, Chapter 8 (Belief Manipulation and Epistemic Power).
Computes and prints every number used in solutions/sol08.tex.
Run from the repository root:  python3 code/solutions/sol08.py      (numpy only)"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def perron(W):
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmin(np.abs(vals - 1))
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def eigs(W):
    lam = np.linalg.eigvals(W)
    return lam[np.argsort(-np.abs(lam))]


def gap(W):
    return 1 - np.abs(eigs(W))[1]


def neff(p):
    return 1 / np.sum(p ** 2)


def halflife(l2):
    return np.log(2) / np.log(1 / l2) if 0 < l2 < 1 else (0.0 if l2 == 0 else np.inf)


def fundamental(W):
    p = perron(W)
    return np.linalg.inv(np.eye(len(W)) - W + np.outer(np.ones(len(W)), p))


def norm_inf(M):
    return np.abs(M).sum(1).max()


# ===================================================================== 8.1
print("=" * 72); print("Exercise 8.1: DeGroot by hand")
W1 = np.array([[.5, .3, .2], [.1, .6, .3], [.3, .2, .5]])
b0 = np.array([1.0, 0.0, 0.5])
b = b0.copy()
for t in range(1, 4):
    print(f"   products W_ij b_j({t-1}):", np.round(W1 * b[None, :], 4).tolist())
    b = W1 @ b
    print(f"b({t}) = {np.round(b, 4)}   spread max-min = {b.max() - b.min():.4f}")
p1 = perron(W1)
print("col sums", W1.sum(0), " pi =", np.round(p1, 4), " exact pi = (14, 19, 17)/50 =", np.array([14, 19, 17]) / 50)
assert np.allclose(p1, np.array([14, 19, 17]) / 50)
print("check first column: pi_j W_j1 =", np.round(p1 * W1[:, 0], 4), " sum =", round(p1 @ W1[:, 0], 4))
bstar = p1 @ b0
print("b* = pi'b0 =", round(bstar, 6), " = 0.28 + 0.5*0.34 =", 0.28 + 0.17)
b = b0.copy()
for t in range(1, 31):
    b = W1 @ b
    if t in (5, 10, 20, 30):
        print(f"  b({t}) = {np.round(b, 6)}   max|b - b*| = {np.abs(b - bstar).max():.2e}")
lam = eigs(W1)
print("eigenvalues:", np.round(lam, 4), " |lambda_2| =", round(abs(lam[1]), 4), " Gamma =", round(gap(W1), 4),
      " half-life =", round(halflife(abs(lam[1])), 4))
tr, det = np.trace(W1), np.linalg.det(W1)
print("trace", round(tr, 4), "det", round(det, 4), " => l2 + l3 =", round(tr - 1, 4), ", l2*l3 =", round(det, 4))

# ===================================================================== 8.2
print("=" * 72); print("Exercise 8.2: effective sample size, N = 100")
N = 100
for p1_ in (0.01, 0.5, 0.9):
    p = np.r_[p1_, np.full(N - 1, (1 - p1_) / (N - 1))]
    print(f"pi_1 = {p1_}: N_eff = {neff(p):.4f}")
# closed form threshold: p^2 + (1-p)^2/99 = 1/10  <=>  100 p^2 - 2 p - 8.9 = 0
pstar = (2 + np.sqrt(4 + 4 * 100 * 8.9)) / 200
print("N_eff = 10 at pi_1 =", round(pstar, 5), " check", round(1 / (pstar ** 2 + (1 - pstar) ** 2 / 99), 6))
assert np.isclose(1 / (pstar ** 2 + (1 - pstar) ** 2 / 99), 10)
print("N_eff at pi_1 = 0 :", round(1 / (1 / 99), 2), "(max N_eff = 100 at pi_1 = 0.01)")
print("TikZ coords (pi_1, N_eff):", " ".join("(%.3f,%.2f)" % (q, 1 / (q ** 2 + (1 - q) ** 2 / 99))
                                           for q in np.r_[0, 0.01, 0.02, 0.03, 0.05, 0.075, np.arange(0.1, 1.0001, 0.05)]))

# ===================================================================== 8.3
print("=" * 72); print("Exercise 8.3: manipulation cost on the network of 8.1")
Z = fundamental(W1)
print("Z =\n", np.round(Z, 4), "\n||Z||_inf =", round(norm_inf(Z), 4), " Z 1 =", np.round(Z @ np.ones(3), 6))
bT = 0.8
db = bT - bstar
print("b* =", round(bstar, 4), " target", bT, " required d b =", round(db, 4))
# minimal-norm shift in pi: dpi = db (b0 - mean b0) / ||b0 - mean b0||^2
c = b0 - b0.mean()
dpi = db * c / (c @ c)
print("min-Euclidean d pi =", np.round(dpi, 4), " new pi =", np.round(p1 + dpi, 4), " check pi'b0 =", round((p1 + dpi) @ b0, 4))
assert np.all(p1 + dpi >= 0) and np.isclose((p1 + dpi) @ b0, bT)
print("reachable range of targets: [", b0.min(), ",", b0.max(), "]")
# gradient of b* wrt an elementary move E = e_j (e_i - e_k)'
z = Z @ (b0 - bstar)
print("z = Z (b0 - b* 1) =", np.round(z, 4))
moves = []
for j in range(3):
    for i in range(3):
        for k in range(3):
            if i != k:
                moves.append((p1[j] * (z[i] - z[k]), j, i, k))
moves.sort(reverse=True)
for g, j, i, k in moves[:4]:
    print(f"  row j={j+1}: move weight from k={k+1} to i={i+1}: db*/dt = {g:.4f}  (capacity W_jk = {W1[j, k]})")
g, j, i, k = moves[0]
t_lin = db / g
print("first-order t for best edge:", round(t_lin, 4), " capacity", W1[j, k])
G1 = np.outer(p1, z); G1t = G1 - G1.mean(1, keepdims=True)
print("first-order min ||dW||_F (no sign constraints) = |db| / ||P(G)||_F =", round(abs(db) / np.linalg.norm(G1t), 4))
# exact along best edge: b*(t)
def bstar_move(W, j, i, k, t):
    W2 = W.copy(); W2[j, i] += t; W2[j, k] -= t
    return perron(W2) @ b0
ts = np.linspace(0, W1[j, k], 2001)
vals = [bstar_move(W1, j, i, k, t) for t in ts]
print(f"exact: moving all of W_{j+1}{k+1} = {W1[j,k]} to W_{j+1}{i+1} gives b* = {vals[-1]:.4f}")
# combined moves: in every row, move all weight from agent 2 to agent 1
W_all = W1.copy(); W_all[:, 0] += W_all[:, 1]; W_all[:, 1] = 0
print("every row shifts all weight from agent 2 to agent 1 (W_.2 -> 0): pi =", np.round(perron(W_all), 4),
      " b* =", round(perron(W_all) @ b0, 4), " ||dW||_F =", round(np.linalg.norm(W_all - W1), 4))
# min-Frobenius steepest-ascent path with nonnegativity (projected gradient)
W = W1.copy(); path = 0.0; h = 1e-3; steps = 0
while perron(W) @ b0 < bT and steps < 100000:
    p = perron(W); Zt = fundamental(W); zt = Zt @ (b0 - p @ b0)
    G = np.outer(p, zt)                         # d b*/d W_jk = pi_j z_k
    G = G - G.mean(1, keepdims=True)            # keep rows summing to one
    G[(W <= 1e-12) & (G < 0)] = 0               # cannot push zero entries negative
    G = G - np.where(W > 1e-12, 1, 0) * (G.sum(1, keepdims=True) / np.maximum((W > 1e-12).sum(1, keepdims=True), 1))
    step = h * G / np.linalg.norm(G)
    Wn = np.clip(W + step, 0, None); Wn = Wn / Wn.sum(1, keepdims=True)
    path += np.linalg.norm(Wn - W); W = Wn; steps += 1
print(f"projected steepest-ascent path to b* = {perron(W) @ b0:.4f}: steps {steps}, path length {path:.4f}, "
      f"endpoint ||dW||_F = {np.linalg.norm(W - W1):.4f}")
print("   endpoint W =\n", np.round(W, 4), "\n   pi =", np.round(perron(W), 4))
W_steep = W.copy()
# first-order lower bound of Prop ch08:thm:manipulation
bound = abs(db) / (np.abs(b0 - bstar).max() * norm_inf(Z))
print("first-order bound: t||E||_inf >=", round(bound, 4), " (max|b0-b*| =", round(np.abs(b0 - bstar).max(), 4), ")")
# 1/|1-lambda_k|
print("max_k 1/|1-lambda_k| =", round(max(1 / abs(1 - l) for l in eigs(W1)[1:]), 4))
# lazier networks W_alpha
P1 = np.outer(np.ones(3), p1)
print("W_alpha = alpha W + (1-alpha) I:")
for a in (1.0, 0.5, 0.25, 0.1):
    Wa = a * W1 + (1 - a) * np.eye(3)
    Za = fundamental(Wa)
    assert np.allclose(perron(Wa), p1)
    assert np.allclose(Za, P1 + (Z - P1) / a)
    ba = abs(db) / (np.abs(b0 - bstar).max() * norm_inf(Za))
    print(f"  alpha={a:4}: Gamma={gap(Wa):.4f} ||Z_a||_inf={norm_inf(Za):.4f} bound={ba:.4f}  "
          f"best-edge db*/dt={max(p1[j_]*((Za@(b0-bstar))[i_]-(Za@(b0-bstar))[k_]) for j_ in range(3) for i_ in range(3) for k_ in range(3) if i_!=k_):.4f}")

# ===================================================================== 8.4
print("=" * 72); print("Exercise 8.4: the algorithm as W^I designer, N = 8")
b8 = np.r_[np.full(4, .8), np.full(4, .2)]
Wc = np.ones((8, 8)) / 8
print("(a) complete: b* =", round(perron(Wc) @ b8, 4), " Gamma =", round(gap(Wc), 4))
Wb = np.zeros((8, 8)); Wb[:4, :4] = Wb[4:, 4:] = 0.25
lim = np.linalg.matrix_power(Wb, 50) @ b8
print("(b) block-diagonal: limit beliefs =", np.round(lim, 4), " eigenvalues", np.round(np.sort(eigs(Wb).real)[::-1], 4),
      " Gamma =", round(gap(Wb), 4))
print("(c)/(d)")
for beta in (0, 1, 5, 20):
    K = np.exp(-beta * np.abs(b8[:, None] - b8[None, :]))
    W = K / K.sum(1, keepdims=True)
    q = np.exp(-0.6 * beta)
    assert np.allclose(W, W.T) and np.allclose(W.sum(0), 1)
    l2 = abs(eigs(W)[1])
    assert np.isclose(l2, np.tanh(0.3 * beta)) and np.isclose(l2, (1 - q) / (1 + q))
    print(f"  beta={beta:2d}: within weight {W[0,0]:.6f}, across {W[0,4]:.3e}; b* = {perron(W) @ b8:.4f}; "
          f"|l2| = tanh(0.3 beta) = {l2:.8f}; Gamma = {gap(W):.3e}; half-life = {halflife(l2):.4g} rounds")
for beta in (5, 10, 20):
    print(f"  asymptotic half-life (ln2/2) e^(0.6 beta) at beta={beta}: {np.log(2)/2*np.exp(0.6*beta):.4g}")
print("TikZ coords (beta, Gamma):", " ".join("(%.1f,%.4f)" % (bb, 1 - np.tanh(0.3 * bb)) for bb in np.arange(0, 20.01, 1)))

# ===================================================================== 8.5
print("=" * 72); print("Exercise 8.5: Bolivia's epistemic structure (illustrative)")
cities = ["LP", "SC", "CB"]
blocks = {"LP": np.array([[.8, .0, .2], [.1, .6, .3], [.5, .2, .3]]),
          "SC": np.array([[.6, .1, .3], [.0, .8, .2], [.15, .55, .3]]),
          "CB": np.array([[.7, .1, .2], [.1, .7, .2], [.35, .35, .3]])}
for c in cities:
    B = blocks[c]; assert np.allclose(B.sum(1), 1)
    q = perron(B)
    print(f"  {c} block alone: pi (State, Private, Citizens) = {np.round(q, 4)}, N_eff = {neff(q):.3f}, Gamma = {gap(B):.4f}")
eta = 0.1
n = 9
Wn = np.zeros((n, n))
for a, c in enumerate(cities):
    Wn[3 * a:3 * a + 3, 3 * a:3 * a + 3] = (1 - eta) * blocks[c]
    for b_, c2 in enumerate(cities):
        if b_ != a:
            for r in range(3):
                Wn[3 * a + r, 3 * b_ + r] = eta / 2      # same type in the other two cities
assert np.allclose(Wn.sum(1), 1)
pn = perron(Wn)
labels = [f"{t}_{c}" for c in cities for t in ("S", "P", "C")]
print("  national pi:", dict(zip(labels, np.round(pn, 4))))
print("  city masses:", np.round(pn.reshape(3, 3).sum(1), 4), " type masses (S, P, C):", np.round(pn.reshape(3, 3).sum(0), 4))
print("  N_eff =", round(neff(pn), 3), " of 9;  Gamma =", round(gap(Wn), 4), " eigenvalues |.|:", np.round(np.abs(eigs(Wn)), 4))
PiS, PiP, PiC = pn.reshape(3, 3).sum(0)
# consensus as a linear functional
bS, bP = 0.9, 0.2
for bc in (0.5, 0.3, 0.7):
    bb = np.tile([bS, bP, bc], 3)
    print(f"  b0: state {bS}, private {bP}, citizens {bc}: b* = {pn @ bb:.4f}")
# which citizen prior makes b* = 0.5 with state 0.9 and private 0.2
c_half = (0.5 - bS * PiS - bP * PiP) / PiC
print("  with state 0.9, private 0.2: b* > 0.5 iff citizens' prior >", round(c_half, 4))
# (d) cost from 0.6 to 0.4: choose citizens' prior so that b* = 0.6
c06 = (0.6 - bS * PiS - bP * PiP) / PiC
b5 = np.tile([bS, bP, c06], 3)
print("  calibration: citizens' prior =", round(c06, 4), " gives b* =", round(pn @ b5, 6))
assert 0 <= c06 <= 1 and np.isclose(pn @ b5, 0.6)
Z5 = fundamental(Wn)
z5 = Z5 @ (b5 - 0.6)
G = np.outer(pn, z5); Gt = G - G.mean(1, keepdims=True)
print("  first-order min ||dW||_F = |db| / ||P(G)||_F =", round(0.2 / np.linalg.norm(Gt), 4),
      " Prop. bound t||E||_inf >=", round(0.2 / (np.abs(b5 - 0.6).max() * norm_inf(Z5)), 4),
      " ||Z||_inf =", round(norm_inf(Z5), 4))
# projected steepest-descent path
W = Wn.copy(); path = 0.0; h = 1e-3; steps = 0
while perron(W) @ b5 > 0.4 and steps < 200000:
    p = perron(W); zt = fundamental(W) @ (b5 - p @ b5)
    G = -np.outer(p, zt)
    G = G - G.mean(1, keepdims=True)
    G[(W <= 1e-12) & (G < 0)] = 0
    G = G - np.where(W > 1e-12, 1, 0) * (G.sum(1, keepdims=True) / np.maximum((W > 1e-12).sum(1, keepdims=True), 1))
    Wn2 = np.clip(W + h * G / np.linalg.norm(G), 0, None); Wn2 /= Wn2.sum(1, keepdims=True)
    path += np.linalg.norm(Wn2 - W); W = Wn2; steps += 1
print(f"  steepest-descent path to b* = {perron(W) @ b5:.4f}: path length {path:.4f}, endpoint ||dW||_F = {np.linalg.norm(W - Wn):.4f}")
dW = W - Wn
rows = np.linalg.norm(dW, axis=1)
print("  row norms of dW:", dict(zip(labels, np.round(rows, 4))))
print("  new pi:", dict(zip(labels, np.round(perron(W), 4))))
# a simple targeted policy: citizens of La Paz move weight x from state to private outlet in LP
for x in (0.1, 0.2, 0.3):
    W2 = Wn.copy(); W2[2, 0] -= x; W2[2, 1] += x
    print(f"  targeted: LP citizens move {x} from S_LP to P_LP: b* = {perron(W2) @ b5:.4f}, ||dW||_F = {np.sqrt(2)*x:.4f}")
# all citizens move y from state to private outlet of their city
for y in (0.1, 0.15, 0.2):
    W2 = Wn.copy()
    for a in range(3):
        W2[3 * a + 2, 3 * a] -= y; W2[3 * a + 2, 3 * a + 1] += y
    ok = np.all(W2 >= 0)
    print(f"  all citizens move {y} from S to P: b* = {perron(W2) @ b5:.4f}, ||dW||_F = {np.sqrt(6)*y:.4f}, feasible {ok}")
# outlets: state outlets move weight z from self to private outlet of the city
for zz in (0.1, 0.2):
    W2 = Wn.copy()
    for a in range(3):
        W2[3 * a, 3 * a] -= zz; W2[3 * a, 3 * a + 1] += zz
    print(f"  state outlets move {zz} from self to the city's private outlet: b* = {perron(W2) @ b5:.4f}, ||dW||_F = {np.sqrt(6)*zz:.4f}")

# ===================================================================== 8.6
print("=" * 72); print("Exercise 8.6: the belief-action-signal loop")
rM, rT = 0.15, 0.05
print("(a) block beliefs: r_hat_M =", rM, " r_hat_T =", rT)
rhat_new = 0.5 * rT + 0.5 * rM
print("(b) r_hat' =", rhat_new, " investment a = theta r_hat rises by theta *", rhat_new - rT, " i.e. +", 100 * (rhat_new / rT - 1), "%")
nM, nT, f = 50, 100, 0.1
nt = int(f * nT)
W = np.zeros((nM + nT, nM + nT))
W[:nM, :nM] = 1 / nM
W[nM:, nM:] = 1 / nT
for r in range(nM, nM + nt):                     # treated traditional agents
    W[r, nM:] = 0.5 / nT; W[r, :nM] = 0.5 / nM
assert np.allclose(W.sum(1), 1)
r0 = np.r_[np.full(nM, rM), np.full(nT, rT)]
rh = W @ r0
print("(c) one-step average traditional belief =", round(rh[nM:].mean(), 4), " = 0.9*0.05 + 0.1*0.10 =", 0.9 * 0.05 + 0.1 * 0.10)
print("    treated belief", round(rh[nM], 4), " untreated", round(rh[-1], 4))
lam = eigs(W)
print("    |lambda_2| =", round(abs(lam[1]), 4), " = 1 - f/2 =", 1 - f / 2, " Gamma =", round(gap(W), 4),
      " (before: Gamma = 0, two closed classes)")
assert np.isclose(abs(lam[1]), 1 - f / 2)
b = r0.copy()
for t in range(1, 101):
    b = W @ b
    if t in (1, 14, 50, 100):
        print(f"    iterated DeGroot t={t}: mean T belief = {b[nM:].mean():.5f}")
print("    half-life of the T-M belief gap =", round(np.log(2) / np.log(1 / (1 - f / 2)), 2), "rounds")
# (d) loop with linear returns g_T(a) = r0 + gamma a, theta = 1
th, gam = 1.0, 0.4
r0T = rT * (1 - gam * th)
print("(d) calibration theta=1, gamma=0.4, r0 =", round(r0T, 4), " => block fixed point r_T* =", round(r0T / (1 - gam * th), 4),
      " stable since gamma*theta =", gam * th, "< 1")
rTs = (r0T + gam * th * (f / 2) * rM) / (1 - gam * th * (1 - f / 2))
print("    with smartphone share f=0.1 (sector-average signal): r_T* =", round(rTs, 5), " slope", gam * th * (1 - f / 2))
# simulate
rt = rT
for t in range(60):
    rt = r0T + gam * th * ((1 - f / 2) * rt + (f / 2) * rM)
print("    simulated limit", round(rt, 5))
assert np.isclose(rt, rTs)

# ===================================================================== 8.7
print("=" * 72); print("Exercise 8.7: Quijano in matrix form (illustrative)")
NE, NI = 2, 8
N7 = NE + NI


def coloniality(eps_EI, eps_IE):
    """E rows: (1-eps_EI) on E (equal), eps_EI on I (equal).
    I rows: eps_IE on E (equal); (1-eps_IE) on I with own share 1/2 of the I block, rest equal."""
    W = np.zeros((N7, N7))
    W[:NE, :NE] = (1 - eps_EI) / NE
    W[:NE, NE:] = eps_EI / NI
    for r in range(NE, N7):
        W[r, :NE] = eps_IE / NE
        W[r, NE:] = (1 - eps_IE) * 0.5 / (NI - 1)
        W[r, r] = (1 - eps_IE) * 0.5
    return W


W7 = coloniality(0.1, 0.6)
assert np.allclose(W7.sum(1), 1)
p7 = perron(W7)
print("W rows: E", np.round(W7[0], 4), "\n        I", np.round(W7[2], 4))
print("pi_E (each) =", round(p7[0], 4), " pi_I (each) =", round(p7[2], 4), " E mass =", round(p7[:NE].sum(), 4),
      " = eps_IE/(eps_EI+eps_IE) =", round(0.6 / 0.7, 4), " N_eff =", round(neff(p7), 4), " Gamma =", round(gap(W7), 4))
bE, bI = 0.9, 0.3
b7 = np.r_[np.full(NE, bE), np.full(NI, bI)]
print("    I mass =", round(p7[NE:].sum(), 4), "; SC citizens weight on state media in the composite =", round(0.9 * 0.15, 4))
print("(b) consensus b* =", round(p7 @ b7, 4), " (population mean of priors =", round(b7.mean(), 4), ")")
print("    consensus error variance / sigma^2 = 1/N_eff =", round(1 / neff(p7), 4), " vs 1/N =", 1 / N7)
# (c) formation illustration: theta* = theta_bar + kappa (signal - r*), Example ch09 values kappa = 1.5
kap9, thbar, rstar = 0.3 / 0.2, 1.0, 0.5
print("(c) theta* for the valuation of formal schooling: colonial consensus", round(thbar + kap9 * (p7 @ b7 - rstar), 4),
      " equal-voice consensus", round(thbar + kap9 * (b7.mean() - rstar), 4))
# (d) decoloniality: equal per-capita power requires E mass = NE/N = 0.2, i.e. eps_EI = 4 eps_IE
best = None
for eIE in np.linspace(0.001, 0.25, 2500):
    Wd = coloniality(4 * eIE, eIE)
    cost = np.linalg.norm(Wd - W7)
    if best is None or cost < best[0]:
        best = (cost, eIE)
cost, eIE = best
Wd = coloniality(4 * eIE, eIE)
pd_ = perron(Wd)
print(f"(d) min ||dW||_F with uniform pi (two-parameter family): eps_IE = {eIE:.4f}, eps_EI = {4*eIE:.4f}, "
      f"||dW||_F = {cost:.4f}, pi = {np.round(pd_, 4)}")
assert np.allclose(pd_, 1 / N7, atol=1e-6)
print("    new consensus b* =", round(pd_ @ b7, 4), " N_eff =", round(neff(pd_), 4), " Gamma =", round(gap(Wd), 4))
for eEI, eIE_ in ((0.1, 0.025), (0.8, 0.6), (1.0, 0.25), (0.8, 0.2)):
    Wx = coloniality(eEI, eIE_)
    px = perron(Wx)
    print(f"    eps_EI={eEI}, eps_IE={eIE_}: E mass {px[:NE].sum():.4f}, b* {px @ b7:.4f}, ||dW||_F {np.linalg.norm(Wx - W7):.4f}")
print("    E listening alone cannot equalize: needs eps_EI = 4*0.6 = 2.4 > 1")
print("all assertions passed")
