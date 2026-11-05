"""Solutions manual, Chapter 6 (Centrality, Power, and the Perron Vector).
Computes and prints every number used in solutions/sol06.tex.
Run from the repository root:  python3 code/solutions/sol06.py
numpy only, except Exercise 6.7(e), which reads data/processed with pandas."""
import pathlib
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)
ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "processed"


def perron(W):
    """Left Perron vector of an irreducible row-stochastic W (sum 1)."""
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmin(np.abs(vals - 1))
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def eigs_by_modulus(W):
    lam = np.linalg.eigvals(W)
    return lam[np.argsort(-np.abs(lam))]


def gap(W):
    return 1 - np.abs(eigs_by_modulus(W))[1]


def neff(p):
    return 1 / np.sum(p ** 2)


def gini(p):
    p = np.sort(p); N = len(p)
    return np.sum((2 * np.arange(1, N + 1) - N - 1) * p) / (N * p.sum())


def rownorm(A):
    return A / A.sum(1, keepdims=True)


def lazy(W):
    return 0.5 * (np.eye(len(W)) + W)


def rank(v):
    """Agents (1-based) from highest to lowest."""
    return (np.argsort(-v, kind="stable") + 1).tolist()


# ===================================================================== 6.1
print("=" * 72); print("Exercise 6.1: five-agent economy (eq. ch06:eq:W5)")
W5 = np.array([[.4, .4, .1, .1, 0],
               [.3, .3, .2, .2, 0],
               [.1, .2, .3, .2, .2],
               [0, .1, .4, .3, .2],
               [0, .1, .1, .3, .5]])
I5 = np.eye(5); one5 = np.ones(5)
assert np.allclose(W5.sum(1), 1)
indeg = W5.sum(0)
outdeg = W5.sum(1)
katz = (np.linalg.inv(I5 - 0.5 * W5.T) - I5) @ one5
bon = np.linalg.inv(I5 - 0.3 * W5.T) @ W5.T @ one5
pi5 = perron(W5)
print("in-degree        ", np.round(indeg, 2), " out-degree", outdeg)
print("Katz (delta=0.5) ", np.round(katz, 4), " sum", round(katz.sum(), 4))
print("Bonacich (b=0.3) ", np.round(bon, 4), " sum", round(bon.sum(), 4))
print("pi               ", np.round(pi5, 4), " sum", pi5.sum())
print("N_eff =", round(neff(pi5), 3), " Gamma =", round(gap(W5), 4))
# sanity: sum of Katz = N (1/(1-d) - 1) since W'1 summed = N
assert np.isclose(katz.sum(), 5 * (1 / (1 - 0.5) - 1))
assert np.isclose(bon.sum(), 5 / (1 - 0.3))
assert np.isclose(pi5.sum(), 1)
print("rankings: Katz", rank(katz), " Bonacich", rank(bon), " pi", rank(pi5))
assert rank(katz) == rank(bon) == rank(pi5) == [3, 4, 2, 5, 1]
# Bonacich at its own beta = Katz(beta)/beta
katz03 = (np.linalg.inv(I5 - 0.3 * W5.T) - I5) @ one5
assert np.allclose(bon, katz03 / 0.3)
print("Katz(0.3) =", np.round(katz03, 4), " = 0.3 * Bonacich(0.3)")
# second-order check of why agent 2 falls: weighted inflows pi_j W_j2
print("pi_j W_j2 terms:", np.round(pi5 * W5[:, 1], 4), " pi_j W_j3 terms:", np.round(pi5 * W5[:, 2], 4))
# normalized Katz path
for d in (0.5, 0.9, 0.99, 0.999):
    kz = np.linalg.inv(I5 - d * W5.T) @ one5
    print(f"  normalized (1'(I-dW)^-1) at delta={d}:", np.round(kz / kz.sum(), 4))
# (c)
W5c = W5.copy(); W5c[2, 0] += 0.1; W5c[2, 2] -= 0.1
assert np.allclose(W5c.sum(1), 1)
pic = perron(W5c)
print("(c) new row 3:", W5c[2], " new pi:", np.round(pic, 4))
print("    d pi =", np.round(pic - pi5, 4), " d pi_1 =", round(pic[0] - pi5[0], 4),
      " relative", round(100 * (pic[0] / pi5[0] - 1), 1), "%")
print("    new ranking", rank(pic), " N_eff", round(neff(pic), 3))
assert pic[0] > pi5[0]
# first-order prediction from Lemma ch06:lem:perturbation: d pi' = pi' E Z, E = e3 (e1 - e3)'
Z = np.linalg.inv(I5 - W5 + np.outer(one5, pi5))
E = np.zeros((5, 5)); E[2, 0] = 1; E[2, 2] = -1
lin = 0.1 * pi5 @ E @ Z
print("    first-order (Lemma) prediction of d pi:", np.round(lin, 4))
# Kac check: pi_1 = 1/m_1 before and after
for name, W in (("before", W5), ("after", W5c)):
    Q = W[1:, 1:]; h = np.linalg.solve(np.eye(4) - Q, np.ones(4)); m1 = 1 + W[0, 1:] @ h
    print(f"    Kac {name}: hitting times of 1 from 2..5 = {np.round(h, 3)}, m_1 = {m1:.4f}, 1/m_1 = {1/m1:.4f}")
    assert np.isclose(1 / m1, perron(W)[0])

# ===================================================================== 6.2
print("=" * 72); print("Exercise 6.2: doubly stochastic or not")
mats = {"a": np.ones((3, 3)) / 3,
        "b": np.array([[.5, .5], [.8, .2]]),
        "c": np.array([[.5, .5, 0], [0, .5, .5], [.5, 0, .5]]),
        "d": np.array([[.4, .3, .3], [.3, .4, .3], [.3, .3, .4]])}
for k, W in mats.items():
    lam = eigs_by_modulus(W)
    print(f"({k}) col sums {W.sum(0)}, DS={np.allclose(W.sum(0), 1)}, pi={np.round(perron(W), 4)}, "
          f"eigs={np.round(lam, 4)}, |eigs|={np.round(np.abs(lam), 4)}, Gamma={gap(W):.4f}")
assert np.allclose(perron(mats["b"]), [8 / 13, 5 / 13])
assert np.isclose(gap(mats["a"]), 1) and np.isclose(gap(mats["c"]), 0.5) and np.isclose(gap(mats["d"]), 0.9)
assert np.isclose(gap(mats["b"]), 0.7)
# (c): powers of the cyclic matrix
print("(c) W^2 =\n", np.linalg.matrix_power(mats["c"], 2), "\n    W^3 =\n", np.linalg.matrix_power(mats["c"], 3))

# ===================================================================== 6.3
print("=" * 72); print("Exercise 6.3: spectral gap and convergence, N = 8")
N = 8
S = np.roll(np.eye(N), 1, 1)
St = np.zeros((N, N)); St[1:, 0] = 1; St[0, 1:] = 1 / (N - 1)
eps = 0.01; h = N // 2
WB = np.full((N, N), eps / h); WB[:h, :h] = WB[h:, h:] = (1 - eps) / h
topo = {"complete": np.ones((N, N)) / N,
        "star (lazy)": lazy(St),
        "ring (lazy)": lazy(0.5 * (S + S.T)),
        "weak blocks": WB}
b0 = np.zeros(N); b0[0] = 1
b0s = np.zeros(N); b0s[1] = 1          # a spoke / non-hub start, used for the star
dist = {}
runs = [(k, W, b0) for k, W in topo.items()] + [("star (lazy), start at spoke", topo["star (lazy)"], b0s)]
for k, W, bb0 in runs:
    assert np.allclose(W.sum(1), 1)
    p = perron(W); lam = eigs_by_modulus(W)
    l2 = abs(lam[1]); l2 = 0.0 if l2 < 1e-12 else l2
    bstar = p @ bb0
    b = bb0.copy(); d = []
    for t in range(1, 51):
        b = W @ b; d.append(np.linalg.norm(b - bstar))
    d = np.array(d); d[d < 1e-13] = 0.0; dist[k] = d
    print(f"{k:12s} Gamma={gap(W):.4f} |l2|={l2:.4f} N_eff={neff(p):.3f} b*={bstar:.4f} "
          f"pi_1={p[0]:.4f} half-life={np.log(2)/-np.log(l2) if 0 < l2 < 1 else 0:.3f}")
    print("      eigenvalues:", np.round(np.sort(lam.real)[::-1], 4))
    print("      ||b(t)-b*1|| t=0,1,2,5,10,20,50:", "%.4f" % np.linalg.norm(bb0 - bstar),
          ["%.3e" % d[t - 1] for t in (1, 2, 5, 10, 20, 50)])
    if np.all(d[10:] > 0):
        ratios = d[1:] / d[:-1]
        print("      ratio d(t+1)/d(t) at t=1,5,10,20,49:", np.round(ratios[[0, 4, 9, 19, 48]], 4))
        tt = np.arange(11, 51); sl = np.polyfit(tt, np.log(d[10:]), 1)[0]
        print("      fitted rate exp(slope) over t=11..50:", round(np.exp(sl), 4))
        assert abs(np.exp(sl) - l2) < 1e-3
    else:
        t0 = int(np.argmax(d == 0)) + 1
        print("      distance is exactly 0 from t =", t0)
        if t0 > 20:
            tt = np.arange(11, 31); sl = np.polyfit(tt, np.log(d[10:30]), 1)[0]
            print("      fitted rate exp(slope) over t=11..30:", round(np.exp(sl), 4))
            assert abs(np.exp(sl) - l2) < 1e-3
gg = {k: gap(W) for k, W in topo.items()}
assert np.isclose(gg["complete"], 1) and np.isclose(gg["star (lazy)"], 0.5)
assert np.isclose(gg["ring (lazy)"], 1 - (0.5 + 0.5 * np.cos(2 * np.pi / 8)))
assert np.isclose(gg["weak blocks"], 0.02)
# star: closed form of the distance from consensus starting at the hub
# without self-loops
ring0 = 0.5 * (S + S.T)
for k, W in (("ring (no loops)", ring0), ("star (no loops)", St)):
    lam = eigs_by_modulus(W)
    print(f"{k:16s} Gamma={gap(W):.4f}; eigenvalues of modulus 1:", np.round(lam[np.isclose(np.abs(lam), 1)], 4),
          " pi=", np.round(perron(W), 4))
    b = b0.copy(); traj = []
    for t in range(1, 7):
        b = W @ b; traj.append(np.round(b, 4))
    print("      b(1), b(2) =", traj[0], traj[1], "; b(5), b(6) first entries", traj[4][0], traj[5][0])
    assert np.isclose(gap(W), 0)
# TikZ coordinates (log10 distance) for the plot in the manual
print("TikZ coords log10 distance (t, log10 d) every 2 steps:")
for k in ("star (lazy), start at spoke", "ring (lazy)", "weak blocks"):
    d = dist[k]
    print("  ", k, " ".join("(%d,%.3f)" % (t, np.log10(d[t - 1])) for t in range(1, 51, 1) if d[t - 1] > 0 and t % 2 == 1))
print("  complete: d(1) =", dist["complete"][0])

# ===================================================================== 6.4
print("=" * 72); print("Exercise 6.4: Bonacich with negative beta")
Wbs = np.array([[0, 0, .5, .5], [0, 0, .5, .5], [.5, .5, 0, 0], [.5, .5, 0, 0]])
I4 = np.eye(4)
for beta in (0.3, -0.3):
    c = np.linalg.inv(I4 - beta * Wbs.T) @ Wbs.T @ np.ones(4)
    print(f"beta={beta:+.1f}: c^Bon =", np.round(c, 4), " 1/(1-beta) =", round(1 / (1 - beta), 4))
    assert np.allclose(c, 1 / (1 - beta))
print("pi =", perron(Wbs), " eigenvalues", np.round(eigs_by_modulus(Wbs), 4), " Gamma =", round(gap(Wbs), 4))
print("lazy version Gamma =", round(gap(lazy(Wbs)), 4), " pi =", perron(lazy(Wbs)))
# asymmetric variant: buyer 1 -- seller 3 only; buyer 2 -- sellers 3,4; seller 4 -- buyer 2 only
Wv = np.array([[0, 0, 1, 0], [0, 0, .5, .5], [.5, .5, 0, 0], [0, 1, 0, 0]])
print("variant W row sums", Wv.sum(1), " col sums (in-degree)", Wv.sum(0))
for beta in (0.3, -0.3):
    c = np.linalg.inv(I4 - beta * Wv.T) @ Wv.T @ np.ones(4)
    print(f"  variant beta={beta:+.1f}: c^Bon =", np.round(c, 4))
cpos = np.linalg.inv(I4 - 0.3 * Wv.T) @ Wv.T @ np.ones(4)
cneg = np.linalg.inv(I4 + 0.3 * Wv.T) @ Wv.T @ np.ones(4)
print("  variant pi =", np.round(perron(lazy(Wv)), 4), "(lazy, since Wv is bipartite)", " Gamma(Wv) =", round(gap(Wv), 4))
print("  variant ratio c_2/c_1: beta=+0.3", round(cpos[1] / cpos[0], 4), " beta=-0.3", round(cneg[1] / cneg[0], 4))
assert cneg[1] / cneg[0] > cpos[1] / cpos[0]

# ===================================================================== 6.5
print("=" * 72); print("Exercise 6.5: power iteration by hand")
W65 = np.array([[.6, .4], [.2, .8]])
v = np.array([1., 1.])
pi65 = np.array([1 / 3, 2 / 3])
assert np.allclose(perron(W65), pi65)
errs = []
for t in range(1, 9):
    vt = W65.T @ v
    v = vt / vt.sum()
    e = np.abs(v - pi65).max(); errs.append(e)
    if t == 5: print("   v(5) = (%.5f, %.5f)" % tuple(v))
    print(f"t={t}: W'v(t-1) = {np.round(vt, 5)}  v(t) = {np.round(v, 5)}  max|v-pi| = {e:.5f}  rounded {np.round(v, 2)}")
errs = np.array(errs)
print("error ratios:", np.round(errs[1:] / errs[:-1], 6))
assert np.allclose(errs[1:] / errs[:-1], 0.4)
lam = np.sort(np.linalg.eigvals(W65).real)
print("eigenvalues", lam, " Gamma =", round(gap(W65), 4), " t with error < 0.005:", int(np.argmax(errs < 0.005)) + 1)

# ===================================================================== 6.6
print("=" * 72); print("Exercise 6.6: President, General, Businesswoman, Farmer")
lab = ["President", "General", "Businesswoman", "Farmer"]
W66 = np.array([[.30, .20, .20, .30],     # President: army, business taxes, food
                [.30, .30, .10, .30],     # General: orders from President, own corps, food rations
                [.10, .10, .40, .40],     # Businesswoman: licences, security, food processing inputs
                [.05, .05, .20, .70]])    # Farmer: own land and family; sells to the businesswoman
assert np.allclose(W66.sum(1), 1)
p66 = perron(W66)
print("col sums (in-degree):", W66.sum(0))
print("pi =", dict(zip(lab, np.round(p66, 4))), " N_eff =", round(neff(p66), 3), " Gamma =", round(gap(W66), 4))
assert np.argmax(p66) == 3
# check the balance equation for the Farmer
print("Farmer balance: sum_j pi_j W_j4 =", round(p66 @ W66[:, 3], 6), " pi_4 =", round(p66[3], 6))
# What it takes: the Farmer leads as long as others' food dependence is large enough.
for f in (0.30, 0.20, 0.10):
    W = W66.copy()
    for i in range(3):
        diff = W[i, 3] - f
        W[i, 3] = f; W[i, i] += diff
    pp = perron(W)
    print(f"  food weight of P,G,B = {f:.2f}: pi = {np.round(pp, 4)}, top = {lab[int(np.argmax(pp))]}")

# ===================================================================== 6.7
print("=" * 72); print("Exercise 6.7: stylized dual economy (Example ch03:ex:bolivia)")
A = np.array([[.15, .03], [.02, .20]])
D = np.diag(A.sum(1)); WP = np.linalg.inv(D) @ A
p = perron(WP)
e1, e2 = WP[0, 1], WP[1, 0]
print("row sums of A:", A.sum(1), "\nW^P =\n", np.round(WP, 4))
print("eps1 (M->T) =", round(e1, 4), " eps2 (T->M) =", round(e2, 4))
print("pi =", np.round(p, 4), " closed form e2/(e1+e2) =", round(e2 / (e1 + e2), 4),
      " Gamma =", round(gap(WP), 4), " = e1+e2 =", round(e1 + e2, 4), " N_eff =", round(neff(p), 4))
assert np.isclose(p[0], e2 / (e1 + e2)) and np.isclose(gap(WP), e1 + e2)
# (b) transfer: pi does not depend on endowments
# (c) plant: a_TM 0.02 -> 0.15
A2 = A.copy(); A2[1, 0] = 0.15
WP2 = np.diag(1 / A2.sum(1)) @ A2; p2 = perron(WP2)
f1, f2 = WP2[0, 1], WP2[1, 0]
print("(c) W^P new =\n", np.round(WP2, 4))
print("    pi new =", np.round(p2, 4), " Gamma =", round(gap(WP2), 4), " |lambda_2| =", round(abs(eigs_by_modulus(WP2)[1]), 4),
      " N_eff =", round(neff(p2), 4))
print("    change in pi_T =", round(p2[1] - p[1], 4))
assert p2[1] < p[1]
# the reverse link: modern buys more from traditional (a_MT 0.03 -> 0.15)
A3 = A.copy(); A3[0, 1] = 0.15
p3 = perron(np.diag(1 / A3.sum(1)) @ A3)
print("    contrast: a_MT 0.03 -> 0.15 gives pi =", np.round(p3, 4))
# Leontief multipliers (computed from A, not W^P) for the plant
L = np.linalg.inv(np.eye(2) - A.T); L2 = np.linalg.inv(np.eye(2) - A2.T)
print("    Leontief L before =\n", np.round(L, 4), "\n    after =\n", np.round(L2, 4))
print("    column sums (output multipliers) before", np.round(L.sum(0), 4), " after", np.round(L2.sum(0), 4))

# (e) data: INE Bolivia 2019 (base 2017), domestic flows, book12
print("-" * 72); print("6.7(e) data: INE supply-use table 2019, domestic, book12 (data/processed/bolivia_Z)")
Zdf = pd.read_csv(DATA / "bolivia_Z" / "Z_2017_2019_book12.csv", index_col=0)
Zdf = Zdf.drop(index="Domestic services", columns="Domestic services")   # no intermediate flows
labels = list(Zdf.index)
Z = Zdf.to_numpy(float)


def wp_from_Z(Z):
    """W^P[s,s'] = Z[s',s] / sum_k Z[k,s]  (input-cost shares, rows = buyers)."""
    return Z.T / Z.sum(0)[:, None]


W12 = wp_from_Z(Z)
p12 = perron(W12)
summ = pd.read_csv(DATA / "bolivia_io_summary.csv")
ref = summ[(summ.base == 2017) & (summ.year == 2019) & (summ.aggregation == "book12") &
           (summ.flows == "domestic") & (summ.normalization == "input_shares")].iloc[0]
print(f"11-sector check: Gamma = {gap(W12):.4f} (summary {ref.Gamma:.4f}), N_eff = {neff(p12):.4f} (summary {ref.Neff:.4f})")
assert np.isclose(gap(W12), ref.Gamma, atol=1e-6) and np.isclose(neff(p12), ref.Neff, atol=1e-6)
for s, v in sorted(zip(labels, p12), key=lambda z: -z[1]):
    print(f"   pi {s:22s} {v:.4f}")
variants = {
    "narrow modern (Hydrocarbons, Mining, Financial)": ["Hydrocarbons", "Mining", "Financial services"],
    "modern + Manufacturing": ["Hydrocarbons", "Mining", "Financial services", "Manufacturing"],
    "broad modern (+ Public adm., Transport, Other services)":
        ["Hydrocarbons", "Mining", "Financial services", "Manufacturing", "Public administration",
         "Transport", "Other services"],
}
agg_res = {}
for name, modern in variants.items():
    m = np.array([s in modern for s in labels])
    M = np.c_[m, ~m].astype(float)
    Z2 = M.T @ Z @ M                    # aggregate flows: rows suppliers, cols buyers
    W2 = wp_from_Z(Z2)
    q = perron(W2)
    agg_res[name] = (W2, q)
    trad = [s for s in labels if s not in modern]
    print(f"{name}\n   traditional = {trad}")
    print("   Z2 (Bs mn, rows supplier, cols buyer) =", np.round(Z2, 0).tolist())
    print(f"   W^P = {np.round(W2, 4).tolist()}  eps_M={W2[0,1]:.4f} eps_T={W2[1,0]:.4f}")
    print(f"   pi = {np.round(q, 4)}  N_eff = {neff(q):.4f} of 2  Gamma = {gap(W2):.4f}")
    assert np.isclose(gap(W2), 1 - abs(1 - W2[0, 1] - W2[1, 0]))
    print(f"   lambda_2 = 1 - eps_M - eps_T = {1 - W2[0,1] - W2[1,0]:.4f}")
print("stylized: pi =", np.round(p, 4), " Gamma =", round(gap(WP), 4), " N_eff =", round(neff(p), 4), " lambda_2 =", round(1 - e1 - e2, 4))
print("11-sector lambda_2 =", np.round(eigs_by_modulus(W12)[1], 4), " N_eff/N =", round(neff(p12) / 11, 4))
for name, (W2, q) in agg_res.items():
    print(f"   N_eff/N for {name}: {neff(q) / 2:.4f}")
print("all assertions passed")
