"""Solutions manual, Chapter 1: every number used in solutions/sol01.tex.
Run from the repository root:  python3 code/solutions/sol01.py
numpy throughout; pandas only to read data/processed (Exercise 1.7)."""
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)


def hdr(s):
    print("\n== " + s)


def U(x, k):
    x = np.asarray(x, float)
    return x - k / 2 * x ** 2


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True)
    return 1 - e[1]


def gini(x):
    x = np.asarray(x, float)
    N = len(x)
    return np.abs(x[:, None] - x[None, :]).sum() / (2 * N ** 2 * x.mean())


# ---------------------------------------------------------------- 1.1
hdr("1.1 Reading the allocation matrix")
X = np.array([[5, 3, 1], [2, 4, 2], [1, 1, 5], [4, 2, 2.]])
print("agent 3 bundle:", X[2])
print("column sums 1'X:", X.sum(0), " housing total:", X[:, 2].sum(),
      " argmax housing (1-based):", X[:, 2].argmax() + 1)
assert np.allclose(X.sum(0), [12, 10, 10])
r = np.linalg.matrix_rank(X)
sv = np.linalg.svd(X, compute_uv=False)
print("rank X =", r, " singular values:", sv.round(4))
print("share of ||X||_F^2 in best rank-1 approx: %.4f" % (sv[0] ** 2 / np.sum(sv ** 2)))
print("expenditure-share rows (each row / row sum):\n", (X / X.sum(1, keepdims=True)).round(4))
M3 = X[[0, 1, 2]]
print("det of rows 1-3:", round(np.linalg.det(M3), 4))
assert r == 3

# ---------------------------------------------------------------- 1.2
hdr("1.2 Cobb-Douglas U = x1^(1/3) x2^(2/3)")
Om = np.array([[4, 1], [1, 3], [3, 2.]])
tot = Om.sum(0)
# market 1: sum_i w_i/3 = 8 with w_i = omega_i1 + p2 omega_i2  ->  (8 + 6 p2)/3 = 8
p2 = (3 * tot[0] - tot[0]) / tot[1]
p = np.array([1, p2])
w = Om @ p
Xs = np.column_stack([w / 3, 2 * w / (3 * p2)])
print("totals", tot, " p2* =", p2, " (8/3 =", 8 / 3, ")")
print("wealth w =", w.round(4))
print("X* =\n", Xs.round(4), "\ncolumn sums", Xs.sum(0))
assert np.isclose(p2, 8 / 3) and np.allclose(Xs.sum(0), tot)
# Walras: value of excess demand is zero at any p
for q in [0.5, 1.0, 4.0]:
    ww = Om @ np.array([1, q]); z = np.array([ww.sum() / 3 - tot[0], 2 * ww.sum() / (3 * q) - tot[1]])
    print("p2=%.1f excess demands" % q, z.round(4), " p'z =", round(z @ [1, q], 10))
    assert abs(z @ [1, q]) < 1e-12
print("expenditure shares p'x: good1 share =", (Xs[:, 0] / w).round(4), " good2 share =", (p2 * Xs[:, 1] / w).round(4))
print("rank X* =", np.linalg.matrix_rank(Xs))

# ---------------------------------------------------------------- 1.3
hdr("1.3 Welfare loss L = (kappa N/2) Var(omega), N=5, kappa=0.1, Omega=5")
k, N = 0.1, 5
cases = {"a": [1, 1, 1, 1, 1], "b": [3, 1, .5, .25, .25], "c": [5, 0, 0, 0, 0]}
for key, om in cases.items():
    om = np.array(om, float)
    V = np.mean(om ** 2) - np.mean(om) ** 2
    L = k * N / 2 * V
    Ldirect = N * U(1, k) - U(om, k).sum()
    print("(%s) omega=%s  sum=%.2f  max kappa*omega=%.2f  mean(om^2)=%.4f  Var=%.4f  L=%.5f  direct=%.5f  Gini=%.4f"
          % (key, om, om.sum(), (k * om).max(), np.mean(om ** 2), V, L, Ldirect, gini(om)))
    assert np.isclose(L, Ldirect) and (k * om).max() < 1
print("W_FB =", N * U(1, k), " share lost in (c) = %.4f" % (1.0 / (N * U(1, k))))
# (d) L as a function of the Gini: one-parameter family omega(t) = (1-t) 1 + t * 5 e1
for t in [0, .25, .5, .75, 1]:
    om = (1 - t) * np.ones(5) + t * np.array([5, 0, 0, 0, 0.])
    V = np.var(om)
    print("family t=%.2f  Gini=%.3f  Var=%.4f  L=%.4f  1.5625*G^2=%.4f" % (t, gini(om), V, k * N / 2 * V, 1.5625 * gini(om) ** 2))
    assert np.isclose(k * N / 2 * V, 1.5625 * gini(om) ** 2)
# two vectors with the same Gini and different L
o1 = np.array([2, 2, 1, 0, 0.]); o2 = np.array([3.4, 0.4, 0.4, 0.4, 0.4])
print("same-Gini pair: G1=%.4f L1=%.4f   G2=%.4f L2=%.4f" % (gini(o1), k * N / 2 * np.var(o1), gini(o2), k * N / 2 * np.var(o2)))
assert np.isclose(gini(o1), gini(o2))
# random scatter: range of L in Gini bins, and the Cauchy-Schwarz lower bound L >= kappa N G^2 (mean 1)
rng = np.random.default_rng(1)
S = rng.dirichlet(np.ones(5) * 0.6, size=200000) * 5
G = np.array([gini(s) for s in S[:20000]])
Ls = k * N / 2 * S[:20000].var(1)
print("Lower bound check (L >= kappa*N*G^2 = 0.5 G^2):", bool(np.all(Ls >= 0.5 * G ** 2 - 1e-12)))
assert np.all(Ls >= 0.5 * G ** 2 - 1e-12)
for lo in [0.2, 0.4, 0.6]:
    sel = (G >= lo - 0.01) & (G < lo + 0.01)
    print("Gini in [%.2f,%.2f): n=%d  L min=%.4f  L max=%.4f  bound 0.5G^2=%.4f  family 1.5625G^2=%.4f"
          % (lo - .01, lo + .01, sel.sum(), Ls[sel].min(), Ls[sel].max(), 0.5 * lo ** 2, 1.5625 * lo ** 2))

# ---------------------------------------------------------------- 1.4
hdr("1.4 Bargaining, N=2, kappa=0.5, Omega=2, d=(0.6,0.2)")
k = 0.5; Om2 = 2.0; d = np.array([.6, .2]); Phi = Om2 - d.sum()
print("surplus Phi =", Phi)
xfb = np.array([1., 1.]); WFB = U(xfb, k).sum()
print("(a) FB x =", xfb, " U =", U(xfb, k), " W =", WFB)
for lam in [np.array([.5, .5]), np.array([.8, .2])]:
    x = d + lam * Phi
    Wv = U(x, k).sum()
    print("lam", lam, " x =", x.round(4), " kappa x =", (k * x).round(3), " U =", U(x, k).round(4),
          " W =", round(Wv, 4), " L vs FB =", round(WFB - Wv, 4), " (k/2)sum z^2 =", round(k / 2 * np.sum((x - 1) ** 2), 4))
    assert np.isclose(WFB - Wv, k / 2 * np.sum((x - 1) ** 2)) and np.all(k * x < 1)
print("U at d:", U(d, k).round(4), " W(d) =", U(d, k).sum().round(4))
xe, xa = d + .5 * Phi, d + np.array([.8, .2]) * Phi
print("(c) squared deviations from 1:", ((xa - 1) ** 2).round(4))
print("(c) utility change vs equal power:", (U(xa, k) - U(xe, k)).round(4))
lamstar = (1 - d) / Phi
print("(d) lambda* = (1-d)/Phi =", lamstar.round(4))
assert np.allclose(d + lamstar * Phi, 1)

# ---------------------------------------------------------------- 1.5
hdr("1.5 Cascade, A = [[.2,.1],[.3,.1]]")
A = np.array([[.2, .1], [.3, .1]]); At = A.T
P = [np.linalg.matrix_power(At, j) for j in range(6)]
for j in range(5):
    print("(A')^%d =\n" % j, P[j].round(5))
S4 = sum(P[:5]); Lx = np.linalg.inv(np.eye(2) - At)
print("S4 =\n", S4.round(5), "\nL exact =\n", Lx.round(5))
err = Lx - S4
print("L - S4 =\n", err.round(5), "\n relative error per entry:\n", (err / Lx).round(5))
print("remainder identity (A')^5 L =\n", (P[5] @ Lx).round(5))
assert np.allclose(err, P[5] @ Lx)
print("det(I-A') =", round(np.linalg.det(np.eye(2) - At), 4))
ev = np.linalg.eigvals(A); print("eigenvalues of A:", np.sort(ev)[::-1].round(5), " rho =", round(max(abs(ev)), 5))
print("trace, det:", A.trace(), round(np.linalg.det(A), 4), " sqrt(0.13) =", round(np.sqrt(0.13), 5), " lambda_2 = %.5f" % ((0.3 - np.sqrt(0.13)) / 2))
print("row sums of A:", A.sum(1), " max row sum bound on rho:", A.sum(1).max())
print("multipliers (col sums of L):", Lx.sum(0).round(4), " forward (row sums):", Lx.sum(1).round(4), " mean:", Lx.sum(0).mean().round(4))
assert max(abs(ev)) < 1

# ---------------------------------------------------------------- 1.6
hdr("1.6 W^CE and the 3x3 W")
Wce = np.ones((3, 3)) / 3
assert np.allclose(Wce @ Wce, Wce)
print("null-space basis of W^CE (N=3): (1,-1,0), (1,0,-1); check:", Wce @ [1, -1, 0], Wce @ [1, 0, -1])
W = np.array([[1, 1, 1], [2, .5, .5], [0, 1.5, 1.5]]) / 3
print("W =\n", W.round(4), "\nrow sums", W.sum(1), " col sums", W.sum(0))
xv = np.array([3, 6, 9.])
print("W x =", W @ xv, " W^CE x =", Wce @ xv, " sums:", (W @ xv).sum(), (Wce @ xv).sum())
print("rank W =", np.linalg.matrix_rank(W), " eigenvalues:", np.round(np.linalg.eigvals(W), 6))
print("pi(W) =", perron(W).round(4), " Gamma(W) =", round(gap(W), 6))
W2 = W @ W
print("W^2 =\n", W2.round(6), " equals W^CE:", np.allclose(W2, Wce))
print("(W - W^CE) =\n", (W - Wce).round(4), "\n(W - W^CE)^2 = 0:", np.allclose((W - Wce) @ (W - Wce), 0))
print("W^2 x =", W2 @ xv)
assert np.allclose(W2, Wce) and np.allclose(W.sum(0), 1)

# ---------------------------------------------------------------- 1.7
hdr("1.7 Bolivia input-output table (INE SUT 2019, base 2017, domestic flows, basic prices)")
Z = pd.read_csv("data/processed/bolivia_Z/Z_2017_2019_book12.csv", index_col=0)
f = pd.read_csv("data/processed/bolivia_io_sector_facts.csv")
f = f[(f.base == 2017) & (f.year == 2019)].set_index("book12")
lab = list(Z.index)
x = f.output.reindex(lab).values
agg = {"Hydrocarbons": "Hydrocarbons", "Mining": "Mining", "Agriculture": "Agriculture",
       "Manufacturing": "Manufacturing", "Construction": "Construction"}
groups = ["Hydrocarbons", "Mining", "Agriculture", "Manufacturing", "Construction", "Services"]
Gm = np.zeros((6, len(lab)))
for j, l in enumerate(lab):
    Gm[groups.index(agg.get(l, "Services")), j] = 1
Zg = Gm @ Z.values @ Gm.T          # rows = sellers, columns = buyers (millions of Bs)
xg = Gm @ x
ex = Gm @ f.exports.reindex(lab).values
print("gross output (million Bs):", dict(zip(groups, xg.round(0))))
print("intermediate sales matrix Z (rows sell to columns), million Bs:\n", Zg.round(0))
Ab = Zg.T / xg[:, None]            # book convention a_ij = purchases by i from j / x_i
print("A (book convention, row i = recipe of i):\n", Ab.round(4))
print("row sums of A (domestic intermediate inputs per unit output):", Ab.sum(1).round(4))
rho = max(abs(np.linalg.eigvals(Ab))); print("rho(A) =", round(rho, 4))
Lb = np.linalg.inv(np.eye(6) - Ab.T)
print("L = (I - A')^-1:\n", Lb.round(4))
BL = Lb.sum(0); FL = Lb.sum(1)
print("output multipliers (col sums):", dict(zip(groups, BL.round(3))))
print("forward (row sums):", dict(zip(groups, FL.round(3))), " mean BL = mean FL:", BL.mean().round(4), FL.mean().round(4))
assert np.isclose(BL.mean(), FL.mean()) and rho < 1
own = np.diag(Lb) / BL
print("share of each multiplier staying in own sector:", dict(zip(groups, own.round(3))))
print("spillover to other sectors (BL - L_jj):", dict(zip(groups, (BL - np.diag(Lb)).round(3))))
print("export share of output:", dict(zip(groups, (ex / xg).round(3))))
# Mining/Gas merged: 5-sector version
G5 = ["Mining/Gas", "Agriculture", "Manufacturing", "Construction", "Services"]
M5 = np.zeros((5, 6)); M5[0, 0] = M5[0, 1] = 1
for i in range(2, 6): M5[i - 1, i] = 1
Z5 = M5 @ Zg @ M5.T; x5 = M5 @ xg
A5 = Z5.T / x5[:, None]; L5 = np.linalg.inv(np.eye(5) - A5.T)
print("5-sector multipliers:", dict(zip(G5, L5.sum(0).round(3))))
# block structure: extractive {Hydrocarbons, Mining} vs rest
blk = np.array([1, 1, 0, 0, 0, 0], bool)
cross = Ab[np.ix_(blk, ~blk)].sum() + Ab[np.ix_(~blk, blk)].sum()
chi = cross / Ab.sum()
print("cross-block mass chi(A), extractive vs rest = %.4f" % chi)
# extractive block purchases: share of its domestic inputs that come from the rest
print("share of extractive-block domestic inputs bought from rest: %.4f"
      % (Zg[np.ix_(~blk, blk)].sum() / Zg[:, blk].sum()))
print("share of rest's domestic inputs bought from extractive block: %.4f"
      % (Zg[np.ix_(blk, ~blk)].sum() / Zg[:, ~blk].sum()))
WP = Ab / Ab.sum(1, keepdims=True)
pi6 = perron(WP)
print("W^P (input shares):\n", WP.round(4))
print("pi(W^P) =", dict(zip(groups, pi6.round(4))), " Gamma(W^P) = %.4f" % gap(WP), " N_eff = %.3f of 6" % (1 / np.sum(pi6 ** 2)))
# hydrocarbons: unit final demand -> where does output go
e = np.zeros(6); e[0] = 1
print("L e_hydro =", (Lb @ e).round(4))
print("hydrocarbons imported share of intermediate inputs (sector facts): %.4f" % f.loc["Hydrocarbons", "imported_share_of_int_inputs"])
print("hydrocarbons compensation of employees / VA: %.4f" % f.loc["Hydrocarbons", "comp_share_of_va"])

# ---------------------------------------------------------------- 1.8
hdr("1.8 Dona Maria's dependency graph")
names = ["Maria", "Wholesaler", "SugarShop", "Cinnamon", "Charcoal"]
G = np.zeros((5, 5)); G[0, 1:] = 1
print("adjacency G (G_ij=1 if i uses j's output):\n", G)
print("symmetric?", np.allclose(G, G.T), " G^2 = 0:", np.allclose(G @ G, 0))
print("in-degree (col sums):", G.sum(0), " out-degree (row sums):", G.sum(1))
Lg = np.linalg.inv(np.eye(5) - 0.25 * G.T)
print("(I - 0.25 G')^-1 =\n", Lg.round(4), "  equals I + 0.25 G':", np.allclose(Lg, np.eye(5) + 0.25 * G.T))
G2 = G.copy(); G2[1, :] = 0; G2[:, 1] = 0
print("after removing the wholesaler, G =\n", G2)
d = np.array([1, 0, 0, 0, 0.])
print("quantity model output y = (I - 0.25 G')^-1 d, before:", (Lg @ d).round(3),
      " after:", (np.linalg.inv(np.eye(5) - 0.25 * G2.T) @ d).round(3))
print("W^CE (N=5) entries 1/5; G has %d nonzero entries of 25" % G.sum())

# ---------------------------------------------------------------- 1.9
hdr("1.9 Countries as sectors")
A9 = np.array([[.05, .30, .30], [.01, .05, .05], [.01, .05, .05]])
L9 = np.linalg.inv(np.eye(3) - A9.T)
print("L =\n", L9.round(4), "\nBL (col sums):", L9.sum(0).round(4), " FL (row sums):", L9.sum(1).round(4))
print("rho(A) =", round(max(abs(np.linalg.eigvals(A9))), 4))
print("(b) L e2 =", (L9 @ [0, 1, 0]).round(4), " total =", (L9 @ [0, 1, 0]).sum().round(4))
print("    L e1 =", (L9 @ [1, 0, 0]).round(4), " total =", (L9 @ [1, 0, 0]).sum().round(4))
Am = A9.T.copy()
Lm = np.linalg.inv(np.eye(3) - Am.T)
print("(c) modified A (= A') =\n", Am, "\nL_mod =\n", Lm.round(4), "\nBL_mod:", Lm.sum(0).round(4), " FL_mod:", Lm.sum(1).round(4))
print("    L_mod = L' :", np.allclose(Lm, L9.T))
WP9 = A9 / A9.sum(1, keepdims=True); WPm = Am / Am.sum(1, keepdims=True)
print("W^P original:\n", WP9.round(4), " col sums", WP9.sum(0).round(4), " pi", perron(WP9).round(4))
print("W^P modified:\n", WPm.round(4), " col sums", WPm.sum(0).round(4), " pi", perron(WPm).round(4))
# (d) symmetric A => L symmetric => BL = FL
Asym = (A9 + A9.T) / 2; Ls = np.linalg.inv(np.eye(3) - Asym.T)
print("(d) symmetric part S=(A+A')/2: BL", Ls.sum(0).round(4), " FL", Ls.sum(1).round(4))
assert np.allclose(Ls.sum(0), Ls.sum(1))
print("BL - FL original:", (L9.sum(0) - L9.sum(1)).round(4))
print("\nAll Chapter 1 solution checks passed.")
