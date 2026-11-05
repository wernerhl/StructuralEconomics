"""Chapter 17: reproduces every number in the chapter's (illustrative) Bolivia
examples, tables, heatmaps and dashboard, and checks the identities used.
Run from the repo root:  python3 code/ch17.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


# ---------------------------------------------------------------------
print("== Illustrative 12-sector Bolivia W^P ==")
sectors = ["Hydrocarbons", "Mining", "Public administration",
           "Manufacturing", "Construction", "Commerce", "Transport", "Financial services",
           "Agriculture", "Other services", "Personal services", "Domestic services"]
blocks = [[0, 1, 2], [3, 4, 5, 6, 7], [8, 9, 10, 11]]
att = np.array([3, 2, 5, 3, 1, 2, 1.5, 1, 2, 1, 0.7, 0.2])      # within-block supplier weights
x = np.array([8, 2, 5, 3, 1, 1.5, 1, 1, 1.5, 0.6, 0.4, 0.1])   # cross-block supplier weights
c = [0.05, 0.09, 0.17]                                          # cross-block share of each block's inputs
N = 12
W = np.zeros((N, N))
for b, bl in enumerate(blocks):
    others = [j for j in range(N) if j not in bl]
    for i in bl:
        W[i, bl] = att[bl] / att[bl].sum() * (1 - c[b])
        W[i, others] = x[others] / x[others].sum() * c[b]
assert np.allclose(W.sum(1), 1)
ev = eig_sorted(W)
p = perron(W)
print("lambda_1..3 =", np.round(np.real(ev[:3]), 3), " Gamma =", round(gap(W), 3))
print("pi =", dict(zip(sectors, np.round(p, 3))))
print(f"power Gini G(pi) = {gini(p):.3f}, Neff = {1/np.sum(p**2):.2f}, Var(pi) = {np.mean((p-1/N)**2):.4f}")
assert abs(gap(W) - 0.127) < 1e-3 and abs(gini(p) - 0.52) < 5e-3

# technical coefficients A = diag(a) W^P (row i: i's recipe), linkages from L = (I - A')^-1
a = np.array([0.25, 0.35, 0.20, 0.55, 0.50, 0.30, 0.45, 0.25, 0.35, 0.30, 0.15, 0.05])
A = np.diag(a) @ W
L = np.linalg.inv(np.eye(N) - A.T)
BL, FL = L.sum(0), L.sum(1)
assert BL.min() >= 1 and FL.min() >= 1 and abs(BL.mean() - FL.mean()) < 1e-12
print("sector: BL, FL")
for s_, bl_, fl_ in zip(sectors, BL, FL):
    print(f"  {s_:22s} {bl_:.2f} {fl_:.2f}")
print(f"mean BL = {BL.mean():.3f} = mean FL = {FL.mean():.3f}")

# alternative normalization: allocation shares over intermediate sales (rows = sellers)
Z = A * 1.0                          # flows per unit output with unit outputs: purchases of i from j
Bm = Z.T / Z.T.sum(1)[:, None]       # seller j's sales shares
print(f"Gamma under input-share normalization = {gap(W):.3f}; under sales-share normalization = {gap(Bm):.3f}")

# employment shares (illustrative) and power-employment ratios
emp = np.array([0.01, 0.04, 0.08, 0.10, 0.07, 0.15, 0.06, 0.02, 0.27, 0.08, 0.09, 0.03])
assert abs(emp.sum() - 1) < 1e-12
print("pi/emp:", dict(zip(sectors, np.round(p / emp, 2))))

# heatmap brightness (normalized to the matrix maximum, rows/cols in block order)
br = np.round(100 * W / W.max()).astype(int)
print("Bolivia heatmap brightness:", br.tolist())

# ---------------------------------------------------------------------
print("\n== Illustrative South Korea comparison matrix ==")
Kraw = np.array([[45 + ((i * 7 + j * 3) % 15) for j in range(1, 13)] for i in range(1, 13)], float)
K = Kraw / Kraw.sum(1)[:, None]
mu2 = abs(eig_sorted(K)[1])
# W_K = (1-gam) I + gam K: own-sector inputs on the diagonal plus uniform-ish cross inputs.
lo, hi = 0.0, 1.0                    # choose gam by bisection so that Gamma(W_K) = 0.45
for _ in range(100):
    gam = 0.5 * (lo + hi)
    if gap((1 - gam) * np.eye(N) + gam * K) < 0.45:
        lo = gam
    else:
        hi = gam
WK = (1 - gam) * np.eye(N) + gam * K
print(f"|mu_2(K)| = {mu2:.3f}; gamma = {gam:.3f}; Gamma(W_K) = {gap(WK):.3f}")
assert abs(gap(WK) - 0.45) < 1e-6
brK = np.round(100 * WK / WK.max()).astype(int)
print("Korea heatmap brightness:", brK.tolist())

# ---------------------------------------------------------------------
print("\n== Between-block variance (illustrative inputs) ==")
y = np.array([18500, 5200, 1800]); w = np.array([0.08, 0.35, 0.57])
m = w @ y; VB = w @ (y - m) ** 2
print(f"mean = {m:.0f}; V_B = {VB/1e6:.2f}M; sigma_B = V_B / 32M = {VB/32e6:.3f}; ratio y1/y3 = {y[0]/y[2]:.1f}")

print("\n== N_eff of four media ecosystems (illustrative) ==")
s = np.array([0.34, 0.37, 0.12, 0.17])
print(f"sum of squares = {np.sum(s**2):.4f}; Neff = {1/np.sum(s**2):.2f} (maximum with 4 nodes = 4); sum of shares = {s.sum():.2f}")
assert abs(1/np.sum(s**2) - 3.4) < 0.05
# audience shares = pi^I only if W^I = 1 s' (rank one), whose Gamma = 1
WI = np.outer(np.ones(4), s)
print(f"W^I = 1 s': pi = {perron(WI)}, Gamma = {1 - np.sort(np.abs(np.linalg.eigvals(WI)))[-2]:.3f}")

print("\n== Psi (illustrative) ==")
print(f"Psi = 0.028/0.018 = {0.028/0.018:.2f}; enforcement share 1/(1+Psi) = {1/(1+0.028/0.018):.2f}")

print("\n== Dashboard bars ==")
for name, G_, sB, ne, ps in [("Bolivia", gap(W), 0.62, 3.4, 1.56), ("Korea", 0.45, 0.20, 12, 0.6)]:
    print(f"{name}: Gamma = {G_:.2f}, 1-sigma_B = {1-sB:.2f}, 1-G(pi) = {1-(gini(p) if name=='Bolivia' else 0.28):.2f}, "
          f"Neff/20 = {ne/20:.2f}, 1/(1+Psi) = {1/(1+ps):.2f}")
print(f"ratios: Gamma {0.45/0.13:.1f}x, sigma_B {0.62/0.20:.1f}x, G {0.52/0.28:.2f}x, Neff {12/3.4:.2f}x, Psi {1.56/0.6:.1f}x")

# ---------------------------------------------------------------------
print("\n== Proof check: Prop. between (exact quadratic decomposition) and Exercise (refinement) ==")
rng_p = np.random.default_rng(17)
for _ in range(200):
    n = 60
    y = rng_p.lognormal(0, 0.6, n)
    kap = 0.9 / y.max()                               # kappa * y < 1
    U = lambda z: z - kap / 2 * z ** 2
    lab = rng_p.integers(0, 3, n)
    def VB_VW(lab):
        ybar = y.mean(); vb = vw = 0.0
        for b in np.unique(lab):
            yb = y[lab == b]; w = len(yb) / n
            vb += w * (yb.mean() - ybar) ** 2; vw += w * yb.var()
        return vb, vw
    vb, vw = VB_VW(lab)
    lhs = np.sum(U(y.mean()) - U(y))
    assert abs(lhs - kap / 2 * n * (vb + vw)) < 1e-9
    assert abs(vb + vw - y.var()) < 1e-12
    fine = lab * 2 + rng_p.integers(0, 2, n)          # refine each block into two
    assert VB_VW(fine)[0] >= vb - 1e-12               # refining never lowers V_B (hence sigma_B)
print("sum[U(ybar) - U(y_i)] = (kappa/2) N (V_B + V_W) and refinement monotonicity verified")
