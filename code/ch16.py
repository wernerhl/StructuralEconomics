"""Chapter 16: reproduces the numbers in the chapter's schematic panels and
checks the matrix facts used in the text.
Run from the repo root:  python3 code/ch16.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def gap(W):
    return 1 - np.sort(np.abs(np.linalg.eigvals(W)))[-2]


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


print("== Schematic 10x10 global W^P (core 1-5, semi-periphery 6-7, periphery 8-10) ==")
N = 10
H = np.zeros((N, N))
for i in range(1, N + 1):
    for j in range(1, N + 1):
        if i <= 5 and j <= 5: v = 65
        elif i <= 5 and 5 < j <= 7: v = 25
        elif i <= 5 and j > 7: v = 5
        elif 5 < i <= 7 and j <= 5: v = 45
        elif 5 < i <= 7 and 5 < j <= 7: v = 35
        elif 5 < i <= 7 and j > 7: v = 8
        elif i > 7 and j <= 5: v = 50
        elif i > 7 and 5 < j <= 7: v = 10
        elif i > 7 and j > 7 and i == j: v = 30
        else: v = 3
        H[i - 1, j - 1] = v
W = H / H.sum(1)[:, None]                  # W = D^-1 K
p = perron(W)
print(f"Gamma = {gap(W):.2f}; pi = {p}; Neff = {1/np.sum(p**2):.2f} (N = {N})")
print("column sums:", W.sum(0), " mean =", W.sum(0).mean())
assert abs(W.sum(0).mean() - 1) < 1e-12   # columns of a row-stochastic matrix average 1
assert gap(W) > 0.85

print("\n== Frictionless gravity: W = 1 s' ==")
s = np.array([0.18, 0.15, 0.08, 0.06, 0.05, 0.04, 0.03, 0.0005, 0.0004, 0.0003])
s_full = np.append(s, (1 - s.sum()) / 190 * np.ones(190))   # 200 countries
Wg = np.outer(np.ones(200), s_full)
print(f"Gamma(1 s') = {gap(Wg):.3f}, pi = s (first 3: {perron(Wg)[:3]}), Neff = {1/np.sum(s_full**2):.1f}")
print(f"pi_US / (1/N) = {0.18 * 200:.0f}")

print("\n== Country-similarity matrix M~ of the ECI ==")
rng = np.random.default_rng(5)
M = (rng.uniform(size=(30, 80)) < np.linspace(0.1, 0.7, 30)[:, None]).astype(float)
M = M[M.sum(1) > 0][:, M.sum(0) > 0]
Mt = (M / M.sum(1)[:, None]) @ (M / M.sum(0)).T
assert np.allclose(Mt.sum(1), 1)
ev = np.linalg.eigvals(Mt)
pM = perron(Mt)
assert np.allclose(pM, M.sum(1) / M.sum())       # left Perron vector = diversity shares
assert np.abs(ev.imag).max() < 1e-10 and ev.real.min() > -1e-10 and ev.real.max() < 1 + 1e-10
Dc, Dp = M.sum(1), M.sum(0)
S = np.diag(Dc ** -0.5) @ M @ np.diag(Dp ** -0.5)
sv = np.linalg.svd(S, compute_uv=False)
assert np.allclose(np.sort(sv ** 2)[::-1][:len(ev)], np.sort(ev.real)[::-1])
print("checked: M~ row-stochastic, eigenvalues real in [0,1] (= squared singular values of D_c^-1/2 M D_p^-1/2),")
print("         left Perron vector = diversity / total diversity")

print("\n== Linkage table (illustrative): internal consistency ==")
BL = {"United States": 2.45, "China": 2.68, "Germany": 2.51, "Japan": 2.33, "Bolivia": 1.28, "Ethiopia": 1.15}
FL = {"United States": 1.98, "China": 2.31, "Germany": 2.12, "Japan": 1.85, "Bolivia": 1.12, "Ethiopia": 1.08}
assert all(v >= 1 for v in BL.values()) and all(v >= 1 for v in FL.values())
print("all BL, FL >= 1 (Leontief-inverse sums include the identity)")
# mean column sum = mean row sum for any square matrix
A = rng.uniform(0, 0.2, (8, 8)); L = np.linalg.inv(np.eye(8) - A.T)
assert abs(L.sum(0).mean() - L.sum(1).mean()) < 1e-12 and L.sum(0).min() >= 1 and L.sum(1).min() >= 1
print("checked: mean backward linkage = mean forward linkage; all >= 1")
g = 0.5; Wr = A / A.sum(1)[:, None]
print("row-normalized W: row sums of (I - g W)^-1 =", np.linalg.inv(np.eye(8) - g * Wr).sum(1)[:3], "... all = 1/(1-g)")

print("\n== Schematic ECI-Gini scatter: least-squares line ==")
pts = np.array([(5.8, 0.9), (5.2, 1.0), (4.5, 1.3), (4.8, 1.9), (4.0, 1.1), (3.8, 1.6), (3.5, 2.0), (3.0, 2.6),
                (2.5, 2.3), (2.8, 3.4), (2.0, 2.9), (1.5, 2.6), (1.2, 1.9), (1.0, 2.2), (0.8, 3.0), (2.2, 1.7)])
b, a = np.polyfit(pts[:, 0], pts[:, 1], 1)
print(f"fitted line y = {a:.3f} + ({b:.3f}) x; endpoints (0.5, {a+0.5*b:.2f}), (6.3, {a+6.3*b:.2f})")
