"""Chapter 24: viz-panel diagnostics and the counter-hegemony model.
Run from the repo root: python3 code/ch24.py  (numpy only)."""
import numpy as np


def stats(K):
    W = K / K.sum(1, keepdims=True)          # row-normalize, W = D^{-1} K
    ev, vl = np.linalg.eig(W.T)
    order = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, order[0]]); pi /= pi.sum()
    return 1 - abs(ev[order[1]]), pi, 1 / np.sum(pi ** 2)


def heg(i, j):      # heatmap cell values (blue!xx) of the left panel
    if i == 1 and j == 1: return 80
    if i == 1: return 50
    if j == 1: return 40
    if i == j: return 30
    return 5


def cnt(i, j):      # right panel
    if i == 1 and j == 1: return 45
    if i == 1: return 20
    if j == 1: return 15
    if i == j: return 35
    if i >= 4 and j >= 4: return 30
    return 12


print("== Viz panel: row-normalized heatmaps ==")
for name, f in (("hegemonic", heg), ("counter-hegemonic", cnt)):
    K = np.array([[f(i, j) for j in range(1, 7)] for i in range(1, 7)], float)
    G, pi, ne = stats(K)
    print(f"{name:18s} Gamma = {G:.2f}  Neff = {ne:.2f}  pi_dom = {pi[0]:.3f}  pi = {np.round(pi,3)}")
    assert np.isclose(pi.sum(), 1)

print("\n== Counter-hegemony model: theta' = rho*theta + (1-rho)*(g(theta) + s) ==")
rho, beta = 0.85, 10.0
g = lambda th: 1 / (1 + np.exp(-beta * (th - 0.5)))
dg = lambda th: beta * g(th) * (1 - g(th))


def fixed_points(s):
    grid = np.linspace(-0.5, 1.5, 200000)
    f = g(grid) + s - grid
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    return grid[idx]


for s in (0.0, 0.10, 0.18, 0.19):
    fps = fixed_points(s)
    slopes = [rho + (1 - rho) * dg(t) for t in fps]
    print(f"s = {s:.2f}: fixed points {np.round(fps,3)}, map slopes {np.round(slopes,3)}")
fps0 = fixed_points(0.0)
assert len(fps0) == 3 and abs(fps0[1] - 0.5) < 1e-3
# saddle-node: g'(theta) = 1 and theta = g(theta) + s
gstar = (1 - np.sqrt(1 - 4 / beta)) / 2
th_star = 0.5 + np.log(gstar / (1 - gstar)) / beta
s_star = th_star - gstar
print(f"saddle-node at s* = {s_star:.3f} (theta = {th_star:.3f})")
assert len(fixed_points(s_star - 0.002)) == 3 and len(fixed_points(s_star + 0.002)) == 1
# one-shot impulse at s = 0: crossing the separatrix 0.5 decides the outcome
for jump in (0.45, 0.55):
    th = fps0[0] + jump
    for _ in range(400):
        th = rho * th + (1 - rho) * g(th)
    print(f"impulse to theta = {fps0[0]+jump:.3f}: converges to {th:.3f}")
