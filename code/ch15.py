"""Chapter 15: reproduces every number in the chapter's worked examples and
appendix, and checks the identities used in the text.
Run from the repo root:  python3 code/ch15.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def eigs(W):
    lam, V = np.linalg.eig(W)
    o = np.argsort(-np.abs(lam))
    lam, V = lam[o], V[:, o]
    U = np.linalg.inv(V)          # rows: left eigenvectors, biorthogonal (U V = I)
    return lam, V, U


def gap(W):
    ev = np.linalg.eigvals(W)
    return 1 - np.sort(np.abs(ev))[-2]


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def var_pi(p):
    return np.mean((p - 1 / len(p)) ** 2)


# ---------------------------------------------------------------------
print("== Proposition: the pi-preserving path to Gamma = 1 ==")
rng = np.random.default_rng(3)
W = rng.uniform(0.05, 1, (4, 4)) ** 2
W /= W.sum(1)[:, None]
p = perron(W)
lam2 = np.sort(np.abs(np.linalg.eigvals(W)))[-2]
for t in [0, 0.25, 0.5, 0.75]:
    Wt = (1 - t) * W + t * np.outer(np.ones(4), p)
    assert np.allclose(perron(Wt), p)
    assert abs(gap(Wt) - (1 - (1 - t) * lam2)) < 1e-9
    print(f"t = {t:.2f}: Gamma = {gap(Wt):.3f}, Var(pi) = {var_pi(perron(Wt)):.4f}")
print("pi held fixed at", p, "; t = 1 gives W = 1 pi', Gamma = 1")

# ---------------------------------------------------------------------
print("\n== Appendix 15A: complex-safe derivative of |lambda_2| ==")
for s in range(500):
    W = rng.uniform(0, 1, (4, 4)) ** 3
    W /= W.sum(1)[:, None]
    lam, V, U = eigs(W)
    if abs(np.imag(lam[1])) > 1e-3:
        break
i, j, l = 0, 1, 3                      # move weight from (i,l) to (i,j)
E = np.zeros((4, 4)); E[i, j] = 1; E[i, l] = -1
h = 1e-7
fd = (gap(W + h * E) - gap(W)) / h
dl = U[1] @ E @ V[:, 1]
correct = -np.real(np.conj(lam[1]) * dl) / abs(lam[1])
wrong = -np.sign(np.real(lam[1])) * np.real(dl)
print(f"lambda_2 = {lam[1]:.3f}: finite difference = {fd:.4f}, Re(conj(l) dl)/|l| formula = {correct:.4f}, "
      f"sign(lambda) formula = {wrong:.4f}")
assert abs(fd - correct) < 1e-4 and abs(fd - wrong) > 1e-3

# ---------------------------------------------------------------------
print("\n== Edge selection: Frobenius vs l1 budget ==")
rng2 = np.random.default_rng(3)
W = rng2.uniform(0.1, 1, (4, 4)); W /= W.sum(1)[:, None]
for _ in range(1):  # replicate the review's draw order
    pass
lam, V, U = eigs(W)
G = np.zeros((4, 4))
for a in range(4):
    for b in range(4):
        G[a, b] = -np.real(np.conj(lam[1]) * U[1, a] * V[b, 1]) / abs(lam[1])
Gc = G - G.mean(1)[:, None]            # projection onto {E : E 1 = 0}
B = 0.05
frob = gap(W + B * Gc / np.linalg.norm(Gc)) - gap(W)
best_edge = max(gap(W + (B / np.sqrt(2)) * (np.eye(4)[a][:, None] * (np.eye(4)[b] - np.eye(4)[c]))) - gap(W)
                for a in range(4) for b in range(4) for c in range(4) if b != c)
print(f"same Frobenius norm B = {B}: projected-gradient step dGamma = {frob:.3f}, best single move = {best_edge:.3f}")
assert frob > best_edge

# ---------------------------------------------------------------------
print("\n== Intuition: connect the cores, not the bridge ==")
A6 = np.array([
    [0.6, 0.38, 0.01, 0.004, 0.003, 0.003],
    [0.38, 0.6, 0.01, 0.004, 0.003, 0.003],
    [0.25, 0.25, 0.0, 0.17, 0.17, 0.16],
    [0.003, 0.003, 0.03, 0.6, 0.182, 0.182],
    [0.003, 0.003, 0.03, 0.182, 0.6, 0.182],
    [0.003, 0.003, 0.03, 0.182, 0.182, 0.6]])
W6 = A6 / A6.sum(1)[:, None]
lam, V, U = eigs(W6)
u2, v2 = np.real(U[1]), np.real(V[:, 1])
sc = 1 / np.abs(u2).max()
print(f"lambda_2 = {np.real(lam[1]):.3f}, Gamma = {gap(W6):.3f}")
print("u2 (scaled to max |.| = 1):", u2 * sc)
print("v2 (scaled to max |.| = 1):", v2 / np.abs(v2).max())
Gd = np.full((6, 6), np.nan)
for a in range(6):
    for b in range(6):
        if a == 2:
            continue          # bridge row has a zero diagonal: uniform compensation infeasible
        Ed = np.zeros((6, 6)); Ed[a, b] = 1; Ed[a, [c for c in range(6) if c != b]] -= 1 / 5
        Gd[a, b] = (gap(W6 + h * Ed) - gap(W6)) / h
k = np.nanargmax(Gd); a, b = divmod(k, 6)
print(f"argmax directional derivative: edge ({a+1},{b+1}) with dGamma/dt = {Gd[a,b]:.3f}; "
      f"best edge into the bridge node: {np.nanmax(Gd[:, 2]):.3f}")
assert a in (0, 1) and b in (3, 4, 5)

# ---------------------------------------------------------------------
print("\n== Example: ISI that works vs ISI that doesn't ==")
A = np.array([[0.33, 0.25, 0.01, 0.01], [0.25, 0.33, 0.01, 0.01],
              [0.012, 0.012, 0.15, 0.126], [0.012, 0.012, 0.126, 0.15]])
assert np.all(A.sum(1) < 1)
WP = A / A.sum(1)[:, None]


def total_output(A, d=np.ones(4)):
    return np.ones(4) @ np.linalg.inv(np.eye(4) - A.T) @ d     # y = A'y + d


print("A =\n", A)
print("W^P = D^-1 A =\n", np.round(WP, 3))
p0 = perron(WP)
print(f"base: Gamma = {gap(WP):.3f}, pi = {p0}, Var(pi) = {var_pi(p0):.4f}, total output 1'L d = {total_output(A):.3f}")
rows = {}
for name, jj in [("ISI-A (within): a_M1,M2 += 0.05", 1), ("ISI-B (cross):  a_M1,T1 += 0.05", 2)]:
    A2 = A.copy(); A2[0, jj] += 0.05
    W2 = A2 / A2.sum(1)[:, None]
    p2 = perron(W2)
    rows[name] = (gap(W2) - gap(WP), var_pi(p2) - var_pi(p0), total_output(A2) - total_output(A))
    print(f"{name}: dGamma = {rows[name][0]:+.3f}, dVar(pi) = {rows[name][1]:+.4f}, d total output = {rows[name][2]:+.3f}, pi = {p2}")
rA, rB = rows.values()
assert rB[0] > rA[0] and rB[1] < rA[1] and rA[2] > rB[2]

print("\n== Row-stochastic matrices kill row-sum linkages ==")
g = 0.5
L = np.linalg.inv(np.eye(4) - g * WP)
print("row sums of (I - g W)^-1:", L.sum(1), "= 1/(1-g) =", 1 / (1 - g))
assert np.allclose(L.sum(1), 1 / (1 - g))
LA = np.linalg.inv(np.eye(4) - A.T)
print("from A: backward linkages (column sums of (I-A')^-1):", LA.sum(0), " forward (row sums):", LA.sum(1))

# ---------------------------------------------------------------------
print("\n== Example: smartphone operator (illustrative) ==")
thT0, thM, th1, thinf = 0.92, 1.09, 0.96, 1.01
print(f"gap closed after one generation: {(th1-thT0)/(thM-thT0):.0%}; long run: {(thinf-thT0)/(thM-thT0):.0%}")
print(f"long-run Delta theta_T = {thinf - thT0:.2f}")

# ---------------------------------------------------------------------
print("\n== Proof checks: Props. optimal_edge and isi_gamma ==")
rng3 = np.random.default_rng(15)
Gr = rng3.normal(size=(5, 5))
Bb = 0.3
best_l1 = Bb / 2 * max(Gr[i, j] - Gr[i, l] for i in range(5) for j in range(5) for l in range(5))
Gt = Gr - Gr.mean(1)[:, None]
best_fr = Bb * np.linalg.norm(Gt)
for _ in range(5000):
    E = rng3.normal(size=(5, 5)) * (rng3.uniform(size=(5, 5)) < 0.5)
    E -= E.mean(1)[:, None]                          # E 1 = 0
    E1 = Bb * E / np.abs(E).sum()                    # ||E||_1 = B
    EF = Bb * E / np.linalg.norm(E)                  # ||E||_F = B
    assert np.sum(Gr * E1) <= best_l1 + 1e-12
    assert np.sum(Gr * EF) <= best_fr + 1e-12
print(f"l1 bound (B/2) max(G_ij - G_il) = {best_l1:.4f} and Frobenius bound B||G~||_F = {best_fr:.4f} never exceeded")

# domestic-Gamma direction  E^ = diag(K1)^-1 (delta - diag(delta 1) K^)  vs finite difference
K = A.copy()
Kh = K / K.sum(1)[:, None]
lam, V, U = eigs(Kh)
Gh = -np.real(np.conj(lam[1]) * np.outer(U[1], V[:, 1])) / abs(lam[1])
for jj in (1, 2):                                    # ISI-A and ISI-B directions
    dl = np.zeros((4, 4)); dl[0, jj] = 1.0
    Eh = np.diag(1 / K.sum(1)) @ (dl - np.diag(dl.sum(1)) @ Kh)
    assert np.allclose(Eh.sum(1), 0)
    s_ = 1e-7
    K2 = K + s_ * dl
    fd = (gap(K2 / K2.sum(1)[:, None]) - gap(Kh)) / s_
    assert abs(fd - np.sum(Gh * Eh)) < 1e-4
    print(f"ISI direction a_M1,{jj+1}: first-order dGamma = {np.sum(Gh*Eh):+.4f}, finite difference = {fd:+.4f}")
