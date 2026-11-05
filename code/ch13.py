"""Chapter 13: reproduces every number in the chapter's worked examples,
viz panels and the exercises that were revised.
Run from the repo root:  python3 code/ch13.py   (numpy only)."""
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


# ---------------------------------------------------------------------
print("== Example: two structural equilibria (reduced co-evolution model) ==")
N = 6
s = np.array([0.25, 0.25, 0.15, 0.15, 0.10, 0.10])   # partner attractiveness (gas, urban, rural)
r = np.array([0.20, 0.20, 0.12, 0.12, 0.06, 0.06])   # block returns
B = np.zeros((N, N))
for b in range(3):
    B[2 * b:2 * b + 2, 2 * b:2 * b + 2] = 0.5          # within-block averaging
K = np.outer(np.ones(N), s)                           # open (market) links 1 s'
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.10
M = mu / (1 - rho)                                    # formation multiplier
beta, c0, a, c_lo, c_hi, eta = 12.0, 0.30, 1.0, 0.05, 0.65, 0.05
Omega, Phi_b, eta_th = 6.0, 1.0, 1.0                  # allocation rule parameters


def Wc(c):
    return (1 - c) * B + c * K


def theta_star(c):
    # Ch10 formation steady state with W^I = W(c): theta* = th0 + M (W r - r*)
    return th0 + M * (Wc(c) @ r - rstar)


def h(c):
    th = theta_star(c)
    pi = perron(Wc(c))
    return Omega / N + (pi - 1 / N) * Phi_b + eta_th * (th - th.mean())


def F(c):
    z = beta * (c - c0) + beta * a * (theta_star(c)[4] - th0)
    return c_lo + (c_hi - c_lo) / (1 + np.exp(-z))


def T(c):  # mean-field evolution map on the family W(c)
    return (1 - eta) * c + eta * F(c)


# fixed points by bisection on F(c) - c
def bisect(f, lo, hi):
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if np.sign(f(mid)) == np.sign(f(lo)):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


g = lambda c: F(c) - c
grid = np.linspace(0, 1, 2001)
vals = np.array([g(c) for c in grid])
idx = np.where(np.diff(np.sign(vals)))[0]
roots = [bisect(g, grid[i], grid[i + 1]) for i in idx]
assert len(roots) == 3
cA, chat, cB = roots
dF = lambda c, e=1e-6: (F(c + e) - F(c - e)) / (2 * e)
for name, c in [("A (fragmented)", cA), ("threshold", chat), ("B (integrated)", cB)]:
    W = Wc(c)
    assert np.allclose(W.sum(1), 1) and W.min() > 0
    th = theta_star(c)
    x = h(c)
    stable = abs(1 - eta + eta * dF(c)) < 1
    print(f"{name}: c = {c:.3f}, F'(c) = {dF(c):.2f}, stable = {stable}")
    print(f"   Gamma = {gap(W):.3f}, pi = {perron(W)}, Neff = {1/np.sum(perron(W)**2):.2f}")
    print(f"   theta* = {th}, x* = {x}, max gap (xmax-xmin)/mean = {(x.max()-x.min())/x.mean():.3f}")
    assert abs(gap(W) - c) < 1e-9          # Gamma(W(c)) = c exactly
    assert np.allclose(perron(W), s)       # power does not depend on c
print("W*_A =\n", np.round(Wc(cA), 3))
print("W*_B =\n", np.round(Wc(cB), 3))
eps_star = (chat - cA) * np.linalg.norm(K - B)
print(f"chat - cA = {chat-cA:.3f}; ||K - B||_F = {np.linalg.norm(K-B):.3f};  push eps* = {eps_star:.3f}")

# heatmap intensities (BOOK.md: normalize to the matrix's own max)
for name, c in [("A", cA), ("B", cB)]:
    W = Wc(c)
    br = np.round(100 * W / W.max()).astype(int)
    print(f"heatmap {name} brightness rows:", br.tolist())

# phase-diagram coordinates: net force F(c) - c  (Delta Gamma = eta (F - c))
print("phase diagram (c, F(c)-c):")
pts = [(c, g(c)) for c in np.linspace(0, 0.8, 33)]
print(" ".join(f"({8*c:.2f},{10*v:.3f})" for c, v in pts))

# ---------------------------------------------------------------------
print("\n== Prop. integration: counterexample to the old 'iff' ==")
e = 0.02
def blk(w):
    return np.array([[1 - w - 2 * e, w, e, e], [w, 1 - w - 2 * e, e, e],
                     [e, e, 1 - w - 2 * e, w], [e, e, w, 1 - w - 2 * e]])
for w in [0.005, 0.02]:
    print(f"within-block weight {w}: Gamma = {gap(blk(w)):.3f}, cross-block weight per row = {2*e}")
assert gap(blk(0.02)) > gap(blk(0.005))

# Cheeger check for the lazy reversible block matrix
def conductance(W):
    p = perron(W); n = len(W); best = np.inf
    for m in range(1, 2 ** n - 1):
        S = [i for i in range(n) if m >> i & 1]
        pS = p[S].sum()
        if pS > 0.5 + 1e-12:
            continue
        Sc = [j for j in range(n) if j not in S]
        best = min(best, sum(p[i] * W[i, j] for i in S for j in Sc) / pS)
    return best
Wl = blk(0.3)
Ph = conductance(Wl)
print(f"lazy block matrix: conductance = {Ph:.3f}, Gamma = {gap(Wl):.3f}, Cheeger: {Ph**2/2:.4f} <= Gamma <= {2*Ph:.3f}")
assert Ph ** 2 / 2 <= gap(Wl) <= 2 * Ph + 1e-12

# ---------------------------------------------------------------------
print("\n== Exercise 13.1 matrix ==")
W1 = np.array([[0.6, 0.38, 0.01, 0.01], [0.38, 0.6, 0.01, 0.01], [0.01, 0.01, 0.6, 0.38], [0.01, 0.01, 0.38, 0.6]])
assert np.allclose(W1.sum(1), 1)
print(f"Gamma = {gap(W1):.3f}")

print("\n== Exercise 13.4: ring relatedness ==")
phi = np.zeros((6, 6))
for p in range(6):
    for q in range(6):
        d = min((p - q) % 6, (q - p) % 6)
        phi[p, q] = {0: 0, 1: 0.3, 2: 0.05, 3: 0}[d]
Mv = np.array([1, 1, 0, 0, 0, 0])
print("omega =", phi @ Mv / phi.sum(1))

print("\n== Exercise 13.5: Sunkel ==")
P = np.array([[0.6, 0.2, 0.2], [0.2, 0.6, 0.2], [0.2, 0.2, 0.6]])
Q = P.copy(); Q[0, 1] += 0.15; Q[0, 0] -= 0.15; Q[1, 2] -= 0.05; Q[1, 0] += 0.05
def dom(W):
    S = W[1:, 1:]; S = S / S.sum(1)[:, None]; return gap(S)
print("W_after =\n", Q)
print(f"global Gamma {gap(P):.3f} -> {gap(Q):.3f}; domestic Gamma {dom(P):.3f} -> {dom(Q):.3f}")
rng = np.random.default_rng(0); up = 0; n = 0
for _ in range(2000):
    d = rng.uniform(0.4, 0.8, 3)
    P0 = np.array([[d[0], (1-d[0])/2, (1-d[0])/2], [(1-d[1])/2, d[1], (1-d[1])/2], [(1-d[2])/2, (1-d[2])/2, d[2]]])
    Q0 = P0.copy(); Q0[0, 1] += 0.15; Q0[0, 0] -= 0.15; Q0[1, 2] -= 0.05; Q0[1, 0] += 0.05
    n += 1; up += gap(Q0) > gap(P0)
print(f"global Gamma rises in {up/n:.0%} of 2000 random diagonal-dominant starts")

# ---------------------------------------------------------------------
print("\n== Proof checks (Props. multiplicity, integration, phase) ==")
# Prop. phase: spectrum of W(c) is {1, 1-c, 1-c, 0, 0, 0}; F is increasing (F' >= 0)
for c in np.linspace(0.05, 0.65, 13):
    ev = np.sort(np.real(np.linalg.eigvals(Wc(c))))[::-1]
    assert np.allclose(ev, [1, 1 - c, 1 - c, 0, 0, 0], atol=1e-8)
assert all(dF(c) >= 0 for c in np.linspace(0.05, 0.65, 61))
assert dF(cA) < 1 and dF(cB) < 1 and dF(chat) > 1
# Prop. integration: Cheeger phi^2/2 <= 1 - lambda_2 <= 2 phi for random lazy reversible chains,
# and eigenvalues of lazy reversible W lie in [0, 1]
rng = np.random.default_rng(13)
for _ in range(200):
    n = 5
    Sym = rng.uniform(0, 1, (n, n)) * (rng.uniform(size=(n, n)) < 0.6)
    Sym = Sym + Sym.T + 1e-3
    P = Sym / Sym.sum(1)[:, None]                    # reversible w.r.t. degree
    Wl = 0.5 * (np.eye(n) + P)                       # lazy
    ev = np.real(np.linalg.eigvals(Wl))
    assert ev.min() > -1e-10
    Ph = conductance(Wl)
    G = gap(Wl)
    assert Ph ** 2 / 2 <= G + 1e-12 and G <= 2 * Ph + 1e-12
print("Cheeger bounds and phase-diagram spectrum verified")
