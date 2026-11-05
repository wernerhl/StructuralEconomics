"""Solutions to the Chapter 13 exercises (How Networks Evolve).
Every number quoted in solutions/sol13.tex is printed here.
Run from the repository root:  python3 code/solutions/sol13.py   (numpy only).
Exercise 13.7 runs 1,300 short simulations; it takes about a minute."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


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


def v(z, kappa=0.5):
    return z - 0.5 * kappa * z ** 2


# ===========================================================================
print("== Exercise 13.1: edge creation by hand")
W1 = np.array([[0.6, 0.38, 0.01, 0.01], [0.38, 0.6, 0.01, 0.01],
               [0.01, 0.01, 0.6, 0.38], [0.01, 0.01, 0.38, 0.6]])
print(f"Gamma = {gap(W1):.4f}; eigenvalues {np.round(np.sort(np.linalg.eigvals(W1).real)[::-1], 4)}")
L = np.where(W1 > 0, 1 / W1, np.inf)
np.fill_diagonal(L, 0)
D = L.copy()
for k in range(4):
    D = np.minimum(D, D[:, [k]] + D[[k], :])
print("(a) weighted geodesic distance (edge length 1/W_ij):\n", D)
Wo = W1 - np.diag(np.diag(W1))
Lap = np.diag(Wo.sum(1)) - Wo
Lp = np.linalg.pinv(Lap)
R = np.array([[(np.eye(4)[i] - np.eye(4)[j]) @ Lp @ (np.eye(4)[i] - np.eye(4)[j]) for j in range(4)] for i in range(4)])
print("    effective resistance distance:\n", R)
pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
cross = lambda i, j: (i < 2) != (j < 2)
for sub in [1.0, 2.0]:
    w = np.array([(sub if cross(i, j) else 1.0) / D[i, j] for i, j in pairs])
    pr = w / w.sum()
    print(f"(b/c) subsidy factor {sub}: " + ", ".join(f"P({i+1},{j+1}) = {p:.4f}" for (i, j), p in zip(pairs, pr))
          + f"; total cross = {sum(p for (i,j),p in zip(pairs,pr) if cross(i,j)):.4f}")


def mf_step(W, sub, delta):
    Dd = W.copy()
    Ll = np.where(Dd > 0, 1 / Dd, np.inf)
    np.fill_diagonal(Ll, 0)
    for k in range(len(W)):
        Ll = np.minimum(Ll, Ll[:, [k]] + Ll[[k], :])
    P = np.zeros_like(W)
    tot = 0.0
    for i, j in pairs:
        w = (sub if cross(i, j) else 1.0) / Ll[i, j]
        P[i, j] = P[j, i] = w
        tot += w
    Wn = W + delta * P / tot
    return Wn / Wn.sum(1)[:, None]


delta = 0.05
paths = {}
for sub in [1.0, 2.0]:
    W = W1.copy()
    g = [gap(W)]
    for t in range(20):
        W = mf_step(W, sub, delta)
        g.append(gap(W))
    paths[sub] = np.array(g)
    print(f"(d) mean-field, subsidy {sub}: Gamma_t at t=0,5,10,15,20 = {np.round(paths[sub][[0,5,10,15,20]], 4)}; "
          f"W_20 diag {W[0,0]:.4f}, within {W[0,1]:.4f}, cross {W[0,2]:.4f}")
print("    plot coordinates (t, Gamma no subsidy, Gamma subsidy):",
      " ".join(f"({t},{paths[1.0][t]:.4f},{paths[2.0][t]:.4f})" for t in range(0, 21, 2)))
assert paths[2.0][-1] > paths[1.0][-1] > paths[1.0][0]
# stochastic version: one new edge of weight delta per period, symmetric, rows renormalized
rng = np.random.default_rng(131)
for sub in [1.0, 2.0]:
    finals, ncross = [], []
    for rep in range(2000):
        W = W1.copy()
        nc = 0
        for t in range(20):
            Ll = np.where(W > 0, 1 / W, np.inf)
            np.fill_diagonal(Ll, 0)
            for k in range(4):
                Ll = np.minimum(Ll, Ll[:, [k]] + Ll[[k], :])
            w = np.array([(sub if cross(i, j) else 1.0) / Ll[i, j] for i, j in pairs])
            i, j = pairs[rng.choice(6, p=w / w.sum())]
            nc += cross(i, j)
            W[i, j] += delta
            W[j, i] += delta
            W = W / W.sum(1)[:, None]
        finals.append(gap(W))
        ncross.append(nc)
    finals = np.array(finals)
    print(f"    stochastic, subsidy {sub}: mean Gamma_20 = {finals.mean():.4f}, 10-90% [{np.quantile(finals,0.1):.4f}, "
          f"{np.quantile(finals,0.9):.4f}], median {np.median(finals):.4f}, mean cross edges formed {np.mean(ncross):.2f} of 20, "
          f"share with no cross edge {np.mean(np.array(ncross) == 0):.3f} (Gamma_20 then {np.mean(finals[np.array(ncross) == 0]):.4f}), "
          f"share with Gamma_20 < 0.04: {np.mean(finals < 0.04):.3f}")

# ===========================================================================
print("\n== Exercise 13.2: structural equilibria (Example ch13:ex:bolivia_two)")
N = 6
s = np.array([0.25, 0.25, 0.15, 0.15, 0.10, 0.10])
r = np.array([0.20, 0.20, 0.12, 0.12, 0.06, 0.06])
B = np.zeros((N, N))
for b in range(3):
    B[2 * b:2 * b + 2, 2 * b:2 * b + 2] = 0.5
K = np.outer(np.ones(N), s)
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.10
M = mu / (1 - rho)
Omega, Phi_b = 6.0, 1.0


def Wc(c):
    return (1 - c) * B + c * K


def theta_star(c):
    return th0 + M * (Wc(c) @ r - rstar)


def h(c):
    th = theta_star(c)
    return Omega / N + (perron(Wc(c)) - 1 / N) * Phi_b + (th - th.mean())


def make_F(prox, surp=12.0):
    def F(c):
        z = prox * (c - 0.30) + surp * (theta_star(c)[4] - th0)
        return 0.05 + 0.60 / (1 + np.exp(-z))
    return F


def roots_of(F):
    g = lambda c: F(c) - c
    grid = np.linspace(0, 1, 4001)
    vals = np.array([g(c) for c in grid])
    out = []
    for i in np.where(np.diff(np.sign(vals)))[0]:
        lo, hi = grid[i], grid[i + 1]
        for _ in range(100):
            m = 0.5 * (lo + hi)
            if np.sign(g(m)) == np.sign(g(lo)):
                lo = m
            else:
                hi = m
        out.append(0.5 * (lo + hi))
    return out


F12 = make_F(12.0)
cA, chat, cB = roots_of(F12)
dF = lambda F, c, e=1e-6: (F(c + e) - F(c - e)) / (2 * e)
nKB = np.linalg.norm(K - B)
for nm, c in [("A", cA), ("threshold", chat), ("B", cB)]:
    W = Wc(c)
    th = theta_star(c)
    Tth = rho * th + (1 - rho) * th0 + mu * (W @ r - rstar)
    x = h(c)
    xh = Omega / N + (perron(W) - 1 / N) * Phi_b + (th - th.mean())
    print(f"(a) {nm}: c = {c:.4f}, F(c) - c = {F12(c)-c:.1e}, |theta - T(theta)| = {np.abs(th-Tth).max():.1e}, "
          f"|x - h| = {np.abs(x-xh).max():.1e}, Gamma = {gap(W):.4f}, pi = {perron(W)}, F'(c) = {dF(F12, c):.3f}, "
          f"g'(c) = 1 - eta + eta F' = {1-0.05+0.05*dF(F12,c):.4f}")
    print(f"     theta* = {np.round(th, 3)}, x* = {np.round(x, 3)}, sum x = {x.sum():.4f}")
print("    s'B =", s @ B, " (equals s')")
for c in [0.05, 0.3, 0.65]:
    ev = np.sort(np.real(np.linalg.eigvals(Wc(c))))[::-1]
    assert np.allclose(ev, [1, 1 - c, 1 - c, 0, 0, 0], atol=1e-8)
    assert np.allclose(perron(Wc(c)), s)
print(f"(b) c_hat = {chat:.4f}; c_hat - c_A = {chat-cA:.4f}; ||1s' - B||_F = {nKB:.4f}; eps* = {(chat-cA)*nKB:.4f}")
assert abs(chat - 0.282) < 5e-4 and abs((chat - cA) * nKB - 0.307) < 5e-4
F8 = make_F(8.0)
rts8 = roots_of(F8)
print(f"(c) proximity weight 8: fixed points {np.round(rts8, 4)}; slopes {[round(dF(F8, c), 3) for c in rts8]}")
if len(rts8) == 3:
    print(f"    c_A = {rts8[0]:.4f}, c_hat = {rts8[1]:.4f}, c_B = {rts8[2]:.4f}, eps* = {(rts8[1]-rts8[0])*nKB:.4f}")
for prox in [12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0]:
    rr = roots_of(make_F(prox))
    print(f"    proximity weight {prox}: fixed points {np.round(rr, 4)}")
# (d) welfare comparison
Ma, kap = 0.8, 0.5
print("(d) welfare at A and B (kappa = 0.5, aspiration multiplier 0.8, W^I = W(c))")
Wd = {}
for nm, c in [("A", cA), ("B", cB)]:
    W = Wc(c)
    x = h(c)
    a1 = Ma * W @ x
    z1 = x - a1
    ahat = Ma * x.mean()
    assert np.all(kap * z1 < 1) and np.all(kap * (x - ahat) < 1) and np.all(kap * x < 1)
    W1_ = v(z1, kap).sum()
    W2_ = v(x - ahat, kap).sum()
    W3_ = v(x, kap).sum()
    pi = perron(W)
    e = (Omega - Phi_b) / N + pi * Phi_b
    Wd[nm] = (W1_, W2_, W3_, e.min(), 1 / np.sum(pi ** 2), gap(W))
    print(f"   {nm}: W1 (endogenous) = {W1_:.4f}, W2 (laundered, a_hat = {ahat:.3f}) = {W2_:.4f}, "
          f"W3 (a = 0) = {W3_:.4f}, Var(x) = {np.var(x):.5f}, min e = {e.min():.4f}, Neff(pi) = {1/np.sum(pi**2):.3f}, "
          f"Neff(pi^I) = {1/np.sum(perron(W)**2):.3f}, Gamma(W^I) = {gap(W):.3f}, rural x = {x[4]:.4f}, rural s1 = {v(z1,kap)[4]:.4f}")
print("   differences B - A:", np.round(np.array(Wd["B"][:4]) - np.array(Wd["A"][:4]), 4))

# ===========================================================================
print("\n== Exercise 13.3: gas nationalization (MNC, Gas company, Government, Traditional)")
Wb = np.array([[0.70, 0.25, 0.03, 0.02],
               [0.55, 0.30, 0.10, 0.05],
               [0.05, 0.15, 0.50, 0.30],
               [0.02, 0.03, 0.15, 0.80]])
Wa = np.array([[0.85, 0.10, 0.03, 0.02],
               [0.10, 0.40, 0.40, 0.10],
               [0.02, 0.35, 0.40, 0.23],
               [0.02, 0.03, 0.15, 0.80]])
for nm, W in [("before", Wb), ("after", Wa)]:
    assert np.allclose(W.sum(1), 1)
    ev = eig_sorted(W)
    p = perron(W)
    print(f"{nm}: eigenvalues {np.round(ev, 4)}, Gamma = {gap(W):.4f}, pi = {np.round(p, 4)}, Neff = {1/np.sum(p**2):.3f}, "
          f"column sums {np.round(W.sum(0), 2)}")
print("   Delta pi =", np.round(perron(Wa) - perron(Wb), 4))
Eab = Wa - Wb
pb = perron(Wb)
Zb = np.linalg.inv(np.eye(4) - Wb + np.outer(np.ones(4), pb))
print("   first-order Delta pi' = pi' E Z =", np.round(pb @ Eab @ Zb, 4), "(large perturbation; exact above)")
# sourced values quoted from the book (INE, base 1990, book12, domestic, input shares)
import csv
import os
DATA = os.path.join("data", "processed")
with open(os.path.join(DATA, "bolivia_io_summary.csv"), newline="", encoding="utf-8") as f:
    gam = {r["year"]: float(r["Gamma"]) for r in csv.DictReader(f)
           if r["base"] == "1990" and r["aggregation"] == "book12" and r["flows"] == "domestic"
           and r["normalization"] == "input_shares"}
with open(os.path.join(DATA, "bolivia_io_pi.csv"), newline="", encoding="utf-8") as f:
    pih = {r["year"]: float(r["pi"]) for r in csv.DictReader(f)
           if r["base"] == "1990" and r["aggregation"] == "book12" and r["flows"] == "domestic"
           and r["normalization"] == "input_shares" and r["sector"] == "Hydrocarbons"}
for y in ["2003", "2005", "2006", "2009", "2013"]:
    print(f"   INE {y}: Gamma(W^P) = {gam[y]:.3f}, pi_hydrocarbons = {pih[y]:.3f}")

# ===========================================================================
print("\n== Exercise 13.4: ring relatedness")
phi = np.zeros((6, 6))
for p in range(6):
    for q in range(6):
        d = min((p - q) % 6, (q - p) % 6)
        phi[p, q] = {0: 0, 1: 0.3, 2: 0.05, 3: 0}[d]
def density(basket):
    Mv = np.zeros(6)
    Mv[list(basket)] = 1
    return phi @ Mv / phi.sum(1)
om = density([0, 1])
print("(a) omega for products 1..6 =", np.round(om, 4))
assert abs(om[2] - 0.5) < 1e-12 and abs(om[5] - 0.5) < 1e-12
for basket in [[0, 1, 5], [0, 1, 2]]:
    print(f"(c) basket {[b+1 for b in basket]}: omega =", np.round(density(basket), 4))
with open(os.path.join(DATA, "bolivia_facts.csv"), newline="", encoding="utf-8") as f:
    facts = {r["fact_id"]: float(r["value"]) for r in csv.DictReader(f)}
print(f"(e) WDI 2019 (data/processed/bolivia_facts.csv): ores and metals {facts['BOL_ores_metals_exports_share_2019']:.1f}%, "
      f"fuels {facts['BOL_fuel_exports_share_2019']:.1f}% of merchandise exports")
for p1 in [2.0, 5.0]:
    print(f"(d) logit gap to make product 5 as likely as product 6: p1*(0.5 - {om[4]:.4f}) = {p1*(0.5-om[4]):.3f} (p1 = {p1})")

# ===========================================================================
print("\n== Exercise 13.5: Sunkel")
def dom(W):
    S = W[1:, 1:]
    S = S / S.sum(1)[:, None]
    return gap(S)
def reroute(P):
    Q = P.copy()
    Q[0, 1] += 0.15; Q[0, 0] -= 0.15; Q[1, 2] -= 0.05; Q[1, 0] += 0.05
    return Q
P0 = np.array([[0.6, 0.2, 0.2], [0.2, 0.6, 0.2], [0.2, 0.2, 0.6]])
Q0 = reroute(P0)
print("(a) W_after =\n", Q0, "\n    row sums", Q0.sum(1))
print(f"(b) domestic Gamma {dom(P0):.4f} -> {dom(Q0):.4f}")
print(f"(c) global Gamma {gap(P0):.4f} -> {gap(Q0):.4f}; eigenvalues after {np.round(eig_sorted(Q0), 4)}")
P1 = np.array([[0.90, 0.05, 0.05], [0.05, 0.60, 0.35], [0.05, 0.35, 0.60]])
Q1 = reroute(P1)
print("(d) W_before' =\n", P1, "\n    eigenvalues", np.round(eig_sorted(P1), 4))
print(f"    global Gamma {gap(P1):.4f} -> {gap(Q1):.4f}; domestic Gamma {dom(P1):.4f} -> {dom(Q1):.4f}")
assert gap(Q1) > gap(P1) and dom(Q1) < dom(P1)
# first-order explanation: d lambda2 = u2' E v2 (symmetric W: u2 = v2, unit norm)
E = Q1 - P1
for nm, P in [("original", P0), ("alternative", P1)]:
    ev, V = np.linalg.eigh(P)
    o = np.argsort(-ev)
    print(f"    {nm}: lambda2 eigenvector(s):")
    lam2 = ev[o[1]]
    idx = [k for k in o[1:] if abs(ev[k] - lam2) < 1e-9]
    Vs = V[:, idx]
    Er = reroute(P) - P
    first = np.sort(np.linalg.eigvals(Vs.T @ Er @ Vs).real)[::-1]
    print(f"       lambda2 = {lam2:.4f} (multiplicity {len(idx)}); first-order shifts of lambda2: {np.round(first, 4)}"
          f" -> predicted lambda2 {lam2 + first.max():.4f}, exact {abs(eig_sorted(reroute(P))[1]):.4f}")
    if len(idx) == 1:
        print(f"       eigenvector v2 = {np.round(Vs[:, 0], 4)}")

# ===========================================================================
print("\n== Exercise 13.6: WWII as structural change (block-level formation model)")
rho6, mu6, th06, rs6, rbar = 0.6, 0.5, 0.5, 0.16, 0.2
M6 = mu6 / (1 - rho6)
def signal(k, l):
    rF = rbar * k               # return to market work observed among women
    return (1 - l) * rF + l * rbar
th_pre = th06 + M6 * (signal(0.0, 0.0) - rs6)
th_M = th06 + M6 * (rbar - rs6)
print(f"M = {M6}; pre-war theta_F* = {th_pre:.3f}; men's theta_M* = {th_M:.3f}")
assert abs(th_pre - 0.3) < 1e-12
th = th_pre
th_war = rho6 * th + (1 - rho6) * th06 + mu6 * (signal(0.5, 0.0) - rs6)
print(f"(b) war generation: within-generation change in theta = 0 (predetermined); "
      f"signal {signal(0,0):.3f} -> {signal(0.5,0):.3f}; next generation theta = {th_war:.4f} (+{th_war-th_pre:.4f})")
traj = [th_pre, th_war]
th = th_war
for g in range(3):
    th = rho6 * th + (1 - rho6) * th06 + mu6 * (signal(0.2, 0.3) - rs6)
    traj.append(th)
th_new = th06 + M6 * (signal(0.2, 0.3) - rs6)
print(f"    post-war signal s = 0.7*0.2*0.2 + 0.3*0.2 = {signal(0.2, 0.3):.3f}")
print("(c) post-war k = 0.2, l = 0.3: theta path (pre, after war, +1, +2, +3) =", np.round(traj, 4),
      f"; new steady state {th_new:.4f}")
th = th_war
tr2 = []
for g in range(3):
    th = rho6 * th + (1 - rho6) * th06 + mu6 * (signal(0.0, 0.0) - rs6)
    tr2.append(th)
print("    counterfactual: full reversal k = l = 0:", np.round(tr2, 4), " -> returns to", th_pre)
th_k0 = th06 + M6 * (signal(0.0, 0.3) - rs6)
print(f"(d) if only W^I persists (k = 0, l = 0.3): steady state {th_k0:.4f}; if only W^P persists (k = 0.2, l = 0): "
      f"{th06 + M6*(signal(0.2,0.0)-rs6):.4f}")

# ===========================================================================
print("\n== Exercise 13.7: network evolution simulation (N = 8)")
N7 = 8
off = ~np.eye(N7, dtype=bool)
SELF = 0.2          # fixed self-weight: rows are SELF on the diagonal + (1 - SELF) spread over links


def evolve(W, p1, T=400, p0=-2.0, p2=2.0, wnew=0.05, d0=0.05, d1=1.0, eta=0.02, Phi7=4.0, s0=0.2, lam=1.0):
    change = np.inf
    for t in range(T):
        p = perron(W)
        x = 1 + (p - 1 / N7) * Phi7                       # allocation rule, Omega/N = 1, Phi = 4
        S = s0 - lam * np.abs(x[:, None] - x[None, :])     # surplus of link (i, j)
        A = W @ W
        A = A / A[off].mean()                             # proximity: two-step reach, mean 1
        C = 1 / (1 + np.exp(-(p0 + p1 * (A - 1) + p2 * S)))   # creation probability (logistic)
        Dm = np.clip(d0 + d1 * np.maximum(0, -S), 0, 1)        # destruction probability
        Wn = (1 - Dm) * W + wnew * C + eta * S * W             # survivors + new edges + reweighting
        Wn[~off] = 0
        Wn = np.maximum(Wn, 0)                                 # operator N: clip ...
        Wn = (1 - SELF) * Wn / Wn.sum(1)[:, None]              # ... and renormalize the links
        Wn[~off] = SELF
        change = np.abs(Wn - W).max()
        W = Wn
    return W, change


def start(W0):
    W = (1 - SELF) * (W0 * off) / (W0 * off).sum(1)[:, None]
    W[~off] = SELF
    return W


rng = np.random.default_rng(137)
W0s = [start(rng.dirichlet(np.ones(N7) * 0.5, size=N7)) for _ in range(100)]
print(f"initial W0: median Gamma {np.median([gap(W) for W in W0s]):.3f}, range "
      f"[{min(gap(W) for W in W0s):.3f}, {max(gap(W) for W in W0s):.3f}]")
grid_p1 = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0]
frac = {}
for p1 in grid_p1:
    gs, ch = [], []
    for W0 in W0s:
        W, c = evolve(W0.copy(), p1)
        gs.append(gap(W))
        ch.append(c)
    gs, ch = np.array(gs), np.array(ch)
    frac[p1] = (np.mean(gs > 0.3), np.mean(gs < 0.1))
    print(f"p1 = {p1:4.2f}: integrated {frac[p1][0]:.2f}, fragmented {frac[p1][1]:.2f}, "
          f"neither {1-frac[p1][0]-frac[p1][1]:.2f}, median Gamma {np.median(gs):.3f}, "
          f"settled (last change < 1e-4) {np.mean(ch < 1e-4):.2f}")
for a_, b_ in zip(grid_p1[:-1], grid_p1[1:]):
    if frac[a_][1] <= 0.5 < frac[b_][1]:
        pstar = a_ + (0.5 - frac[a_][1]) / (frac[b_][1] - frac[a_][1]) * (b_ - a_)
        print(f"fragmented share crosses 50% between p1 = {a_} and {b_}; linear interpolation p1* = {pstar:.3f}")
print("robustness: horizon T = 800")
for p1 in [1.0, 1.25, 1.5]:
    gs = np.array([gap(evolve(W0.copy(), p1, T=800)[0]) for W0 in W0s])
    print(f"   p1 = {p1}: integrated {np.mean(gs > 0.3):.2f}, fragmented {np.mean(gs < 0.1):.2f}")
W, _ = evolve(W0s[0].copy(), 3.0)
print("example final W (p1 = 3, draw 0):\n", np.round(W, 2), "\n Gamma =", round(gap(W), 4),
      " eigenvalue moduli", np.round(np.abs(eig_sorted(W)), 3))
W, _ = evolve(W0s[0].copy(), 0.0)
print("same draw, p1 = 0: Gamma =", round(gap(W), 4), "; off-diagonal entries in", np.round([W[off].min(), W[off].max()], 4))
