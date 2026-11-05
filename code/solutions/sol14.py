"""Solutions to the Chapter 14 exercises (Coalitions and Class Alliances).
Every number quoted in solutions/sol14.tex is printed here.
Run from the repository root:  python3 code/solutions/sol14.py   (numpy only; csv for data/processed)."""
import csv
import itertools
import os
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def gap(W):
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


names = "FBMW"
pi0 = np.array([0.35, 0.30, 0.20, 0.15])
DP = {"A1 {F,B}": np.array([+0.05, +0.03, -0.04, -0.04]),
      "A2 {B,M,W}": np.array([-0.10, +0.02, +0.05, +0.03]),
      "A3 {F,B,M}": np.array([+0.03, +0.04, +0.02, -0.09]),
      "A4 {M,W}": np.array([-0.10, -0.02, +0.05, +0.07])}
MEM = {"A1 {F,B}": [0, 1], "A2 {B,M,W}": [1, 2, 3], "A3 {F,B,M}": [0, 1, 2], "A4 {M,W}": [2, 3]}

# ===========================================================================
print("== Exercise 14.1: coalition payoffs")
Omega, Phi, N, kap = 4.0, 2.0, 4, 0.5
x = Omega / N + (pi0 - 1 / N) * Phi
print("(a) x* =", x, " kappa x =", kap * x, " sum =", x.sum())
assert np.allclose(x, [1.2, 1.1, 0.9, 0.8]) and np.all(kap * x < 1)
u = np.full(4, 0.25)
v = Phi * (u - pi0)
print("(b) v_i = Phi (pi_new - pi) =", v, " -> v_3, v_4 =", v[2], v[3])
U = lambda z: z - 0.5 * kap * z ** 2
print("    utility changes with U = x - (kappa/2) x^2:", np.round(U(x + v) - U(x), 4),
      " sum", round(float(np.sum(U(x + v) - U(x))), 4), " L_power before =", 0.5 * kap * N * Phi ** 2 * np.mean((pi0 - 0.25) ** 2))
# (c) minimum-norm Delta W that makes u the Perron vector: Delta W = 1 r' / ||u||^2 * u ... with r' = u'(I - W)
for gam_ in [1.0, 0.5, 0.2]:
    W = gam_ * np.outer(np.ones(4), pi0) + (1 - gam_) * np.eye(4)
    r = u @ (np.eye(4) - W)
    dW = np.outer(u, r) / (u @ u)
    assert np.allclose(perron(W + dW), u) and np.allclose(dW.sum(1), 0) and np.all(W + dW >= -1e-12)
    print(f"(c) W = gamma 1 pi' + (1-gamma) I with gamma = {gam_}: min ||Delta W||_F = {np.linalg.norm(dW):.4f} "
          f"(= gamma*2*||u - pi|| = {gam_*2*np.linalg.norm(u-pi0):.4f})")
print(f"    ||u - pi||_2 = {np.linalg.norm(u - pi0):.4f}; feasibility with phi(p) = p: bound {0.35:.2f}; "
      f"with phi(p) = p^2: bound {0.35**2:.4f} -> gamma <= {0.35**2/(2*np.linalg.norm(u-pi0)):.3f}")
print("(d) agent 2: Delta pi_2 in A1 {1,2} =", DP["A1 {F,B}"][1], "; in A2 {2,3,4} =", DP["A2 {B,M,W}"][1])

# ===========================================================================
print("\n== Exercise 14.2: alliance stability")
for k, d in DP.items():
    assert abs(d.sum()) < 1e-12
    gains = {names[i]: d[i] for i in MEM[k]}
    print(f"{k}: members' Delta pi = {gains}, all positive: {all(val > 0 for val in gains.values())}")
stable = {}
for k in DP:
    issues = []
    for i in MEM[k]:
        alt = {k2: DP[k2][i] for k2 in DP if k2 != k}
        best = max(alt, key=alt.get)
        if alt[best] > DP[k][i] + 1e-12:
            issues.append(f"{names[i]} prefers {best} ({alt[best]:+.2f} > {DP[k][i]:+.2f})")
        elif abs(alt[best] - DP[k][i]) < 1e-12:
            issues.append(f"{names[i]} indifferent with {best} ({alt[best]:+.2f})")
    stable[k] = not any("prefers" in s for s in issues)
    print(f"   {k}: internally stable = {stable[k]}; " + "; ".join(issues))
assert [k for k in stable if stable[k]] == ["A4 {M,W}"]
# conditions for each to hold: required Delta pi for members to resist best alternatives
for k in DP:
    req = {names[i]: max(DP[k2][i] for k2 in DP if k2 != k) for i in MEM[k]}   # best rival offer
    outs = [i for i in range(4) if i not in MEM[k]]
    print(f"   {k}: to deter all defections members need >= {req}; sum = {sum(req.values()):.2f} "
          f"(program delivers {DP[k][MEM[k]].sum():.2f}); outsiders would have to lose {sum(req.values()):.2f} "
          f"(their current power {pi0[outs].sum():.2f})")
print("   B's payoff from A4 =", DP["A4 {M,W}"][1], "; compensation needed (1+Psi) Phi * 0.02 =", Phi * 0.02, "(Psi = 0)")

# ===========================================================================
print("\n== Exercise 14.3: the Latin American pendulum (simulation)")
AL = {"A1": [0, 1], "A2": [1, 2, 3], "A3": [0, 1, 2], "A4": [2, 3]}
W0 = 0.5 * np.eye(4) + 0.5 * np.outer(np.ones(4), pi0)
assert np.allclose(perron(W0), pi0)
PHI0, DELTA, FLOOR = 0.3, 0.2, 0.02


def program(W, S):
    """Coalition S's power-maximizing Delta W, first-order (Lemma ch06:lem:perturbation) on the
    Frobenius ball ||Delta W|| <= PHI0 * pi(S), rows kept stochastic; returns new W and exact Delta pi."""
    p = perron(W)
    Z = np.linalg.inv(np.eye(4) - W + np.outer(np.ones(4), p))
    z = Z @ np.isin(np.arange(4), S).astype(float)
    G = np.outer(p, z)
    G = G - G.mean(1)[:, None]                 # project onto Delta W 1 = 0
    eps = PHI0 * p[S].sum()
    Wn = np.maximum(W + eps * G / np.linalg.norm(G), FLOOR)
    Wn = Wn / Wn.sum(1)[:, None]
    return Wn, perron(Wn) - p, eps


W = W0.copy()
S = "A2"
rows = []
for t in range(8):
    p_before = perron(W)
    Wn, dp, eps = program(W, AL[S])
    W = (1 - DELTA) * Wn + DELTA * W0           # institutional inertia toward the initial structure
    p = perron(W)
    prog = {k: program(W, v_)[1] for k, v_ in AL.items()}
    viable = sorted(k for k, v_ in AL.items() if all(prog[k][i] > 0 for i in v_))
    best = None
    for i in AL[S]:
        for k in viable:
            if k != S and i in AL[k] and prog[k][i] > prog[S][i] + 1e-4:
                g_ = prog[k][i] - prog[S][i]
                if best is None or g_ > best[0]:
                    best = (g_, i, k)
    rows.append((t, S, p, best))
    print(f"t={t}: alliance {S} pushes ||dW|| = {eps:.3f}; pi_(t+1) = {np.round(p, 3)}, Gamma = {gap(W):.3f}; "
          f"viable next: {viable}; " + (f"{names[best[1]]} defects to {best[2]} (gain {best[0]:+.3f} vs staying)"
                                         if best else "no defection"))
    print("      prospective Delta pi by alliance:", {k: np.round(v_, 3).tolist() for k, v_ in prog.items()})
    if best:
        S = best[2]
seq = [r_[1] for r_ in rows]
print("alliance sequence:", seq)
print("pi_B path:", np.round([r_[2][1] for r_ in rows], 3), " pi_F path:", np.round([r_[2][0] for r_ in rows], 3))
assert seq[2:] == ["A1", "A2", "A1", "A2", "A1", "A2"]      # a two-alliance cycle from t = 2 on
assert np.all(np.diff([r_[2][1] for r_ in rows]) > -0.01)   # B's power ratchets up
print("variant without inertia (delta = 0):")
DELTA0 = DELTA
DELTA = 0.0
W = W0.copy(); S = "A2"; seq0 = []; pB = []
for t in range(8):
    Wn, dp, eps = program(W, AL[S])
    W = Wn
    p = perron(W)
    prog = {k: program(W, v_)[1] for k, v_ in AL.items()}
    viable = [k for k, v_ in AL.items() if all(prog[k][i] > 0 for i in v_)]
    best = None
    for i in AL[S]:
        for k in viable:
            if k != S and i in AL[k] and prog[k][i] > prog[S][i] + 1e-4:
                if best is None or prog[k][i] - prog[S][i] > best[0]:
                    best = (prog[k][i] - prog[S][i], i, k)
    seq0.append(S); pB.append(p[1])
    if best:
        S = best[2]
print("   sequence", seq0, " pi_B", np.round(pB, 3))
DELTA = DELTA0

# ===========================================================================
print("\n== Exercise 14.4: Bolivia (sourced values quoted)")
DATA = os.path.join("data", "processed")
with open(os.path.join(DATA, "bolivia_facts.csv"), newline="", encoding="utf-8") as f:
    facts = {r["fact_id"]: float(r["value"]) for r in csv.DictReader(f)}
for k in ["BOL_fuel_exports_share_2006", "BOL_fuel_exports_share_2013", "BOL_fuel_exports_share_2019",
          "BOL_natgas_rents_pct_gdp_2006", "BOL_natgas_rents_pct_gdp_2013", "BOL_natgas_rents_pct_gdp_2019"]:
    print(f"   {k}: {facts[k]:.2f}")
with open(os.path.join(DATA, "bolivia_io_summary.csv"), newline="", encoding="utf-8") as f:
    gam = {r["year"]: float(r["Gamma"]) for r in csv.DictReader(f)
           if r["base"] == "1990" and r["aggregation"] == "book12" and r["flows"] == "domestic"
           and r["normalization"] == "input_shares"}
with open(os.path.join(DATA, "bolivia_io_pi.csv"), newline="", encoding="utf-8") as f:
    pih = {r["year"]: float(r["pi"]) for r in csv.DictReader(f)
           if r["base"] == "1990" and r["aggregation"] == "book12" and r["flows"] == "domestic"
           and r["normalization"] == "input_shares" and r["sector"] == "Hydrocarbons"}
for y in ["1990", "2000", "2006", "2013"]:
    print(f"   INE base 1990, {y}: Gamma(W^P) = {gam[y]:.3f}, pi_hydrocarbons = {pih[y]:.3f}")

# ===========================================================================
print("\n== Exercise 14.5: transfers and the Perron vector")
W = W0.copy()
p = perron(W)
om = np.full(4, 0.5)
dom_ = np.array([-0.1, 0.0, 0.0, 0.1])
print("within the period: pi(W) depends on W only; Delta omega =", dom_, "-> Delta pi = 0")
ETA = 0.5


def E_map(W, x):
    """One step of an evolution map: each agent reweights partners by their allocation (attractiveness)."""
    Wn = W * (1 + ETA * x[None, :])
    return Wn / Wn.sum(1)[:, None]


x0 = om + p * Phi
x1 = om + dom_ + p * Phi
Wa, Wb = E_map(W, x0), E_map(W, x1)
pa, pb = perron(Wa), perron(Wb)
print("x without / with transfer:", x0, x1)
print("pi_(t+1) without transfer:", np.round(pa, 5), " with:", np.round(pb, 5), " difference:", np.round(pb - pa, 5))
# first-order: E = D_x E [Delta x]; Delta pi' = pi_a' E Z_a
h = 1e-6
Ed = (E_map(W, x0 + h * dom_) - Wa) / h
Za = np.linalg.inv(np.eye(4) - Wa + np.outer(np.ones(4), pa))
lin = pa @ Ed @ Za
print("first-order Delta pi' = pi' (D_x E . Delta x) Z:", np.round(lin, 5))
assert np.allclose(lin, pb - pa, atol=2e-4)
# compare with a direct Delta W of the same Frobenius size as the induced change
dW_ind = Wb - Wa
print(f"   size of induced Delta W: {np.linalg.norm(dW_ind):.5f}; |Delta pi_4| per unit of transfer = {(pb-pa)[3]/0.1:.4f}")
up = np.full(4, 0.25)
r_ = up @ (np.eye(4) - W)
print(f"   for comparison, the minimum-norm Delta W that equalizes power has norm {np.linalg.norm(np.outer(up, r_)/(up@up)):.4f}")

# ===========================================================================
print("\n== Exercise 14.6: formation amplifier in politics")
dpi, Ph = 0.05, 2.0
for nm, Ps in [("A (W^P)", 1.0), ("B (W^I)", 3.0)]:
    print(f"   {nm}: Psi = {Ps}, v_i = (1+Psi) Phi dpi = {(1+Ps)*Ph*dpi:.3f}")
amp = lambda Ps, lam, d: 1 + Ps * d * (1 - lam) / (1 - d * lam)
for d in [0.9, 0.5]:
    for lamA, lamB in [(0.8, 0.8), (0.6, 0.9)]:
        aA, aB = amp(1.0, lamA, d), amp(3.0, lamB, d)
        print(f"   discount {d}, lambda_A = {lamA}, lambda_B = {lamB}: effective amplifier A = {aA:.3f}, B = {aB:.3f}")

# ===========================================================================
print("\n== Exercise 14.7: empty core")
gam7, Phi7, n7 = 0.5, 2.0, 3
phi_f = lambda P: 0.6 * max(0.0, P - 0.5)       # minority coalitions cannot act
pi7 = np.full(3, 1 / 3)


def v_pow(S):
    """max Phi * sum_{i in S} Delta pi_i over Delta W = gamma 1 Delta pi' with ||Delta W||_F <= phi(pi(S))."""
    if len(S) in (0, n7):
        return 0.0
    rS = phi_f(pi7[list(S)].sum()) / (gam7 * np.sqrt(n7))      # bound on ||Delta pi||_2
    proj = np.sqrt(len(S) * (n7 - len(S)) / n7)                # norm of 1_S projected on sum-zero
    return Phi7 * rS * proj


vals = {S: v_pow(S) for r_ in range(1, 4) for S in itertools.combinations(range(3), r_)}
print("v(S) with Delta Omega = 0:", {tuple(i + 1 for i in S): round(val, 5) for S, val in vals.items()})
g = vals[(0, 1)]
print(f"   pair value g = {g:.6f} = 0.4 sqrt(2)/3 = {0.4*np.sqrt(2)/3:.6f}")
S12 = [0, 1]
rS = phi_f(2 / 3) / (gam7 * np.sqrt(3))
d = np.array([1, 1, -2]) / 3
d = rS * d / np.linalg.norm(d)
W7 = gam7 * np.ones((3, 3)) / 3 + (1 - gam7) * np.eye(3)
W7n = W7 + gam7 * np.outer(np.ones(3), d)
print("   optimal Delta pi for {1,2}:", np.round(d, 5), " new pi:", np.round(perron(W7n), 5),
      f" ||Delta W||_F = {np.linalg.norm(W7n - W7):.4f} = phi(2/3) = {phi_f(2/3):.4f}")
assert np.allclose(perron(W7n), pi7 + d)
thr = 1.5 * g
print(f"   core nonempty iff (1+Psi) Delta Omega >= 1.5 g = {thr:.5f} (= 0.2 sqrt 2 = {0.2*np.sqrt(2):.5f})")


def core_nonempty(D, g):
    # brute-force search on a grid of the simplex x1+x2+x3 = D
    grid = np.linspace(0, D, 301)
    for a in grid:
        for b in grid:
            c = D - a - b
            if c < -1e-12:
                continue
            if a >= -1e-12 and b >= -1e-12 and a + b >= g - 1e-9 and a + c >= g - 1e-9 and b + c >= g - 1e-9:
                return True, (a, b, c)
    return False, None


for D in [0.0, 0.2, 0.28, thr, 0.35]:
    ok, xx = core_nonempty(D, g)
    print(f"   (1+Psi) Delta Omega = {D:.5f}: core nonempty = {ok}" + (f", e.g. x = {np.round(xx, 4)}" if ok else ""))
