"""Solutions, Chapter 24 (Hegemony and Counter-Hegemony).
Prints every number used in solutions/sol24.tex.  Run from the repo root:
    python3 code/solutions/sol24.py      (numpy only)
All parameter values are illustrative unless stated otherwise."""
import itertools
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi = pi / pi.sum()
    return pi, 1 - abs(ev[o[1]])


def g(th, beta=10.0):
    return 1 / (1 + np.exp(-beta * (th - 0.5)))


def dg(th, beta=10.0):
    return beta * g(th, beta) * (1 - g(th, beta))


def fixed_points(s, beta=10.0):
    grid = np.linspace(-0.5, 1.5, 400001)
    f = g(grid, beta) + s - grid
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    out = []
    for i in idx:                      # refine by bisection
        lo, hi = grid[i], grid[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if np.sign(g(lo, beta) + s - lo) == np.sign(g(mid, beta) + s - mid):
                lo = mid
            else:
                hi = mid
        r = 0.5 * (lo + hi)
        if not out or abs(r - out[-1]) > 1e-6:
            out.append(r)
    return np.array(out)


def s_star(beta):
    gs = (1 - np.sqrt(1 - 4 / beta)) / 2
    a = 0.5 + np.log(gs / (1 - gs)) / beta
    return a - gs, a


rho, beta = 0.85, 10.0
SSTAR, A_FOLD = s_star(beta)

# ---------------------------------------------------------------------------
print("=== 24.1 Hegemony as a fixed point (two groups E, M) ===")
WP = np.array([[0.9, 0.1], [0.4, 0.6]])
piP, GP = perron(WP)
neffP = 1 / np.sum(piP ** 2)
varP = np.mean((piP - 0.5) ** 2)
dpi = 0.5 * np.abs(piP - 0.5).sum()
print(f"W^P = {WP.tolist()}: pi^P = {piP}, Gamma = {GP:.3f}, Neff = {neffP:.4f} of 2, "
      f"Var = {varP:.4f} (max 0.25), d_pi = {dpi:.3f} (max 0.5)")
assert np.allclose(piP, [0.8, 0.2]) and abs(varP - 0.09) < 1e-12 and np.isclose(piP.sum(), 1)
Omega, Phi = 2.0, 1.0
x = Omega / 2 + (piP - 0.5) * Phi
print(f"allocation x = Omega/N + (pi - 1/N) Phi = {x}  (sum {x.sum():.1f})")
assert np.allclose(x, [1.3, 0.7])
q = 0.7
WI = np.array([[0.9, 0.1], [q, 1 - q]])
piI, GI = perron(WI)
print(f"W^I (q = {q}) = {WI.tolist()}: pi^I = {piI}, Neff^I = {1/np.sum(piI**2):.3f}, Gamma^I = {GI:.3f}")
assert np.allclose(piI, [0.875, 0.125])
sigma = 0.5
s_heg = sigma * (1 - q) * (x[0] - x[1])
print(f"counter-signal s = sigma (1-q)(x_E - x_M) = {sigma}*{1-q:.1f}*{x[0]-x[1]:.1f} = {s_heg:.3f}; s* = {SSTAR:.4f}")
fp = fixed_points(s_heg)
slopes = rho + (1 - rho) * dg(fp)
print(f"fixed points of theta at s = {s_heg:.2f}: {np.round(fp, 4)}, map slopes {np.round(slopes, 4)}")
th_h = fp[0]
assert len(fp) == 3 and slopes[0] < 1 and slopes[1] > 1 and slopes[2] < 1
c_mob = 0.15
gain = 1.0 - x[1]                       # x_M under W^ctr minus x_M under W^heg
thbar = c_mob / gain
print(f"mobilization: contest iff theta * {gain:.1f} > c = {c_mob} iff theta > {thbar:.2f}")
print(f"consent at theta_h = {th_h:.4f}: perceived gain {th_h*gain:.4f} < c = {c_mob}")
assert th_h * gain < c_mob and abs(thbar - 0.5) < 1e-12
# counter-hegemonic structure: doubly stochastic W^P, q = 0.3
WPc = np.array([[0.75, 0.25], [0.25, 0.75]])
piPc, _ = perron(WPc)
fpc = fixed_points(0.0)
print(f"W^ctr: pi^P = {piPc}, x = (1, 1), s = 0; fixed points {np.round(fpc, 4)}; theta_c = {fpc[-1]:.4f} > 1/2")
assert np.allclose(piPc, 0.5) and fpc[-1] > 0.5
qstar = 1 - SSTAR / (sigma * (x[0] - x[1]))
print(f"hegemonic equilibrium exists iff q > q* = 1 - s*/(sigma*Dx) = {qstar:.3f}")
assert 0.39 < qstar < 0.40


def simulate(q_, th0, T=200):
    th, regime, switch = th0, "heg", None
    for t in range(T):
        if regime == "heg":
            s = sigma * (1 - q_) * (x[0] - x[1])
        else:
            s = 0.0
        th = rho * th + (1 - rho) * (g(th) + s)
        if regime == "heg" and th > thbar:
            regime, switch = "ctr", t + 1
    return th, regime, switch


for q_ in (0.7, 0.35):
    th, reg, sw = simulate(q_, 0.0)
    print(f"simulate from theta_0 = 0 with q = {q_}: theta_200 = {th:.4f}, regime {reg}, switch at t = {sw}")
th, reg, sw = simulate(0.7, 0.6)
print(f"simulate from theta_0 = 0.6 with q = 0.7: theta_200 = {th:.4f}, regime {reg}, switch at t = {sw}")
assert simulate(0.7, 0.0)[1] == "heg" and simulate(0.35, 0.0)[1] == "ctr" and simulate(0.7, 0.6)[1] == "ctr"
for bb in (4.0, 6.0, 10.0, 20.0):
    ss, aa = s_star(bb) if bb > 4 else (0.0, 0.5)
    print(f"beta = {bb:4.1f}: s* = {ss:.4f}, theta_c(s=0) = {fixed_points(0.0, bb)[-1]:.4f}")

# ---------------------------------------------------------------------------
print("\n=== 24.2 Maintenance cost: threshold Neff with one controlled outlet, n outlets ===")
n = 4
pi_state = 0.34
neff_now = 1 / (pi_state ** 2 + (1 - pi_state) ** 2 / (n - 1))
print(f"n = 4, controlled share 0.34, independents symmetric: Neff = {neff_now:.3f}")
for kappa in (0.30, 0.50, 0.66):
    ne = 1 / ((1 - kappa) ** 2 + kappa ** 2 / (n - 1))
    print(f"kappa = WTP/(mu_H delta (1-lambda2)) = {kappa:.2f}: cost exceeds WTP once independent share > {kappa:.2f}, "
          f"i.e. Neff > {ne:.3f}")
# Neff is increasing in the independent share on (0, (n-1)/n)
sh = np.linspace(0.01, (n - 1) / n, 200)
ne = 1 / ((1 - sh) ** 2 + sh ** 2 / (n - 1))
assert np.all(np.diff(ne) > 0) and abs(ne[-1] - n) < 1e-9

# ---------------------------------------------------------------------------
print("\n=== 24.3 Counter-hegemony simulation: lambda^I_t = 0.02 t, 50 periods ===")
fp0 = fixed_points(0.0)
print(f"start: theta_h(s=0) = {fp0[0]:.4f}; s* = {SSTAR:.4f} (fold at theta = {A_FOLD:.4f})")


def run(rho_, conv, T=50):
    th = fixed_points(0.0)[0]
    t_exist, t_cross = None, None
    for t in range(1, T + 1):
        lam = 0.02 * t
        s = 0.5 * lam if conv == "A" else 0.075 * lam / (1 - rho_)
        if t_exist is None and s > SSTAR:
            t_exist = t
        th = rho_ * th + (1 - rho_) * (g(th) + s)
        if t_cross is None and th > 0.5:
            t_cross = t
    return t_exist, t_cross, th


for conv in ("A", "B"):
    for r_ in (0.50, 0.70, 0.85, 0.95):
        te, tc, th = run(r_, conv)
        print(f"convention {conv}, rho = {r_:.2f}: equilibrium lost at t = {te}, theta crosses 1/2 at t = {tc}, theta_50 = {th:.3f}")
assert run(0.85, "A")[0] == 19 and run(0.85, "B")[0] == 19
assert run(0.50, "A")[1] < run(0.85, "A")[1]
assert run(0.95, "B")[0] < run(0.70, "B")[0] and run(0.50, "B")[0] is None

# ---------------------------------------------------------------------------
print("\n=== 24.4 Class-alliance vectors of Ch14 (illustrative), internal-stability sign check ===")
alliances = {"Oligarchic {F,B}": ((0.05, 0.03, -0.04, -0.04), "FB"),
             "National-populist {B,M,W}": ((-0.10, 0.02, 0.05, 0.03), "BMW"),
             "Techno-bureaucratic {F,B,M}": ((0.03, 0.04, 0.02, -0.09), "FBM"),
             "Left-populist {M,W}": ((-0.10, -0.02, 0.05, 0.07), "MW")}
idx = dict(F=0, B=1, M=2, W=3)
for k, (dp, mem) in alliances.items():
    assert abs(sum(dp)) < 1e-12 and all(dp[idx[m]] > 0 for m in mem)
    print(f"{k}: dpi = {dp}, sum = {sum(dp):+.2f}, all members gain")
periods = [("Tin oligarchy", None, 1952, 1.5), ("MNR", 1952, 1964, 2.0), ("Military", 1964, 1982, 1.8),
           ("Neoliberal", 1985, 2003, 2.5), ("MAS", 2006, 2019, 3.2)]
for nm, a, b, ne in periods:
    print(f"  {nm}: {a}-{b}, duration {None if a is None else b - a}, Neff^I {ne}")

# ---------------------------------------------------------------------------
print("\n=== 24.5 Shock needed to cross the separatrix, as a function of s ===")
need = {}
for s in (0.0, 0.05, 0.10, 0.15, 0.17, 0.18):
    f = fixed_points(s)
    need[s] = f[1] - f[0]
    print(f"s = {s:.2f}: theta_h = {f[0]:.4f}, theta_u = {f[1]:.4f}, required shock = {need[s]:.4f}")
vals = list(need.values())
assert all(np.diff(vals) < 0)

# ---------------------------------------------------------------------------
print("\n=== 24.6 Duration vs Neff^I (illustrative inputs) ===")


def corr_perm(d, ne):
    d, ne = np.asarray(d, float), np.asarray(ne, float)
    r = np.corrcoef(d, ne)[0, 1]
    perms = list(itertools.permutations(range(len(d))))
    rs = np.array([np.corrcoef(d, ne[list(p)])[0, 1] for p in perms])
    p_one = np.mean(rs <= r + 1e-12)
    b, a = np.polyfit(ne, d, 1)
    return r, p_one, b, a, len(perms)


LA_d = [40, 34, 21, 18, 13]          # oligarchic 1890-1930, nat-pop 1930-64, BA 1964-85, neolib 1985-2003, left-pop 2006-19
LA_n = [1.5, 2.0, 1.8, 2.5, 3.0]
r, p, b, a, npm = corr_perm(LA_d, LA_n)
print(f"Latin America (5 epochs): r = {r:.3f}, one-sided permutation p = {p:.3f} ({npm} permutations), "
      f"OLS duration = {a:.2f} + ({b:.2f}) Neff; at Neff = 3.5: {a + 3.5*b:.1f} years")
LAr = (r, p, b, a)
BO_d = [12, 18, 18, 13]
BO_n = [2.0, 1.8, 2.5, 3.2]
r2, p2, b2, a2, npm2 = corr_perm(BO_d, BO_n)
print(f"Bolivia (4 closed periods): r = {r2:.3f}, one-sided permutation p = {p2:.3f} ({npm2} permutations), "
      f"slope {b2:.2f}; at Neff = 3.5: {a2 + 3.5*b2:.1f} years")
yrs = np.array([1930, 1964, 1985, 2003, 2019])
print(f"Latin America: corr(duration, end year) = {np.corrcoef(LA_d, yrs)[0,1]:.3f}; corr(Neff, end year) = {np.corrcoef(LA_n, yrs)[0,1]:.3f}")
assert LAr[0] < -0.8 and r2 > -0.5

# ---------------------------------------------------------------------------
print("\n=== 24.7 Consent armoured by coercion: collapse time = max(T_consent, T_coercion) ===")
s0, ds = 0.05, 0.01
K0, Kbar = 1.0, 0.4
Tcons = (SSTAR - s0) / ds
for dK in (0.01, 0.02, 0.05):
    Tco = (K0 - Kbar) / dK
    dom = "Skocpol (coercion fails last)" if Tco > Tcons else "Gramsci (consent fails last)"
    print(f"s0 = {s0}, ds = {ds}/yr: T_consent = {Tcons:.1f}; K0 = {K0}, Kbar = {Kbar}, dK = {dK}/yr: "
          f"T_coercion = {Tco:.1f}; collapse at {max(Tcons, Tco):.1f} -> {dom}")
print("all chapter-24 solution checks passed")
