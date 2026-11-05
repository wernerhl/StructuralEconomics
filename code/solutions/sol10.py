"""Solutions to the Chapter 10 exercises (The Formal Formation Model).
Every number quoted in solutions/sol10.tex is printed here.
Run from the repository root:  python3 code/solutions/sol10.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def tau1(W):
    n = len(W)
    return 0.5 * max(np.abs(W[i] - W[j]).sum() for i in range(n) for j in range(n))


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def gini(v):
    v = np.sort(np.asarray(v, float))
    n = len(v)
    return np.sum((2 * np.arange(1, n + 1) - n - 1) * v) / (n * v.sum())


# ===========================================================================
print("== Exercise 10.1: affine formation map")
rho, mu, th0, rs = 0.8, 0.3, 1.0, 0.1
a1, b1 = 0.05, 0.1
slope = rho + mu * b1                       # 0.83
icpt = (1 - rho) * th0 + mu * (a1 - rs)     # 0.185
th_star = icpt / (1 - slope)
print(f"T(theta) = {slope:.2f} theta + {icpt:.3f};  theta* = {th_star:.6f};  J = {slope:.2f}")
assert abs(slope - 0.83) < 1e-12 and abs(icpt - 0.185) < 1e-12
assert abs(th_star - 0.185 / 0.17) < 1e-12
print(f"r(theta*) = {a1 + b1 * th_star:.6f}  (signal gap r - r* = {a1 + b1*th_star - rs:.6f})")
print(f"half-life of deviations: ln2/(-ln 0.83) = {np.log(2)/-np.log(slope):.3f} generations")
th = 0.5
path = [th]
for t in range(20):
    th = slope * th + icpt
    path.append(th)
path = np.array(path)
closed = th_star + slope ** np.arange(21) * (0.5 - th_star)
assert np.allclose(path, closed)
for t in [0, 1, 2, 3, 5, 10, 15, 20]:
    print(f"   t={t:2d}: theta_t = {path[t]:.4f}   gap to theta* = {th_star - path[t]:.4f}")
print(f"   share of initial gap left at t=20: 0.83^20 = {slope**20:.4f}")

# ===========================================================================
print("\n== Exercise 10.2: convex (quadratic) returns")
rho, mu, th0, rs = 0.85, 0.25, 1.0, 0.1
c0, c1, c2 = 0.05, 0.1, 0.05
# T(theta) = rho th + (1-rho) th0 + mu (c0 + c1 th + c2 th^2 - r*)
A2 = mu * c2
B2 = rho + mu * c1
C2 = (1 - rho) * th0 + mu * (c0 - rs)
print(f"T(theta) = {A2:.4f} theta^2 + {B2:.3f} theta + {C2:.4f}")
print(f"T(theta) - theta = {A2:.4f} theta^2 + ({B2-1:.3f}) theta + {C2:.4f} = 0.0125 (theta^2 - 10 theta + 11)")
roots = np.sort(np.roots([A2, B2 - 1, C2]).real)
print("fixed points:", roots, " (5 -+ sqrt(14) =", 5 - np.sqrt(14), 5 + np.sqrt(14), ")")
assert np.allclose(roots, [5 - np.sqrt(14), 5 + np.sqrt(14)])
Tp = lambda t: B2 + 2 * A2 * t
print("slopes T'(theta*) =", [round(Tp(t), 4) for t in roots])
assert Tp(roots[0]) < 1 < Tp(roots[1])
Tq = lambda t: A2 * t ** 2 + B2 * t + C2
for t0 in [0.9, 1.1, 5.0, 9.0]:
    t = t0
    traj = [t]
    for k in range(400):
        t = Tq(t)
        traj.append(t)
        if t > 1e6:
            break
    lab = f"-> {t:.4f}" if t < 1e6 else f"diverges (theta exceeds 1e6 after {k+1} periods)"
    print(f"   theta0 = {t0}: theta_1 = {traj[1]:.4f}, theta_5 = {traj[5]:.4f} ... {lab}")
# lower edge of the basin of the low fixed point (T is not monotone below its vertex)
vertex = -B2 / (2 * A2)
other = -B2 / A2 - roots[1]      # second solution of T(theta) = theta_H (sum of roots = -B2/A2)
assert abs(Tq(other) - roots[1]) < 1e-9
print(f"T has its minimum at theta = {vertex:.1f}; T(theta) = theta_H also at theta = {other:.3f};"
      f" basin of theta_L = ({other:.2f}, {roots[1]:.4f})")

# ===========================================================================
print("\n== Exercise 10.3: Jacobian of a 3-agent economy (reduced-form returns)")
WI3 = np.array([[0.5, 0.5, 0], [0.5, 0.5, 0], [0, 0, 1.0]])
rho3, mu3, b3 = 0.8, 0.2, 0.1
J3 = rho3 * np.eye(3) + mu3 * b3 * WI3
print("J =\n", J3)
ev = np.sort(np.linalg.eigvals(J3).real)[::-1]
print("eigenvalues:", ev, " rho(J) =", ev.max())
assert np.allclose(ev, [0.82, 0.82, 0.80])
J3b = rho3 * np.eye(3) + 2 * mu3 * b3 * WI3
ev2 = np.sort(np.linalg.eigvals(J3b).real)[::-1]
print("mu doubled (0.4): eigenvalues", ev2)
mu_crit = (1 - rho3) / b3
print("stability lost when rho + mu*b = 1: mu =", mu_crit)
assert abs(mu_crit - 2.0) < 1e-12
# comparison: the full model (allocation depends on theta_i - theta_bar, phi = 1, beta = 0.1, gamma = 0)
P3 = np.eye(3) - np.ones((3, 3)) / 3
Jfull = rho3 * np.eye(3) + mu3 * 1.0 * WI3 @ (b3 * np.eye(3)) @ P3
evf = np.linalg.eigvals(Jfull)
print("full-model Jacobian rho I + mu phi W^I C P, eigenvalues:", np.sort(evf.real)[::-1])
evs = np.sort(np.linalg.eigvals(WI3 @ P3).real)[::-1]
print("eigenvalues of W^I P:", evs, " -> full-model instability at mu =", (1 - rho3) / (b3 * evs[0]))

# ===========================================================================
print("\n== Exercise 10.4: three policies in the stylized Bolivia (Example ch10:ex:policies)")
N = 9
blk = np.repeat([0, 1, 2], 3)


def block_matrix(within_total, cross_each):
    return np.where(blk[:, None] == blk[None, :], within_total / 3, cross_each).astype(float)


WI = block_matrix(0.85, 0.025)
WP = block_matrix(0.94, 0.01)
alpha = np.repeat([0.06, 0.04, 0.03], 3)
beta = np.repeat([0.10, 0.07, 0.04], 3)
gam = 0.02
pi = np.repeat([0.15, 0.11, 0.22 / 3], 3)
omega = np.full(N, 0.5)
Phi, phi = 4.5, 1.0
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.1
P = np.eye(N) - np.ones((N, N)) / N


def model(WI=WI, WP=WP, omega=omega):
    C = np.diag(beta) + gam * WP
    J = rho * np.eye(N) + mu * phi * WI @ C @ P
    c = (1 - rho) * th0 + mu * (WI @ (alpha + C @ (omega + Phi * pi)) - rstar)
    return C, J, c


def outcomes(th, WP=WP, omega=omega):
    C = np.diag(beta) + gam * WP
    x = omega + Phi * pi + phi * P @ th
    return x, alpha + C @ x


C, J, c = model()
th = np.linalg.solve(np.eye(N) - J, c)
x, r = outcomes(th)
gap0 = th[0] - th[6]
xgap0 = x[0] - x[6]
print("(a) theta* by block", th[::3], " x* by block", x[::3], " r* by block", r[::3])
print(f"    rho(J) = {max(abs(np.linalg.eigvals(J))):.4f}; gap theta1-theta3 = {gap0:.4f}; "
      f"allocation gap x1-x3 = {xgap0:.4f}; sum x = {x.sum():.6f}")
assert abs(x.sum() - 9) < 1e-12 and np.allclose(J @ np.ones(N), rho * np.ones(N))
assert np.allclose(th, T := (rho * th + (1 - rho) * th0 + mu * (WI @ r - rstar)))  # fixed point
M = mu / (1 - rho)
assert np.allclose(th, th0 + M * (WI @ r - rstar))
print(f"    check theta_i = theta0 + M (W^I r - r*): M = {M:.4f}, max error {np.abs(th - th0 - M*(WI@r - rstar)).max():.1e}")
print(f"    pi sums to {pi.sum():.6f}; N_eff = {1/np.sum(pi**2):.3f}; tau1(W^I) = {tau1(WI):.3f}")

om_A = omega.copy()
om_A[blk == 0] -= 0.1
om_A[blk == 2] += 0.1
WP_B = WP.copy()
for i in range(N):
    for j in range(N):
        if {blk[i], blk[j]} == {0, 2}:
            WP_B[i, j] = 0.08
        elif blk[i] == blk[j] and blk[i] in (0, 2):
            WP_B[i, j] = 0.73 / 3
WI_C = WI.copy()
for i in np.where(blk == 2)[0]:
    WI_C[i] = np.where(blk == 0, 0.075, np.where(blk == 2, 0.70 / 3, 0.025))
assert np.allclose(WP_B.sum(1), 1) and np.allclose(WI_C.sum(1), 1) and abs(om_A.sum() - omega.sum()) < 1e-12

policies = {"A": (WI, WP, om_A), "B": (WI, WP_B, omega), "C": (WI_C, WP, omega),
            "ABC": (WI_C, WP_B, om_A)}
res = {}
for name, (wi, wp, om) in policies.items():
    xE, rE = outcomes(th, wp, om)                    # new structure, old preferences
    Cn, Jn, cn = model(wi, wp, om)
    thn = np.linalg.solve(np.eye(N) - Jn, cn)
    xn, rn = outcomes(thn, wp, om)
    th1 = Jn @ th + cn
    x1, r1 = outcomes(th1, wp, om)
    E, F, F1 = rE - r, rn - rE, r1 - rE
    Ex, Fx = xE - x, xn - xE
    sig = Jn @ th + cn - th                          # signal shock mu Delta s
    assert np.allclose(thn - th, np.linalg.solve(np.eye(N) - Jn, sig))
    dgap = (thn[0] - thn[6]) - gap0
    dxgap = (xn[0] - xn[6]) - xgap0
    res[name] = dict(dth=thn - th, dgap=dgap, dxgap=dxgap, E=E, F=F, F1=F1, Ex=Ex, Fx=Fx,
                     rhoJ=max(abs(np.linalg.eigvals(Jn))), sig=sig)
    sE = np.abs(E).sum()
    psi = np.abs(F).sum() / sE if sE > 1e-12 else np.nan
    psi1 = np.abs(F1).sum() / sE if sE > 1e-12 else np.nan
    res[name]["psi"], res[name]["psi1"] = psi, psi1
    print(f"(b) Policy {name}: rho(J_new) = {res[name]['rhoJ']:.4f}")
    print(f"    signal shock mu*Delta s by block {np.round(sig[::3], 5)}")
    print(f"    returns: E by block {np.round(E[::3], 4)}, F(inf) {np.round(F[::3], 4)}, F(1) {np.round(F1[::3], 5)}")
    print(f"    allocations: E by block {np.round(Ex[::3], 4)}, F(inf) {np.round(Fx[::3], 4)}")
    print(f"    sum|E| = {sE:.4f}, sum|F(inf)| = {np.abs(F).sum():.4f}, sum|F(1)| = {np.abs(F1).sum():.5f}, "
          f"Psi_bar(inf) = {psi:.3f}, Psi_bar(1) = {psi1:.3f}")
    print(f"    dtheta* by block {np.round(thn[::3]-th[::3], 4)}; d(theta1-theta3) = {dgap:.4f} "
          f"({-100*dgap/gap0:.1f}% of {gap0:.4f}); d(x1-x3) = {dxgap:.4f} ({-100*dxgap/xgap0:.1f}% of {xgap0:.4f})")
s3 = res["A"]["dgap"] + res["B"]["dgap"] + res["C"]["dgap"]
print(f"(d) sum of individual d(gap) = {s3:.5f}; combined = {res['ABC']['dgap']:.5f}; "
      f"interaction = {res['ABC']['dgap'] - s3:.6f}")
sx3 = res["A"]["dxgap"] + res["B"]["dxgap"] + res["C"]["dxgap"]
print(f"    sum of individual d(x-gap) = {sx3:.5f}; combined = {res['ABC']['dxgap']:.5f}; "
      f"interaction = {res['ABC']['dxgap'] - sx3:.6f}")
assert res["C"]["dgap"] < res["A"]["dgap"] < res["B"]["dgap"] < 0
assert res["A"]["dxgap"] < res["C"]["dxgap"]
assert abs(res["C"]["E"]).max() < 1e-14

# ===========================================================================
print("\n== Exercise 10.5: formation Leontief inverse")
J5 = np.array([[0.85, 0.05], [0.03, 0.90]])
ev5 = np.sort(np.linalg.eigvals(J5).real)[::-1]
print("trace", np.trace(J5), " det", np.linalg.det(J5), " eigenvalues", ev5,
      " (1.75 -+ sqrt(0.0085))/2 =", (1.75 + np.sqrt(0.0085)) / 2, (1.75 - np.sqrt(0.0085)) / 2)
assert ev5.max() < 1
print("||J||_inf =", np.abs(J5).sum(1).max())
inv = np.linalg.inv(np.eye(2) - J5)
print("det(I-J) =", np.linalg.det(np.eye(2) - J5), " (I-J)^{-1} =\n", inv)
assert np.allclose(inv, np.array([[200, 100], [60, 300]]) / 27)
S10 = sum(np.linalg.matrix_power(J5, k) for k in range(11))
print("I + J + ... + J^10 =\n", S10)
print("ratio partial/exact entrywise =\n", S10 / inv)
R = inv - S10
assert np.allclose(R, np.linalg.matrix_power(J5, 11) @ inv)
print("remainder J^11 (I-J)^{-1} =\n", R, "\n rho(J)^11 =", ev5.max() ** 11)
kneed = int(np.ceil(np.log(0.01) / np.log(ev5.max())))
print("terms needed for the dominant mode to fall below 1%:", kneed)
ds = np.array([0.01, 0.0])
dth = inv @ (0.2 * ds)
print("Delta theta* = (I-J)^{-1} mu Delta s =", dth, "; impact (round 0) =", 0.2 * ds)
assert np.allclose(dth, [0.2 * 0.01 * 200 / 27, 0.2 * 0.01 * 60 / 27])

# ===========================================================================
print("\n== Exercise 10.6: basin boundaries and policy design (Example ch10:ex:trap)")
rL, rH, cen, sl = 0.04, 0.22, 1.0, 25.0
rho_t, mu_t, th0_t, rs_t = 0.85, 0.25, 1.0, 0.12


def g(t):
    return rL + (rH - rL) / (1 + np.exp(-sl * (t - cen)))


def gp(t):
    e = np.exp(-sl * (t - cen))
    return (rH - rL) * sl * e / (1 + e) ** 2


def Tmap(t, da=0.0):
    return rho_t * t + (1 - rho_t) * th0_t + mu_t * (g(t) + da - rs_t)


def fixed_points(da=0.0):
    grid = np.linspace(0.0, 3.0, 300007)
    f = Tmap(grid, da) - grid
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    out = []
    for i in idx:
        lo, hi = grid[i], grid[i + 1]
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if np.sign(Tmap(lo, da) - lo) == np.sign(Tmap(mid, da) - mid):
                lo = mid
            else:
                hi = mid
        out.append(0.5 * (lo + hi))
    return out


fL, fU, fH = fixed_points()
print(f"(a) theta_L = {fL:.4f}, theta_U = {fU:.4f}, theta_H = {fH:.4f}; minimum delta = {fU - fL:.4f}")
print(f"    T'(theta_L) = {rho_t + mu_t*gp(fL):.4f}, T'(theta_U) = {rho_t + mu_t*gp(fU):.4f}; theta_H - theta_L = {fH - fL:.4f}")
for dl in [fU - fL - 0.002, fU - fL + 0.002]:
    t = fL + dl
    for _ in range(3000):
        t = Tmap(t)
    print(f"    push delta = {dl:.4f} -> long run {t:.4f}")
# (b) threshold derivative and finite shifts
dthU = -mu_t / (rho_t + mu_t * gp(fU) - 1)
dthL = mu_t / (1 - rho_t - mu_t * gp(fL))
print(f"(b) d theta_U / d alpha = -mu/(T'(theta_U)-1) = {dthU:.3f};  d theta_L / d alpha = {dthL:.3f};"
      f"  d(gap)/d alpha = {dthU - dthL:.3f}")
for da in [0.005, 0.01]:
    fp = fixed_points(da)
    print(f"    Delta alpha = {da}: fixed points {np.round(fp, 4)}; theta_U moves {fp[1]-fU:+.4f}; "
          f"new gap {fp[1]-fp[0]:.4f}")
# critical alpha shift at which theta_L and theta_U merge (saddle-node)
lo, hi = 0.0, 0.1
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if len(fixed_points(mid)) >= 3:
        lo = mid
    else:
        hi = mid
da_c = hi
print(f"    the low steady state disappears for Delta alpha > {da_c:.4f}")
# (c) costs: one-time push vs permanent vs temporary return subsidy
one_period = (fU - fL) / mu_t
print(f"(c) one-period return subsidy needed: (theta_U - theta_L)/mu = {one_period:.4f}")


def escapes(da, n):
    t = fL
    for k in range(n):
        t = Tmap(t, da)
    for _ in range(3000):
        t = Tmap(t)
    return t > fU


best = []
for n in [1, 2, 3, 5, 10, 20, 40, 80]:
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if escapes(mid, n):
            hi = mid
        else:
            lo = mid
    best.append((n, hi, n * hi))
    print(f"    temporary subsidy for n = {n:2d} periods: minimum Delta alpha = {hi:.4f}, total cost n*Delta alpha = {n*hi:.4f}")
assert best[0][2] == min(b[2] for b in best)
for dlt in [0.9, 0.95]:
    print(f"    permanent subsidy Delta alpha_c = {da_c:.4f}: PV at discount factor {dlt} = {da_c/(1-dlt):.4f}")

# ===========================================================================
print("\n== Exercise 10.7: Cardoso-Faletto situations (four agents: F, B, M, W)")
N4 = 4
P4 = np.eye(N4) - np.ones((N4, N4)) / N4
sit = {
    "Enclave": (
        np.array([[0.85, 0.05, 0.05, 0.05],
                  [0.60, 0.30, 0.05, 0.05],
                  [0.50, 0.10, 0.30, 0.10],
                  [0.60, 0.05, 0.05, 0.30]]),
        np.array([[0.70, 0.20, 0.05, 0.05],
                  [0.40, 0.50, 0.05, 0.05],
                  [0.05, 0.05, 0.60, 0.30],
                  [0.02, 0.03, 0.30, 0.65]])),
    "National populism": (
        np.array([[0.60, 0.20, 0.10, 0.10],
                  [0.05, 0.45, 0.20, 0.30],
                  [0.05, 0.30, 0.40, 0.25],
                  [0.05, 0.35, 0.20, 0.40]]),
        np.array([[0.70, 0.10, 0.10, 0.10],
                  [0.05, 0.35, 0.30, 0.30],
                  [0.05, 0.30, 0.35, 0.30],
                  [0.05, 0.30, 0.30, 0.35]])),
    "Associated dependent": (
        np.array([[0.50, 0.40, 0.05, 0.05],
                  [0.40, 0.45, 0.10, 0.05],
                  [0.20, 0.30, 0.40, 0.10],
                  [0.15, 0.30, 0.15, 0.40]]),
        np.array([[0.60, 0.35, 0.025, 0.025],
                  [0.40, 0.50, 0.05, 0.05],
                  [0.20, 0.30, 0.40, 0.10],
                  [0.02, 0.03, 0.15, 0.80]])),
}
al4 = np.array([0.06, 0.05, 0.04, 0.03])
be4 = np.array([0.10, 0.08, 0.06, 0.04])
ga4, om4, Phi4, phi4 = 0.02, 0.5, 2.0, 1.0
rho4, mu4, th04, rs4 = 0.85, 0.2, 1.0, 0.1
out = {}
for name, (wp, wi) in sit.items():
    assert np.allclose(wp.sum(1), 1) and np.allclose(wi.sum(1), 1)
    Wbar = 0.5 * wp + 0.5 * wi
    p = perron(Wbar)
    C4 = np.diag(be4) + ga4 * wp
    J4 = rho4 * np.eye(N4) + mu4 * phi4 * wi @ C4 @ P4
    c4 = (1 - rho4) * th04 + mu4 * (wi @ (al4 + C4 @ (om4 + Phi4 * p)) - rs4)
    t4 = np.linalg.solve(np.eye(N4) - J4, c4)
    x4 = om4 + Phi4 * p + phi4 * P4 @ t4
    r4 = al4 + C4 @ x4
    assert abs(x4.sum() - 4) < 1e-12 and np.all(x4 > 0) and max(abs(np.linalg.eigvals(J4))) < 1
    out[name] = dict(pi=p, th=t4, x=x4)
    print(f"{name}: pi = {np.round(p, 3)} (N_eff {1/np.sum(p**2):.2f}); rho(J) = {max(abs(np.linalg.eigvals(J4))):.3f}")
    print(f"   theta* = {np.round(t4, 3)}; r* = {np.round(r4, 3)}")
    print(f"   x* = {np.round(x4, 3)}; x_max/x_min = {x4.max()/x4.min():.2f}; Gini(x) = {gini(x4):.3f}; "
          f"Var(x) = {np.var(x4):.4f}; power part Phi*(pi-1/4) = {np.round(Phi4*(p-0.25),3)}; "
          f"preference part = {np.round(phi4*P4@t4,3)} (range {np.ptp(phi4*P4@t4):.3f})")
order = sorted(out, key=lambda k: gini(out[k]["x"]))
print("most equal -> least equal (Gini of x*):", order)
assert order[0] == "National populism"
