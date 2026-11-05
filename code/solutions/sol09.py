"""Solutions manual, Chapter 9 (Where Preferences Come From).
Computes and prints every number used in solutions/sol09.tex.
Run from the repository root:  python3 code/solutions/sol09.py      (numpy only)"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)

# ===================================================================== 9.1
print("=" * 72); print("Exercise 9.1: formation equation by hand")


def path_91(rho, theta0=0.5, T=60, thbar=1.0, mu=0.2, r=0.15, rstar=0.1):
    th = [theta0]
    for _ in range(T):
        th.append(rho * th[-1] + (1 - rho) * thbar + mu * (r - rstar))
    return np.array(th)


for rho in (0.8, 0.95):
    c = (1 - rho) * 1.0 + 0.2 * 0.05
    ss = c / (1 - rho)
    th = path_91(rho)
    print(f"rho={rho}: theta_(t+1) = {rho} theta_t + {c:.4f}; steady state = {ss:.4f} "
          f"(= 1 + mu/(1-rho)*0.05 = {1 + 0.2 / (1 - rho) * 0.05:.4f})")
    print("   theta_1..10 =", np.round(th[1:11], 4))
    assert np.allclose(th, ss - (ss - 0.5) * rho ** np.arange(len(th)))
    print("   0.9 theta* =", round(0.9 * ss, 4))
    t90lev = int(np.argmax(th >= 0.9 * ss))
    t90gap = int(np.argmax(np.abs(th - ss) <= 0.1 * (ss - 0.5) + 1e-15))
    print(f"   first t with theta_t >= 0.9 theta* : {t90lev} (exact solution of {ss-0.5:.2f} rho^t = {0.1*ss:.3f}: "
          f"t = {np.log(0.1 * ss / (ss - 0.5)) / np.log(rho):.2f})")
    print(f"   first t closing 90% of the gap    : {t90gap} (t = ln 0.1/ln rho = {np.log(0.1) / np.log(rho):.2f})")
    print(f"   half-life = {np.log(2) / np.log(1 / rho):.2f} generations;  theta_10 gap closed = {1 - rho ** 10:.4f}")
    print("   TikZ coords:", " ".join("(%d,%.4f)" % (t, th[t]) for t in range(0, 61, 2)))

# ===================================================================== 9.2
print("=" * 72); print("Exercise 9.2: heterogeneous steady states (Example ch09:ex:two_steady)")
rho, mu, rstar, thbar = 0.8, 0.3, 0.1, 1.0
kap = mu / (1 - rho)
a = np.array([0.05, 0.03]); G = np.diag([0.1, 0.02])


def steady(W):
    M = np.eye(2) - kap * W @ G
    return np.linalg.solve(M, thbar * np.ones(2) + kap * (W @ a - rstar))


def Ww(w):
    return np.array([[1 - w, w], [w, 1 - w]])


th0 = steady(np.eye(2))
print("map slopes rho + mu g'_M =", rho + mu * 0.1, ", rho + mu g'_T =", rho + mu * 0.02)
print("(a) W = I: theta* =", np.round(th0, 4), " exact 0.925/0.85 =", round(0.925 / 0.85, 4), ", 0.895/0.97 =", round(0.895 / 0.97, 4),
      " gap =", round(th0[0] - th0[1], 4))
assert np.allclose(th0, [0.925 / 0.85, 0.895 / 0.97])
for w in (0.3, 0.4, 0.5, 0.6, 1.0):
    th = steady(Ww(w)); r = a + G @ th
    gapf = kap * (1 - 2 * w) * (r[0] - r[1])
    print(f"   w={w}: theta* = {np.round(th, 4)}, gap = {th[0]-th[1]:+.4f}, returns = {np.round(r, 4)}, "
          f"kappa(1-2w)(r_M - r_T) = {gapf:+.4f}")
    assert np.isclose(th[0] - th[1], gapf)
th3 = steady(Ww(0.3))
print("(b) w = 0.3: theta* =", np.round(th3, 4), " signals", np.round(Ww(0.3) @ (a + G @ th3), 4))
J = rho * np.eye(2) + mu * Ww(0.3) @ G
print("    local rate lambda = rho(J) =", round(max(abs(np.linalg.eigvals(J))), 4))
print("    mean theta*: W=I", round(th0.mean(), 4), " w=0.3", round(th3.mean(), 4), " w=0.5", round(steady(Ww(0.5)).mean(), 4))

# ===================================================================== 9.3
print("=" * 72); print("Exercise 9.3: persistence in female labor supply")
s3 = 0.85 + 0.1 * 0.3
c3 = 0.15 * 0.5 - 0.1 * 0.1
ss3 = c3 / (1 - s3)
print(f"theta' = {s3:.2f} theta + {c3:.3f};  steady state = {c3:.3f}/{1 - s3:.2f} = {ss3:.5f};  half-life = {np.log(2)/np.log(1/s3):.3f} generations")
base = [0.2]
for t in range(50):
    base.append(s3 * base[-1] + c3)
base = np.array(base)
war = base.copy(); war[10] = 0.6
for t in range(10, 50):
    war[t + 1] = s3 * war[t] + c3
print("closed forms: no war theta_t = %.4f - %.4f*0.88^t; war theta_t = %.4f + %.4f*0.88^(t-10)" % (ss3, ss3 - 0.2, ss3, 0.6 - ss3))
print("anchor term 0.15*0.5 =", 0.15 * 0.5, "; rho without the network effect: half-life", round(np.log(2) / np.log(1 / 0.85), 3))
print("no war: theta_10 =", round(base[10], 4), " theta_20 =", round(base[20], 4), " theta_50 =", round(base[50], 4))
print("war   : theta_10 = 0.6, theta_11 =", round(war[11], 4), " theta_20 =", round(war[20], 4), " theta_30 =", round(war[30], 4),
      " theta_50 =", round(war[50], 4))
print("difference war - no war at t=10,15,20,30,50:", np.round((war - base)[[10, 15, 20, 30, 50]], 5),
      " = (0.6 - theta_10) * 0.88^(t-10)")
assert np.allclose((war - base)[10:], (0.6 - base[10]) * s3 ** np.arange(41))
tback = int(np.argmax((np.arange(51) > 10) & (war <= base[10] + 1e-12)))
print("does the post-war path fall back to theta_10 (pre-war level 0.4465)? min post-war =", round(war[10:].min(), 4))
print("TikZ coords base:", " ".join("(%d,%.4f)" % (t, base[t]) for t in range(0, 51, 2)))
print("TikZ coords war :", " ".join("(%d,%.4f)" % (t, war[t]) for t in [10] + list(range(11, 51, 2))))

# ===================================================================== 9.4
print("=" * 72); print("Exercise 9.4: aspiration and the smartphone")
thM, thT = th0
rM, rT = a + G @ th0
print("sector steady states (W = I): theta_M =", round(thM, 4), " theta_T =", round(thT, 4), " r_M =", round(rM, 5), " r_T =", round(rT, 5))
# individual: observes sector returns, which her own theta does not move
th_old = thbar + kap * (rT - rstar)
th_new = thbar + kap * (0.7 * rT + 0.3 * rM - rstar)
print("(a) individual: old steady state", round(th_old, 4), " new", round(th_new, 5), " gain", round(th_new - th_old, 5))
print("    share of the M-T aspiration gap closed in the long run =", round((th_new - th_old) / (thM - thT), 4), "(= w = 0.3)")
assert np.isclose((th_new - th_old) / (thM - thT), 0.3)
traj = [0.8]
for t in range(30):
    sig = rT if t < 5 else 0.7 * rT + 0.3 * rM
    traj.append(rho * traj[-1] + (1 - rho) * thbar + mu * (sig - rstar))
traj = np.array(traj)
print("(b) theta_t, t=0..12:", np.round(traj[:13], 4))
gap5 = th_new - traj[5]
half = int(np.argmax((np.arange(len(traj)) >= 5) & (th_new - traj <= gap5 / 2)))
print(f"    gap at t=5: {gap5:.4f}; half closed by t = {half} ({half-5} generations after the phone); "
      f"continuous half-life ln2/ln(1/0.8) = {np.log(2)/np.log(1/rho):.3f}")
assert np.allclose(th_new - traj[5:], gap5 * rho ** np.arange(len(traj) - 5))
# gap measured from the M-T gap: never half with w = 0.3; w needed = 0.5
print("    weight needed to close half of the M-T gap: w = 0.5")
# (c) many agents: share s of the T sector treated, sector return depends on sector mean theta
rMf = rM


def sector(s):
    """unknowns (theta_tr, theta_un); r_T = 0.03 + 0.02 (s th_tr + (1-s) th_un)."""
    A = np.array([[1 - kap * 0.7 * 0.02 * s, -kap * 0.7 * 0.02 * (1 - s)],
                  [-kap * 0.02 * s, 1 - kap * 0.02 * (1 - s)]])
    b = np.array([thbar + kap * (0.7 * 0.03 + 0.3 * rMf - rstar),
                  thbar + kap * (0.03 - rstar)])
    return np.linalg.solve(A, b)


ind = th_new - th_old
print("(c) individual effect (negligible share):", round(ind, 5))
for s in (0.001, 0.1, 0.5, 1.0):
    ttr, tun = sector(s)
    total = s * (ttr - thT) + (1 - s) * (tun - thT)
    print(f"    share treated s={s}: treated theta {ttr:.5f}, untreated {tun:.5f}; aggregate effect per treated agent "
          f"= {total / s:.5f} = {total / s / ind:.4f} x individual")
ttr1, _ = sector(1.0)
lamT = rho + mu * 0.7 * 0.02
print("    whole sector treated: theta_T =", round(ttr1, 5), " multiplier 1/(1 - kappa*0.7*0.02) =", round(1 / (1 - kap * 0.7 * 0.02), 4),
      " local rate =", round(lamT, 4), " half-life", round(np.log(2) / np.log(1 / lamT), 3))
# crowding: if modern returns observed fall as more traditional agents aim at them, e.g. r_M_obs = r_M - 0.05 s
for s in (0.1, 0.5, 1.0):
    rMf = rM - 0.05 * s
    ttr, tun = sector(s)
    total = s * (ttr - thT) + (1 - s) * (tun - thT)
    print(f"    with crowding r_M,obs = r_M - 0.05 s, s={s}: per-treated effect {total / s:.5f} = {total / s / ind:.4f} x individual")
rMf = rM

# ===================================================================== 9.5
print("=" * 72); print("Exercise 9.5: Bernard et al. as a W^I perturbation")
# (c) time aggregation: n = 50 six-month sub-periods per 25-year generation
n = 50
for r6 in (1.0, 0.99, 0.95, 0.9):
    factor = n if r6 == 1.0 else (1 - r6 ** n) / (1 - r6)
    print(f"   rho_6m = {r6}: mu_gen / mu_6m = {factor:.2f};  rho_gen = rho_6m^50 = {r6 ** n:.4f}")
print("   kappa invariant: mu_gen/(1-rho_gen) = mu_6m/(1-rho_6m)")
for r6 in (0.99, 0.95, 0.9):
    m6 = 1.0
    mg = m6 * (1 - r6 ** n) / (1 - r6); rg = r6 ** n
    assert np.isclose(mg / (1 - rg), m6 / (1 - r6))
# numeric illustration of (c): implied mu_6m = 0.1 / (w * dr) in SD units of theta
for w, dr in ((0.05, 1.0), (0.1, 1.0), (0.1, 2.0)):
    print(f"   w = {w}, (r_film - rbar) = {dr} SD of theta: implied mu_6m = {0.1 / (w * dr):.3f}")
# (d) convex savings: elasticity >= 1
for th0_, dth in ((1.0, 0.1),):
    for nm, f in (("theta^2", lambda t: t ** 2), ("exp(theta)-1", lambda t: np.exp(t) - 1), ("linear 2 theta", lambda t: 2 * t)):
        print(f"   f = {nm}: 10% rise in theta raises savings by {100 * (f(th0_ * 1.1) / f(th0_) - 1):.1f}%")

# ===================================================================== 9.6
print("=" * 72); print("Exercise 9.6: the Psi ratio for education")
c, g = 0.10, 0.05
E = 1 / (1 - c) - 1
F = (1 + g) / (1 - c) - 1 / (1 - c)
print(f"(a) E/x0 = 1/0.9 - 1 = {E:.4f}  (first-order approx 0.10)")
print(f"(b) F/x0 = 1.05/0.9 - 1/0.9 = {F:.4f}")
print(f"(c) Psi = F/E = {F / E:.4f} = g/c = {g / c:.4f};  total (E+F)/x0 = {E + F:.4f} = 1.05/0.9 - 1;  "
      f"share missed by ignoring formation = Psi/(1+Psi) = {F / E / (1 + F / E):.4f}")
assert np.isclose(F / E, g / c) and np.isclose(F / E, 0.5)
lam = 0.8
ginf = g / (1 - lam)
print(f"(d) if theta follows lambda = {lam} from a steady state: g(T) = g_inf (1 - lambda^T), g(1) = 0.05 => g_inf = {ginf:.2f};")
for T in (1, 2, 3, 5, 10):
    print(f"    T={T}: g(T) = {ginf * (1 - lam ** T):.4f}, Psi = {ginf * (1 - lam ** T) / c:.4f}")
Tstar = np.log(1 - c / ginf) / np.log(lam)
print(f"    Psi > 1 once g(T) > c, i.e. T > ln(1 - c/g_inf)/ln(lambda) = {Tstar:.3f} generations; Psi_inf = {ginf / c:.2f}")

# ===================================================================== 9.7
print("=" * 72); print("Exercise 9.7: Furtado's demonstration effect")
rho7, mu7, thbar7, rstar7 = 0.9, 0.4, 0.70, 0.70
kap7 = mu7 / (1 - rho7)
cC, cP = 0.85, 0.70
for w in (0.4,):
    sig = w * cC + (1 - w) * cP
    th7 = thbar7 + kap7 * (sig - rstar7)
    print(f"(b) signal = {w}*0.85 + {1-w}*0.70 = {sig:.4f}; theta* = 0.70 + {kap7:.0f}*({sig:.2f} - 0.70) = {th7:.4f}")
    assert np.isclose(th7, 0.94)
wstar = (1 - thbar7) / (kap7 * (cC - cP))
print(f"(c) theta*(w) = 0.70 + {kap7 * (cC - cP):.2f} w > 1  iff  w > {wstar:.4f}")
print(f"    half-life of adjustment = {np.log(2) / np.log(1 / rho7):.2f} generations (slope rho = 0.9, signal fixed)")
traj = [0.70]
for t in range(60):
    traj.append(rho7 * traj[-1] + (1 - rho7) * thbar7 + mu7 * (0.4 * cC + 0.6 * cP - rstar7))
print("    path from 0.70: t=1,5,10,20,60:", np.round(np.array(traj)[[1, 5, 10, 20, 60]], 4))
# variant: periphery signal = elite's own consumption ratio (feedback)
slope = rho7 + mu7 * 0.6
ss_fb = (thbar7 * (1 - rho7) + mu7 * (0.4 * cC - rstar7)) / (1 - slope)
print(f"    variant with own-ratio feedback: slope = {slope:.2f} > 1; fixed point {ss_fb:.4f} is unstable")
print("all assertions passed")
