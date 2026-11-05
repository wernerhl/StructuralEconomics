"""Chapter 10: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch10.py
numpy only."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def tau1(W):
    """Dobrushin coefficient 1/2 max_{i,j} sum_k |W_ik - W_jk|."""
    n = len(W)
    return 0.5 * max(np.abs(W[i] - W[j]).sum() for i in range(n) for j in range(n))


def neff(p):
    return 1 / np.sum(p ** 2)


# ---------------------------------------------------------------------------
# Example ch10:ex:bolivia_three  (all parameters)
N = 9
blk = np.repeat([0, 1, 2], 3)


def block_matrix(within_total, cross_each):
    W = np.where(blk[:, None] == blk[None, :], within_total / 3, cross_each)
    return W.astype(float)


WI = block_matrix(0.85, 0.025)        # information network
WP = block_matrix(0.94, 0.01)         # production input shares
alpha = np.repeat([0.06, 0.04, 0.03], 3)
beta = np.repeat([0.10, 0.07, 0.04], 3)
gam = 0.02
pi = np.repeat([0.15, 0.11, 0.22 / 3], 3)     # power, per agent; sums to 1
omega = np.full(N, 0.5)
Phi = 4.5                                      # bargaining surplus
phi = 1.0                                      # preference-to-allocation sensitivity
rho, mu, th0, rstar = 0.85, 0.2, 1.0, 0.1
P = np.eye(N) - np.ones((N, N)) / N
assert np.allclose(WI.sum(1), 1) and np.allclose(WP.sum(1), 1) and abs(pi.sum() - 1) < 1e-12


def model(WI=WI, WP=WP, omega=omega):
    C = np.diag(beta) + gam * WP
    J = rho * np.eye(N) + mu * phi * WI @ C @ P
    c = (1 - rho) * th0 + mu * (WI @ (alpha + C @ (omega + Phi * pi)) - rstar)
    return C, J, c


def outcomes(th, WP=WP, omega=omega):
    C = np.diag(beta) + gam * WP
    x = omega + Phi * pi + phi * P @ th
    r = alpha + C @ x
    return x, r


C, J, c = model()
th = np.linalg.solve(np.eye(N) - J, c)
x, r = outcomes(th)
print("== Example bolivia_three")
print("theta* by block:", th[::3])
print("x* by block    :", x[::3], " sum x =", x.sum(), " Omega/N =", x.mean())
print("r* by block    :", r[::3])
print("pi per agent   :", pi[::3], " N_eff(pi) =", neff(pi), " Var(pi) =", np.mean((pi - 1 / N) ** 2))
print("x/(Omega/N)    :", x[::3] / x.mean())
gap0 = th[0] - th[6]
print("preference gap gas-rural:", gap0)
print("rho(J) =", max(abs(np.linalg.eigvals(J))))
# J 1 = rho 1  (uniform shifts are not amplified)
assert np.allclose(J @ np.ones(N), rho * np.ones(N))
# norm bound for stability
bound = rho + mu * phi * 1 * max(beta + gam) * 2 * (N - 1) / N
print("norm bound rho + mu*phi*max(beta+gamma)*2(N-1)/N =", bound)
assert max(abs(np.linalg.eigvals(J))) <= np.abs(J).sum(1).max() <= bound + 1e-12
# Proposition two_types bound |theta_i - theta_j| <= M tau1(W^I)(max r - min r)
Mf = mu / (1 - rho)
print("M =", Mf, " tau1(W^I) =", tau1(WI), " bound =", Mf * tau1(WI) * (r.max() - r.min()))
assert abs(th.max() - th.min()) <= Mf * tau1(WI) * (r.max() - r.min()) + 1e-12
# what the old table values were: pure block diagonal W^I with the OLD returns
r_old = np.repeat([0.182, 0.112, 0.063], 3)
print("(review check) old table theta under pure block W^I:", (1 + (0.2 / 0.15) * (block_matrix(1, 0) @ r_old - 0.1))[::3],
      " under stated W^I:", (1 + (0.2 / 0.15) * (WI @ r_old - 0.1))[::3])
# heatmap colour for W^I panel (normalized to matrix max)
print("W^I heatmap colours: within", round(100 * WI.max() / WI.max()), " cross", round(100 * 0.025 / WI.max()))

# ---------------------------------------------------------------------------
# Example ch10:ex:policies  (outcome y = returns r; Definition ch09:def:psi)
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
assert np.allclose(WP_B.sum(1), 1) and np.allclose(WI_C.sum(1), 1)

print("\n== Example policies (outcome = returns)")
for name, (wi, wp, om) in {"A transfer": (WI, WP, om_A), "B structural": (WI, WP_B, omega),
                            "C information": (WI_C, WP, omega)}.items():
    _, rE = outcomes(th, wp, om)
    E = rE - r
    Cn, Jn, cn = model(wi, wp, om)
    thn = np.linalg.solve(np.eye(N) - Jn, cn)
    _, rn = outcomes(thn, wp, om)
    F = rn - rE
    # horizon T = 1
    th1 = Jn @ th + cn
    _, r1 = outcomes(th1, wp, om)
    F1 = r1 - rE
    # formation Leontief inverse check: dtheta = (I - J_new)^{-1} (T_new(theta0) - theta0)
    assert np.allclose(thn - th, np.linalg.solve(np.eye(N) - Jn, Jn @ th + cn - th))
    dgap = (thn[0] - thn[6]) - gap0
    sE, sF = np.abs(E).sum(), np.abs(F).sum()
    print(f"{name}: dtheta by block {np.round((thn-th)[::3],4)}, d(gap) = {dgap:.4f} "
          f"({-100*dgap/gap0:.1f}% of gap)")
    print(f"   E by block {np.round(E[::3],4)}  F(inf) by block {np.round(F[::3],4)}")
    if sE > 1e-12:
        print(f"   aggregate Psi(inf) = {sF/sE:.3f},  Psi(T=1) = {np.abs(F1).sum()/sE:.3f}")
    else:
        print(f"   E = 0: Psi undefined;  sum|F(inf)| = {sF:.4f}, sum|F(1)| = {np.abs(F1).sum():.4f}")
    print(f"   rho(J_new) = {max(abs(np.linalg.eigvals(Jn))):.4f}")

# ---------------------------------------------------------------------------
# Proposition psi_onelop: Psi(T) = k (1 - lam^T)/(1 - lam),  k = mu beta eta, lam = rho + k
print("\n== One-loop Psi formula check (simulation)")
for (rh, m, be, et) in [(0.85, 0.2, 0.1, 1.0), (0.62, 0.18, 0.5, 1.0)]:
    k = m * be * et
    lam = rh + k
    th_ = 0.0
    de = 1.0
    for T in range(1, 6):
        th_ = lam * th_ + m * be * de      # deviation dynamics after the shock
        assert abs(et * th_ / de - k * (1 - lam ** T) / (1 - lam)) < 1e-12
    print(f"rho={rh}, mu={m}, beta*eta={be*et}: Psi(inf) = {k/(1-lam):.4f}, diverges at rho = {1-k:.3f}")

# ---------------------------------------------------------------------------
# Example ch10:ex:trap  (S-shaped returns)
rL, rH, cen, s = 0.04, 0.22, 1.0, 25.0
rho_t, mu_t, th0_t, rs_t = 0.85, 0.25, 1.0, 0.12


def g(t):
    return rL + (rH - rL) / (1 + np.exp(-s * (t - cen)))


def gp(t):
    e = np.exp(-s * (t - cen))
    return (rH - rL) * s * e / (1 + e) ** 2


def Tmap(t, m=mu_t, rs=rs_t):
    return rho_t * t + (1 - rho_t) * th0_t + m * (g(t) - rs)


def fixed_points(m=mu_t, rs=rs_t):
    grid = np.linspace(0.0, 3.0, 300007)
    f = Tmap(grid, m, rs) - grid
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    out = []
    for i in idx:
        lo, hi = grid[i], grid[i + 1]
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if np.sign(Tmap(lo, m, rs) - lo) == np.sign(Tmap(mid, m, rs) - mid):
                lo = mid
            else:
                hi = mid
        out.append(0.5 * (lo + hi))
    return out


fp = fixed_points()
print("\n== Example trap")
print("fixed points:", np.round(fp, 4), " slopes:", np.round([rho_t + mu_t * gp(t) for t in fp], 4))
assert len(fp) == 3
for t0 in [0.8, 0.95, 0.97, 0.99, 1.05, 1.3]:
    t = t0
    for _ in range(2000):
        t = Tmap(t)
    print(f"  start {t0} -> {t:.4f}")
print("aspiration gap theta_U - theta_L =", fp[1] - fp[0])
print("returns at steady states:", np.round([g(t) for t in fp], 4))
print("contraction (uniqueness) bound on mu: (1-rho)/max g' =", (1 - rho_t) / ((rH - rL) * s / 4))
# saddle-node: smallest mu with three fixed points
lo, hi = 0.13, 0.25
for _ in range(60):
    mid = 0.5 * (lo + hi)
    if len(fixed_points(mid)) >= 3:
        hi = mid
    else:
        lo = mid
print("saddle-node threshold mu_c ~", round(hi, 4), " i.e. M_c =", round(hi / (1 - rho_t), 4))
# Exercise 10.6(b): raising alpha (shift all returns up by d) moves the threshold down
d = 0.005
fp2 = fixed_points(rs=rs_t - d)
print("10.6(b) alpha + 0.005: fixed points", np.round(fp2, 4), " threshold moves by", fp2[1] - fp[1])

# ---------------------------------------------------------------------------
print("\n== Exercises")
# 10.1
print("10.1 Phi(theta) = 0.83 theta + 0.185; theta* =", 0.185 / 0.17, " J = 0.83")
# 10.2 (quadratic returns: at most two fixed points)
cq = [0.25 * 0.05, 0.85 + 0.25 * 0.1 - 1, 0.15 + 0.25 * (0.05 - 0.1)]
roots = np.sort(np.roots(cq).real)
print("10.2 fixed points", roots, " slopes", [0.85 + 0.25 * (0.1 + 0.1 * t) for t in roots])
# 10.3
WI3 = np.array([[0.5, 0.5, 0], [0.5, 0.5, 0], [0, 0, 1.0]])
J3 = 0.8 * np.eye(3) + 0.2 * 0.1 * WI3
print("10.3 eigenvalues", np.sort(np.linalg.eigvals(J3).real)[::-1], " instability at mu =", (1 - 0.8) / 0.1)
# 10.5
J5 = np.array([[0.85, 0.05], [0.03, 0.90]])
inv = np.linalg.inv(np.eye(2) - J5)
S = sum(np.linalg.matrix_power(J5, k) for k in range(11))
print("10.5 eig", np.linalg.eigvals(J5), "\n inverse", inv, "\n partial sum", S, "\n dtheta", inv @ (0.2 * np.array([0.01, 0])))
# Appendix 10B pitchfork variant: r* = 0.13 (logistic midpoint) -> symmetric about theta = 1
mu_pf = (1 - rho_t) / ((rH - rL) * s / 4)
print("\n10B pitchfork variant: mu where T'(1)=1:", mu_pf,
      " fixed points at mu=0.12:", np.round(fixed_points(0.12, 0.13), 4),
      " at mu=0.16:", np.round(fixed_points(0.16, 0.13), 4))
