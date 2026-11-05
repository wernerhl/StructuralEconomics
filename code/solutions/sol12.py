"""Solutions to the Chapter 12 exercises (Welfare Under Endogenous Preferences).
Every number quoted in solutions/sol12.tex is printed here.
Run from the repository root:  python3 code/solutions/sol12.py   (numpy only)."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def v(z, kappa=0.5):
    return z - 0.5 * kappa * z ** 2


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def gap(W):
    ev = np.sort(np.abs(np.linalg.eigvals(W)))[::-1]
    return 1 - ev[1]


N = 9
blk = np.repeat([0, 1, 2], 3)
xb = np.array([1.2486, 0.9852, 0.7662])      # Example ch10:ex:bolivia_three (code/ch10.py)
x = np.repeat(xb, 3)
WI = np.where(blk[:, None] == blk[None, :], 0.85 / 3, 0.025)
Wce = np.ones((N, N)) / N
kap = 0.5

# ===========================================================================
print("== Exercise 12.1: satisfaction paradox")
gaps = {}
for M in [0.0, 0.4, 0.8]:
    a = M * WI @ x
    z = x - a
    assert np.all(kap * z < 1)
    s = v(z, kap)
    ai = M * Wce @ x
    si = v(x - ai, kap)
    gaps[M] = (s[0] - s[6], si[0] - si[6])
    print(f"M = {M}: a_b = {np.round(a[::3],4)}, s_b = {np.round(s[::3],4)}, gap = {s[0]-s[6]:.4f}; "
          f"integrated s_b = {np.round(si[::3],4)}, gap = {si[0]-si[6]:.4f}; ratio = {(s[0]-s[6])/(si[0]-si[6]):.4f}")
print("gap curve (M, block-W^I gap, integrated gap):",
      " ".join(f"({M:.1f},{v(x-M*WI@x)[0]-v(x-M*WI@x)[6]:.4f},{v(x-M*Wce@x)[0]-v(x-M*Wce@x)[6]:.4f})"
               for M in np.arange(0, 0.91, 0.1)))
# slope of the gap in M (block W^I) is linear? check second differences
gm = np.array([v(x - M * WI @ x)[0] - v(x - M * WI @ x)[6] for M in np.linspace(0, 0.9, 10)])
print("second differences of gap(M):", np.round(np.diff(gm, 2), 6))
cq = np.polyfit(np.linspace(0, 0.9, 10), gm, 2)
print(f"gap(M) is quadratic; vertex at M = {-cq[1]/(2*cq[0]):.3f}")
# two-block check of Prop ch12:thm:paradox
for xM, xT, M, k in [(1.2, 0.8, 0.4, 0.5), (1.2, 0.8, 0.8, 0.5), (1.5, 0.5, 0.6, 0.6)]:
    sb = v((1 - M) * xM, k) - v((1 - M) * xT, k)
    xbar = 0.5 * (xM + xT)
    si = v(xM - M * xbar, k) - v(xT - M * xbar, k)
    print(f"two-block x=({xM},{xT}), M={M}, kappa={k}: block gap {sb:.5f}, integrated gap {si:.5f}, "
          f"ratio {sb/si:.5f} (1-M = {1-M:.2f}); midpoints {(1-M)*xbar:.4f} = {xbar - M*xbar:.4f}")
    assert abs(sb / si - (1 - M)) < 1e-12
# non-quadratic v: ratio is not exactly 1-M
vl = lambda z: np.log(1 + z)
xM, xT, M = 1.2, 0.8, 0.8
r_log = (vl((1 - M) * xM) - vl((1 - M) * xT)) / (vl(xM - M * 1.0) - vl(xT - M * 1.0))
print(f"log(1+z) utility, x=(1.2,0.8), M=0.8: ratio = {r_log:.4f} (vs 1-M = 0.2)")

# ===========================================================================
print("\n== Exercise 12.2: laundered aspirations")
M = 0.8
ahat = M * x.mean()
print(f"mean x = {x.mean():.4f}; laundered aspiration a_hat = {ahat:.4f}")
print("laundered satisfaction by block:", np.round(v(xb - ahat), 4))
assert v(xb[2] - ahat) < 0
xeq = np.ones(N)
crit = {}
for name, xx in [("current", x), ("equal", xeq)]:
    endo = v(xx - M * WI @ xx).sum()
    laund = v(xx - ahat).sum()
    crit[name] = (endo, laund)
    print(f"{name}: endogenous W1 = {endo:.4f}, laundered W2 = {laund:.4f}, "
          f"linear part sum(x - a) endogenous = {np.sum(xx - M*WI@xx):.4f}")
print(f"differences equal - current: endogenous {crit['equal'][0]-crit['current'][0]:.4f}, "
      f"laundered {crit['equal'][1]-crit['current'][1]:.4f}")
print("column sums of W^I:", np.round(WI.sum(0), 4))
# with the media-hub W^I of Chapter 12 (30% of each row moved to the gas block)
Wrich = 0.7 * WI + 0.3 * np.tile(np.where(blk == 0, 1 / 3, 0.0), (N, 1))
for name, xx in [("current", x), ("equal", xeq)]:
    print(f"media-hub W^I, {name}: endogenous W1 = {v(xx - M * Wrich @ xx).sum():.4f}")

# ===========================================================================
print("\n== Exercise 12.3: capability")
pi2 = np.array([0.7, 0.3])
e = 5 + pi2 * 2
print("entitlements e =", e, " sum =", e.sum())
for c2 in [1.0, 0.8]:
    cap = np.array([1.0, c2]) * e
    print(f"c2 = {c2}: capability bounds {cap}, difference {cap[0]-cap[1]:.2f}, ratio {cap[0]/cap[1]:.4f}")

# ===========================================================================
print("\n== Exercise 12.4: Gamma is not a welfare criterion")
def W_of(g, p):
    return g * np.outer(np.ones(len(p)), p) + (1 - g) * np.eye(len(p))
cases = {"A": (0.9, np.array([0.55, 0.25, 0.15, 0.05])),
         "B": (0.3, np.array([0.30, 0.25, 0.25, 0.20])),
         "C": (0.3, np.array([0.55, 0.25, 0.15, 0.05]))}
Om, Ph, k4, n = 4.0, 2.0, 0.5, 4
for nm, (g, p) in cases.items():
    W = W_of(g, p)
    assert np.allclose(W.sum(1), 1) and np.allclose(perron(W), p)
    Vp = np.mean((p - 1 / n) ** 2)
    xx = Om / n + (p - 1 / n) * Ph
    assert np.all(k4 * xx < 1)
    U = v(xx, k4).sum()
    U0 = n * v(Om / n, k4)
    loss = 0.5 * k4 * n * Ph ** 2 * Vp
    assert abs((U0 - U) - loss) < 1e-12
    print(f"W_{nm}: Gamma = {gap(W):.3f}, pi = {p}, Var(pi) = {Vp:.6f}, N_eff = {1/np.sum(p**2):.3f}, "
          f"x = {xx}, sum v(x) = {U:.4f}, loss = {loss:.4f}")
print("benchmark N v(Omega/N) =", n * v(Om / n, k4))
print("W_A =\n", W_of(*cases["A"]))

# ===========================================================================
print("\n== Exercise 12.5: Easterlin pattern")
C, m = 10, 5
ybar = np.linspace(0.5, 2.0, C)
dev = np.array([-0.4, -0.2, 0.0, 0.2, 0.4])
X = (ybar[:, None] * (1 + dev[None, :])).ravel()
cty = np.repeat(np.arange(C), m)
Nn = C * m
Wnat = np.where(cty[:, None] == cty[None, :], 1 / m, 0.0)
Wglob = np.ones((Nn, Nn)) / Nn
Ma = 0.8


def slopes(WI_, kappa):
    a = Ma * WI_ @ X
    z = X - a
    assert np.all(kappa * z < 1)
    s = z - 0.5 * kappa * z ** 2
    sm = np.array([s[cty == c].mean() for c in range(C)])
    xm = np.array([X[cty == c].mean() for c in range(C)])
    cross = np.polyfit(xm, sm, 1)[0]
    sd = s - sm[cty]
    xd = X - xm[cty]
    within = np.sum(sd * xd) / np.sum(xd ** 2)
    return cross, within


for kappa in [0.0, 0.2]:
    for name, Wi in [("national", Wnat), ("partial eps=0.5", 0.5 * Wnat + 0.5 * Wglob),
                     ("global", Wglob)]:
        cr, wi = slopes(Wi, kappa)
        print(f"kappa = {kappa}, {name:16s}: cross-country slope = {cr:.4f}, within-country slope = {wi:.4f}")
print("analytic (linear v): cross slope = 1 - M(1-eps): national", 1 - Ma, " partial", 1 - Ma * 0.5, " global", 1.0)

# ===========================================================================
print("\n== Exercise 12.6: evaluation and preference vintage")
WIc = WI.copy()
for i in np.where(blk == 2)[0]:
    WIc[i] = np.where(blk == 0, 0.075, np.where(blk == 2, 0.70 / 3, 0.025))
a_old = Ma * WI @ x
a_new = Ma * WIc @ x
s0 = v(x[6] - a_old[6])
print(f"reproduce: rural aspiration {a_old[6]:.4f} -> {a_new[6]:.4f}; rural satisfaction {s0:.4f} -> {v(x[6]-a_new[6]):.4f}")
assert abs(a_old[6] - 0.655) < 5e-4 and abs(a_new[6] - 0.713) < 5e-4
first = {}
for rho_ in [0.5, 0.75, 0.9]:
    mu_ = Ma * (1 - rho_)
    a = a_old.copy()
    rows = []
    fi = {"old": None, "current": None, "new_ss": None}
    for t in range(0, 21):
        xt = x.copy()
        xt[blk == 2] += 0.02 * t
        s_old = v(xt[6] - a_old[6])
        s_cur = v(xt[6] - a[6])
        s_new = v(xt[6] - (Ma * WIc @ xt)[6])
        for key, val in [("old", s_old), ("current", s_cur), ("new_ss", s_new)]:
            if fi[key] is None and val > s0 + 1e-12:
                fi[key] = t
        rows.append((t, s_old, s_cur, s_new))
        a = rho_ * a + mu_ * WIc @ xt            # a_{t+1} = rho a_t + mu W_new x_t
    first[rho_] = fi
    if rho_ == 0.75:
        for t, so, sc, sn in rows:
            if t in (0, 1, 2, 3, 5, 7, 10, 15, 20):
                print(f"   rho=0.75 t={t:2d}: old {so:.4f}, current {sc:.4f}, new-ss {sn:.4f}")
    print(f"rho = {rho_} (mu = {mu_:.3f}): first t with rural s > {s0:.4f}: {fi}; current-vintage minimum "
          f"{min(r[2] for r in rows):.4f} at t={min(rows, key=lambda r: r[2])[0]}")
# closed form for the new-steady-state criterion
w_own = 0.70
print(f"new-ss: z_t = z_0 + 0.02 t (1 - M*0.70) ; z_0 = {x[6]-a_new[6]:.4f}, slope {0.02*(1-Ma*w_own):.4f}, "
      f"target z = {x[6]-a_old[6]:.4f}")

# ===========================================================================
print("\n== Exercise 12.7: Sen vs Bentham")
k7 = 0.1
aL, aW = 4.0, 1.5
a7 = np.array([aL, aW])
res = {}
for lab, p in [("unequal", np.array([0.8, 0.2])), ("equal", np.array([0.5, 0.5]))]:
    x7 = 5 + (p - 0.5) * 4
    e7 = 3 + p * 4
    z = x7 - a7
    assert np.all(k7 * z < 1) and np.all(k7 * (x7 - 3) < 1)
    s = v(z, k7)
    sl = v(x7 - 3.0, k7)
    res[lab] = dict(x=x7, s=s, sl=sl, e=e7)
    print(f"{lab}: x = {x7}, kappa*z = {k7*z}, s(current) = {np.round(s,4)}, Bentham = {s.sum():.4f}, "
          f"laundered s = {np.round(sl,4)}, sum = {sl.sum():.4f}, Rawls = {s.min():.4f}, e = {e7}, min e = {e7.min():.2f}")
# (d) adaptation: aspirations proportional to own allocation (a_i / x_i fixed at its current ratio)
ratio = a7 / res["unequal"]["x"]
xe = res["equal"]["x"]
a_ad = ratio * xe
s_ad = v(xe - a_ad, k7)
print(f"(d) adapted aspirations a' = {np.round(a_ad,4)} (ratios {np.round(ratio,4)}), s' = {np.round(s_ad,4)}, "
      f"Bentham {s_ad.sum():.4f}, Rawls {s_ad.min():.4f}")
print(f"    landlord loss: fixed aspirations {res['unequal']['s'][0]-res['equal']['s'][0]:.4f}, "
      f"adapted {res['unequal']['s'][0]-s_ad[0]:.4f}; laundered {res['unequal']['sl'][0]-res['equal']['sl'][0]:.4f}; "
      f"entitlement {res['unequal']['e'][0]-res['equal']['e'][0]:.2f}")
