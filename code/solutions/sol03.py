"""Solutions manual, Chapter 3: every number used in solutions/sol03.tex.
Run from the repository root:  python3 code/solutions/sol03.py
numpy throughout; pandas only to read data/processed (Exercise 3.5)."""
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)


def hdr(s):
    print("\n== " + s)


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eigs(W):
    e = np.linalg.eigvals(W)
    return e[np.argsort(-abs(e))]


def gap(W):
    return 1 - abs(eigs(W)[1])


def leontief(A, g=1.0):
    return np.linalg.inv(np.eye(len(A)) - g * A.T)


def rownorm(A):
    return A / A.sum(1, keepdims=True)


# ---------------------------------------------------------------- 3.1
hdr("3.1 Solow with two groups (alpha=0.3, theta_M=2, theta_F=0.5)")
a = 0.3
sig = lambda th: th / (1 + th)
sM, sF = sig(2), sig(.5); sbar = (sM + sF) / 2
Kss = lambda s: (s * (1 - a)) ** (1 / (1 - a))
print("sigma_M=%.4f sigma_F=%.4f sbar=%.4f  s(1-a) for s=2/3: %.4f, for sbar: %.4f" % (sM, sF, sbar, sig(2) * (1 - a), sbar * (1 - a)))
K1, K2 = Kss(sbar), Kss(sig(2))
print("K*(sbar)=%.4f  K*(both theta=2, s=%.4f)=%.4f  ratio=%.4f  gap=%.2f%%" % (K1, sig(2), K2, K1 / K2, 100 * (1 - K1 / K2)))
print("ratio formula (0.5/(2/3))^(1/0.7) = %.4f ; output ratio (K ratio)^alpha = %.4f, gap %.2f%%"
      % ((sbar / sig(2)) ** (1 / (1 - a)), (K1 / K2) ** a, 100 * (1 - (K1 / K2) ** a)))
# check stability: K_{t+1} = s(1-a)K^a converges
K = 1.0
for t in range(200):
    K = sbar * (1 - a) * K ** a
assert np.isclose(K, K1)
print("wage w=(1-a)K^a at K1: %.4f, at K2: %.4f" % ((1 - a) * K1 ** a, (1 - a) * K2 ** a))
print("consumption check: savings s^g = sigma w; young M saves %.4f, F saves %.4f (at K1)" % (sM * (1 - a) * K1 ** a, sF * (1 - a) * K1 ** a))
# golden rule in this OLG: steady-state consumption c = (1-a)K^a - K ... report marginal product
print("marginal product of capital a K^(a-1): at K1 %.4f, at K2 %.4f" % (a * K1 ** (a - 1), a * K2 ** (a - 1)))
Kgr = a ** (1 / (1 - a))
C = lambda K: K ** a - K
print("golden rule K_gr = a^(1/(1-a)) = %.4f ; steady-state consumption C=K^a-K: C(K1)=%.4f C(K2)=%.4f C(Kgr)=%.4f"
      % (Kgr, C(K1), C(K2), C(Kgr)))
print("output Y=K^a: Y(K1)=%.4f Y(K2)=%.4f" % (K1 ** a, K2 ** a))

# ---------------------------------------------------------------- 3.2
hdr("3.2 Propagation in a 3-sector economy")
A = np.array([[.10, .20, .05], [.15, .10, .15], [.05, .10, .10]])
L = leontief(A)
print("A row sums:", A.sum(1), " col sums:", A.sum(0))
print("L =\n", L.round(4))
print("(b) L e2 =", L[:, 1].round(4), " total =", L[:, 1].sum().round(4))
print("(c) multipliers (col sums):", L.sum(0).round(4), " forward (row sums):", L.sum(1).round(4))
rho = max(abs(np.linalg.eigvals(A))); print("(d) rho(A) = %.4f" % rho)
WP = rownorm(A); p = perron(WP)
print("W^P =\n", WP.round(4), "\n eig:", np.round(eigs(WP), 4), " Gamma(W^P)=%.4f  pi=%s  N_eff=%.4f" % (gap(WP), p.round(4), 1 / np.sum(p ** 2)))
off = L - np.diag(np.diag(L))
print("mean off-diagonal L entry %.4f ; norm bound gamma||A||_2/(1-gamma||A||_2) = %.4f (||A||_2 = %.4f)"
      % (off.sum() / 6, np.linalg.norm(A, 2) / (1 - np.linalg.norm(A, 2)), np.linalg.norm(A, 2)))
blk = np.array([1, 1, 0], bool)
chi = (A[np.ix_(blk, ~blk)].sum() + A[np.ix_(~blk, blk)].sum()) / A.sum()
print("cross-block mass chi(A) for {1,2} vs {3}: %.4f" % chi)
L_conn = L

# ---------------------------------------------------------------- 3.3
hdr("3.3 The enclave (a31=a32=a13=a23=0.01)")
Ae = A.copy(); Ae[2, 0] = Ae[2, 1] = Ae[0, 2] = Ae[1, 2] = .01
Le = leontief(Ae)
print("A_enclave =\n", Ae, "\nL_enclave =\n", Le.round(4))
print("(a) shock to sector 3: L e3 =", Le[:, 2].round(4), " spillover to 1+2 = %.4f" % Le[:2, 2].sum())
print("(b) connected: L e3 =", L_conn[:, 2].round(4), " spillover to 1+2 = %.4f" % L_conn[:2, 2].sum())
print("    factors: to sector 1 %.2f, to sector 2 %.2f, total %.2f"
      % (L_conn[0, 2] / Le[0, 2], L_conn[1, 2] / Le[1, 2], L_conn[:2, 2].sum() / Le[:2, 2].sum()))
print("    reverse direction (sectors 1,2 -> 3): connected", L_conn[2, :2].round(4), " enclave", Le[2, :2].round(4))
print("    multipliers enclave:", Le.sum(0).round(4), " connected:", L_conn.sum(0).round(4))
WPe = rownorm(Ae); pe = perron(WPe)
chie = (Ae[np.ix_(blk, ~blk)].sum() + Ae[np.ix_(~blk, blk)].sum()) / Ae.sum()
print("(c) W^P_enclave =\n", WPe.round(4), "\n eig:", np.round(eigs(WPe), 4), " Gamma=%.4f (connected %.4f)  pi=%s  chi=%.4f (connected %.4f)"
      % (gap(WPe), gap(WP), pe.round(4), chie, chi))
print("    rho(A_enclave) = %.4f" % max(abs(np.linalg.eigvals(Ae))))
print("    sector-3 cross-block input share in W^P: %.4f (connected %.4f); sectors 1,2 share bought from 3: %s (connected %s)"
      % (WPe[2, :2].sum(), WP[2, :2].sum(), WPe[:2, 2].round(4), WP[:2, 2].round(4)))
# lumped 2-block chain bound: Gamma <= eps1+eps2 when exactly lumpable; here approximate
pib = pe[:2].sum()
print("    block {1,2} Perron mass %.4f (connected %.4f)" % (pib, p[:2].sum()))

# ---------------------------------------------------------------- 3.4
hdr("3.4 Transfer vs structural policy (stylized dual economy)")
aM, aT, s_, d_ = .4, .3, .2, .05
Zp = np.array([2., 1.]); al = np.array([aM, aT]); Lab = np.array([1., 1.])
Adual = np.array([[.15, .03], [.02, .20]]); Alink = np.array([[.15, .10], [.10, .20]])


def f(K):
    return Zp * K ** al * Lab ** (1 - al)


def steady(Am, tau=0.0, iters=5000):
    K = np.ones(2)
    for _ in range(iters):
        Y = leontief(Am) @ f(K)
        inc = np.array([(1 - tau) * Y[0], Y[1] + tau * Y[0]])
        K = s_ * inc + (1 - d_) * K
    return K


def simulate(Am, K0, T, mode, tau=0.1):
    K = K0.copy(); out = []
    for t in range(T + 1):
        Y = leontief(Am) @ f(K)
        if mode == "transfer":
            inc = np.array([(1 - tau) * Y[0], Y[1] + tau * Y[0]])
        elif mode == "linkage":
            inc = np.array([(1 - tau) * Y[0], Y[1]])
        else:
            inc = Y.copy()
        out.append((t, K.copy(), Y.copy(), inc.copy()))
        K = s_ * inc + (1 - d_) * K
    return out


K0 = steady(Adual)
Y0 = leontief(Adual) @ f(K0)
print("marginal products of capital in f at baseline steady state: M %.4f T %.4f"
      % (aM * Zp[0] * K0[0] ** (aM - 1), aT * Zp[1] * K0[1] ** (aT - 1)))
print("Ldual =\n", leontief(Adual).round(4), "\nLlink =\n", leontief(Alink).round(4))
print("baseline steady state K* =", K0.round(4), " Y* =", Y0.round(4), " Y_M/Y_T = %.4f" % (Y0[0] / Y0[1]))
print("(a) impact of transfer: I_T from %.4f to %.4f (+%.4f = 0.1 Y_M), I_M from %.4f to %.4f; Y unchanged on impact"
      % (Y0[1], Y0[1] + .1 * Y0[0], .1 * Y0[0], Y0[0], .9 * Y0[0]))
Lg = leontief(Alink)
print("(b) cross-sector propagation: dual L_TM=%.4f L_MT=%.4f ; linked L_TM=%.4f L_MT=%.4f ; ratio TM %.2f"
      % (leontief(Adual)[1, 0], leontief(Adual)[0, 1], Lg[1, 0], Lg[0, 1], Lg[1, 0] / leontief(Adual)[1, 0]))
print("    multipliers dual", leontief(Adual).sum(0).round(4), " linked", Lg.sum(0).round(4))
res = {}
for name, Am, mode in [("none", Adual, "none"), ("transfer", Adual, "transfer"), ("linkage", Alink, "linkage")]:
    out = simulate(Am, K0, 10, mode)
    res[name] = out
    t, K, Y, inc = out[-1]
    VA = Y * (1 - Am.sum(1))
    print("(c) %-8s t=10  K=%s  Y=%s  Y_M/Y_T=%.4f  income=%s  I_M/I_T=%.4f  total income=%.4f  VA=%s VA_M/VA_T=%.4f"
          % (name, K.round(4), Y.round(4), Y[0] / Y[1], inc.round(4), inc[0] / inc[1], inc.sum(), VA.round(4), VA[0] / VA[1]))
for name in res:
    print("   %-8s Y_M/Y_T at t=0,1,5,10:" % name, [round(float(res[name][t][2][0] / res[name][t][2][1]), 4) for t in [0, 1, 5, 10]],
          " I_M/I_T:", [round(float(res[name][t][3][0] / res[name][t][3][1]), 4) for t in [0, 1, 5, 10]])
for name, Am, tau in [("transfer", Adual, .1), ("linkage", Alink, 0.0)]:
    Ks = steady(Am, tau) if name == "transfer" else None
    if name == "transfer":
        Ys = leontief(Am) @ f(Ks)
        print("   long-run (steady state) transfer: Y_M/Y_T = %.4f" % (Ys[0] / Ys[1]))
# long run linkage policy with continuing cost
K = K0.copy()
for _ in range(5000):
    Y = leontief(Alink) @ f(K); inc = np.array([.9 * Y[0], Y[1]]); K = s_ * inc + (1 - d_) * K
Y = leontief(Alink) @ f(K)
print("   long-run linkage policy (cost continues): Y_M/Y_T = %.4f ; Y_T = %.4f" % (Y[0] / Y[1], Y[1]))
# value-added productivity: VA_i = Y_i (1 - sum_j a_ij) per worker
for name, Am in [("dual", Adual), ("link", Alink)]:
    print("   value added share of gross output, %s: %s" % (name, (1 - Am.sum(1)).round(3)))

# ---------------------------------------------------------------- 3.5
hdr("3.5 Prebisch in the global table: trade-based proxy, Growth Lab HS92, 2019")
tr = pd.read_csv("data/processed/atlas_bilateral_selected_years.csv.gz")
tr = tr[tr.year == 2019]
nodes = ["USA", "DEU", "CHN", "BOL", "NGA"]
tr = tr.assign(e=tr.exporter.where(tr.exporter.isin(nodes), "ROW"), m=tr.importer.where(tr.importer.isin(nodes), "ROW"))
F = tr.groupby(["m", "e"]).trade_value.sum().astype(float).unstack(fill_value=0.0)   # rows importers, cols exporters
order = nodes + ["ROW"]
F = F.reindex(index=order, columns=order, fill_value=0.0)
Wt = F.values / F.values.sum(1, keepdims=True)
print("rows = importer, columns = exporter (share of importer's goods imports):")
print(pd.DataFrame(Wt, index=order, columns=order).round(4).to_string())
cs = Wt.sum(0)
print("(a) column sums:", dict(zip(order, cs.round(4))))
pt = perron(Wt)
print("(e) pi (6-node aggregation):", dict(zip(order, pt.round(5))), " Gamma=%.4f" % gap(Wt))
print("    eigenvalue moduli:", np.round(abs(eigs(Wt)), 4))
print("(c) exposure to Bolivian supply: W[USA,BOL]=%.6f  W[NGA,BOL]=%.6f  ratio %.1f" % (Wt[0, 3], Wt[4, 3], Wt[0, 3] / Wt[4, 3]))
print("    Bolivia's imports by source:", dict(zip(order, Wt[3].round(4))))
print("    Bolivia's exports to USA: %.0f, to NGA: %.0f (USD)" % (F.loc["USA", "BOL"], F.loc["NGA", "BOL"]))
# 5 named countries only (rows renormalised over the five)
F5 = F.loc[nodes, nodes].values
W5 = F5 / F5.sum(1, keepdims=True); p5 = perron(W5)
print("5-country closed system: column sums", dict(zip(nodes, W5.sum(0).round(4))), " pi", dict(zip(nodes, p5.round(4))), " Gamma=%.4f" % gap(W5))
# cross-block mass, center {USA,DEU,CHN} vs periphery {BOL,NGA}, in the 5-country system (flows)
cen = np.array([1, 1, 1, 0, 0], bool)
chi5 = (F5[np.ix_(cen, ~cen)].sum() + F5[np.ix_(~cen, cen)].sum()) / F5.sum()
print("(d) cross-block mass of flows center/periphery (5-country system): %.4f" % chi5)
print("    periphery rows' weight on center (5-country): BOL %.4f NGA %.4f ; center rows' weight on periphery: %s"
      % (W5[3, cen].sum(), W5[4, cen].sum(), W5[cen][:, ~cen].sum(1).round(5)))
g = pd.read_csv("data/processed/global_wp_pi.csv")
g = g[(g.year == 2019) & g.country.isin(nodes)].set_index("country")
print("full 230-economy W (data/processed/global_wp_pi.csv), 2019:")
print(g[["pi", "colsum_W", "world_export_share"]].reindex(nodes).to_string())
print("pi_USA/pi_BOL: 6-node %.0f, 230-economy %.0f" % (pt[0] / pt[3], g.loc["USA", "pi"] / g.loc["BOL", "pi"]))

# ---------------------------------------------------------------- 3.6
hdr("3.6 Kuznets revisited (T=50)")


def run(c, T=50, Z=np.array([2., 1.]), alph=np.array([.4, .3]), K0=np.ones(2)):
    Am = np.array([[.15, c], [c, .20]])
    Lm = leontief(Am); K = K0.copy(); path = []
    for t in range(T + 1):
        Y = Lm @ (Z * K ** alph)
        path.append((K.copy(), Y.copy()))
        K = s_ * Y + (1 - d_) * K
    return path


for name, c in [("A connected", .15), ("B disconnected", .02)]:
    path = run(c)
    ratio = np.array([Y[0] / Y[1] for K, Y in path])
    VB = np.array([0.25 * np.log(Y[0] / Y[1]) ** 2 for K, Y in path])
    print("%s: Y_M/Y_T at t=0,1,2,5,10,20,30,40,50:" % name, [round(float(ratio[t]), 4) for t in [0, 1, 2, 5, 10, 20, 30, 40, 50]])
    print("   between-sector var of log output/worker (equal labour shares):", [round(float(VB[t]), 4) for t in [0, 1, 2, 5, 10, 20, 30, 40, 50]])
    print("   monotone ratio path:", bool(np.all(np.diff(ratio) > 0) or np.all(np.diff(ratio) < 0)), " argmax VB:", int(VB.argmax()))
    long = run(c, T=3000)[-1][1]
    print("   long-run Y_M/Y_T = %.4f" % (long[0] / long[1]))
# Kuznets mechanism: weight w on modern sector, V_B = w(1-w) (dlog y)^2
print("long-run ratio A/B = %.4f" % (run(.15, T=3000)[-1][1][0] / run(.15, T=3000)[-1][1][1] / (run(.02, T=3000)[-1][1][0] / run(.02, T=3000)[-1][1][1])))
dl = np.log(run(.02, T=3000)[-1][1][0] / run(.02, T=3000)[-1][1][1])
print("Case B long-run log gap %.4f; V_B = w(1-w) gap^2 at w=0.1,0.3,0.5,0.7,0.9:" % dl,
      [round(float(w * (1 - w) * dl ** 2), 4) for w in [.1, .3, .5, .7, .9]])
print("(d) long-run ratio by cross-linkage c:")
for c in [0, .02, .05, .10, .15, .20, .30, .40]:
    Y = run(c, T=4000)[-1][1]
    rho = max(abs(np.linalg.eigvals(np.array([[.15, c], [c, .20]]))))
    print("   c=%.2f rho(A)=%.4f  Y_M/Y_T=%.4f  Y=%s" % (c, rho, Y[0] / Y[1], Y.round(4)))
# AK case (alpha = 1): ratio diverges iff decoupled with different growth factors
for c in [0, .02]:
    path = run(c, T=200, alph=np.array([1., 1.]), Z=np.array([.3, .2]))
    r = [path[t][1][0] / path[t][1][1] for t in [0, 50, 100, 200]]
    Am = np.array([[.15, c], [c, .20]]); G = s_ * leontief(Am) @ np.diag([.3, .2]) + (1 - d_) * np.eye(2)
    ev, V = np.linalg.eig(G); v = np.abs(V[:, np.argmax(ev.real)].real); Ylim = leontief(Am) @ (np.array([.3, .2]) * v)
    print("   AK (Z=(0.3,0.2), alpha=1) c=%.2f: Y_M/Y_T at t=0,50,100,200: %s ; growth factors (eig G) %s ; Perron-direction limit ratio %s"
          % (c, np.round(r, 4), np.round(np.sort(ev.real)[::-1], 4), "inf" if c == 0 else round(float(Ylim[0] / Ylim[1]), 2)))

# ---------------------------------------------------------------- 3.7
hdr("3.7 The formation connection")
rho_, mu, thb, rstar, r = .8, .1, 2., .1, .1
thM = 2.


def path(lam_fn, T=20, th0=.5):
    th = [th0]
    for t in range(T):
        th.append(rho_ * th[-1] + (1 - rho_) * thb + mu * (lam_fn(t) * r - rstar))
    return np.array(th)


th = path(lambda t: .7)
gapv = thM - th
print("fixed point theta*(lambda=0.7) = %.4f ; long-run gap %.4f ; M = mu/(1-rho) = %.2f" % ((0.4 - .003) / .2, thM - (0.4 - .003) / .2, mu / (1 - rho_)))
print("theta_t, t=0..20:", th.round(4), " initial deviation from theta*: %.3f" % (th[0] - (0.4 - .003) / .2 + 0))
print("gap_t:", gapv.round(4))
th_half = np.argmax(gapv <= gapv[0] / 2)
print("gap halves (<= %.4f) first at generation %d ; continuous solution t = %.3f" % (gapv[0] / 2, th_half, np.log((0.75 - 0.015) / 1.485) / np.log(.8)))
th2 = path(lambda t: .7 if t < 10 else 1.0)
g2 = thM - th2
print("switch lambda=1 at t=10: gap_10 = %.4f ; gaps t=10..20:" % g2[10], g2[10:].round(4))
k = np.argmax(g2[10:] <= g2[10] / 2)
print("additional generations to halve from t=10: %d ; continuous ln2/-ln(0.8) = %.3f" % (k, np.log(2) / -np.log(.8)))
assert th_half == 4 and k == 4
print("\nAll Chapter 3 solution checks passed.")
