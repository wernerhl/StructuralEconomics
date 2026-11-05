"""Solutions manual, Chapter 4: every number used in solutions/sol04.tex.
Run from the repository root:  python3 code/solutions/sol04.py
numpy throughout; pandas only to read data/processed (Exercise 4.4)."""
import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)


def hdr(s):
    print("\n== " + s)


def U(x, k):
    x = np.asarray(x, float)
    return x - k / 2 * x ** 2


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eigs(W):
    e = np.linalg.eigvals(W)
    return e[np.argsort(-abs(e))]


def gap(W):
    return 1 - abs(eigs(W)[1])


def neff(p):
    return 1 / np.sum(p ** 2)


def var(p):
    return np.mean((p - 1 / len(p)) ** 2)


def dpi(p):
    return 0.5 * np.abs(p - 1 / len(p)).sum()


def gini(p):
    p = np.sort(np.asarray(p, float)); N = len(p)
    return float(np.sum((2 * np.arange(1, N + 1) - N - 1) * p) / (N * p.sum()))


# ---------------------------------------------------------------- 4.2
hdr("4.2 Measuring departure")
W = np.array([[.5, .3, .2], [.1, .7, .2], [.2, .1, .7]])
cs = W.sum(0)
print("row sums", W.sum(1), " column sums", cs, " delta_col = %.4f (sqrt(0.06) = %.4f)" % (np.linalg.norm(cs - 1), np.sqrt(.06)))
p = perron(W)
print("pi =", p.round(4), " exact check pi'W - pi' =", (p @ W - p).round(12))
print("Var(pi)=%.6f  N_eff=%.4f  N_eff/N=%.4f  d_pi=%.4f  max Var=(N-1)/N^2=%.4f  Gamma=%.4f  eig=%s"
      % (var(p), neff(p), neff(p) / 3, dpi(p), 2 / 9, gap(W), np.round(eigs(W), 4)))
# exact pi by solving
M = np.vstack([(W.T - np.eye(3))[:2], np.ones(3)]); pex = np.linalg.solve(M, [0, 0, 1])
print("pi exact (solve) =", pex, " = (", np.round(pex * 30, 6), ")/30")
P1 = np.outer(np.ones(3), p)
print("distance from anonymity ||W - 1 pi'||_F = %.4f ; ||W - W^CE||_F = %.4f" % (np.linalg.norm(W - P1), np.linalg.norm(W - np.ones((3, 3)) / 3)))

# ---------------------------------------------------------------- 4.3
hdr("4.3 The hierarchy in action (N=4, Omega=4, kappa=0.5, Phi=1, gamma_U=0.5)")
N, Om, kap, Phi, gam = 4, 4., .5, 1., .5
m = Om / N
WU = np.array([[0, .5, .25, .25], [.6, 0, .2, .2], [.6, .2, 0, .2], [.6, .2, .2, 0]])
c = WU.sum(0)
Wce = np.ones((4, 4)) / 4


def star(h, a=.8, s=.8):
    S = np.zeros((4, 4))
    for i in range(4):
        if i == h:
            S[i] = (1 - s) / 3; S[i, i] = s
        else:
            S[i, h] = a; S[i, i] = 1 - a
    return S


def losses(p, cvec, g):
    alpha = 1 + g * cvec
    H = N / np.sum(1 / alpha)
    Lext = (1 - kap * m) ** 2 / (2 * kap) * np.sum((alpha - H) ** 2 / alpha)
    Lpow = kap * N / 2 * Phi ** 2 * var(p)
    cov = np.mean((p - p.mean()) * (cvec - cvec.mean()))
    Lint = kap * g * Phi ** 2 / 2 * np.sum(cvec * (p - 1 / N) ** 2) - g * Phi * N * (1 - kap * m) * cov
    xs = m + (1 - kap * m) * (alpha - H) / (kap * alpha)
    xb = m + (p - 1 / N) * Phi
    Ltot = np.sum(alpha * U(xs, kap)) - np.sum(alpha * U(xb, kap))
    assert np.isclose(Ltot, Lext + Lpow + Lint) and kap * max(xs.max(), xb.max()) < 1
    return Lext, Lpow, Lint, Ltot, cov, xb


print("W^U (zero diagonal, row-stochastic):\n", WU, "\n c (column sums) =", c, " Var(c) = %.4f" % np.var(c))
sc = {}
aU, aM = .2, .8
Wa = aU * WU + aM * Wce
sc["a"] = (Wa, c, gam)
Wb = star(3)
sc["b"] = (Wb, np.zeros(4), 0.0)
sc["c1"] = (aU * WU + aM * star(0), c, gam)
sc["c2"] = (aU * WU + aM * star(3), c, gam)
res = {}
for k_, (Wm, cv, g) in sc.items():
    p = perron(Wm)
    Le, Lp, Li, Lt, cov, xb = losses(p, cv, g)
    res[k_] = (Le, Lp, Li, Lt)
    print("(%s) col sums %s  pi %s  N_eff %.3f  Gamma %.4f  Cov(pi,c) %.5f  x %s\n     L_ext %.5f  L_power %.5f  L_int %.5f  L %.5f"
          % (k_, Wm.sum(0).round(3), p.round(4), neff(p), gap(Wm), cov, xb.round(4), Le, Lp, Li, Lt))
# A1 fails with doubly stochastic W^U: pi uniform
WUds = np.array([[0, .5, 0, .5], [.5, 0, .5, 0], [0, .5, 0, .5], [.5, 0, .5, 0]])
pds = perron(aU * WUds + aM * Wce)
print("(a') doubly stochastic W^U (ring of care): pi =", pds.round(4), " c =", WUds.sum(0))
for k_ in ["c1", "c2"]:
    s_ = res["a"][3] + res["b"][3]
    print("(%s) L = %.5f vs L(a)+L(b) = %.5f + %.5f = %.5f ; difference %.5f" % (k_, res[k_][3], res["a"][3], res["b"][3], s_, res[k_][3] - s_))
print("star hub pi_h formula a/(a+1-s) = %.4f" % (.8 / (.8 + 1 - .8)))

# ---------------------------------------------------------------- 4.4
hdr("4.4 Bolivia's six assumptions: measured departure measures")
f = pd.read_csv("data/processed/bolivia_io_sector_facts.csv")


def A_of(base, year):
    Z = pd.read_csv("data/processed/bolivia_Z/Z_%d_%d_book12.csv" % (base, year), index_col=0)
    ff = f[(f.base == base) & (f.year == year)].set_index("book12")
    x = ff.output.reindex(Z.index).values.astype(float)
    keep = x > 0
    A = (Z.values.T / x[:, None])[np.ix_(keep, keep)]
    return A, [l for l, k in zip(Z.index, keep) if k]


for base, y0, y1 in [(1990, 2003, 2013), (1990, 1990, 2000), (1990, 2000, 2009), (2017, 2017, 2024)]:
    A0, l0 = A_of(base, y0); A1, l1 = A_of(base, y1)
    assert l0 == l1
    ch = np.linalg.norm(A1 - A0) / np.linalg.norm(A0)
    print("base %d, %d->%d (%d sectors): ||dA||_F/||A||_F = %.4f ; per decade (x10/%d) = %.4f"
          % (base, y0, y1, len(l0), ch, len(l0) and (y1 - y0), ch * 10 / (y1 - y0)))
s = pd.read_csv("data/processed/bolivia_io_summary.csv")
row = s[(s.base == 2017) & (s.year == 2019) & (s.aggregation == "book12") & (s.flows == "domestic") & (s.normalization == "input_shares")].iloc[0]
print("W^P 2019 book12 domestic input shares: Gamma=%.4f  N_eff=%.4f of N=%d  N_eff/N=%.4f  power Gini=%.4f  d_pi=%.4f  pi_max=%.4f (%s)"
      % (row.Gamma, row.Neff, row.N, row.Neff / row.N, row.power_gini, row.d_pi, row.pi_max, row.pi_argmax))
fa = pd.read_csv("data/processed/bolivia_facts.csv").set_index("fact_id")
for k_ in ["BOL_INE_gini_2005", "BOL_INE_gini_2019"]:
    print("fact:", k_, fa.loc[k_, "value"], fa.loc[k_, "unit"])
# stylized 2x2 W^P of Chapter 3
Ad = np.array([[.15, .03], [.02, .20]]); WPd = Ad / Ad.sum(1, keepdims=True); pd2 = perron(WPd)
print("stylized dual W^P =", WPd.round(3).tolist(), " Gamma=%.4f pi=%s N_eff/N=%.4f" % (gap(WPd), pd2.round(4), neff(pd2) / 2))

# ---------------------------------------------------------------- 4.5
hdr("4.5 Sunkel: one rewiring, five observables")
lab = ["C", "E", "T1", "T2"]
W0 = np.array([[.90, .10, 0, 0], [.20, .40, .20, .20], [.05, .35, .40, .20], [.05, .25, .20, .50]])
W1 = np.array([[.80, .20, 0, 0], [.55, .40, .025, .025], [.05, .10, .65, .20], [.05, .05, .20, .70]])
E = W1 - W0
assert np.allclose(E.sum(1), 0)
OmS, PhiS, NS = 4., 2., 4
for t in [0, .5, 1]:
    Wt = (1 - t) * W0 + t * W1
    p = perron(Wt)
    x = OmS / NS + (p - 1 / NS) * PhiS
    dom = p[1:].sum()
    dep = np.sum(p[1:] * Wt[1:, 0]) / dom                 # domestic attention to the center, pi-weighted
    trans = (Wt[0, 1] + Wt[1, 0])                          # center-elite cross weight
    natl = (Wt[1, 2] + Wt[1, 3] + Wt[2, 1] + Wt[3, 1])     # elite-traditional cross weights
    Wd = Wt[1:, 1:] / Wt[1:, 1:].sum(1, keepdims=True)     # domestic network, renormalised
    print("t=%.1f pi=%s  x=%s  dependency=%.4f  transnational=%.3f  national=%.3f  Gamma(W)=%.4f  Gamma(domestic)=%.4f  x_E-x_Tbar=%.4f  pi_T=%.4f"
          % (t, p.round(4), x.round(4), dep, trans, natl, gap(Wt), gap(Wd), x[1] - x[2:].mean(), p[2:].sum()))
p0 = perron(W0)
Zf = np.linalg.inv(np.eye(4) - W0 + np.outer(np.ones(4), p0))
d_formula = p0 @ E @ Zf
h = 1e-6
d_num = (perron(W0 + h * E) - p0) / h
print("d pi/dt at t=0: formula pi'EZ = %s ; numerical = %s" % (d_formula.round(4), d_num.round(4)))
assert np.allclose(d_formula, d_num, atol=1e-4)

# ---------------------------------------------------------------- 4.6
hdr("4.6 A6 is the deepest: a transfer that lowers long-run welfare")
kap6 = .05
yE, base_t, beta = 10., 4., 2.
thbar, rho6, mu6 = 1., .9, .5
rM, rT, rstar, w0, om = .25, .05, .15, .5, .2
M6 = mu6 / (1 - rho6)


def w_of(tau):
    return w0 + om * max(tau, 0.)


def rhat(tau):
    return w_of(tau) * rT + (1 - w_of(tau)) * rM


def theta_ss(tau):
    return thbar + M6 * (rhat(tau) - rstar)


def welfare(tau, theta):
    return U(yE - tau, kap6) + U(base_t + beta * theta + tau, kap6)


yt0 = base_t + beta * theta_ss(0)
tauA = (yE - yt0) / 2
print("M = %.2f ; theta*(0) = %.4f ; y_t(0) = %.4f ; Economist A's tau* = %.4f ; w(tau*) = %.2f" % (M6, theta_ss(0), yt0, tauA, w_of(tauA)))
phi = M6 * om * (rM - rT)
print("d theta*/d tau = -%.4f ; d y_t/d tau (long run) = 1 - beta*phi = %.4f" % (phi, 1 - beta * phi))
print("A's welfare (theta fixed): W(0)=%.4f  W(tau*)=%.4f" % (welfare(0, theta_ss(0)), welfare(tauA, theta_ss(0))))
print("B's long-run welfare: W(0)=%.4f  W(tau*)=%.4f  theta*(tau*)=%.4f  y_t=%.4f"
      % (welfare(0, theta_ss(0)), welfare(tauA, theta_ss(tauA)), theta_ss(tauA), base_t + beta * theta_ss(tauA) + tauA))
grid = np.linspace(0, 3, 301)
WB = np.array([welfare(t, theta_ss(t)) for t in grid])
print("B's long-run optimum over tau in [0,3]: tau=%.2f W=%.4f ; dW_B/dtau at 0+ = %.4f" % (grid[WB.argmax()], WB.max(), -(1 - kap6 * yE) + (1 - beta * phi) * (1 - kap6 * yt0)))
th = theta_ss(0); path = []
for t in range(31):
    path.append(welfare(tauA, th))
    th = rho6 * th + (1 - rho6) * thbar + mu6 * (rhat(tauA) - rstar)
path = np.array(path)
cross = int(np.argmax(path < welfare(0, theta_ss(0))))
print("welfare path under A's policy, generations 0,1,2,5,10,20,30:", [round(float(path[t]), 4) for t in [0, 1, 2, 5, 10, 20, 30]])
print("half-life of formation deviations ln2/-ln(rho) = %.3f generations" % (np.log(2) / -np.log(rho6)))
print("falls below no-transfer welfare %.4f in generation %d" % (welfare(0, theta_ss(0)), cross))
# an unconditional transfer (omega = 0) or an education-conditioned one (omega < 0)
for om_alt in [0., -.1]:
    phi_alt = M6 * om_alt * (rM - rT)
    th_alt = thbar - phi_alt * tauA
    print("transfer with omega=%.1f: long-run theta=%.4f welfare=%.4f" % (om_alt, th_alt, welfare(tauA, th_alt)))

# ---------------------------------------------------------------- 4.7
hdr("4.7 The spectral gap of the assumptions (N=8, Omega=8, kappa=0.5, Phi=1)")
N = 8; kap7 = .5; Phi7 = 1.; Om7 = 8.
Wce8 = np.ones((N, N)) / N


def layer_cluster(members):
    F = np.ones((N, N)) / N
    for i in members:
        F[i] = 0; F[i, members] = 1 / len(members)
    return F


def layer_blocks():
    F = np.zeros((N, N))
    for i in range(N):
        same = (np.arange(N) < 4) == (i < 4)
        F[i, same] = .95 / 4; F[i, ~same] = .05 / 4
    return F


def layer_hub(h, w=.5):
    F = np.ones((N, N)) * (1 - w) / N
    F[:, h] += w
    return F


def layer_core(core):
    F = np.zeros((N, N)); per = [j for j in range(N) if j not in core]
    for i in range(N):
        if i in core:
            F[i, core] = .9 / len(core); F[i, per] = .1 / len(per)
        else:
            F[i, core] = .7 / len(core); F[i, per] = .3 / len(per)
    return F


aligned = [layer_cluster([0, 1, 2]), layer_blocks(), layer_hub(0), layer_hub(1, .4), layer_core([0, 1])]
misaligned = [layer_cluster([0, 1, 2]), layer_blocks(), layer_hub(7), layer_hub(5, .4), layer_core([6, 7])]
names = ["W^U cluster {1,2,3}", "W^P blocks {1-4},{5-8}", "W^I hub", "market hub", "core-periphery"]
for tag, feats in [("aligned", aligned), ("misaligned", misaligned)]:
    print("--", tag)
    out = []
    for k in range(6):
        Wk = (Wce8 + sum(feats[:k])) / (k + 1) if k else Wce8
        assert np.allclose(Wk.sum(1), 1)
        p = perron(Wk)
        x = Om7 / N + (p - 1 / N) * Phi7
        Lp = kap7 * N / 2 * Phi7 ** 2 * var(p)
        out.append((gap(Wk), neff(p), Lp))
        print("k=%d %-24s Gamma=%.4f  N_eff=%.4f  L_power=%.5f  max kappa x=%.4f  pi_max=%.4f"
              % (k, names[k - 1] if k else "W^CE", gap(Wk), neff(p), Lp, kap7 * x.max(), p.max()))
        assert kap7 * x.max() < 1
    Ls = np.array([o[2] for o in out]); res47 = res47 if 'res47' in dir() else {}; res47[tag] = Ls[-1]
    print("   increments of L_power:", np.diff(Ls).round(5), " second differences:", np.diff(Ls, 2).round(5))
print("aligned/misaligned final L ratio = %.3f" % (res47["aligned"] / res47["misaligned"]))
# rank-one hub layers: adding 1v' to row-stochastic W scales non-unit eigenvalues by k/(k+1)
W2 = (Wce8 + aligned[0] + aligned[1]) / 3; W3 = (3 * W2 + aligned[2]) / 4
ev2 = np.sort(abs(np.linalg.eigvals(W2)))[::-1][1:]; ev3 = np.sort(abs(np.linalg.eigvals(W3)))[::-1][1:]
print("non-unit |eigenvalues| of W3 = 0.75 * those of W2:", np.allclose(ev3, .75 * ev2), " Gamma(W3) = 1 - 0.75*(1-Gamma(W2)) = %.4f" % (1 - .75 * (1 - gap(W2))))
print("\nAll Chapter 4 solution checks passed.")
