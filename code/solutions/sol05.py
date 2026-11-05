"""Solutions manual, Chapter 5: every number used in solutions/sol05.tex.
Run from the repository root:  python3 code/solutions/sol05.py
numpy throughout; pandas only to read data/processed (Exercise 5.8), through the
repository's own SUT -> IOT conversion (data/compute/compute_bolivia_io.py)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True)


def hdr(s):
    print("\n== " + s)


def U(x, k=1.0):
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


def rownorm(K, selfloop=True):
    K = np.asarray(K, float).copy()
    for i in range(len(K)):
        if K[i].sum() == 0 and selfloop:
            K[i, i] = 1.0
    return K / K.sum(1, keepdims=True)


# ---------------------------------------------------------------- 5.1
hdr("5.1 Utility externalities (Example 5.2)")
WU = np.array([[0, .30, .02], [.25, 0, .05], [.02, .05, 0]])
g = .5; x = np.array([.6, .4, .8]); u = U(x)
V = u + g * WU @ u
print("U(x) =", u, " V =", V.round(4), " (V1 check 0.473)")
print("externality terms V-U:", (V - u).round(4))
print("V2 = 0.32 + 0.5(0.25*0.42 + 0.05*0.48) = %.4f ; V3 = 0.48 + 0.5(0.02*0.42 + 0.05*0.32) = %.4f" % (V[1], V[2]))
c = WU.sum(0)
for gg in [.5, 1.0, .1, 0.0]:
    al = 1 + gg * c
    nu = (3 - 1.8) / np.sum(1 / al)
    xs = 1 - nu / al
    H = 3 / np.sum(1 / al); m = .6
    Lext = (1 - m) ** 2 / 2 * np.sum((al - H) ** 2 / al)
    gain = np.sum(al * U(xs)) - np.sum(al * U(np.full(3, .6)))
    print("gamma=%.1f c=%s alpha=%s nu=%.5f x_social=%s sum=%.4f  gain over equal split=%.6f (formula %.6f)"
          % (gg, c, al.round(4), nu, xs.round(4), xs.sum(), gain, Lext))
    assert np.isclose(xs.sum(), 1.8) and np.isclose(gain, Lext)
print("sum U at x_social vs equal split (gamma=0.5):", np.sum(U(1 - ((3 - 1.8) / np.sum(1 / (1 + .5 * c))) / (1 + .5 * c))).round(6), np.sum(U(np.full(3, .6))).round(6))

# ---------------------------------------------------------------- 5.2
hdr("5.2 A factory town (F, W, S, M)")
lab = ["Owner", "Worker", "Shopkeeper", "Mayor"]
KU = np.array([[0, .1, .1, .3], [.1, 0, .4, .1], [.1, .3, 0, .1], [.3, .2, .2, 0]])
KP = np.array([[.1, .5, .05, .1], [.6, .2, .1, .1], [.3, 0, .3, .1], [.5, .1, .2, .2]])
KI = np.array([[.5, .1, .1, .3], [.3, .2, .3, .2], [.2, .3, .3, .2], [.4, .1, .2, .3]])
WU2, WP2, WI2 = rownorm(KU), rownorm(KP), rownorm(KI)
for n_, M in [("W^U", WU2), ("W^P", WP2), ("W^I", WI2)]:
    print(n_, "=\n", M.round(4), " pi(layer) =", perron(M).round(4))
Wc = (WU2 + WP2 + WI2) / 3
pc = perron(Wc)
print("composite W =\n", Wc.round(4), "\n pi =", dict(zip(lab, pc.round(4))), " N_eff=%.4f  Gamma=%.4f" % (neff(pc), gap(Wc)))
print("column sums of composite:", Wc.sum(0).round(4))
KP2 = KP.copy(); KP2[:, 0] = 0; KP2[0, :] = 0
WP3 = rownorm(KP2)
Wc2 = (WU2 + WP3 + WI2) / 3; pc2 = perron(Wc2)
print("factory closes: W^P =\n", WP3.round(4), "\n pi =", dict(zip(lab, pc2.round(4))), " change", (pc2 - pc).round(4),
      " N_eff=%.4f Gamma=%.4f" % (neff(pc2), gap(Wc2)))

# ---------------------------------------------------------------- 5.3
hdr("5.3 DeGroot learning")
WI = np.array([[.5, .3, .2], [.1, .6, .3], [.2, .2, .6]])
print("row sums:", WI.sum(1))
for t in [2, 5]:
    P = np.linalg.matrix_power(WI, t)
    print("(W^I)^%d =\n" % t, P.round(4), " max column spread:", (P.max(0) - P.min(0)).round(4))
p = perron(WI); print("pi =", p.round(4), " eig:", np.round(eigs(WI), 4), " |lambda2| = %.4f  Gamma = %.4f" % (abs(eigs(WI)[1]), gap(WI)))
M = np.vstack([(WI.T - np.eye(3))[:2], np.ones(3)]); pex = np.linalg.solve(M, [0, 0, 1])
print("pi exact:", pex, " x 43 =", (pex * 43).round(6))
b0 = np.array([0, .5, 1.])
bstar = p @ b0
for t in [1, 2, 5, 10, 50]:
    bt = np.linalg.matrix_power(WI, t) @ b0
    print("t=%2d b=%s  max|b-b*|=%.2e" % (t, bt.round(4), abs(bt - bstar).max()))
print("consensus pi'b0 = %.4f" % bstar)
# power iteration of Appendix 5B
v = np.ones(3) / 3
for t in range(1, 6):
    v = v @ WI
    if t in (1, 2, 5):
        print("power iteration t=%d" % t, v.round(3))

# ---------------------------------------------------------------- 5.4
hdr("5.4 The phone as structural change (Mother, Friend, Cousin, TikTok, Teacher)")
lab5 = ["Mother", "Friend", "Cousin", "TikTok", "Teacher"]
pre = np.array([[.6, .2, 0, 0, .2], [.3, .5, 0, 0, .2], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [.2, .1, 0, 0, .7]])
post = np.array([[.55, .15, .05, .05, .2], [.2, .4, .05, .25, .1], [.1, 0, .8, .1, 0], [0, .05, .05, .85, .05], [.2, .1, .05, .05, .6]])
assert np.allclose(pre.sum(1), 1) and np.allclose(post.sum(1), 1)
b0 = np.array([.3, .3, .8, .6, .7])
Lpre = np.array([.5, .3, 0, 0, .2]); Lpost = np.array([.35, .25, .15, .15, .10])
ev = eigs(pre)
print("pre-phone: eigenvalues", np.round(ev, 4), " multiplicity of 1:", int(np.sum(np.isclose(ev, 1))))
Pt = np.linalg.matrix_power(pre, 2000)
print("pre-phone limit beliefs (W^t b0, t->inf):", (Pt @ b0).round(4))
sub = pre[np.ix_([0, 1, 4], [0, 1, 4])]; ps = perron(sub)
print("Potosi class {Mother, Friend, Teacher} Perron vector:", ps.round(4), " consensus = %.4f" % (ps @ b0[[0, 1, 4]]))
ppost = perron(post)
print("post-phone: pi =", dict(zip(lab5, ppost.round(4))), " N_eff=%.4f Gamma=%.4f  consensus=%.4f" % (neff(ppost), gap(post), ppost @ b0))
lim_pre = Pt @ b0
print("Luisa post row weight on Potosi class: %.2f" % Lpost[[0, 1, 4]].sum())
print("Luisa (outside, listens with her row): pre network & pre row: %.4f ; pre network & post row: %.4f ; post network & post row: %.4f ; post network & pre row: %.4f"
      % (Lpre @ lim_pre, Lpost @ lim_pre, Lpost @ (np.linalg.matrix_power(post, 2000) @ b0), Lpre @ (np.linalg.matrix_power(post, 2000) @ b0)))
# Luisa as a 6th agent who nobody listens to: same as outside listener
W6 = np.zeros((6, 6)); W6[:5, :5] = post; W6[5, :5] = Lpost
print("Luisa as 6th agent (nobody listens to her): pi =", perron(W6).round(4))

# ---------------------------------------------------------------- 5.5
hdr("5.5 When networks align (village)")
Wv = np.array([[.4, .3, .2, .1], [.3, .4, .1, .2], [.2, .1, .4, .3], [.1, .2, .3, .4]])
print("village W row sums", Wv.sum(1), " col sums", Wv.sum(0), " pi", perron(Wv).round(4), " Gamma %.4f" % gap(Wv))
for al in [(1 / 3, 1 / 3, 1 / 3), (.2, .5, .3)]:
    assert np.allclose(al[0] * Wv + al[1] * Wv + al[2] * Wv, Wv)
KPv = Wv.copy(); imports = .15; KPv[0, 1] -= imports   # agent 1 replaces half of its purchases from agent 2
WPv = rownorm(KPv)
print("perturbed W^P row 1 (domestic shares) =", WPv[0].round(4))
for al in [(1 / 3, 1 / 3, 1 / 3), (.2, .5, .3)]:
    Wn = al[0] * Wv + al[1] * WPv + al[2] * Wv
    pn = perron(Wn)
    print("alpha=%s composite row 1 = %s  pi = %s  change = %s  N_eff=%.4f Gamma=%.4f"
          % (np.round(al, 3), Wn[0].round(4), pn.round(4), (pn - .25).round(4), neff(pn), gap(Wn)))
# Leontief view: technical coefficients A = 0.5 * K (domestic), imports leak
A0 = .5 * Wv; A1 = .5 * KPv
L0 = np.linalg.inv(np.eye(4) - A0.T); L1 = np.linalg.inv(np.eye(4) - A1.T)
print("multipliers (col sums) before:", L0.sum(0).round(4), " after:", L1.sum(0).round(4))
print("forward (row sums) before:", L0.sum(1).round(4), " after:", L1.sum(1).round(4))

# ---------------------------------------------------------------- 5.6
hdr("5.6 Prebisch in W^I (6x6; center 1-2, periphery 3-6)")


def build(e1, e2):
    W = np.zeros((6, 6))
    W[:2, :2] = (1 - e1) / 2; W[:2, 2:] = e1 / 4
    W[2:, :2] = e2 / 2; W[2:, 2:] = (1 - e2) / 4
    return W


b0 = np.array([0, 0, 1, 1, 1, 1.])
for e1, e2 in [(.1, .7), (.1, .3), (.4, .3), (.7, .1)]:
    W = build(e1, e2); p = perron(W)
    print("eps1=%.2f eps2=%.2f  pi=%s  center mass=%.4f (formula %.4f)  consensus=%.4f  Gamma=%.4f  N_eff=%.4f"
          % (e1, e2, p.round(4), p[:2].sum(), e2 / (e1 + e2), p @ b0, gap(W), neff(p)))
    assert np.isclose(p[:2].sum(), e2 / (e1 + e2))
W = build(.1, .7)
print("check rows of baseline:\n", W.round(4))
bt = np.linalg.matrix_power(W, 30) @ b0
print("beliefs after 30 rounds:", bt.round(4))

# ---------------------------------------------------------------- 5.7
hdr("5.7 Composite power in the La Paz model")
KU3 = np.array([[0, .3, 0], [.2, 0, 0], [0, 0, 0.]])
KP3 = np.array([[.1, .4, .1], [.3, .2, 0], [0, 0, .3]])
WI3 = np.array([[.3, .5, .2], [.4, .3, .3], [.2, .1, .7]])
WU3, WP3_ = rownorm(KU3), rownorm(KP3)
Wbase = .2 * WU3 + .5 * WP3_ + .3 * WI3
print("baseline pi =", perron(Wbase).round(4), " N_eff = %.4f" % neff(perron(Wbase)))
grid = [0, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95, .99, 1.0]
pis = []
for aP in grid:
    aU, aI = .4 * (1 - aP), .6 * (1 - aP)
    W = aU * WU3 + aP * WP3_ + aI * WI3
    ev = eigs(W)
    p = perron(W); pis.append(p[2])
    print("alpha_P=%.2f (alpha_U=%.2f alpha_I=%.2f) pi=%s  pi_3=%.4f  Gamma=%.4f" % (aP, aU, aI, p.round(4), p[2], gap(W)))
fine = np.linspace(0, 1, 1001)
vals = np.array([perron(.4 * (1 - a) * WU3 + a * WP3_ + .6 * (1 - a) * WI3)[2] for a in fine])
print("pi_3 monotone increasing on the 1001-point grid:", bool(np.all(np.diff(vals) > 0)))
print("min pi_3 over alpha_P = %.4f at alpha_P = %.3f ; max %.4f at %.3f" % (vals.min(), fine[vals.argmin()], vals.max(), fine[vals.argmax()]))
print("W^P alone: classes -> vendor absorbing; pi(W^P) =", perron(WP3_).round(4))
# (c) add Luisa (student, agent 4)
KU4 = np.array([[0, .3, 0, 0], [.2, 0, 0, 0], [0, 0, 0, .6], [0, 0, .5, 0]])
KP4 = np.array([[.1, .4, .1, 0], [.3, .2, 0, 0], [0, 0, .3, .1], [0, 0, 0, 0]])
KI4 = np.array([[.3, .5, .2, 0], [.4, .3, .3, 0], [.2, .1, .6, .1], [.1, .2, .3, .4]])
W4 = .2 * rownorm(KU4) + .5 * rownorm(KP4) + .3 * rownorm(KI4)
p4 = perron(W4)
print("with Luisa: W^U4 =\n", rownorm(KU4).round(4), "\n W^P4 =\n", rownorm(KP4).round(4), "\n W^I4 =\n", rownorm(KI4).round(4))
print(" composite W4 =\n", W4.round(4), "\n pi =", p4.round(4), " N_eff=%.4f  vendor+Luisa=%.4f" % (neff(p4), p4[2:].sum()))
# (d) vendor's utility self-loop -> uniform row
WU3u = WU3.copy(); WU3u[2] = 1 / 3
Wd = .2 * WU3u + .5 * WP3_ + .3 * WI3; pd_ = perron(Wd)
print("(d) uniform utility row for vendor: pi =", pd_.round(4), " pi_3 change = %.4f (from %.4f)" % (pd_[2] - perron(Wbase)[2], perron(Wbase)[2]))
WP3u = WP3_.copy(); WP3u[2] = 1 / 3
pdd = perron(.2 * WU3u + .5 * WP3u + .3 * WI3)
print("    and also uniform production row: pi =", pdd.round(4))

# ---------------------------------------------------------------- 5.8
hdr("5.8 From gas well to school (INE SUT 2019, base 2017, detail 39 industries)")
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path("data/compute").resolve()))
import compute_bolivia_io as cb  # noqa: E402

d = pd.read_csv("data/processed/ine_sut_long.csv", dtype={"row_code": str, "col_code": str})
t = cb.load_b17(d, 2019)
io = cb.sut_to_iot(t)
codes = t["inds"]; xi = io["x"]; Zd, Zt = io["Zd"], io["Zt"]
idx = {c: k for k, c in enumerate(codes)}
name = {"05": "gas extraction", "17": "refining", "23": "electricity", "09": "milling & bakery",
        "06": "mining", "35": "non-market education"}


def link(buyer, seller):
    j, i = idx[buyer], idx[seller]
    a_tot = Zt[i, j] / xi[j]; a_dom = Zd[i, j] / xi[j]
    wp = Zd[i, j] / Zd[:, j].sum()
    return a_tot, a_dom, (1 - a_dom / a_tot) if a_tot > 0 else np.nan, wp


for buyer, seller in [("17", "05"), ("23", "05"), ("23", "17"), ("09", "23"), ("09", "09"), ("09", "01"), ("06", "09"), ("06", "17"), ("06", "23")]:
    a_tot, a_dom, imp, wp = link(buyer, seller)
    print("a[%s <- %s]: %s buys from %s: total %.4f  domestic %.4f  imported share %.4f  W^P (domestic share of buyer's inputs) %.4f"
          % (buyer, seller, name.get(buyer, buyer), name.get(seller, seller), a_tot, a_dom, imp, wp))
for c_ in ["05", "17", "23", "09", "06"]:
    j = idx[c_]
    print("industry %s %-22s output %.0f  dom int/x %.4f  imported share of int inputs %.4f  export share %.4f"
          % (c_, name[c_], xi[j], Zd[:, j].sum() / xi[j], io["imp_int"][j] / (io["imp_int"][j] + Zd[:, j].sum()), io["exports"][j] / xi[j]))
plab = d[(d.base == 2017) & (d.year == 2019) & (d.table == "supply")].drop_duplicates("row_code").set_index("row_code").row_label
for pc_ in ["01", "05", "10", "18", "22"]:
    k = t["prods"].index(pc_)
    print("product %s %-45s imports/(output+imports) = %.4f" % (pc_, plab[pc_][:45], t["m"][k] / (t["q"][k] + t["m"][k])))
# household consumption of bakery products and education
F = d[(d.base == 2017) & (d.year == 2019) & (d.table == "use") & d.col_code.isna() & (d.row_code == "10")]
print("product 10 (molineria/panaderia) final uses:", dict(zip(F.col_label.str[:28], F.value.round(1))))
# (c) two blocks from the 12-sector domestic Z
Z12 = pd.read_csv("data/processed/bolivia_Z/Z_2017_2019_book12.csv", index_col=0)
Z12 = Z12.drop(index="Domestic services", columns="Domestic services")
labs = list(Z12.index); Zv = Z12.values
WPd = rownorm(Zv.T)
p11 = perron(WPd)
for blockname, members in [("gas = {Hydrocarbons}", ["Hydrocarbons"]), ("extractive = {Hydrocarbons, Mining}", ["Hydrocarbons", "Mining"])]:
    bm = np.array([l in members for l in labs])
    # aggregated 2-block W^P from flows (buyer rows)
    K = Zv.T
    eps1 = K[np.ix_(bm, ~bm)].sum() / K[bm].sum()
    eps2 = K[np.ix_(~bm, bm)].sum() / K[~bm].sum()
    rowcross = WPd[~bm][:, bm].sum(1)
    print("%s: eps1 (block buys from rest) = %.4f ; eps2 (rest buys from block) = %.4f ; predicted block power eps2/(eps1+eps2) = %.4f ; actual pi mass (11 sectors) = %.4f"
          % (blockname, eps1, eps2, eps2 / (eps1 + eps2), p11[bm].sum()))
    print("    per-sector weight of the rest on the block: min %.4f max %.4f (Prop 6.12 needs these equal)" % (rowcross.min(), rowcross.max()))
print("pi (11 sectors):", dict(zip(labs, p11.round(4))))
print("\nAll Chapter 5 solution checks passed.")
