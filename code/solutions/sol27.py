"""Solutions, Chapter 27 (Algorithmic W^I).
Prints every number used in solutions/sol27.tex.  Run from the repo root:
    python3 code/solutions/sol27.py      (numpy + pandas; reads data/processed for WDI internet use)
All network parameters are illustrative."""
import itertools
import pathlib
import numpy as np
import pandas as pd

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "processed"
np.set_printoptions(precision=4, suppress=True)


def stats(W):
    W = W / W.sum(1, keepdims=True)
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi /= pi.sum()
    return 1 - abs(ev[o[1]]), 1 / np.sum(pi ** 2), pi


def two_group(lam, m=3):
    W = np.full((2 * m, 2 * m), lam / m)
    W[:m, :m] = (1 - lam) / m
    W[m:, m:] = (1 - lam) / m
    return W


def sinkhorn(K, iters=5000):
    W = K.copy()
    for _ in range(iters):
        W = W / W.sum(1, keepdims=True)
        W = W / W.sum(0, keepdims=True)
    return W / W.sum(1, keepdims=True)


# ---------------------------------------------------------------------------
print("=== 27.1 Four groups (UM, UW, RM, RW): gender factor lg = 0.5, geography factor lr = 0.3 ===")
lg, lr = 0.5, 0.3
gender = np.array([0, 1, 0, 1]); geo = np.array([0, 0, 1, 1])
K = np.array([[(lg if gender[i] != gender[j] else 1) * (lr if geo[i] != geo[j] else 1) for j in range(4)] for i in range(4)])
W0 = K / K.sum(1, keepdims=True)
crossg = lambda W: np.mean([(W[i] * (gender != gender[i])).sum() for i in range(4)])
crossr = lambda W: np.mean([(W[i] * (geo != geo[i])).sum() for i in range(4)])
G0, N0, p0 = stats(W0)
print("W0 =\n", W0)
print(f"baseline: cross-gender share {crossg(W0):.4f}, cross-geography share {crossr(W0):.4f}, Gamma {G0:.4f}, Neff {N0:.4f}")
print(f"eigenvalues: {np.round(np.sort(np.real(np.linalg.eigvals(W0)))[::-1], 4)}")
assert abs(G0 - (1 - (1 - lr) / (1 + lr))) < 1e-12 and abs(N0 - 4) < 1e-9
M = np.where(np.eye(4) > 0, 1.2, 0.85)
W1 = W0 * M; W1 /= W1.sum(1, keepdims=True)
G1, N1, p1 = stats(W1)
print("W1 =\n", W1)
print(f"after (+20% own group, -15% others): cross-gender {crossg(W1):.4f} (d = {crossg(W1)-crossg(W0):+.4f}), "
      f"cross-geography {crossr(W1):.4f} (d = {crossr(W1)-crossr(W0):+.4f}), Gamma {G1:.4f} (d = {G1-G0:+.4f}), "
      f"Neff {N1:.4f} (d = {N1-N0:+.4f})")
print(f"eigenvalues after: {np.round(np.sort(np.real(np.linalg.eigvals(W1)))[::-1], 4)}")
assert np.allclose(W1, W1.T) and abs(N1 - 4) < 1e-9 and G1 < G0
# unequal group sizes: attention proportional to source-group size
pop = np.array([0.3, 0.3, 0.2, 0.2])
Wp0 = K * pop; Wp0 /= Wp0.sum(1, keepdims=True)
Wp1 = Wp0 * M; Wp1 /= Wp1.sum(1, keepdims=True)
Gp0, Np0, pp0 = stats(Wp0); Gp1, Np1, pp1 = stats(Wp1)
print(f"unequal sizes p = {pop}: Neff {Np0:.4f} -> {Np1:.4f}; pi {pp0} -> {pp1}")
assert abs(Np1 - Np0) > 1e-3

# ---------------------------------------------------------------------------
print("\n=== 27.2 (a) user-specific engagement: weighted-mean condition, random check ===")
rng = np.random.default_rng(272)
viol = 0
for _ in range(20000):
    n = int(rng.integers(3, 9))
    grp = rng.integers(0, 2, n); grp[0] = 0; grp[1] = 1
    w = rng.dirichlet(np.ones(n)); e = rng.uniform(0.1, 3, n)
    same = grp == grp[0]
    lam = w[~same].sum()
    es = (w[same] * e[same]).sum() / (1 - lam); ec = (w[~same] * e[~same]).sum() / lam
    lamp = (w[~same] * e[~same]).sum() / (w * e).sum()
    assert abs(lamp - lam * ec / (lam * ec + (1 - lam) * es)) < 1e-12
    assert (lamp < lam) == (es > ec)
print("checked: lambda' = lambda ebar_c/(lambda ebar_c + (1-lambda) ebar_s) and lambda' < lambda <=> ebar_s > ebar_c (20000 draws)")

print("\n=== 27.2 (b) two attributes (the code/ch27.py example) ===")
types = [("U", "L")] * 8 + [("U", "R")] * 2 + [("R", "L")] * 2 + [("R", "R")] * 8
n = len(types)
K2 = np.array([[0.0 if i == j else (4.0 if types[i][0] == types[j][0] else 1.0) for j in range(n)] for i in range(n)])
V0 = K2 / K2.sum(1, keepdims=True)
boost = np.array([[3.0 if types[i][1] == types[j][1] else 1.0 for j in range(n)] for i in range(n)])
V1 = V0 * boost; V1 /= V1.sum(1, keepdims=True)
cb = lambda W: np.mean([W[i, [j for j in range(n) if types[j][1] != types[i][1]]].sum() for i in range(n)])
sub = lambda W, t: np.mean([W[i, [j for j in range(n) if types[j][0] != t[0]]].sum() for i in range(n) if types[i] == t])
print(f"cross-belief exposure {cb(V0):.3f} -> {cb(V1):.3f}")
for t in (("R", "L"), ("U", "R"), ("R", "R"), ("U", "L")):
    print(f"  cross-geography exposure of {t}: {sub(V0, t):.3f} -> {sub(V1, t):.3f}")
print(f"  analytic (R,L): 10/46 = {10/46:.3f} -> 26/70 = {26/70:.3f}; (R,R): 10/46 -> 14/106 = {14/106:.3f}")
assert round(cb(V0), 3) == 0.426 and round(cb(V1), 3) == 0.218
assert round(sub(V1, ("R", "L")), 3) == 0.371 and round(sub(V0, ("R", "L")), 3) == 0.217

print("\n=== 27.2 (c) top-K ranking ===")
w = np.array([0.45, 0.45, 0.10]); e = np.array([2.0, 0.2, 1.0]); same = np.array([True, True, False])
lam = w[~same].sum()
es = (w[same] * e[same]).sum() / (1 - lam); ec = (w[~same] * e[~same]).sum() / lam
mult = (w * e)[~same].sum() / (w * e).sum()
top = np.argsort(-e)[:2]
wk = np.zeros(3); wk[top] = w[top]; wk /= wk.sum()
print(f"user: weights {w}, engagement {e}; ebar_s = {es:.3f} > ebar_c = {ec:.3f}")
print(f"  multiplicative: lambda {lam:.3f} -> {mult:.4f}; top-2: lambda -> {wk[~same].sum():.4f}")
assert mult < lam < wk[~same].sum()
Wt = two_group(0.22)
Wtop = np.where(np.kron(np.eye(2), np.ones((3, 3))) > 0, Wt, 0.0); Wtop /= Wtop.sum(1, keepdims=True)
ev = np.sort(np.abs(np.linalg.eigvals(Wtop)))[::-1]
print(f"  two equal groups, group-level engagement, K = 3 (= own-group size): |lambda_2| = {ev[1]:.3f}, Gamma = {1-ev[1]:.3f}")

# ---------------------------------------------------------------------------
print("\n=== 27.3 Internet use, % of population (WDI IT.NET.USER.ZS) ===")
wdi = pd.read_csv(DATA / "wdi_long.csv")
net = wdi[(wdi.country_iso3 == "BOL") & (wdi.indicator == "IT.NET.USER.ZS")].set_index("year").value
for y in (2010, 2019, 2023, 2024):
    print(f"  {y}: {net[y]:.1f}%")
print(f"  illustrative upper bound on the digital share of national W^I (omega_d <= 0.50): {0.5*net[2024]:.1f}% (2024)")

# ---------------------------------------------------------------------------
print("\n=== 27.4 Objective design: engagement vs Neff on a 6-user network ===")
m = 3
W0 = two_group(0.22, m)
grp = np.array([0] * m + [1] * m)
Eng = np.where(grp[:, None] == grp[None, :], 2.07, 1.0)
Eng[:, 0] *= 3.0                              # source 1 is viral
eng = lambda W: np.mean((W * Eng).sum(1))
lam_mean = lambda W: np.mean([(W[i] * (grp != grp[i])).sum() for i in range(2 * m)])
WA = W0 * Eng; WA /= WA.sum(1, keepdims=True)
WB = sinkhorn(W0 * Eng)
res = {}
for name, W in (("organic", W0), ("engagement ranking", WA), ("balanced ranking (Sinkhorn)", WB)):
    G_, N_, p_ = stats(W)
    res[name] = (eng(W), N_, G_, lam_mean(W))
    print(f"  {name:28s}: engagement {eng(W):.3f}, Neff {N_:.3f}, Gamma {G_:.3f}, mean cross-group share {lam_mean(W):.3f}")
assert res["balanced ranking (Sinkhorn)"][1] > 5.999 and res["engagement ranking"][1] < 6
assert res["engagement ranking"][0] > res["balanced ranking (Sinkhorn)"][0] > res["organic"][0]
best = max(itertools.permutations(range(6)), key=lambda s: sum(Eng[i, s[i]] for i in range(6)))
Pm = np.zeros((6, 6)); Pm[range(6), best] = 1
print(f"  best permutation matrix (unregularized max over doubly stochastic W): engagement {eng(Pm):.3f}, "
      f"Neff {1/np.sum(np.full(6, 1/6)**2):.0f}, max |eigenvalue| other than 1 = 1 -> Gamma = 0")
for tau in (0.5, 1.0, 2.0, 4.0):
    Wt_ = sinkhorn(W0 * Eng ** tau)
    G_, N_, _ = stats(Wt_)
    print(f"    tau = {tau}: balanced W0 o e^tau: engagement {eng(Wt_):.3f}, Neff {N_:.3f}, Gamma {G_:.3f}")

# ---------------------------------------------------------------------------
print("\n=== 27.5 AI mentoring: rural cross-group share 0.22 -> 0.32 ===")
Wm = two_group(0.22, m)
Wm[m:, :m] = 0.32 / m; Wm[m:, m:] = 0.68 / m
G_b, N_b, _ = stats(two_group(0.22, m)); G_m, N_m, p_m = stats(Wm)
print(f"  uniform mentors: Gamma {G_b:.3f} -> {G_m:.3f}, Neff {N_b:.3f} -> {N_m:.3f}, pi_urban total {p_m[:m].sum():.4f}")
Wh = two_group(0.22, m)
Wh[m:, :m] = 0.22 / m; Wh[m:, 0] += 0.10; Wh[m:, m:] = 0.68 / m
G_h, N_h, p_h = stats(Wh)
print(f"  one star mentor: Gamma {G_h:.3f}, Neff {N_h:.3f}, pi_star {p_h[0]:.4f}")
assert G_m > G_b and N_m < N_b and N_h < N_m
rho_, mudr = 0.85, 0.09
d1 = mudr * 0.10
print(f"  formation: one generation d(theta) = mu*dlam*dr = {d1:.4f}; steady state = {d1/(1-rho_):.4f}")
from math import erf, sqrt


def zq(p):
    lo, hi = -10.0, 10.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.5 * (1 + erf(mid / sqrt(2))) < p:
            lo = mid
        else:
            hi = mid
    return mid


za, zb = zq(1 - 0.05 / 2), zq(0.80)
icc, nper, mde, P = 0.05, 20, 0.10, 0.5
Deff = icc + (1 - icc) / nper
J = (za + zb) ** 2 * Deff / (P * (1 - P) * mde ** 2)
Jr = int(2 * np.ceil(J / 2))
print(f"  power: z = {za:.4f} + {zb:.4f} = {za+zb:.4f}; design term {Deff:.4f}; J = {J:.1f} -> {Jr} communities "
      f"({Jr//2} per arm), {Jr*nper} entrepreneurs")
assert 300 < J < 310

# ---------------------------------------------------------------------------
print("\n=== 27.7 Four regulatory regimes on the 6-user network ===")
R = W0 * np.where(grp[:, None] == grp[None, :], 2.07, 1.0); R /= R.sum(1, keepdims=True)
a = 0.3
Wa = (1 - a) * R; Wa[:, 0] += a
zeta = 0.25
Wb = (1 - zeta) * Wa + zeta * W0
Wc = Wa.copy()
for i in range(2 * m):
    cross = grp != grp[i]
    lam_i = Wc[i, cross].sum()
    if lam_i < 0.3:
        Wc[i, cross] *= 0.3 / lam_i
        Wc[i, ~cross] *= 0.7 / (1 - lam_i)
Wd = sinkhorn(Wa)
Wcd = sinkhorn(Wc)
out = {}
for name, W in (("organic", W0), ("(a) no regulation", Wa), ("(b) transparency", Wb),
                ("(c) diversity mandate 0.3", Wc), ("(d) public-interest", Wd), ("(c)+(d) combined", Wcd)):
    G_, N_, p_ = stats(W)
    out[name] = (N_, G_)
    print(f"  {name:26s}: Neff {N_:.3f}, Gamma {G_:.3f}, mean cross-group share {lam_mean(W):.3f}, "
          f"min cross share {min((W[i]*(grp!=grp[i])).sum() for i in range(2*m)):.3f}, pi_1 {p_[0]:.3f}, engagement {eng(W):.3f}")
assert out["(d) public-interest"][0] > 5.999 and out["(a) no regulation"][0] < out["(b) transparency"][0]
print("all chapter-27 solution checks passed")
