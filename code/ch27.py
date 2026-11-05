"""Chapter 27: engagement reweighting, Gamma vs Neff, bounded confidence,
and the sourced internet-use shares (data/processed/wdi_long.csv).
Run from the repo root: python3 code/ch27.py  (numpy + csv)."""
import csv
import pathlib
import numpy as np

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "processed"


def stats(W):
    W = W / W.sum(1, keepdims=True)
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi /= pi.sum()
    return 1 - abs(ev[o[1]]), 1 / np.sum(pi ** 2)


def two_group(lam, m=3):
    W = np.full((2 * m, 2 * m), lam / m)
    W[:m, :m] = (1 - lam) / m
    W[m:, m:] = (1 - lam) / m
    return W


print("== Engagement reweighting: lambda' = lambda e_c / (lambda e_c + (1-lambda) e_s) ==")
lam = 0.22
for target in (0.17, 0.12):
    r = (lam / target - lam) / (1 - lam)        # required engagement premium e_s/e_c
    lp = lam / (lam + (1 - lam) * r)
    print(f"premium e_s/e_c = {r:.2f} -> lambda' = {lp:.3f} (drop {lam-lp:.2f} = {100*(lam-lp)/lam:.0f}%)")
    assert abs(lp - target) < 1e-12
# formation gap in steady state: mu/(1-rho) * dlam * dr, with mu*dr = 0.09, rho = 0.85
for dl in (0.05, 0.10):
    print(f"steady-state aspiration gap widening for dlambda = {dl}: {1/(1-0.85)*dl*0.09:.3f}")

print("\n== Gamma vs Neff: echo chamber vs viral hub (6 nodes, 2 groups) ==")
base = two_group(0.22)
echo = two_group(0.12)
mats = {"baseline": base, "echo": echo}
for a in (0.3, 0.6):
    V = (1 - a) * base.copy()
    V[:, 0] += a
    mats[f"hub a={a}"] = V
for k, M in mats.items():
    G, ne = stats(M)
    print(f"{k:10s}: Gamma = {G:.3f}, Neff = {ne:.2f}")
assert abs(stats(base)[0] - 0.44) < 1e-9 and abs(stats(echo)[0] - 0.24) < 1e-9
print("two equal groups: Gamma = 2*lambda (lambda <= 1/2), Neff = N (doubly stochastic)")

print("\n== Heatmap shading (blue!xx), xx = round(100 W_ij / max over panels) ==")
show = {"baseline": base, "echo": echo, "hub a=0.6": mats["hub a=0.6"]}
mx = max(M.max() for M in show.values())
for k, M in show.items():
    vals = sorted(set(np.round(100 * M / mx).astype(int).ravel()))
    print(f"{k}: distinct cells {np.unique(np.round(M,4))} -> shading {vals}")

print("\n== Two attributes: belief homophily can RAISE cross-geography exposure ==")
# types: (geo, belief). Urban mostly L, rural mostly R.
types = [("U", "L")] * 8 + [("U", "R")] * 2 + [("R", "L")] * 2 + [("R", "R")] * 8
n = len(types)
K = np.zeros((n, n))
for i, (gi, bi) in enumerate(types):
    for j, (gj, bj) in enumerate(types):
        if i != j:
            K[i, j] = 4.0 if gi == gj else 1.0   # organic feed: geographic homophily only
W0 = K / K.sum(1, keepdims=True)
boost = np.array([[3.0 if types[i][1] == types[j][1] else 1.0 for j in range(n)] for i in range(n)])
W1 = W0 * boost; W1 /= W1.sum(1, keepdims=True)
cross = lambda W: np.mean([W[i, [j for j in range(n) if types[j][0] != types[i][0]]].sum() for i in range(n)])
cross_belief = lambda W: np.mean([W[i, [j for j in range(n) if types[j][1] != types[i][1]]].sum() for i in range(n)])
print(f"cross-belief exposure: {cross_belief(W0):.3f} -> {cross_belief(W1):.3f} (falls)")
print(f"average cross-geography exposure: {cross(W0):.3f} -> {cross(W1):.3f}")
sub = lambda W, g, bel: np.mean([W[i, [j for j in range(n) if types[j][0] != g]].sum()
                                 for i in range(n) if types[i] == (g, bel)])
for g, bel in (("R", "L"), ("U", "R")):
    print(f"cross-geography exposure of ({g},{bel}) users: {sub(W0, g, bel):.3f} -> {sub(W1, g, bel):.3f} (rises)")
    assert sub(W1, g, bel) > sub(W0, g, bel)
assert cross_belief(W1) < cross_belief(W0)

print("\n== Primitive W(b*) forces consensus; bounded confidence gives clusters ==")
rng = np.random.default_rng(27)
b0 = np.sort(rng.uniform(0, 1, 60))
for eps in (0.05, 0.15, 0.30):
    b = b0.copy()
    for _ in range(500):
        A = (np.abs(b[:, None] - b[None, :]) <= eps).astype(float)
        b = (A / A.sum(1, keepdims=True)) @ b
    clusters = np.unique(np.round(b, 6))
    print(f"Hegselmann-Krause eps = {eps:.2f}: {len(clusters)} cluster(s) {np.round(clusters,3)}")
# a primitive W: fixed point must be consensus
Wp = rng.dirichlet(np.ones(5), size=5)
ev, vr = np.linalg.eig(Wp)
fixed = np.real(vr[:, np.argmin(np.abs(ev - 1))])
print("primitive W: eigenvector for eigenvalue 1 is constant:", np.allclose(fixed / fixed[0], 1))

print("\n== Proof checks for Proposition ch27:thm:echo ==")
rng = np.random.default_rng(27)
for _ in range(500):
    # (i) general row-stochastic W, two groups of random sizes, random engagement
    n = rng.integers(2, 9)
    grp = rng.integers(0, 2, n)
    W = rng.dirichlet(np.ones(n), size=n)
    ec = rng.uniform(0.1, 1.0); es = ec * rng.uniform(1.01, 4.0)
    E = np.where(grp[:, None] == grp[None, :], es, ec)
    Wp = W * E; Wp /= Wp.sum(1, keepdims=True)
    cross = grp[:, None] != grp[None, :]
    lam_i = (W * cross).sum(1); lam_p = (Wp * cross).sum(1)
    assert np.allclose(lam_p, lam_i * ec / (lam_i * ec + (1 - lam_i) * es))
    assert np.all(lam_p[(lam_i > 0) & (lam_i < 1)] < lam_i[(lam_i > 0) & (lam_i < 1)])
    # (ii) two equal uniform groups: same family, Gamma = 1-|1-2 lam|, doubly stochastic
    m = int(rng.integers(1, 6)); lam = rng.uniform(0.01, 0.99)
    T = two_group(lam, m)
    Tp = T * np.where(np.kron(np.eye(2), np.ones((m, m))) > 0, es, ec)
    Tp /= Tp.sum(1, keepdims=True)
    lp = lam * ec / (lam * ec + (1 - lam) * es)
    assert np.allclose(Tp, two_group(lp, m))
    G, ne = stats(T)
    assert abs(G - (1 - abs(1 - 2 * lam))) < 1e-9 and abs(ne - 2 * m) < 1e-6
    if lam <= 0.5:
        assert stats(Tp)[0] < G
print("checked: lambda' formula on random W, same-family reweighting, Gamma = 1-|1-2lambda|, Neff = N")


print("\n== Sourced: internet users, % of population (WDI IT.NET.USER.ZS) ==")
with open(DATA / "wdi_long.csv", newline="") as fh:
    net = {int(r["year"]): float(r["value"]) for r in csv.DictReader(fh)
           if r["country_iso3"] == "BOL" and r["indicator"] == "IT.NET.USER.ZS" and r["value"]}
for y, v in ((2010, 22.4), (2019, 47.5), (2023, 77.3)):
    assert abs(net[y] - v) <= 0.05 + 1e-9, (y, net[y], v)
    print(f"  {y}: {net[y]:.1f}%")
