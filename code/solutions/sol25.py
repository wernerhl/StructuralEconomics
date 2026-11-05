"""Solutions, Chapter 25 (Truncated Industrialization).
Prints every number used in solutions/sol25.tex.  Run from the repo root:
    python3 code/solutions/sol25.py      (numpy + pandas; reads data/processed)
Sourced values: data/processed (retrieved 2026-11-05; see data/README.md).
Synthetic product spaces below are NOT data: they only demonstrate the method."""
import pathlib
import numpy as np
import pandas as pd

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "processed"


def rownorm(K):
    return K / K.sum(1, keepdims=True)


def gap(W):
    ev = np.linalg.eigvals(W)
    ev = ev[np.argsort(-np.abs(ev))]
    return 1 - abs(ev[1])


def peel(phi):
    """Charikar greedy peeling for max e(S)/|S| on a symmetric weight matrix (zero diagonal)."""
    n = len(phi)
    alive = np.ones(n, bool)
    deg = phi.sum(1).astype(float)
    e = phi.sum() / 2
    best, bestS, k = e / n, alive.copy(), n
    order = []
    while k > 1:
        d = np.where(alive, deg, np.inf)
        v = int(np.argmin(d))
        e -= deg[v]
        alive[v] = False
        deg -= phi[:, v]
        k -= 1
        order.append(v)
        if e / k > best:
            best, bestS = e / k, alive.copy()
    return np.where(bestS)[0], best


# ---------------------------------------------------------------------------
print("=== 25.1 (a) The Goldberg objective on a dense weighted proximity matrix ===")
P, k, phiD, phiAll = 775, 100, 0.35, 0.08
obj_core = (k - 1) * phiD / 2
obj_all = (P - 1) * phiAll / 2
print(f"e(S)/|S| = (|S|-1) phibar_S / 2: core (100, 0.35) = {obj_core:.3f}; whole space (775, 0.08) = {obj_all:.3f}")
assert obj_all > obj_core
npairs = P * (P - 1) // 2
ncore = k * (k - 1) // 2
phi_rest = (phiAll * npairs - phiD * ncore) / (npairs - ncore)
print(f"implied mean proximity of the non-core pairs = {phi_rest:.4f}")
kmin = 2 * obj_all / phiD + 1
print(f"a set with mean proximity 0.35 must have more than {kmin:.0f} products to beat the whole space")

print("\n=== 25.1 (b) Synthetic demonstration (NOT data) ===")
rng = np.random.default_rng(2501)
core = np.arange(k)
A = rng.beta(0.6, 0.6 * (1 - phi_rest) / phi_rest, (P, P))
A[np.ix_(core, core)] = rng.beta(3.5, 6.5, (k, k))
phi = np.triu(A, 1); phi = phi + phi.T
dens = lambda S: phi[np.ix_(S, S)].sum() / (len(S) * (len(S) - 1))
print(f"synthetic means: core {dens(core):.3f}, all {dens(np.arange(P)):.3f}")
S1, v1 = peel(phi)
inter1 = len(np.intersect1d(S1, core))
print(f"weighted peeling: |S| = {len(S1)}, objective {v1:.3f}, contains {inter1} of the 100 core products, "
      f"density ratio {dens(S1)/dens(np.arange(P)):.2f}")
assert len(S1) > 500
thr = 0.25
B = (phi >= thr).astype(float); np.fill_diagonal(B, 0)
S2, v2 = peel(B)
inter2 = len(np.intersect1d(S2, core))
jac = inter2 / len(np.union1d(S2, core))
print(f"thresholded at phi >= {thr}: |S| = {len(S2)}, objective {v2:.3f}, core products found {inter2}, Jaccard {jac:.3f}, "
      f"density ratio (weighted) {dens(S2)/dens(np.arange(P)):.2f}")
print(f"share of pairs above threshold: core {100*B[np.ix_(core, core)].sum()/(k*(k-1)):.1f}%, "
      f"non-core {100*B[k:, k:].sum()/((P-k)*(P-k-1)):.2f}%")
assert jac > 0.9

# ---------------------------------------------------------------------------
eci = pd.read_csv(DATA / "eci_facts.csv")
E = eci.set_index(["country", "year"])
print("\n=== 25.2 ECI (Growth Lab HS92), 2019, five countries ===")
five = ["BOL", "BRA", "CHL", "KOR", "DEU"]
vals = {c: E.loc[(c, 2019), "eci_hs92"] for c in five}
for c in sorted(five, key=lambda c: -vals[c]):
    print(f"  {c}: ECI {vals[c]:+.3f}, rank {int(E.loc[(c, 2019), 'eci_rank_hs92'])}")
assert sorted(five, key=lambda c: -vals[c]) == ["KOR", "DEU", "BRA", "CHL", "BOL"]

# ---------------------------------------------------------------------------
print("\n=== 25.3 Gamma (2019, book12_isic, domestic, input shares) and ECI ===")
gc = pd.read_csv(DATA / "gamma_comparison.csv")
g19 = gc[(gc.year == 2019) & (gc.aggregation == "book12_isic") & (gc.normalization == "input_shares")].copy()
g19["iso"] = g19.country.str[:3]
G = g19.set_index("iso").Gamma
print(f"  {len(G)} economies with 12-sector Gamma: min {G.min():.3f} ({G.idxmin()}), max {G.max():.3f} ({G.idxmax()}), "
      f"mean {G.mean():.3f}, sd {G.std(ddof=1):.3f}")
both = sorted(set(G.index) & set(eci[eci.year == 2019].country))
x = np.array([E.loc[(c, 2019), "eci_hs92"] for c in both])
y = np.array([G[c] for c in both])
r = np.corrcoef(x, y)[0, 1]
b, a = np.polyfit(x, y, 1)
res = y - (a + b * x)
se_b = np.sqrt(res @ res / (len(x) - 2) / ((x - x.mean()) @ (x - x.mean())))
rk = lambda v: np.argsort(np.argsort(v))
rs = np.corrcoef(rk(x), rk(y))[0, 1]
print(f"  {len(both)} with ECI too: {both}")
print(f"  corr(Gamma, ECI) = {r:.3f}, Spearman = {rs:.3f}; OLS Gamma = {a:.3f} + {b:.4f} ECI (se {se_b:.4f})")
for c in both:
    print(f"    {c}: Gamma {G[c]:.3f}, ECI {E.loc[(c, 2019), 'eci_hs92']:+.3f}")
o = pd.read_csv(DATA / "oecd_io_summary.csv")
o19 = o[(o.year == 2019) & (o.flows == "domestic") & (o.normalization == "input_shares") & (o.edition == "ICIO2023")]
g12 = o19[o19.aggregation == "book12_isic"].set_index("country").Gamma
gdet = o19[o19.aggregation.str.startswith("detail")].set_index("country").Gamma
cs = sorted(set(g12.index) & set(gdet.index))
u, v = g12[cs].values, gdet[cs].values
print(f"  ICIO economies with both resolutions: {len(cs)}; Pearson(12-sector, detail) = {np.corrcoef(u, v)[0,1]:.3f}, "
      f"Spearman = {np.corrcoef(rk(u), rk(v))[0,1]:.3f}; mean 12-sector {u.mean():.3f}, mean detail {v.mean():.3f}")
for c in ("KOR", "DEU", "USA", "MEX", "CHL", "CHN"):
    print(f"    {c}: 12-sector {g12[c]:.3f}, detail {gdet[c]:.3f}")
bs = pd.read_csv(DATA / "bolivia_io_summary.csv")
b39 = bs[(bs.base == 2017) & (bs.year == 2019) & (bs.aggregation == "detail39") & (bs.flows == "domestic")
         & (bs.normalization == "input_shares")].Gamma.iloc[0]
print(f"  Bolivia detail39 (INE base 2017, 2019): {b39:.3f}")
lo12 = g12[cs].idxmin(); lodet = gdet[cs].idxmin()
print(f"  lowest 12-sector: {lo12}; lowest detail: {lodet}")

# ---------------------------------------------------------------------------
print("\n=== 25.4 Counterfactual: truncation arithmetic and the two toy product spaces (Ch25 code) ===")
print(f"  tau: 1 - 2/100 = {1-2/100:.2f} -> 1 - 22/100 = {1-22/100:.2f} (illustrative counts)")
nP, nD = 6, 4
n = nP + nD
phiB = np.full((n, n), 0.05); phiB[:nP, :nP] = 0.3; phiB[nP:, nP:] = 0.9; np.fill_diagonal(phiB, 0)
gB, gB_no = gap(rownorm(phiB)), gap(rownorm(phiB[:nP, :nP]))
blocks, bsz = 3, 3
m = blocks * bsz + nD
phiA = np.full((m, m), 0.001)
for kk in range(blocks):
    phiA[kk * bsz:(kk + 1) * bsz, kk * bsz:(kk + 1) * bsz] = 0.5
phiA[blocks * bsz:, :] = 0.3; phiA[:, blocks * bsz:] = 0.3; phiA[blocks * bsz:, blocks * bsz:] = 0.9
np.fill_diagonal(phiA, 0)
gA, gA_no = gap(rownorm(phiA)), gap(rownorm(phiA[:blocks * bsz, :blocks * bsz]))
print(f"  bridge core: Gamma without D {gA_no:.3f} -> with D {gA:.3f}; module core: without {gB_no:.3f} -> with {gB:.3f}")
assert round(gA, 3) == 0.548 and round(gA_no, 3) == 0.009 and round(gB, 3) == 0.218 and round(gB_no, 3) == 0.800
es = pd.read_csv(DATA / "eci_selected.csv")
D = es[es.classification == "hs92"].set_index(["country_iso3", "year"])
divB = D.loc[("BOL", 2019), "diversity"]
print(f"  Bolivia diversity (products with RCA >= 1, HS92 4-digit), 2019: {divB:.0f}; adding 20 products: +{20/divB*100:.1f}%")

# ---------------------------------------------------------------------------
print("\n=== 25.5 ECI (SITC) for the divergence countries, earliest available years ===")
for c in ("BRA", "MEX", "KOR"):
    print(f"  {c}: ECI_SITC 1995 {E.loc[(c, 1995), 'eci_sitc']:+.3f}, 2000 {E.loc[(c, 2000), 'eci_sitc']:+.3f}, "
          f"2019 {E.loc[(c, 2019), 'eci_sitc']:+.3f}")
assert "TWN" not in set(eci.country)

# ---------------------------------------------------------------------------
print("\n=== 25.6 Net change in diversity (products with RCA), HS92, 1995-2019 ===")
for c in ("BOL", "KOR"):
    d0, d1 = D.loc[(c, 1995), "diversity"], D.loc[(c, 2019), "diversity"]
    e0, e1 = E.loc[(c, 1995), "eci_hs92"], E.loc[(c, 2019), "eci_hs92"]
    print(f"  {c}: diversity {d0:.0f} -> {d1:.0f} (net {d1-d0:+.0f}, {(d1-d0)/24:+.2f} per year); ECI {e0:+.2f} -> {e1:+.2f}")
assert D.loc[("BOL", 2019), "diversity"] < D.loc[("BOL", 1995), "diversity"]
assert D.loc[("KOR", 2019), "diversity"] < D.loc[("KOR", 1995), "diversity"]
bol_div = D.loc["BOL"].diversity
print(f"  Bolivia diversity range 1995-2024: {bol_div.min():.0f} ({int(bol_div.idxmin())}) to {bol_div.max():.0f} ({int(bol_div.idxmax())})")
vel = pd.read_csv(DATA / "gamma_velocity.csv")
for _, rr in vel.iterrows():
    if rr.country in ("BOL", "KOR"):
        print(f"  v_Gamma {rr.country} [{rr.vintage}]: {rr.v_endpoint_per_decade:+.3f} per decade")

# ---------------------------------------------------------------------------
print("\n=== 25.7 PCI along the lithium chain (Growth Lab HS92) ===")
pci = pd.read_csv(DATA / "pci_selected_products.csv", dtype={"product_hs92_code": str})
Pq = pci.set_index(["product_hs92_code", "year"]).pci
for code in ("2711", "2836", "2825", "8507"):
    print(f"  HS {code}: PCI 2019 {Pq[(code, 2019)]:+.3f}, 2023 {Pq[(code, 2023)]:+.3f}")
d1 = Pq[("2836", 2019)] - Pq[("2711", 2019)]
d2 = Pq[("8507", 2019)] - Pq[("2836", 2019)]
print(f"  PCI step gas -> carbonate {d1:+.3f}; carbonate -> accumulators {d2:+.3f}; "
      f"accumulators vs illustrative mean PCI of D (1.8): {Pq[('8507', 2019)] - 1.8:+.3f}")
assert round(Pq[("8507", 2019)], 3) == 0.549 and round(Pq[("2836", 2019)], 3) == -0.573
print("all chapter-25 solution checks passed")
