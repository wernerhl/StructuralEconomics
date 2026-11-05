"""Chapter 25: dense subgraph vs bridge, Cheeger bound, truncation numbers,
and the sourced 12-sector test (values read from data/processed/*.csv).
Run from the repo root: python3 code/ch25.py  (numpy + csv)."""
import csv
import itertools
import pathlib
import numpy as np

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "processed"


def rownorm(K):
    return K / K.sum(1, keepdims=True)


def gap(W):
    ev = np.linalg.eigvals(W)
    ev = ev[np.argsort(-np.abs(ev))]
    return 1 - abs(ev[1])


def goldberg_obj(phi, S):
    S = list(S)
    sub = phi[np.ix_(S, S)]
    return np.triu(sub, 1).sum() / len(S)


def charikar_peel(phi):
    """Greedy peeling (Charikar 2000): 1/2-approximation to max e(D)/|D|."""
    alive = list(range(len(phi)))
    best, bestS = -1, None
    while alive:
        val = goldberg_obj(phi, alive)
        if val > best:
            best, bestS = val, list(alive)
        deg = phi[np.ix_(alive, alive)].sum(1)
        alive.pop(int(np.argmin(deg)))
    return bestS, best


def conductance(phi, S):
    S = np.array(S)
    T = np.setdiff1d(np.arange(len(phi)), S)
    d = phi.sum(1)
    return phi[np.ix_(S, T)].sum() / min(d[S].sum(), d[T].sum())


print("== Case B: D is a dense MODULE (6 periphery products well mixed, 4 core) ==")
nP, nD = 6, 4
n = nP + nD
phiB = np.full((n, n), 0.05)
phiB[:nP, :nP] = 0.3
phiB[nP:, nP:] = 0.9
np.fill_diagonal(phiB, 0)
dens = lambda M, S: M[np.ix_(S, S)].sum() / (len(S) * (len(S) - 1))
D = list(range(nP, n))
print(f"density of D = {dens(phiB, D):.3f}, of all = {dens(phiB, list(range(n))):.3f}, ratio = {dens(phiB, D)/dens(phiB, list(range(n))):.2f}")
S, val = charikar_peel(phiB)
# exact Goldberg optimum by brute force (n = 10)
best = max((goldberg_obj(phiB, c), c) for k in range(1, n + 1) for c in itertools.combinations(range(n), k))
print(f"peeling finds {S} (obj {val:.3f}); brute-force optimum {list(best[1])} (obj {best[0]:.3f})")
gB, gB_noD = gap(rownorm(phiB)), gap(rownorm(phiB[:nP, :nP]))
print(f"Gamma with D = {gB:.3f}; without D = {gB_noD:.3f}  (removing the dense module RAISES Gamma)")
assert gB_noD > gB

print("\n== Case A: D is a BRIDGE between three periphery blocks ==")
blocks, b = 3, 3
m = blocks * b + nD
phiA = np.full((m, m), 0.001)
for k in range(blocks):
    phiA[k * b:(k + 1) * b, k * b:(k + 1) * b] = 0.5
phiA[blocks * b:, :] = 0.3
phiA[:, blocks * b:] = 0.3
phiA[blocks * b:, blocks * b:] = 0.9
np.fill_diagonal(phiA, 0)
gA, gA_noD = gap(rownorm(phiA)), gap(rownorm(phiA[:blocks * b, :blocks * b]))
print(f"Gamma with D = {gA:.3f}; without D = {gA_noD:.3f}  (removing the bridge LOWERS Gamma)")
assert gA_noD < gA
# Cheeger upper bound after removal: Gamma <= 2 * conductance of block-1 cut
rest = phiA[:blocks * b, :blocks * b]
phi_cut = conductance(rest, list(range(b)))
print(f"conductance of block-1 cut without D = {phi_cut:.4f}; bound 2*Phi = {2*phi_cut:.4f} >= Gamma {gA_noD:.4f}")
assert gA_noD <= 2 * phi_cut + 1e-12
# with D present the same block-1 cut has conductance
print(f"conductance of block-1 cut with D present = {conductance(phiA, list(range(b))):.3f}")

# random check of Gamma <= 2 Phi(S) for reversible W = D^-1 phi
rng = np.random.default_rng(25)
for _ in range(300):
    k = 6
    P = rng.uniform(0, 1, (k, k)); P = (P + P.T) / 2; np.fill_diagonal(P, rng.uniform(0, 1, k))
    W = rownorm(P)
    ev = np.sort(np.real(np.linalg.eigvals(W)))[::-1]
    for r_ in range(1, k):
        for Sset in itertools.combinations(range(k), r_):
            assert 1 - ev[1] <= 2 * conductance(P, list(Sset)) + 1e-9
            assert gap(W) <= 1 - ev[1] + 1e-9
print("checked Gamma <= 1 - lambda_2 <= 2*Phi(S) for 300 random reversible W, all cuts")

# Proof check for Proposition ch25:prop:bridge_bound: Rayleigh quotient of the test vector
for _ in range(300):
    k = int(rng.integers(3, 9))
    P = rng.uniform(0, 1, (k, k)); P = (P + P.T) / 2; np.fill_diagonal(P, 0.0)
    d = P.sum(1)
    Sset = rng.choice(k, size=int(rng.integers(1, k)), replace=False)
    T = np.setdiff1d(np.arange(k), Sset)
    vS, vT = d[Sset].sum(), d[T].sum()
    f = np.zeros(k); f[Sset] = 1 / vS; f[T] = -1 / vT
    assert abs(f @ d) < 1e-12
    num = 0.5 * sum(P[p_, q_] * (f[p_] - f[q_]) ** 2 for p_ in range(k) for q_ in range(k))
    den = f @ (d * f)
    kap = 1 / vS + 1 / vT
    cut = P[np.ix_(Sset, T)].sum()
    assert abs(num - cut * kap ** 2) < 1e-10 and abs(den - kap) < 1e-10
    lam2 = np.sort(np.real(np.linalg.eigvals(P / d[:, None])))[::-1][1]
    assert 1 - lam2 <= num / den + 1e-10 <= 2 * conductance(P, list(Sset)) + 2e-10
print("checked: test-vector Rayleigh quotient = cut*(1/vol S + 1/vol S^c) >= 1 - lambda_2, <= 2 Phi(S)")

print("\n== Degenerate definition: average including diagonal phi_pp = 1 ==")
phid = phiB.copy(); np.fill_diagonal(phid, 1.0)
single = phid[np.ix_([0], [0])].mean()
print(f"|D| = 1 gives density {single:.2f} = maximum: the old definition picks a single product")

print("\n== Truncation index (illustrative) ==")
for c, k in (("Bolivia", 2), ("South Korea", 82)):
    print(f"{c}: tau = 1 - {k}/100 = {1 - k/100:.2f}")

print("\n== Relaxation-time penalty (1-G)/G, stylized truncated economy G = 0.13 (illustrative) ==")
f = lambda g: (1 - g) / g
print(f"(1-0.13)/0.13 = {f(0.13):.3f}; (1-0.35)/0.35 = {f(0.35):.3f}; difference = {f(0.13)-f(0.35):.3f}")
for g0 in (0.25, 0.30, 0.35, 0.45, 0.52):
    print(f"counterfactual Gamma(0) = {g0:.2f}: share = {(f(0.13)-f(g0))/f(0.13):.3f}")

print("\n== Bolivia's incremental path timing ==")
lo, hi = 4 * 20, 4 * 30
print(f"4 steps x 20-30 years = {lo}-{hi} years; at 7.5x velocity: {lo/7.5:.1f}-{hi/7.5:.1f} years")
print("\n== Korea path (Ch23 linear): Gamma(1973) =", round(0.15 + 0.0075 * 13, 3), " Gamma(1993) =", round(0.15 + 0.0075 * 33, 3))


# ---------------------------------------------------------------------------
# Sourced values (data/processed), asserted against the text's rounding
# ---------------------------------------------------------------------------
def rows(name):
    with open(DATA / name, newline="") as fh:
        return list(csv.DictReader(fh))


def close(x, text, nd):
    """x rounds (half away from zero) to the value printed in the text with nd decimals."""
    assert abs(x - text) <= 0.5 * 10 ** -nd + 1e-9, (x, text)
    return x


print("\n== Sourced: 12-sector Gamma, 2019 (book12_isic, domestic, input shares) ==")
gam = {}
for r in rows("gamma_comparison.csv"):
    if r["year"] == "2019" and r["aggregation"] == "book12_isic" and r["normalization"] == "input_shares":
        gam["BOL" if r["country"].startswith("BOL") else r["country"]] = float(r["Gamma"])
text_gamma = {"BOL": 0.475, "KOR": 0.511, "DEU": 0.544, "USA": 0.463, "CHL": 0.557,
              "CHN": 0.582, "MEX": 0.649, "ARG": 0.522}
for c, v in text_gamma.items():
    close(gam[c], v, 3)
    print(f"  {c}: Gamma = {gam[c]:.3f}")
comp = [gam[c] for c in ("KOR", "DEU", "USA", "CHL", "CHN", "MEX")]
assert min(comp) < gam["BOL"] < max(comp)
print(f"  Bolivia inside comparator range [{min(comp):.3f}, {max(comp):.3f}]; above USA: {gam['BOL'] > gam['USA']}")

print("\n== Sourced: ECI (Growth Lab HS92) ==")
eci = {(r["country"], int(r["year"])): float(r["eci_hs92"]) for r in rows("eci_facts.csv")}
rank = {(r["country"], int(r["year"])): int(r["eci_rank_hs92"]) for r in rows("eci_facts.csv")}
for c, v in {"BOL": -1.30, "KOR": 1.81, "DEU": 1.54, "USA": 1.29, "CHL": -0.06, "CHN": 1.30, "MEX": 0.95}.items():
    close(eci[(c, 2019)], v, 2)
    print(f"  {c} 2019: ECI = {eci[(c, 2019)]:+.3f}")
assert rank[("BOL", 2019)] == 131 and rank[("KOR", 2019)] == 2
close(eci[("KOR", 1995)], 1.10, 2)
close(eci[("BOL", 2006)], -1.21, 2)
close(eci[("BOL", 2014)], -1.34, 2)
print(f"  ranks 2019: BOL {rank[('BOL', 2019)]}, KOR {rank[('KOR', 2019)]}; KOR 1995 {eci[('KOR', 1995)]:+.2f}; "
      f"BOL 2006 {eci[('BOL', 2006)]:+.2f}, 2014 {eci[('BOL', 2014)]:+.2f}")

print("\n== Sourced: Gamma over time ==")
bsum = rows("bolivia_io_summary.csv")
def bol(base, year, agg="book12", norm="input_shares"):
    return [float(r["Gamma"]) for r in bsum if r["base"] == str(base) and r["year"] == str(year)
            and r["aggregation"] == agg and r["flows"] == "domestic" and r["normalization"] == norm][0]
for base, y, v in ((1990, 1990, 0.567), (1990, 2006, 0.518), (1990, 2013, 0.484), (2017, 2019, 0.467)):
    close(bol(base, y), v, 3)
    print(f"  BOL {y} (base {base}): {bol(base, y):.3f}")
osum = rows("oecd_io_summary.csv")
def oecd(c, y, ed, agg="book12_isic"):
    return [float(r["Gamma"]) for r in osum if r["country"] == c and r["year"] == str(y) and r["edition"] == ed
            and r["aggregation"] == agg and r["flows"] == "domestic" and r["normalization"] == "input_shares"][0]
close(oecd("KOR", 1995, "ICIO2021"), 0.626, 3)
close(oecd("KOR", 2019, "ICIO2023"), 0.511, 3)
close(oecd("KOR", 2019, "ICIO2023", "detail45"), 0.455, 3)
close(bol(2017, 2019, "detail39"), 0.410, 3)
print(f"  KOR 1995 {oecd('KOR', 1995, 'ICIO2021'):.3f} -> 2019 {oecd('KOR', 2019, 'ICIO2023'):.3f}; "
      f"detail: BOL39 {bol(2017, 2019, 'detail39'):.3f}, KOR45 {oecd('KOR', 2019, 'ICIO2023', 'detail45'):.3f}")
allg = [float(r["Gamma"]) for r in bsum if r["base"] == "2017" and r["year"] == "2019" and r["flows"] == "domestic"
        and r["aggregation"] in ("book12", "detail39")]
close(min(allg), 0.33, 2); close(max(allg), 0.69, 2)
gapBK = gam["KOR"] - gam["BOL"]
print(f"  normalization range {min(allg):.3f}-{max(allg):.3f} (width {max(allg)-min(allg):.2f}) vs BOL-KOR gap {gapBK:.3f}")
assert max(allg) - min(allg) > 5 * gapBK

vel = {(r["country"], r["vintage"]): float(r["v_endpoint_per_decade"]) for r in rows("gamma_velocity.csv")}
close(vel[("BOL", "INE base 1990, 1990-2013")], -0.036, 3)
close(vel[("KOR", "ICIO2021")], -0.055, 3)
print(f"  v_Gamma: BOL 1990-2013 {vel[('BOL', 'INE base 1990, 1990-2013')]:+.3f}/decade; KOR 1995-2014 {vel[('KOR', 'ICIO2021')]:+.3f}/decade")

print("\n== Sourced: relaxation time, Bolivia vs Korea ==")
fB, fK = f(0.475), f(0.511)
share = (fB - fK) / fB
print(f"  1/0.475-1 = {fB:.2f}; 1/0.511-1 = {fK:.2f}; share with Korea as counterfactual = {share:.3f}")
assert round(fB, 2) == 1.11 and round(fK, 2) == 0.96 and round(100 * share) == 13

print("\n== Sourced: hydrocarbons enclave and multipliers (INE 2019, book12, domestic) ==")
pi_h = [float(r["pi"]) for r in rows("bolivia_io_pi.csv") if r["year"] == "2019" and r["aggregation"] == "book12"
        and r["flows"] == "domestic" and r["normalization"] == "input_shares" and r["sector"] == "Hydrocarbons"][0]
close(pi_h, 0.039, 3)
facts = {r["fact_id"]: float(r["value"]) for r in rows("bolivia_facts.csv")}
close(facts["BOL_hydrocarbons_export_share_of_output_2019_b2017"] / 100, 0.70, 2)
close(facts["BOL_hydrocarbons_dom_inputs_per_output_2019_b2017"], 0.225, 3)
link = {r["sector"]: r for r in rows("bolivia_io_linkages.csv") if r["year"] == "2019" and r["aggregation"] == "book12"
        and r["flows"] == "domestic"}
close(float(link["Manufacturing"]["intermediate_input_share"]), 0.558, 3)
close(float(link["Hydrocarbons"]["BL"]), 1.34, 2)
close(float(link["Manufacturing"]["BL"]), 1.86, 2)
print(f"  pi_hydro = {pi_h:.3f}; exported {facts['BOL_hydrocarbons_export_share_of_output_2019_b2017']:.1f}%; "
      f"dom. inputs/output {facts['BOL_hydrocarbons_dom_inputs_per_output_2019_b2017']:.3f} (manuf "
      f"{float(link['Manufacturing']['intermediate_input_share']):.3f}); BL hydro {float(link['Hydrocarbons']['BL']):.3f}, "
      f"manuf {float(link['Manufacturing']['BL']):.3f}")

print("\n== Sourced: GDP per capita PPP 2019 (WDI), GDP growth 2006-2014, PCI ==")
close(facts["BOL_gdppc_ppp_current_2019"], 11005, 0)
close(facts["KOR_gdppc_ppp_current_2019"], 46511, 0)
close(facts["BOL_gdp_growth_mean_2006_2014"], 5.06, 2)
pci = {(r["product_hs92_code"], int(r["year"])): float(r["pci"]) for r in rows("pci_selected_products.csv")}
close(pci[("8507", 2019)], 0.549, 3); close(pci[("8507", 2023)], 0.811, 3)
close(pci[("2836", 2019)], -0.573, 3); close(pci[("2711", 2019)], -2.226, 3)
print(f"  BOL ${facts['BOL_gdppc_ppp_current_2019']:,.0f}, KOR ${facts['KOR_gdppc_ppp_current_2019']:,.0f}; "
      f"growth 2006-14 {facts['BOL_gdp_growth_mean_2006_2014']:.2f}%; PCI 8507 {pci[('8507', 2019)]:.3f} "
      f"({pci[('8507', 2023)]:.3f} in 2023), 2836 {pci[('2836', 2019)]:.3f}, 2711 {pci[('2711', 2019)]:.3f}")
print("all sourced values match the text")
