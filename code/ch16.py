"""Chapter 16: reproduces the numbers in the chapter's schematic panels,
checks the matrix facts used in the text, and checks the sourced values
(data/processed/*.csv) quoted in the text against the CSVs.
Run from the repo root:  python3 code/ch16.py   (numpy + csv)."""
import csv
import os
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
DATA = os.path.join(ROOT, "data", "processed")


def read_csv(name):
    with open(os.path.join(DATA, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def check(label, value, text, nd):
    """assert that the CSV value rounds to the number printed in the text"""
    print(f"  {label}: {value:.6g} -> text {text}")
    assert abs(value - text) <= 0.5 * 10 ** -nd + 1e-9, (label, value, text)   # half-up tolerance

np.set_printoptions(precision=3, suppress=True)


def gap(W):
    return 1 - np.sort(np.abs(np.linalg.eigvals(W)))[-2]


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


print("== Schematic 10x10 global W^P (core 1-5, semi-periphery 6-7, periphery 8-10) ==")
N = 10
H = np.zeros((N, N))
for i in range(1, N + 1):
    for j in range(1, N + 1):
        if i <= 5 and j <= 5: v = 65
        elif i <= 5 and 5 < j <= 7: v = 25
        elif i <= 5 and j > 7: v = 5
        elif 5 < i <= 7 and j <= 5: v = 45
        elif 5 < i <= 7 and 5 < j <= 7: v = 35
        elif 5 < i <= 7 and j > 7: v = 8
        elif i > 7 and j <= 5: v = 50
        elif i > 7 and 5 < j <= 7: v = 10
        elif i > 7 and j > 7 and i == j: v = 30
        else: v = 3
        H[i - 1, j - 1] = v
W = H / H.sum(1)[:, None]                  # W = D^-1 K
p = perron(W)
print(f"Gamma = {gap(W):.2f}; pi = {p}; Neff = {1/np.sum(p**2):.2f} (N = {N})")
print("column sums:", W.sum(0), " mean =", W.sum(0).mean())
assert abs(W.sum(0).mean() - 1) < 1e-12   # columns of a row-stochastic matrix average 1
assert gap(W) > 0.85

print("\n== Frictionless gravity: W = 1 s' ==")
s = np.array([0.18, 0.15, 0.08, 0.06, 0.05, 0.04, 0.03, 0.0005, 0.0004, 0.0003])   # toy shares
s_full = np.append(s, (1 - s.sum()) / 190 * np.ones(190))   # 200 countries
Wg = np.outer(np.ones(200), s_full)
print(f"Gamma(1 s') = {gap(Wg):.3f}, pi = s (first 3: {perron(Wg)[:3]}), Neff = {1/np.sum(s_full**2):.1f}")
assert abs(gap(Wg) - 1) < 1e-9 and np.allclose(perron(Wg), s_full)

# ---------------------------------------------------------------------------
print("\n== Sourced: 2019 trade W^P (Growth Lab HS92, imports only, 230 economies) ==")
summ = {r["year"]: r for r in read_csv("global_wp_summary.csv")}["2019"]
pi19 = {r["country"]: r for r in read_csv("global_wp_pi.csv") if r["year"] == "2019"}
N19 = int(summ["N"])
assert N19 == 230 and len(pi19) == 230
f = lambda k: float(summ[k])
check("Gamma_global", f("Gamma"), 0.452, 3)
check("pi_CHN", f("pi_CHN"), 0.130, 3)
check("pi_USA", f("pi_USA"), 0.089, 3)
check("pi_DEU", f("pi_DEU"), 0.082, 3)
check("pi_BOL", f("pi_BOL"), 0.00053, 5)
check("top-10 share (%)", 100 * f("top10_share"), 51.9, 1)
assert summ["top10"] == "CHN;USA;DEU;JPN;FRA;KOR;NLD;ITA;GBR;RUS"
check("Neff", f("Neff"), 22.6, 1)
check("power Gini", f("power_gini"), 0.85, 2)
cs = lambda c: float(pi19[c]["colsum_W"])
check("column sum USA", cs("USA"), 24.2, 1)
check("column sum CHN", cs("CHN"), 30.7, 1)
check("column sum DEU", cs("DEU"), 11.1, 1)
check("column sum BOL", cs("BOL"), 0.10, 2)
assert abs(np.mean([cs(c) for c in pi19]) - 1) < 1e-9          # column sums average 1
ws = lambda c: float(pi19[c]["world_export_share"])
check("export share CHN", ws("CHN"), 0.142, 3)
check("export share USA", ws("USA"), 0.091, 3)
check("export share BOL", ws("BOL"), 0.00046, 5)
# viz panel bars (right panel), rounded as printed
bars = {"CHN": (0.130, 3), "USA": (0.089, 3), "DEU": (0.082, 3), "JPN": (0.043, 3), "KOR": (0.032, 3),
        "GBR": (0.026, 3), "BRA": (0.0136, 4), "BOL": (0.00053, 5), "LAO": (0.00046, 5), "ETH": (0.00014, 5)}
for c, (v, nd) in bars.items():
    check(f"bar pi_{c}", float(pi19[c]["pi"]), v, nd)
check("1/N", 1 / N19, 0.00435, 5)
print(f"  pi*N: CHN {f('pi_CHN')*N19:.1f} (text 30), USA {f('pi_USA')*N19:.1f} (text 20)")
assert round(f("pi_CHN") * N19) == 30 and round(f("pi_USA") * N19) == 20
# pi_BOL over the boom (Example ch16:ex:bot_terms)
pib = {r["year"]: float(r["pi"]) for r in read_csv("global_wp_pi.csv") if r["country"] == "BOL"}
for y, v in [("2000", 0.00017), ("2006", 0.00047), ("2014", 0.00071), ("2018", 0.00051), ("2019", 0.00053)]:
    check(f"pi_BOL {y}", pib[y], v, 5)

print("\n== Sourced: ECI (Growth Lab HS92) and PCI ==")
eci = {(r["country"], r["year"]): r for r in read_csv("eci_facts.csv")}
for c, y, v, rank in [("BOL", "2006", -1.21, None), ("BOL", "2018", -1.46, None), ("BOL", "2019", -1.30, 131),
                      ("KOR", "2019", 1.81, 2), ("JPN", "2019", 2.05, 1)]:
    check(f"ECI {c} {y}", float(eci[(c, y)]["eci_hs92"]), v, 2)
    if rank is not None:
        assert int(eci[(c, y)]["eci_rank_hs92"]) == rank
pci = {(r["product_hs92_code"], r["year"]): float(r["pci"]) for r in read_csv("pci_selected_products.csv")}
check("PCI 8507 (2019)", pci[("8507", "2019")], 0.55, 2)
check("PCI 2836 (2019)", pci[("2836", "2019")], -0.57, 2)
check("PCI 2711 (2019)", pci[("2711", "2019")], -2.23, 2)

print("\n== Sourced: domestic W^P, common basis (book12_isic, domestic, input shares, 2019) ==")
bsum = read_csv("bolivia_io_summary.csv")
def bol(year, agg, key, base=None):
    r = [x for x in bsum if x["year"] == year and x["aggregation"] == agg and x["flows"] == "domestic"
         and x["normalization"] == "input_shares" and (base is None or x["base"] == base)]
    assert len(r) == 1
    return float(r[0][key])
check("Gamma BOL (book12_isic)", bol("2019", "book12_isic", "Gamma"), 0.475, 3)
check("power Gini BOL", bol("2019", "book12_isic", "power_gini"), 0.435, 3)
check("Neff BOL", bol("2019", "book12_isic", "Neff"), 6.67, 2)
assert int(float([x for x in bsum if x["year"] == "2019" and x["aggregation"] == "book12_isic"][0]["N"])) == 11
check("Gamma BOL 2006 (book12, base 1990)", bol("2006", "book12", "Gamma", "1990"), 0.518, 3)
check("Gamma BOL 2013 (book12, base 1990)", bol("2013", "book12", "Gamma", "1990"), 0.484, 3)
g19 = [float(x["Gamma"]) for x in bsum if x["year"] == "2019" and x["base"] == "2017" and x["flows"] == "domestic"
       and x["aggregation"] in ("book12", "detail39")]          # 2 aggregations x 4 normalizations
print(f"  Gamma BOL 2019 range over normalizations/aggregations: {min(g19):.2f}-{max(g19):.2f}")
assert round(min(g19), 2) == 0.33 and round(max(g19), 2) == 0.69
osum = read_csv("oecd_io_summary.csv")
def icio(c, key):
    r = [x for x in osum if x["edition"] == "ICIO2023" and x["country"] == c and x["year"] == "2019"
         and x["aggregation"] == "book12_isic" and x["flows"] == "domestic" and x["normalization"] == "input_shares"]
    assert len(r) == 1, (c, len(r))
    return float(r[0][key])
for c, v in [("KOR", 0.511), ("USA", 0.463), ("DEU", 0.544), ("CHL", 0.557), ("CHN", 0.582), ("MEX", 0.649)]:
    check(f"Gamma {c}", icio(c, "Gamma"), v, 3)
check("power Gini KOR", icio("KOR", "power_gini"), 0.658, 3)
check("power Gini DEU", icio("DEU", "power_gini"), 0.611, 3)
check("Neff KOR", icio("KOR", "Neff"), 3.62, 2)
check("Neff DEU", icio("DEU", "Neff"), 3.89, 2)

print("\n== Sourced: Bolivia domestic linkages and hydrocarbons (INE SUT 2019, book12, domestic) ==")
lk = [r for r in read_csv("bolivia_io_linkages.csv") if r["year"] == "2019" and r["aggregation"] == "book12"
      and r["flows"] == "domestic"]
assert len(lk) == 12
BLb = {r["sector"]: float(r["BL"]) for r in lk}
FLb = {r["sector"]: float(r["FL"]) for r in lk}
check("mean BL", np.mean(list(BLb.values())), 1.445, 3)
assert abs(np.mean(list(BLb.values())) - np.mean(list(FLb.values()))) < 1e-9
check("BL hydrocarbons", BLb["Hydrocarbons"], 1.34, 2)
check("FL hydrocarbons", FLb["Hydrocarbons"], 1.20, 2)
sf = {r["book12"]: r for r in read_csv("bolivia_io_sector_facts.csv") if r["base"] == "2017" and r["year"] == "2019"}
check("hydro dom. inputs / output", float(sf["Hydrocarbons"]["dom_int_inputs_per_output"]), 0.225, 3)
check("manuf dom. inputs / output", float(sf["Manufacturing"]["dom_int_inputs_per_output"]), 0.558, 3)
check("hydro export share (%)", 100 * float(sf["Hydrocarbons"]["export_share_of_output"]), 70, 0)

print("\n== Sourced: distribution, growth, fuel exports (INE, WDI via bolivia_facts.csv) ==")
facts = {r["fact_id"]: float(r["value"]) for r in read_csv("bolivia_facts.csv")}
for k, v, nd in [("BOL_INE_gini_2005", 0.600, 3), ("BOL_INE_gini_2006", 0.592, 3), ("BOL_INE_gini_2019", 0.426, 3),
                 ("BOL_INE_extreme_poverty_oldline_2005", 38.2, 1), ("BOL_INE_extreme_poverty_oldline_2018", 15.2, 1),
                 ("BOL_gdp_growth_mean_2006_2019", 4.52, 2), ("BOL_fuel_exports_share_2006", 48.8, 1),
                 ("BOL_fuel_exports_share_2013", 54.5, 1), ("BOL_fuel_exports_share_2019", 32.1, 1)]:
    check(k, facts[k], v, nd)

print("\n== Country-similarity matrix M~ of the ECI ==")
rng = np.random.default_rng(5)
M = (rng.uniform(size=(30, 80)) < np.linspace(0.1, 0.7, 30)[:, None]).astype(float)
M = M[M.sum(1) > 0][:, M.sum(0) > 0]
Mt = (M / M.sum(1)[:, None]) @ (M / M.sum(0)).T
assert np.allclose(Mt.sum(1), 1)
ev = np.linalg.eigvals(Mt)
pM = perron(Mt)
assert np.allclose(pM, M.sum(1) / M.sum())       # left Perron vector = diversity shares
assert np.abs(ev.imag).max() < 1e-10 and ev.real.min() > -1e-10 and ev.real.max() < 1 + 1e-10
Dc, Dp = M.sum(1), M.sum(0)
S = np.diag(Dc ** -0.5) @ M @ np.diag(Dp ** -0.5)
sv = np.linalg.svd(S, compute_uv=False)
assert np.allclose(np.sort(sv ** 2)[::-1][:len(ev)], np.sort(ev.real)[::-1])
print("checked: M~ row-stochastic, eigenvalues real in [0,1] (= squared singular values of D_c^-1/2 M D_p^-1/2),")
print("         left Perron vector = diversity / total diversity")

print("\n== Linkage table (illustrative): internal consistency ==")
BL = {"United States": 2.45, "China": 2.68, "Germany": 2.51, "Japan": 2.33, "Bolivia": 1.28, "Ethiopia": 1.15}
FL = {"United States": 1.98, "China": 2.31, "Germany": 2.12, "Japan": 1.85, "Bolivia": 1.12, "Ethiopia": 1.08}
assert all(v >= 1 for v in BL.values()) and all(v >= 1 for v in FL.values())
print("all BL, FL >= 1 (Leontief-inverse sums include the identity)")
# mean column sum = mean row sum for any square matrix
A = rng.uniform(0, 0.2, (8, 8)); L = np.linalg.inv(np.eye(8) - A.T)
assert abs(L.sum(0).mean() - L.sum(1).mean()) < 1e-12 and L.sum(0).min() >= 1 and L.sum(1).min() >= 1
print("checked: mean backward linkage = mean forward linkage; all >= 1")
g = 0.5; Wr = A / A.sum(1)[:, None]
print("row-normalized W: row sums of (I - g W)^-1 =", np.linalg.inv(np.eye(8) - g * Wr).sum(1)[:3], "... all = 1/(1-g)")

print("\n== Schematic ECI-Gini scatter: least-squares line ==")
pts = np.array([(5.8, 0.9), (5.2, 1.0), (4.5, 1.3), (4.8, 1.9), (4.0, 1.1), (3.8, 1.6), (3.5, 2.0), (3.0, 2.6),
                (2.5, 2.3), (2.8, 3.4), (2.0, 2.9), (1.5, 2.6), (1.2, 1.9), (1.0, 2.2), (0.8, 3.0), (2.2, 1.7)])
b, a = np.polyfit(pts[:, 0], pts[:, 1], 1)
print(f"fitted line y = {a:.3f} + ({b:.3f}) x; endpoints (0.5, {a+0.5*b:.2f}), (6.3, {a+6.3*b:.2f})")
