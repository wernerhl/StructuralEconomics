"""Chapter 26: counterexamples to the old proof, viz-panel classification
(sourced 2019 Gamma and ECI), Korea check, the Bolivia remark, and plan arithmetic.
Sourced values are read from data/processed/*.csv and asserted against the text.
Run: python3 code/ch26.py (numpy + csv)."""
import csv
import pathlib
import numpy as np

DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "processed"


def rows(name):
    with open(DATA / name, newline="") as fh:
        return list(csv.DictReader(fh))


def close(x, text, nd):
    """x rounds (half away from zero) to the value printed in the text with nd decimals."""
    assert abs(x - text) <= 0.5 * 10 ** -nd + 1e-9, (x, text)
    return x


def stats(W):
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi /= pi.sum()
    N = len(W)
    var = np.mean((pi - 1 / N) ** 2)
    return 1 - abs(ev[o[1]]), pi, var, 1 / np.sum(pi ** 2)


print("== Gamma does not bound Var(pi) (Ch6 Prop gap-power-independent) ==")
N = 4
pi = np.array([0.97, 0.01, 0.01, 0.01])
W1 = np.outer(np.ones(N), pi)
G, p, v, ne = stats(W1)
print(f"rank-one 1pi': Gamma = {G:.3f}, Var(pi) = {v:.4f} (max (N-1)/N^2 = {(N-1)/N**2:.4f}), Neff = {ne:.2f}")
assert abs(G - 1) < 1e-9 and v <= (N - 1) / N ** 2
e = 0.01
B = np.full((2, 2), 0.5)
W2 = np.block([[(1 - e) * B, e * B], [e * B, (1 - e) * B]])
G, p, v, ne = stats(W2)
print(f"doubly stochastic weak blocks: Gamma = {G:.3f}, Var(pi) = {v:.2e}, Neff = {ne:.2f}")
assert abs(G - 0.02) < 1e-9 and v < 1e-20

print("\n== Viz panel (sourced, 2019): Gamma book12_isic (INE / OECD ICIO 2023), ECI HS92 (Growth Lab) ==")
gam = {}
for r in rows("gamma_comparison.csv"):
    if r["year"] == "2019" and r["aggregation"] == "book12_isic" and r["normalization"] == "input_shares":
        gam["BOL" if r["country"].startswith("BOL") else r["country"]] = float(r["Gamma"])
eci = {(r["country"], int(r["year"])): float(r["eci_hs92"]) for r in rows("eci_facts.csv")}
rank = {(r["country"], int(r["year"])): int(r["eci_rank_hs92"]) for r in rows("eci_facts.csv")}
# (Gamma, ECI) as printed in the text / plotted
text_pts = {"BOL": (0.475, -1.30), "BRA": (0.496, 0.48), "COL": (0.460, 0.14), "MEX": (0.649, 0.95),
            "KOR": (0.511, 1.81), "JPN": (0.510, 2.05), "DEU": (0.544, 1.54), "CZE": (0.486, 1.48),
            "CHL": (0.557, -0.06), "ARG": (0.522, 0.06), "MYS": (0.381, 1.03), "POL": (0.552, 1.06),
            "CHN": (0.582, 1.30)}
plotted = {"BOL": (4.50, -1.30), "BRA": (4.92, 0.48), "COL": (4.21, 0.14), "MEX": (7.98, 0.95),
           "KOR": (5.22, 1.81), "JPN": (5.20, 2.05), "DEU": (5.87, 1.54), "CZE": (4.72, 1.48),
           "CHL": (6.14, -0.06), "ARG": (5.43, 0.06), "MYS": (2.62, 1.03), "POL": (6.03, 1.06),
           "CHN": (6.63, 1.30)}
for k, (g, c) in text_pts.items():
    close(gam[k], g, 3); close(eci[(k, 2019)], c, 2)
    x = 20 * (gam[k] - 0.25)
    assert abs(x - plotted[k][0]) < 0.006 and plotted[k][1] == c, (k, x)
    print(f"  {k}: Gamma = {gam[k]:.3f}, ECI = {eci[(k, 2019)]:+.3f}, x = 20(Gamma-0.25) = {x:.2f}")
close(eci[("VEN", 2019)], -1.02, 2)
assert "VEN" not in gam
inside = sorted(k for k, (g, c) in text_pts.items() if g > 0.30 and c > 0)
print("inside the zone:", inside)
assert "CHL" not in inside and "CHN" in inside and "BOL" not in inside
assert all(g > 0.30 for g, _ in text_pts.values())  # every country right of Gamma*: panel separates by ECI only
gs = [g for g, _ in text_pts.values()]
close(min(gs), 0.38, 2); close(max(gs), 0.65, 2)
print(f"Gamma range of plotted countries: {min(gs):.3f}-{max(gs):.3f}")
for k, v in (("BRA", 0.496), ("MEX", 0.649), ("KOR", 0.511), ("USA", 0.463), ("DEU", 0.544)):
    close(gam[k], v, 3)

print("\n== Korea, sourced from 1995 (OECD ICIO 2021/2023, book12_isic) and ECI ==")
osum = rows("oecd_io_summary.csv")
def oecd(c, y, ed):
    return [float(r["Gamma"]) for r in osum if r["country"] == c and r["year"] == str(y) and r["edition"] == ed
            and r["aggregation"] == "book12_isic" and r["flows"] == "domestic" and r["normalization"] == "input_shares"][0]
close(oecd("KOR", 1995, "ICIO2021"), 0.626, 3); close(oecd("KOR", 2019, "ICIO2023"), 0.511, 3)
close(eci[("KOR", 1995)], 1.10, 2); close(eci[("KOR", 2019)], 1.81, 2)
print(f"  KOR Gamma 1995 {oecd('KOR', 1995, 'ICIO2021'):.3f} -> 2019 {oecd('KOR', 2019, 'ICIO2023'):.3f}; "
      f"ECI {eci[('KOR', 1995)]:+.2f} -> {eci[('KOR', 2019)]:+.2f}")

print("\n== Korea on Chapter 23's illustrative path, Gamma(t) = 0.15 + 0.0075 (t - 1960) ==")
for t in (1965, 1970, 1975, 1980, 1985):
    g = 0.15 + 0.0075 * (t - 1960)
    print(f"  {t}: Gamma = {g:.3f}  {'< Gamma* = 0.30' if g < 0.30 else '>= Gamma*'}")

print("\n== Remark: Confronting Bolivia (sourced) ==")
facts = {r["fact_id"]: float(r["value"]) for r in rows("bolivia_facts.csv")}
F = lambda k, v, nd: close(facts[k], v, nd)
F("BOL_INE_gini_2005", 0.600, 3); F("BOL_INE_gini_2019", 0.426, 3)
F("BOL_INE_gini_2024", 0.418, 3); F("BOL_INE_gini_2025", 0.426, 3)
F("BOL_INE_extreme_poverty_oldline_2005", 38.2, 1); F("BOL_INE_extreme_poverty_oldline_2018", 15.2, 1)
F("BOL_INE_extreme_poverty_EPF_2019", 13.2, 1); F("BOL_INE_extreme_poverty_EPF_2024", 13.3, 1)
F("BOL_INE_extreme_poverty_EPF_2025", 16.7, 1)
F("BOL_PIP_gini_2005", 0.585, 3); F("BOL_PIP_gini_2019", 0.416, 3)
F("BOL_gdp_growth_mean_2006_2019", 4.52, 2); F("BOL_gdp_growth_mean_2006_2024", 3.46, 2)
F("BOL_fuel_exports_share_2006", 48.8, 1); F("BOL_fuel_exports_share_2013", 54.5, 1)
F("BOL_fuel_exports_share_2019", 32.1, 1); F("BOL_natgas_rents_pct_gdp_2013", 5.56, 2)
F("BOL_hydrocarbons_export_share_of_output_2019_b2017", 70, 0)
F("BOL_hydrocarbons_dom_inputs_per_output_2019_b2017", 0.225, 3)
F("BOL_hydrocarbons_comp_share_of_va_2019_b2017", 6, 0)
print(f"  Gini {facts['BOL_INE_gini_2005']:.3f} -> {facts['BOL_INE_gini_2019']:.3f}; extreme poverty "
      f"{facts['BOL_INE_extreme_poverty_oldline_2005']}% -> {facts['BOL_INE_extreme_poverty_oldline_2018']}% (2018, old lines), "
      f"{facts['BOL_INE_extreme_poverty_EPF_2019']}% (2019, new lines); growth 2006-19 {facts['BOL_gdp_growth_mean_2006_2019']}%")
bsum = rows("bolivia_io_summary.csv")
def bol(base, year, agg="book12"):
    return [float(r["Gamma"]) for r in bsum if r["base"] == str(base) and r["year"] == str(year)
            and r["aggregation"] == agg and r["flows"] == "domestic" and r["normalization"] == "input_shares"][0]
for base, y, v in ((1990, 2005, 0.502), (1990, 2006, 0.518), (1990, 2013, 0.484), (2017, 2019, 0.467)):
    close(bol(base, y), v, 3)
allg = [float(r["Gamma"]) for r in bsum if r["base"] == "2017" and r["year"] == "2019" and r["flows"] == "domestic"
        and r["aggregation"] in ("book12", "detail39")]
close(min(allg), 0.33, 2); close(max(allg), 0.69, 2)
for y, v in ((2005, -1.10), (2006, -1.21), (2018, -1.46), (2019, -1.30)):
    close(eci[("BOL", y)], v, 2)
assert rank[("BOL", 2005)] == 120 and rank[("BOL", 2019)] == 131
d_bol = eci[("BOL", 2019)] - eci[("BOL", 1995)]; d_kor = eci[("KOR", 2019)] - eci[("KOR", 1995)]
close(d_bol, -0.44, 2); close(d_kor, 0.70, 2)
print(f"  Gamma 2005 {bol(1990, 2005):.3f}, 2006 {bol(1990, 2006):.3f}, 2013 {bol(1990, 2013):.3f}, 2019 {bol(2017, 2019):.3f}; "
      f"ECI 2005 {eci[('BOL', 2005)]:+.2f} -> 2019 {eci[('BOL', 2019)]:+.3f}, rank 120 -> 131")
print(f"  Delta ECI 1995-2019: BOL {d_bol:+.2f}, KOR {d_kor:+.2f}; needed for ECI = 0 from 2019: {-eci[('BOL', 2019)]:+.3f} (text: +1.30)")
pi_h = [float(r["pi"]) for r in rows("bolivia_io_pi.csv") if r["year"] == "2019" and r["aggregation"] == "book12"
        and r["flows"] == "domestic" and r["normalization"] == "input_shares" and r["sector"] == "Hydrocarbons"][0]
close(pi_h, 0.039, 3)
print(f"  hydrocarbons pi = {pi_h:.3f}")

print("\n== Plan arithmetic ==")
costs = {"broadband": 1.5, "polytechnics": 2, "corridors": 5, "lithium": 3, "agri": 1, "industrial": 4}
tot = sum(costs.values())
F("BOL_gdp_usd_2019", 49.1, 1); F("BOL_gdp_usd_2024", 54.9, 1)
sh19 = 17 / 30 / facts["BOL_gdp_usd_2019"] * 100
sh24 = 17 / 30 / facts["BOL_gdp_usd_2024"] * 100
print(f"costed items = ${tot}B (rounded to $17B); per year over 30 = ${17/30*1000:.0f}M; "
      f"share of 2019 GDP ${facts['BOL_gdp_usd_2019']}B = {sh19:.2f}%; of 2024 GDP ${facts['BOL_gdp_usd_2024']}B = {sh24:.2f}%")
close(sh19, 1.2, 1); close(sh24, 1.0, 1)
print("Neff after Phase 1: 3.4 + 2.0 (broadband) =", 3.4 + 2.0, "-> text ~5.5 with community journalism (illustrative)")
print("all sourced values match the text")
