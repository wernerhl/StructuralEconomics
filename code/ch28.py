"""Chapter 28: checks behind the Framing Principle, the Bolivia table and the
Bolivia rows of 'The Book in One Page'. Sourced values are read from
data/processed/*.csv and asserted against the text's rounding.
Run from the repo root: python3 code/ch28.py  (numpy + csv)."""
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
    return 1 - abs(ev[o[1]]), pi, 1 / np.sum(pi ** 2)


N = 5
s = np.array([0.4, 0.25, 0.15, 0.12, 0.08])          # wealth shares
G, pi, ne = stats(np.outer(np.ones(N), s))
print(f"W_mkt = 1 s': Gamma = {G:.3f}, pi = s: {np.allclose(pi, s)}, Neff = {ne:.2f}")
G, pi, ne = stats(np.outer(np.ones(N), np.eye(N)[0]))
print(f"dictatorship 1 e1': Gamma = {G:.3f}, Neff = {ne:.2f}  -> Gamma = 1 does not identify the benchmark")
# Gamma and pi independent: W = g 1 pi' + (1-g) I
p = np.array([0.6, 0.1, 0.1, 0.1, 0.1])
for g in (0.1, 0.5, 1.0):
    G, pi, ne = stats(g * np.outer(np.ones(N), p) + (1 - g) * np.eye(N))
    print(f"gamma = {g}: Gamma = {G:.3f}, pi = {np.round(pi, 3)}")
    assert abs(G - g) < 1e-9 and np.allclose(pi, p)
# doubly stochastic but not W_CE: uniform pi, Delta != 0
C = np.array([[0.6, 0.3, 0.1], [0.1, 0.6, 0.3], [0.3, 0.1, 0.6]])
G, pi, ne = stats(C)
print(f"circulant: Gamma = {G:.3f}, pi uniform = {np.allclose(pi, 1/3)}, ||Delta||_F = {np.linalg.norm(C - 1/3):.3f}")

print("\n== Bolivia table and 'Book in One Page' rows 17, 23, 25, 26 (sourced) ==")
bsum = rows("bolivia_io_summary.csv")
def bol(base, year, agg="book12", col="Gamma"):
    return [float(r[col]) for r in bsum if r["base"] == str(base) and r["year"] == str(year)
            and r["aggregation"] == agg and r["flows"] == "domestic" and r["normalization"] == "input_shares"][0]
close(bol(2017, 2019), 0.467, 3)
close(bol(2017, 2019, "book12_isic"), 0.475, 3)
close(bol(2017, 2019, col="power_gini"), 0.445, 3)
close(bol(1990, 1990), 0.567, 3); close(bol(1990, 2013), 0.484, 3)
allg = [float(r["Gamma"]) for r in bsum if r["base"] == "2017" and r["year"] == "2019" and r["flows"] == "domestic"
        and r["aggregation"] in ("book12", "detail39")]
close(min(allg), 0.33, 2); close(max(allg), 0.69, 2)
print(f"  Gamma 2019 {bol(2017, 2019):.3f} (like-for-like {bol(2017, 2019, 'book12_isic'):.3f}; range {min(allg):.3f}-{max(allg):.3f}); "
      f"G(pi) {bol(2017, 2019, col='power_gini'):.3f}; Gamma 1990 {bol(1990, 1990):.3f} -> 2013 {bol(1990, 2013):.3f}")

comp = {}
for r in rows("oecd_io_summary.csv"):
    if r["aggregation"] == "book12_isic" and r["flows"] == "domestic" and r["normalization"] == "input_shares":
        comp[(r["country"], r["edition"], r["year"])] = (float(r["Gamma"]), float(r["power_gini"]))
six = [comp[(c, "ICIO2023", "2019")][0] for c in ("KOR", "DEU", "USA", "CHL", "CHN", "MEX")]
close(min(six), 0.463, 3); close(max(six), 0.649, 3)
close(comp[("KOR", "ICIO2023", "2019")][1], 0.658, 3); close(comp[("DEU", "ICIO2023", "2019")][1], 0.611, 3)
close(comp[("KOR", "ICIO2023", "2019")][0], 0.511, 3); close(comp[("KOR", "ICIO2021", "1995")][0], 0.626, 3)
assert bol(2017, 2019, col="power_gini") < comp[("KOR", "ICIO2023", "2019")][1]
four = [bol(2017, 2019, "book12_isic")] + [comp[(c, "ICIO2023", "2019")][0] for c in ("KOR", "DEU", "USA")]
assert 0.46 <= min(four) and max(four) <= 0.55   # limit 7: "between 0.46 and 0.55"
print(f"  comparators {min(six):.3f}-{max(six):.3f}; G(pi) KOR {comp[('KOR', 'ICIO2023', '2019')][1]:.3f}, "
      f"DEU {comp[('DEU', 'ICIO2023', '2019')][1]:.3f}; KOR Gamma 1995 {comp[('KOR', 'ICIO2021', '1995')][0]:.3f} -> 2019 "
      f"{comp[('KOR', 'ICIO2023', '2019')][0]:.3f}; BOL/KOR/DEU/USA in [{min(four):.3f}, {max(four):.3f}]")

pi_h = [float(r["pi"]) for r in rows("bolivia_io_pi.csv") if r["year"] == "2019" and r["aggregation"] == "book12"
        and r["flows"] == "domestic" and r["normalization"] == "input_shares" and r["sector"] == "Hydrocarbons"][0]
close(pi_h, 0.039, 3)
facts = {r["fact_id"]: float(r["value"]) for r in rows("bolivia_facts.csv")}
F = lambda k, v, nd: close(facts[k], v, nd)
F("BOL_hydrocarbons_export_share_of_output_2019_b2017", 70, 0)
F("BOL_hydrocarbons_dom_inputs_per_output_2019_b2017", 0.225, 3)
F("BOL_hydrocarbons_comp_share_of_va_2019_b2017", 6, 0)
F("BOL_INE_gini_2005", 0.600, 3); F("BOL_INE_gini_2019", 0.426, 3)
F("BOL_INE_extreme_poverty_oldline_2005", 38.2, 1); F("BOL_INE_extreme_poverty_oldline_2018", 15.2, 1)
F("BOL_INE_extreme_poverty_EPF_2019", 13.2, 1)
F("BOL_gdp_usd_2019", 49.1, 1)
close(17 / 30 / facts["BOL_gdp_usd_2019"] * 100, 1.2, 1)
print(f"  pi_hydro {pi_h:.3f}; exported {facts['BOL_hydrocarbons_export_share_of_output_2019_b2017']}%; "
      f"Gini {facts['BOL_INE_gini_2005']} -> {facts['BOL_INE_gini_2019']}; plan cost {17/30/facts['BOL_gdp_usd_2019']*100:.2f}% of 2019 GDP")

eci = {(r["country"], int(r["year"])): (float(r["eci_hs92"]), int(r["eci_rank_hs92"])) for r in rows("eci_facts.csv")}
close(eci[("BOL", 2006)][0], -1.21, 2); close(eci[("BOL", 2019)][0], -1.30, 2); close(eci[("KOR", 2019)][0], 1.81, 2)
assert eci[("BOL", 2019)][1] == 131
pig = [float(r["pi"]) for r in rows("global_wp_pi.csv") if r["year"] == "2019" and r["country"] == "BOL"][0]
close(pig, 0.00053, 5)
vel = {(r["country"], r["vintage"]): float(r["v_endpoint_per_decade"]) for r in rows("gamma_velocity.csv")}
close(vel[("BOL", "INE base 1990, 1990-2013")], -0.036, 3)
print(f"  ECI 2006 {eci[('BOL', 2006)][0]:+.3f} -> 2019 {eci[('BOL', 2019)][0]:+.3f} (rank {eci[('BOL', 2019)][1]}); "
      f"pi_G {pig:.5f}; v_Gamma {vel[('BOL', 'INE base 1990, 1990-2013')]:+.3f}/decade")
# SCI with Gamma_ref = 0.50 is undefined for Bolivia and Korea: Gamma_t0 is above Gamma_ref
print(f"  SCI denominator Gamma_ref - Gamma_1990 = {0.50 - bol(1990, 1990):+.3f} (< 0: undefined)")
assert 0.50 - bol(1990, 1990) < 0 and 0.50 - comp[("KOR", "ICIO2021", "1995")][0] < 0
print("all sourced values match the text")
