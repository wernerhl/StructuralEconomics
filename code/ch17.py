"""Chapter 17: reproduces every number in the chapter's Bolivia examples, tables,
heatmaps, dashboard and trajectory, and checks the identities used.

Sourced values are read from data/processed/*.csv (csv module) and recomputed
where the raw flows are available (INE domestic flow matrices in
data/processed/bolivia_Z/, OECD ICIO domestic blocks in
data/processed/oecd_icio_domestic/).  Every value quoted in the text is
asserted to the rounding used there.  The between-block, information and Psi
examples remain illustrative (inventory classes B and C).

Run from the repo root:  python3 code/ch17.py   (numpy + csv only)."""
import csv
import gzip
import os
from decimal import ROUND_HALF_UP, Decimal

import numpy as np

np.set_printoptions(precision=3, suppress=True)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "data", "processed")


def rd(name):
    with open(os.path.join(PROC, name), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def r(x, d):
    """Round half away from zero on the decimal representation, as quoted in the text
    (so the Growth Lab's -1.295 is quoted as -1.30)."""
    q = Decimal(1).scaleb(-d)
    return float(Decimal(repr(float(x))).quantize(q, rounding=ROUND_HALF_UP))


def check(name, value, text, d):
    assert r(value, d) == text, f"{name}: data {value:.6f} -> {r(value, d)} but text says {text}"


def read_Z(fname):
    with open(os.path.join(PROC, "bolivia_Z", fname), newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    labels = rows[0][1:]
    Z = np.array([[float(v) for v in row[1:]] for row in rows[1:]])
    assert [row[0] for row in rows[1:]] == labels
    return Z, labels


def wp_from_Z(Z, labels, sales=False):
    """W^P = D^-1 A (input shares: row s = Z[:, s] / sum) after dropping, iteratively,
    sectors with no intermediate purchases or sales (data/README.md, section 1)."""
    keep = np.ones(len(labels), bool)
    while True:
        Zs = Z[np.ix_(keep, keep)]
        bad = (Zs.sum(0) <= 0) | (Zs.sum(1) <= 0)
        if not bad.any():
            break
        keep[np.where(keep)[0][bad]] = False
    Zk = Z[np.ix_(keep, keep)]
    K = Zk if sales else Zk.T
    return K / K.sum(1, keepdims=True), [l for l, k in zip(labels, keep) if k]


SUMMARY = rd("bolivia_io_summary.csv")
PI = rd("bolivia_io_pi.csv")
LINK = rd("bolivia_io_linkages.csv")
FACTS = {row["fact_id"]: float(row["value"]) for row in rd("bolivia_facts.csv")}
OECD = rd("oecd_io_summary.csv")


def summ(year, agg="book12", norm="input_shares", flows="domestic"):
    m = [s for s in SUMMARY if s["year"] == str(year) and s["aggregation"] == agg
         and s["normalization"] == norm and s["flows"] == flows]
    assert len(m) == 1
    return m[0]


def oecd(country, year, edition="ICIO2023"):
    m = [s for s in OECD if s["country"] == country and s["year"] == str(year) and s["edition"] == edition
         and s["aggregation"] == "book12_isic" and s["normalization"] == "input_shares"]
    assert len(m) == 1
    return m[0]


def pi_csv(year, agg="book12"):
    return {s["sector"]: float(s["pi"]) for s in PI if s["year"] == str(year) and s["aggregation"] == agg
            and s["flows"] == "domestic" and s["normalization"] == "input_shares"}


# ---------------------------------------------------------------------
print("== Example ch17:ex:bolivia_gamma: INE 2019 (base 2017), domestic, basic prices, book12 ==")
Z12, sectors12 = read_Z("Z_2017_2019_book12.csv")
W, sectors = wp_from_Z(Z12, sectors12)          # 11 sectors: Domestic services has no flows
N = len(sectors)
assert N == 11 and "Domestic services" not in sectors
ev = eig_sorted(W)
p = perron(W)
S19 = summ(2019)
assert abs(gap(W) - float(S19["Gamma"])) < 1e-9 and abs(gini(p) - float(S19["power_gini"])) < 1e-9
assert abs(1 / np.sum(p ** 2) - float(S19["Neff"])) < 1e-9 and S19["dropped_sectors"] == "Domestic services"
print("lambda_1..3 =", np.round(ev[:3], 3), " Gamma =", round(gap(W), 3))
assert abs(ev[1].imag) < 1e-12                         # lambda_2 is real
check("lambda2", ev[1].real, 0.533, 3); check("|lambda3|", abs(ev[2]), 0.313, 3)
check("Gamma 2019 book12", gap(W), 0.467, 3)
Wsales, _ = wp_from_Z(Z12, sectors12, sales=True)
check("Gamma sales shares", gap(Wsales), 0.694, 3)
assert abs(gap(Wsales) - float(summ(2019, norm="sales_shares")["Gamma"])) < 1e-9
print(f"Gamma input shares = {gap(W):.3f}; sales shares = {gap(Wsales):.3f}")
# normalization range: book12 and detail39, domestic, input/sales shares, with/without diagonal
var8 = [float(summ(2019, agg, nm)["Gamma"]) for agg in ("book12", "detail39")
        for nm in ("input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag")]
print(f"Gamma range over the eight variants: {min(var8):.3f}-{max(var8):.3f}")
check("Gamma min", min(var8), 0.33, 2); check("Gamma max", max(var8), 0.69, 2)
check("Gamma detail39", float(summ(2019, "detail39")["Gamma"]), 0.410, 3)
print(f"power Gini G(pi) = {gini(p):.3f}, Neff = {1/np.sum(p**2):.2f} of {N}")
check("G(pi)", gini(p), 0.445, 3); check("Neff", 1 / np.sum(p ** 2), 6.48, 2)

# output shares and Rasmussen linkages from A (never from W^P); x from bolivia_io_linkages.csv
lk = {s["sector"]: s for s in LINK if s["year"] == "2019" and s["base"] == "2017"
      and s["aggregation"] == "book12" and s["flows"] == "domestic"}
x = np.array([float(lk[s]["output"]) for s in sectors12])
A = Z12 / x[None, :]                                   # data convention; L is the book's L
L = np.linalg.inv(np.eye(12) - A)
BL, FL = L.sum(0), L.sum(1)
assert BL.min() >= 1 and FL.min() >= 1 and abs(BL.mean() - FL.mean()) < 1e-12
TABLE_GAMMA = {  # sector: (output share, BL, FL, type) as printed in the chapter
    "Hydrocarbons": (0.047, 1.34, 1.20, "weak"), "Mining": (0.045, 1.36, 1.30, "weak"),
    "Public administration": (0.107, 1.24, 1.02, "weak"), "Manufacturing": (0.208, 1.86, 2.21, "key sector"),
    "Construction": (0.093, 1.57, 1.13, "backward"), "Commerce": (0.078, 1.42, 1.58, "forward"),
    "Transport": (0.063, 1.54, 1.52, "key sector"), "Financial services": (0.041, 1.43, 1.57, "forward"),
    "Agriculture": (0.084, 1.55, 1.70, "key sector"), "Other services": (0.144, 1.36, 2.02, "forward"),
    "Personal services": (0.087, 1.66, 1.10, "backward"), "Domestic services": (0.004, 1.00, 1.00, "weak")}
print("sector: output share, BL, FL, type")
for k, s in enumerate(sectors12):
    sh, bl, fl, ty = TABLE_GAMMA[s]
    assert abs(BL[k] - float(lk[s]["BL"])) < 1e-6 and abs(FL[k] - float(lk[s]["FL"])) < 1e-6
    check(s + " share", x[k] / x.sum(), sh, 3); check(s + " BL", BL[k], bl, 2); check(s + " FL", FL[k], fl, 2)
    assert lk[s]["type"] == ty
    print(f"  {s:22s} {x[k]/x.sum():.3f} {BL[k]:.2f} {FL[k]:.2f} {lk[s]['type']}")
print(f"mean BL = {BL.mean():.3f} = mean FL = {FL.mean():.3f}")
check("mean BL", BL.mean(), 1.44, 2)

print("\n-- comparators: OECD ICIO 2023, 2019, book12_isic, domestic, input shares --")
B19i = summ(2019, "book12_isic")
check("BOL isic Gamma", float(B19i["Gamma"]), 0.475, 3); check("BOL isic G", float(B19i["power_gini"]), 0.435, 3)
COMP = {"KOR": (0.511, 0.658), "DEU": (0.544, 0.611), "USA": (0.463, 0.614), "CHL": (0.557, 0.568),
        "CHN": (0.582, 0.622), "MEX": (0.649, 0.608)}
for c, (g_, G_) in COMP.items():
    o = oecd(c, 2019)
    check(c + " Gamma", float(o["Gamma"]), g_, 3); check(c + " G", float(o["power_gini"]), G_, 3)
    print(f"  {c}: Gamma = {float(o['Gamma']):.3f}, G(pi) = {float(o['power_gini']):.3f}, pi_max = {o['pi_argmax']}")
print(f"  BOL (book12_isic): Gamma = {float(B19i['Gamma']):.3f}, G(pi) = {float(B19i['power_gini']):.3f}")
assert float(B19i["power_gini"]) < min(float(oecd(c, 2019)["power_gini"]) for c in COMP)
assert min(g for g, _ in COMP.values()) < float(B19i["Gamma"]) < max(g for g, _ in COMP.values())

# ---------------------------------------------------------------------
print("\n== Heatmaps: Bolivia (INE) and South Korea (ICIO 2023), 2019, book12_isic, input shares ==")
Zi, lab_i = read_Z("Z_2017_2019_book12_isic.csv")
Wb, sec_b = wp_from_Z(Zi, lab_i)
assert abs(gap(Wb) - float(B19i["Gamma"])) < 1e-9


def icio_group(code):
    """ICIO 2023 industry -> book12_isic group (as data/compute/compute_oecd_io.py)."""
    if code.startswith("A"):
        return "Agriculture"
    if code in ("B05_06", "B09"):
        return "Hydrocarbons"
    for prefix, g_ in (("B", "Mining"), ("C", "Manufacturing"), ("F", "Construction"), ("G", "Commerce"),
                       ("H", "Transport"), ("K", "Financial services"), ("O", "Public administration"),
                       ("T", "Domestic services")):
        if code.startswith(prefix):
            return g_
    if code[0] in "DEJLMNPQ":
        return "Other services"
    if code[0] in "IRS":
        return "Personal services"
    raise KeyError(code)


with gzip.open(os.path.join(PROC, "oecd_icio_domestic", "KOR.csv.gz"), "rt", newline="") as f:
    rows = [row for row in csv.DictReader(f) if row["year"] == "2019"]
gi = {s: k for k, s in enumerate(lab_i)}
Zk = np.zeros((12, 12))
for row in rows:
    Zk[gi[icio_group(row["row"])], gi[icio_group(row["col"])]] += float(row["value"])
Wk, sec_k = wp_from_Z(Zk, lab_i)
assert sec_k == sec_b and len(sec_b) == 11
okor = oecd("KOR", 2019)
assert abs(gap(Wk) - float(okor["Gamma"])) < 1e-6 and abs(gini(perron(Wk)) - float(okor["power_gini"])) < 1e-6
print(f"Bolivia Gamma = {gap(Wb):.3f}; Korea Gamma = {gap(Wk):.3f}")
print(f"largest entries: Bolivia {Wb.max():.3f} ({sec_b[np.unravel_index(Wb.argmax(), Wb.shape)[0]]} row), "
      f"Korea {Wk.max():.3f} ({sec_k[np.unravel_index(Wk.argmax(), Wk.shape)[0]]} row)")
pk = perron(Wk)
print("Korea pi:", {s: round(float(v), 3) for s, v in zip(sec_k, pk)})
check("pi_KOR manuf", pk[sec_k.index("Manufacturing")], 0.434, 3)
print("diagonal share of inputs (mean): Bolivia %.2f, Korea %.2f" % (np.diag(Wb).mean(), np.diag(Wk).mean()))
print("column sums (incoming weight): Bolivia", {s: round(float(v), 2) for s, v in zip(sec_b, Wb.sum(0))})
print("column sums (incoming weight): Korea  ", {s: round(float(v), 2) for s, v in zip(sec_k, Wk.sum(0))})


def tikz(Wm):
    br = np.round(100 * Wm / Wm.max()).astype(int)
    return ",".join(f"{i+1}/{j+1}/{br[i, j]}" for i in range(len(Wm)) for j in range(len(Wm)))


HEAT_BOL, HEAT_KOR = tikz(Wb), tikz(Wk)
check("Gamma BOL/KOR", gap(Wb) / gap(Wk), 0.93, 2)          # "7% below Korea's"
assert np.argmax(Wb.sum(0)) == sec_b.index("Manufacturing") and np.argmax(Wk.sum(0)) == sec_k.index("Manufacturing")
assert sorted(np.argsort(-Wb.sum(0))[:2]) == sorted(np.argsort(-Wk.sum(0))[:2]) == [3, 9]  # manuf, other services
assert np.argmax(Wb[0]) == sec_b.index("Transport")          # hydrocarbons' main supplier
print("Bolivia heatmap triples:", HEAT_BOL)
print("Korea heatmap triples:  ", HEAT_KOR)
check("Wb max", Wb.max(), 0.60, 2); check("Wk max", Wk.max(), 0.66, 2)

# ---------------------------------------------------------------------
print("\n== Example ch17:ex:bolivia_pi: pi (INE 2019, book12) and ILOSTAT 2019 employment ==")
assert all(abs(pi_csv(2019)[s] - v) < 1e-9 for s, v in zip(sectors, p))
PI_TEXT = {"Hydrocarbons": 0.039, "Mining": 0.037, "Public administration": 0.005, "Manufacturing": 0.159,
           "Construction": 0.029, "Commerce": 0.095, "Transport": 0.092, "Financial services": 0.162,
           "Agriculture": 0.091, "Other services": 0.270, "Personal services": 0.021}
for s, v in zip(sectors, p):
    check("pi " + s, v, PI_TEXT[s], 3)
print("pi =", {s: round(float(v), 3) for s, v in zip(sectors, p)})
# Employment rows: ISIC Rev.4 sections mapped to the book's sectors (approximate; ILOSTAT does
# not split hydrocarbons from other mining, nor public from private education and health).
EMP_MAP = {"Hydrocarbons and mining": "B", "Public administration": "O", "Manufacturing": "C",
           "Construction": "F", "Commerce": "G", "Transport": "H", "Financial services": "K",
           "Agriculture": "A", "Other services": "DEJLMNPQ", "Personal services": "IRS",
           "Domestic services": "T"}
ILO = {(row["sex"], row["eco"]): float(row["value_thousands"])
       for row in rd("ilostat_bolivia_emp_eco_sex.csv") if row["year"] == "2019"}
emp_tot = ILO[("SEX_T", "ECO_ISIC4_TOTAL")]
emp_k = {s: sum(ILO[("SEX_T", "ECO_ISIC4_" + c)] for c in cs) for s, cs in EMP_MAP.items()}
fem_k = {s: sum(ILO[("SEX_F", "ECO_ISIC4_" + c)] for c in cs) for s, cs in EMP_MAP.items()}
emp = {s: v / emp_tot for s, v in emp_k.items()}            # employment shares
fem_share = {s: fem_k[s] / emp_k[s] for s in EMP_MAP}        # female share of each sector's employment
print(f"mapped employment covers {sum(emp.values()):.4f} of total (ISIC U, extraterritorial, excluded)")
pd_ = dict(zip(sectors, p))
pi_emp = {s: (pd_["Hydrocarbons"] + pd_["Mining"] if s == "Hydrocarbons and mining" else pd_.get(s, np.nan))
          for s in EMP_MAP}                                   # Domestic services: no W^P row (nan)
ratio = {s: pi_emp[s] / emp[s] for s in EMP_MAP}
EMP_TEXT = {"Hydrocarbons and mining": (0.013, 0.077, 5.7), "Public administration": (0.026, 0.005, 0.19),
            "Manufacturing": (0.112, 0.159, 1.42), "Construction": (0.079, 0.029, 0.37),
            "Commerce": (0.171, 0.095, 0.55), "Transport": (0.071, 0.092, 1.30),
            "Financial services": (0.007, 0.162, 23.4), "Agriculture": (0.274, 0.091, 0.33),
            "Other services": (0.121, 0.270, 2.23), "Personal services": (0.107, 0.021, 0.20),
            "Domestic services": (0.019, None, None)}
print("sector: employment share, pi, pi/share, female share")
for s in EMP_MAP:
    e_, pi_t, ra = EMP_TEXT[s]
    check("emp " + s, emp[s], e_, 3)
    if pi_t is not None:
        check("pi_emp " + s, pi_emp[s], pi_t, 3)
        check("ratio " + s, ratio[s], ra, 1 if ra >= 5 else 2)
    print(f"  {s:24s} {emp[s]:.3f} {pi_emp[s]:.3f} {ratio[s]:6.2f}  F={fem_share[s]:.3f}")
check("ILO B share (README)", FACTS["BOL_emp_share_mining_quarrying_incl_hydrocarbons_2019"] / 100, 0.013, 3)
check("ILO A share (README)", FACTS["BOL_emp_share_agriculture_2019"] / 100, r(emp["Agriculture"], 3), 3)
ranked = sorted((v, s) for s, v in ratio.items() if not np.isnan(v))
print("highest ratio:", ranked[-1][1], "lowest:", ranked[0][1])
assert ranked[-1][1] == "Financial services"
ia = sectors.index("Agriculture")
check("manuf input share from agric", W[sectors.index("Manufacturing"), ia], 0.34, 2)       # "a third"
check("personal input share from agric", W[sectors.index("Personal services"), ia], 0.14, 2)  # "a seventh"
check("1/11", 1 / N, 0.091, 3)
d_ga = pd_["Hydrocarbons"] - pd_["Agriculture"]
print(f"pi_hydro - pi_agr = {d_ga:+.3f}")
check("pi_gas - pi_agr", d_ga, -0.052, 3)
check("agric/hydro+mining workers", emp["Agriculture"] / emp["Hydrocarbons and mining"], 21, 0)

# hydrocarbons in the 2019 SUT (bolivia_facts.csv)
hc = {k: FACTS[f"BOL_hydrocarbons_{k}_2019_b2017"] for k in
      ("export_share_of_output", "dom_inputs_per_output", "comp_share_of_va", "va_share", "output_share")}
print("hydrocarbons 2019:", hc)
check("hydro exported", hc["export_share_of_output"], 70, 0); check("hydro dom inputs", hc["dom_inputs_per_output"], 0.225, 3)
check("hydro comp/VA", hc["comp_share_of_va"], 6, 0); check("hydro VA share", hc["va_share"], 6.1, 1)
sf = {s["book12"]: s for s in rd("bolivia_io_sector_facts.csv") if s["year"] == "2019" and s["base"] == "2017"}
check("manuf dom inputs", float(sf["Manufacturing"]["dom_int_inputs_per_output"]), 0.558, 3)
check("hydro VA share (sector facts)", float(sf["Hydrocarbons"]["va_share"]), 0.061, 3)

# ---------------------------------------------------------------------
print("\n== Between-block variance (illustrative inputs) ==")
y = np.array([18500, 5200, 1800]); w = np.array([0.08, 0.35, 0.57])
m = w @ y; VB = w @ (y - m) ** 2
print(f"mean = {m:.0f}; V_B = {VB/1e6:.2f}M; sigma_B = V_B / 32M = {VB/32e6:.3f}; ratio y1/y3 = {y[0]/y[2]:.1f}")

print("\n== N_eff of four media ecosystems (illustrative) ==")
s = np.array([0.34, 0.37, 0.12, 0.17])
print(f"sum of squares = {np.sum(s**2):.4f}; Neff = {1/np.sum(s**2):.2f} (maximum with 4 nodes = 4); sum of shares = {s.sum():.2f}")
assert abs(1/np.sum(s**2) - 3.4) < 0.05
WI = np.outer(np.ones(4), s)
print(f"W^I = 1 s': pi = {perron(WI)}, Gamma = {1 - np.sort(np.abs(np.linalg.eigvals(WI)))[-2]:.3f}")
check("population 2019 (M)", FACTS["BOL_population_2019_mn"], 11.67, 2)

print("\n== Psi (illustrative) ==")
print(f"Psi = 0.028/0.018 = {0.028/0.018:.2f}; enforcement share 1/(1+Psi) = {1/(1+0.028/0.018):.2f}")
check("self-employment 2019", FACTS["BOL_self_employed_2019"], 68.1, 1)

# ---------------------------------------------------------------------
print("\n== Dashboard 2019 (Gamma, G(pi): book12_isic like-for-like; rest illustrative) ==")
eci = {(row["country"], row["year"]): float(row["eci_hs92"]) for row in rd("eci_facts.csv") if row["eci_hs92"]}
dash = {"Bolivia": (float(B19i["Gamma"]), 0.62, float(B19i["power_gini"]), 3.4, 1.56, eci[("BOL", "2019")]),
        "Korea": (float(okor["Gamma"]), 0.20, float(okor["power_gini"]), 12, 0.6, eci[("KOR", "2019")])}
check("ECI BOL 2019", eci[("BOL", "2019")], -1.30, 2); check("ECI KOR 2019", eci[("KOR", "2019")], 1.81, 2)
check("ECI BOL 2006", eci[("BOL", "2006")], -1.21, 2); check("ECI BOL 2018", eci[("BOL", "2018")], -1.46, 2)
for name, (G_, sB, Gp, ne, ps, ec) in dash.items():
    print(f"{name}: Gamma = {G_:.3f}, 1-sigma_B = {1-sB:.2f}, 1-G(pi) = {1-Gp:.3f}, "
          f"Neff/20 = {ne/20:.2f}, 1/(1+Psi) = {1/(1+ps):.2f}, ECI = {r(ec, 2):+.2f}")
check("bar 1-G BOL", 1 - dash["Bolivia"][2], 0.565, 3); check("bar 1-G KOR", 1 - dash["Korea"][2], 0.342, 3)
print(f"ratios: G(pi) Korea/Bolivia {dash['Korea'][2]/dash['Bolivia'][2]:.2f}x; Gamma Korea/Bolivia "
      f"{dash['Korea'][0]/dash['Bolivia'][0]:.2f}x; sigma_B {0.62/0.20:.1f}x, Neff {12/3.4:.2f}x, Psi {1.56/0.6:.1f}x")
check("G ratio", dash["Korea"][2] / dash["Bolivia"][2], 1.5, 1)
assert eci[("KOR", "2019")] - eci[("BOL", "2019")] > 3          # "more than three standard deviations"

# ---------------------------------------------------------------------
print("\n== Example ch17:ex:trajectory: book12, domestic, input shares (INE base 1990 | base 2017) ==")
TRAJ = {1990: (0.567, 0.420, 10, 0.037), 2000: (0.494, 0.420, 10, 0.036), 2006: (0.518, 0.424, 11, 0.064),
        2011: (0.478, 0.369, 11, 0.088), 2013: (0.484, 0.371, 11, 0.085),
        2017: (0.506, 0.446, 11, 0.045), 2019: (0.467, 0.445, 11, 0.039)}
for yr, (g_, G_, n_, ph) in TRAJ.items():
    base = "1990" if yr < 2017 else "2017"
    Zy, ly = read_Z(f"Z_{base}_{yr}_book12.csv")
    Wy, sy = wp_from_Z(Zy, ly)
    py = perron(Wy)
    sm = summ(yr)
    assert abs(gap(Wy) - float(sm["Gamma"])) < 1e-9 and int(sm["N"]) == len(sy) == n_
    check(f"Gamma {yr}", gap(Wy), g_, 3); check(f"G {yr}", gini(py), G_, 3)
    check(f"pi_hydro {yr}", py[sy.index("Hydrocarbons")], ph, 3)
    print(f"  {yr} (base {base}): N = {len(sy)}, Gamma = {gap(Wy):.3f}, G(pi) = {gini(py):.3f}, "
          f"pi_hydro = {py[sy.index('Hydrocarbons')]:.3f}")
check("Gamma 2005", float(summ(2005)["Gamma"]), 0.502, 3); check("Gamma 2009", float(summ(2009)["Gamma"]), 0.504, 3)
vel = {row["vintage"]: row for row in rd("gamma_velocity.csv") if row["country"] == "BOL"}
v9013 = float(vel["INE base 1990, 1990-2013"]["v_endpoint_per_decade"])
print(f"structural velocity 1990-2013 = {v9013:+.3f} per decade")
check("velocity", v9013, -0.036, 3)
check("KOR 1995", float(oecd("KOR", 1995, "ICIO2021")["Gamma"]), 0.626, 3)
check("KOR 2019", float(okor["Gamma"]), 0.511, 3)

print("\n== Distribution 2005-2019 (INE, WDI) ==")
DIST = {"BOL_INE_extreme_poverty_oldline_2005": 38.2, "BOL_INE_extreme_poverty_oldline_2018": 15.2,
        "BOL_INE_extreme_poverty_EPF_2019": 13.2, "BOL_INE_gini_2005": 0.600, "BOL_INE_gini_2019": 0.426,
        "BOL_gdp_growth_mean_2006_2019": 4.52}
for k, v in DIST.items():
    assert abs(FACTS[k] - v) < 1e-9, k
    print(f"  {k} = {FACTS[k]}")

# ---------------------------------------------------------------------
print("\n== Proof check: Prop. between (exact quadratic decomposition) and Exercise (refinement) ==")
rng_p = np.random.default_rng(17)
for _ in range(200):
    n = 60
    y = rng_p.lognormal(0, 0.6, n)
    kap = 0.9 / y.max()                               # kappa * y < 1
    U = lambda z: z - kap / 2 * z ** 2
    lab = rng_p.integers(0, 3, n)
    def VB_VW(lab):
        ybar = y.mean(); vb = vw = 0.0
        for b in np.unique(lab):
            yb = y[lab == b]; w = len(yb) / n
            vb += w * (yb.mean() - ybar) ** 2; vw += w * yb.var()
        return vb, vw
    vb, vw = VB_VW(lab)
    lhs = np.sum(U(y.mean()) - U(y))
    assert abs(lhs - kap / 2 * n * (vb + vw)) < 1e-9
    assert abs(vb + vw - y.var()) < 1e-12
    fine = lab * 2 + rng_p.integers(0, 2, n)          # refine each block into two
    assert VB_VW(fine)[0] >= vb - 1e-12               # refining never lowers V_B (hence sigma_B)
print("sum[U(ybar) - U(y_i)] = (kappa/2) N (V_B + V_W) and refinement monotonicity verified")
