"""Chapter 23: reproduces every number used in the chapter's examples,
tables, viz panel and the inhomogeneous-products section.
Sourced trajectories are read from data/processed/*.csv (csv module) and
checked against the numbers printed in chapters/chapter23.tex.
Run from the repo root:  python3 code/ch23.py   (numpy only)."""
import csv
import pathlib
import numpy as np

np.set_printoptions(precision=3, suppress=True)


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def tau1(W):
    n = len(W)
    return max(0.5 * np.abs(W[i] - W[j]).sum() for i in range(n) for j in range(n))


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmax(np.real(ev))
    p = np.real(vl[:, k])
    return p / p.sum()


print("== Section: products of time-varying W ==")
X = np.array([[1, 0, 0], [0, 0, 1], [0.85, 0, 0.15]])
Y = np.array([[0, 0.95, 0.05], [0, 1, 0], [1, 0, 0]])
for name, M in [("X", X), ("Y", Y), ("XY", X @ Y), ("YX", Y @ X)]:
    assert np.allclose(M.sum(1), 1)
    print(f"{name}: Gamma = {gap(M):.3f}, tau1 = {tau1(M):.3f}, |lambda2| = {abs(eig_sorted(M)[1]):.3f}")
assert abs(gap(X) - 0.850) < 1e-3 and abs(gap(Y) - 0.776) < 1e-3
assert abs(gap(X @ Y) - 0.042) < 1e-3
# realized-product gap per period
print(f"realized two-period gap 1-|l2(XY)|^(1/2) = {1-abs(eig_sorted(X@Y)[1])**0.5:.3f}")

# tau1 properties on random matrices
rng = np.random.default_rng(23)
for _ in range(2000):
    A = rng.dirichlet(np.ones(4) * 0.5, size=4)
    B = rng.dirichlet(np.ones(4) * 0.5, size=4)
    assert tau1(A @ B) <= tau1(A) * tau1(B) + 1e-12
    assert abs(eig_sorted(A)[1]) <= tau1(A) + 1e-12
print("checked: tau1(AB) <= tau1(A) tau1(B) and |lambda2| <= tau1 on 2000 random pairs")

# Proof checks for Proposition ch23:prop:dobrushin
# Step 1: ||delta' A||_1 <= tau1(A) ||delta||_1 on zero-sum delta, equality at e_i - e_j
n_complex = 0
for _ in range(1000):
    n = int(rng.integers(2, 7))
    A = rng.dirichlet(np.ones(n) * 0.5, size=n)
    d = rng.normal(size=n); d -= d.mean()
    assert np.abs(d @ A).sum() <= tau1(A) * np.abs(d).sum() + 1e-12
    i, j = max(((i, j) for i in range(n) for j in range(n)),
               key=lambda ij: np.abs(A[ij[0]] - A[ij[1]]).sum())
    e = np.zeros(n); e[i], e[j] = 1, -1
    assert abs(np.abs(e @ A).sum() - 2 * tau1(A)) < 1e-12
    # (ii) every eigenvalue other than one copy of 1 (including complex ones) has modulus <= tau1
    ev = eig_sorted(A)
    n_complex += int(np.any(np.abs(np.imag(ev)) > 1e-9))
    assert np.all(np.abs(ev[1:]) <= tau1(A) + 1e-10)
# repeated unit eigenvalue (reducible): tau1 = 1, Gamma = 0
R = np.array([[0.6, 0.4, 0, 0], [0.3, 0.7, 0, 0], [0, 0, 0.5, 0.5], [0, 0, 0.2, 0.8]])
assert abs(tau1(R) - 1) < 1e-12 and abs(gap(R)) < 1e-12
# Hajnal bound: prod tau1 <= exp(-sum(1 - tau1))
taus = rng.uniform(0, 1, 200)
assert np.prod(taus) <= np.exp(-np.sum(1 - taus)) + 1e-300
print(f"checked: tau1 operator-norm identity, eigenvalue bound ({n_complex} cases with complex spectrum), reducible case, Hajnal bound")

# weak ergodicity: product of scrambling matrices collapses rows; absolute probability sequence
T = 60
Ws = [rng.dirichlet(np.ones(4), size=4) for _ in range(T)]
P = np.eye(4)
bound = 1.0
for t in range(T):
    P = Ws[t] @ P           # backward product W_{t} ... W_0
    bound *= tau1(Ws[t])
print(f"tau1 of 60-period backward product = {tau1(P):.2e} <= prod tau1 = {bound:.2e}")
assert tau1(P) <= bound + 1e-15
# absolute probability sequence: pi_s' = row of lim P (any row), pi_{t+1}' W_t = pi_t'
pis = [None] * (T + 1)
pis[T] = np.ones(4) / 4     # any terminal vector; backward recursion
for t in range(T - 1, -1, -1):
    pis[t] = pis[t + 1] @ Ws[t]
b0 = rng.uniform(0, 1, 4)
b = b0.copy()
for t in range(T):
    b = Ws[t] @ b
print(f"consensus after 60 periods = {b.mean():.4f}; pi_0'b_0 = {pis[0] @ b0:.4f}")
assert np.allclose(b, pis[0] @ b0, atol=1e-8)

print("\n== Perron-vector sensitivity (perturbation lemma) ==")


def two_block(c, m=3):
    W = np.full((2 * m, 2 * m), c / m)
    W[:m, :m] = (1 - c) / m
    W[m:, m:] = (1 - c) / m
    return W


sens = {}
for G in (0.13, 0.45):
    W = two_block(G / 2)
    p = perron(W)
    Z = np.linalg.inv(np.eye(6) - W + np.outer(np.ones(6), p))
    E = np.zeros((6, 6))
    E[0, 1] = -1.0
    E[0, 4] = 1.0               # move weight from within block to cross block in one row
    dpi = p @ E @ Z             # perturbation lemma: dpi' = pi' E Z
    h = 1e-6
    fd = (perron(W + h * E) - p) / h
    assert np.allclose(dpi, fd, atol=1e-4)
    sens[G] = dpi[:3].sum()
    print(f"Gamma = {gap(W):.2f}: d(pi_block1)/d(delta) = {sens[G]:+.3f}")
print(f"sensitivity ratio Gamma 0.13 vs 0.45 = {sens[0.13]/sens[0.45]:.2f}")

print("\n== HP filter: gain-1/2 cutoff period for annual data ==")
for lam in (6.25, 100, 400, 1600):
    w = np.arccos(1 - 1 / (2 * np.sqrt(lam)))
    print(f"lambda = {lam:7.2f}: cutoff period = {2*np.pi/w:5.1f} years")



# ============================================================================
# Sourced trajectories (data/processed).  Gamma is that of the 12-sector
# domestic input-share matrix W^P = D^-1 A, diagonal kept: `book12` for
# Bolivia's own series, `book12_isic` (ISIC-based) for cross-country use.
# ============================================================================
ROOT = pathlib.Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
TEX = (ROOT / "chapters" / "chapter23.tex").read_text(encoding="utf-8")


def rows(name):
    with open(PROC / name, encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def rnd(x, nd):                # round half away from zero, as in the text
    return round(x + (1e-9 if x >= 0 else -1e-9), nd)


def r2(x):
    return rnd(x, 2)


def r3(x):
    return rnd(x, 3)


def in_tex(s):
    assert s in TEX, f"not found in chapter23.tex: {s!r}"


BOL = {"1990": {}, "2017": {}}          # base -> year -> Gamma (book12)
BOL_ISIC, BOL_N, BOL_G = {}, {}, {}
for r in rows("bolivia_io_summary.csv"):
    if r["flows"] == "domestic" and r["normalization"] == "input_shares":
        y = int(r["year"])
        if r["aggregation"] == "book12":
            BOL[r["base"]][y] = float(r["Gamma"])
            BOL_N[(r["base"], y)] = int(r["N"])
            BOL_G[(r["base"], y)] = float(r["power_gini"])
        elif r["aggregation"] == "book12_isic":
            BOL_ISIC[y] = float(r["Gamma"])
PIH = {"1990": {}, "2017": {}}
for r in rows("bolivia_io_pi.csv"):
    if (r["aggregation"], r["flows"], r["normalization"], r["sector"]) == \
            ("book12", "domestic", "input_shares", "Hydrocarbons"):
        PIH[r["base"]][int(r["year"])] = float(r["pi"])
ICIO = {}                                # (edition, iso3) -> year -> (Gamma, G(pi))
for r in rows("oecd_io_summary.csv"):
    if (r["aggregation"], r["flows"], r["normalization"]) == ("book12_isic", "domestic", "input_shares"):
        ICIO.setdefault((r["edition"], r["country"]), {})[int(r["year"])] = \
            (float(r["Gamma"]), float(r["power_gini"]))
VEL = {(r["country"], r["vintage"]): r for r in rows("gamma_velocity.csv")}
EPOCH = {r["epoch"]: r for r in rows("eci_epochs_bolivia.csv")}
ECI = {(r["country"], int(r["year"])): r for r in rows("eci_facts.csv")}
FACT = {r["fact_id"]: float(r["value"]) for r in rows("bolivia_facts.csv")}
GPI = {(int(r["year"]), r["country"]): float(r["pi"]) for r in rows("global_wp_pi.csv") if r["country"] == "BOL"}
b90, b17 = BOL["1990"], BOL["2017"]

print("\n== Bolivia: Gamma by year (INE, book12, domestic input shares) ==")
for base in ("1990", "2017"):
    print(f"base {base}:", {y: r3(g) for y, g in sorted(BOL[base].items())})
assert 2010 not in b90 and 2015 not in b90 and 2016 not in b17      # no 2010, 2015, 2016 tables
assert min(b90) == 1988 and max(b90) == 2014 and min(b17) == 2017 and max(b17) == 2024
assert all(abs(BOL_ISIC[y] - b17[y]) < 0.03 for y in b17)
for y, v in {1988: 0.563, 1990: 0.567, 1993: 0.573, 1994: 0.567, 1999: 0.642, 2000: 0.494,
             2005: 0.502, 2006: 0.518, 2009: 0.504, 2011: 0.478, 2013: 0.484, 2014: 0.361}.items():
    assert r3(b90[y]) == v, (y, b90[y])
for y, v in {2017: 0.506, 2019: 0.467, 2023: 0.381, 2024: 0.419}.items():
    assert r3(b17[y]) == v, (y, b17[y])
assert r3(BOL_ISIC[2019]) == 0.475 and min(b17, key=b17.get) == 2023
# number of sectors with a defined W^P row: 10 in diagonal-make years, 11 otherwise (base 1990)
assert BOL_N[("1990", 2006)] == 11 and BOL_N[("1990", 2007)] == 10 and BOL_N[("1990", 2011)] == 11

print("\n== Bolivia: epoch table ==")
ep = [("1. Adjustment", "1990", 1988, 1993, "+0.01"), ("2. Gas discovery", "1990", 1994, 2005, "-0.07"),
      ("3. Nationalization", "1990", 2006, 2009, "-0.01"), ("4a. Redistribution", "1990", 2011, 2013, "+0.01"),
      ("4b. Redistribution", "2017", 2017, 2019, "-0.04"), ("5. Crisis", "2017", 2019, 2024, "-0.05")]
for name, base, y0, y1, txt in ep:
    g0, g1 = BOL[base][y0], BOL[base][y1]
    d = g1 - g0
    print(f"{name:20s} base {base} {y0}-{y1}: {g0:.3f} -> {g1:.3f}, dGamma = {d:+.3f}")
    assert f"{r2(d):+.2f}" == txt, (name, d)
    in_tex(f"{r3(g0):.3f}") ; in_tex(f"{r3(g1):.3f}")
    in_tex(f"${txt}$")
d_eci = {k: float(v["d_eci_hs92"]) for k, v in EPOCH.items()}
print("dECI by epoch (HS92):", {k: r2(v) for k, v in d_eci.items()})
for key, txt in (("2. Gas discovery 1994-2005 (1995->2005)", "-0.24"), ("3. Nationalization 2006-2009", "-0.14"),
                 ("4. Redistribution 2010-2019", "-0.02"), ("2019->2024", "+0.20"), ("Full 1995-2019", "-0.44")):
    assert f"{r2(d_eci[key]):+.2f}" == txt, (key, d_eci[key])
    in_tex(f"${txt}$")
eci = lambda c, y: float(ECI[(c, y)]["eci_hs92"])
for (c, y), txt in {("BOL", 1995): "-0.85", ("BOL", 2010): "-1.28", ("BOL", 2019): "-1.30", ("BOL", 2024): "-1.10",
                    ("KOR", 1995): "1.10", ("KOR", 2014): "1.70", ("KOR", 2019): "1.81",
                    ("CHN", 2010): "1.15", ("CHN", 2019): "1.30", ("CHL", 2019): "-0.06",
                    ("NOR", 1995): "1.24", ("NOR", 2014): "1.01"}.items():
    assert f"{r2(eci(c, y)):.2f}" == txt, (c, y, eci(c, y))
assert abs(eci("CHL", 1995) - 0.1) < 0.006                               # "about 0.1"
assert int(ECI[("BOL", 2010)]["eci_rank_hs92"]) == 127 and int(ECI[("BOL", 2019)]["eci_rank_hs92"]) == 131
assert f"{r2(eci('KOR', 2019) - eci('KOR', 1995)):+.2f}" == "+0.70"

print("\n== Bolivia: hydrocarbons pi in W^P, and sector facts ==")
print("base 1990:", {y: r3(v) for y, v in sorted(PIH["1990"].items())})
print("base 2017:", {y: r3(v) for y, v in sorted(PIH["2017"].items())})
for (base, y), v in {("1990", 1988): 0.025, ("1990", 1994): 0.022, ("1990", 2003): 0.047, ("1990", 2005): 0.074,
                     ("1990", 2006): 0.064, ("1990", 2007): 0.087, ("1990", 2009): 0.077, ("1990", 2013): 0.085,
                     ("1990", 2014): 0.117, ("2017", 2017): 0.045, ("2017", 2019): 0.039,
                     ("2017", 2024): 0.021}.items():
    assert r3(PIH[base][y]) == v, (base, y, PIH[base][y])
for fid, v in {"BOL_hydrocarbons_export_share_of_output_2000_b1990": 22.0,
               "BOL_hydrocarbons_export_share_of_output_2006_b1990": 69.2,
               "BOL_hydrocarbons_export_share_of_output_2019_b2017": 70.1,
               "BOL_hydrocarbons_dom_inputs_per_output_2006_b1990": 0.381,
               "BOL_hydrocarbons_dom_inputs_per_output_2019_b2017": 0.225,
               "BOL_hydrocarbons_va_share_2000_b1990": 3.5, "BOL_hydrocarbons_va_share_2006_b1990": 7.9,
               "BOL_hydrocarbons_va_share_2013_b1990": 10.0, "BOL_hydrocarbons_va_share_2017_b2017": 6.3,
               "BOL_hydrocarbons_va_share_2019_b2017": 6.1, "BOL_hydrocarbons_va_share_2024_b2017": 3.3,
               "BOL_fuel_exports_share_2019": 32.1, "BOL_fuel_exports_share_2024": 18.9,
               "BOL_INE_gini_2005": 0.600, "BOL_INE_gini_2019": 0.426,
               "BOL_INE_extreme_poverty_oldline_2005": 38.2, "BOL_INE_extreme_poverty_oldline_2018": 15.2,
               "BOL_gdp_growth_mean_1990_2019": 4.0}.items():
    nd = len(str(v).split(".")[1])
    assert rnd(FACT[fid], nd) == v, (fid, FACT[fid], v)
assert f"{GPI[(1995, 'BOL')]:.5f}" == "0.00018" and f"{GPI[(2006, 'BOL')]:.5f}" == "0.00047"
assert r3(BOL_G[("2017", 2019)]) == 0.445

print("\n== Break test on Bolivia's base-1990 series (Definition ch23:def:break) ==")
ys = sorted(b90)
dG = {y: b90[y] - b90[p] for p, y in zip(ys, ys[1:]) if y - p == 1}
quiet = [dG[y] for y in dG if y <= 1997]
sd = float(np.std(quiet, ddof=1))
print(f"sd of year-to-year changes 1989-1997 = {sd:.4f}")
assert r2(sd) == 0.01
flag2 = sorted(y for y in dG if abs(dG[y]) > 2 * sd)
flag3 = sorted(y for y in dG if abs(dG[y]) > 3 * sd)
big = sorted(y for y in dG if abs(dG[y]) > 0.05)
print("flagged k=3:", flag3, " moves > 0.05:", {y: r2(dG[y]) for y in big})
assert big == [1999, 2000, 2014] and flag3 == [1999, 2000, 2007, 2012, 2014]
assert [f"{r2(dG[y]):+.2f}" for y in big] == ["+0.10", "-0.15", "-0.12"]
assert f"{r2(dG[2007]):+.2f}" == "-0.03" and f"{r2(dG[2012]):+.2f}" == "+0.03"
# make matrices are published only for 2002, 2003, 2006, 2011-2014: N = 11 in those years, else 10
assert all(BOL_N[("1990", y)] == (11 if y in (2002, 2003, 2006, 2011, 2012, 2013, 2014) else 10) for y in b90)

print("\n== Korea and comparators (OECD ICIO, book12_isic) ==")
kor21, kor23 = ICIO[("ICIO2021", "KOR")], ICIO[("ICIO2023", "KOR")]
for (ed, y), v in {("ICIO2021", 1995): 0.626, ("ICIO2021", 2000): 0.517, ("ICIO2021", 2007): 0.397,
                   ("ICIO2021", 2014): 0.522, ("ICIO2021", 2011): 0.508, ("ICIO2023", 2011): 0.491,
                   ("ICIO2023", 2019): 0.511, ("ICIO2023", 2020): 0.541}.items():
    assert r3(ICIO[(ed, "KOR")][y][0]) == v, (ed, y)
assert min(kor21, key=lambda y: kor21[y][0]) == 2007
assert r3(kor21[1995][1]) == 0.618 and r3(kor23[2019][1]) == 0.658
gam19 = {c: ICIO[("ICIO2023", c)][2019][0] for c in ("KOR", "DEU", "USA", "CHL", "CHN", "MEX")}
assert [r3(gam19[c]) for c in ("KOR", "DEU", "USA", "CHL", "CHN", "MEX")] == [0.511, 0.544, 0.463, 0.557, 0.582, 0.649]
assert r2(min(BOL_ISIC[2019], gam19["USA"], gam19["KOR"], gam19["DEU"])) == 0.46 and r2(gam19["DEU"]) == 0.54
assert r2(gam19["KOR"] - BOL_ISIC[2019]) == 0.04
nor = ICIO[("ICIO2021", "NOR")]
assert r3(nor[1995][0]) == 0.644 and r3(nor[2014][0]) == 0.313

print("\n== Structural velocities (gamma_velocity.csv; endpoints per decade) ==")
vel_txt = {("KOR", "ICIO2021"): "-0.055", ("KOR", "ICIO2023"): "+0.055", ("CHL", "ICIO2021"): "+0.091",
           ("CHL", "ICIO2023"): "+0.107", ("ARG", "ICIO2021"): "+0.037", ("ARG", "ICIO2023"): "-0.050",
           ("DEU", "ICIO2021"): "-0.002", ("DEU", "ICIO2023"): "+0.043", ("USA", "ICIO2021"): "-0.010",
           ("USA", "ICIO2023"): "-0.016", ("CHN", "ICIO2023"): "-0.038", ("BRA", "ICIO2021"): "+0.020",
           ("COL", "ICIO2021"): "-0.036", ("PER", "ICIO2021"): "-0.060", ("PER", "ICIO2023"): "-0.047",
           ("BOL", "INE base 1990, 1990-2013"): "-0.036"}
for k, txt in vel_txt.items():
    v = float(VEL[k]["v_endpoint_per_decade"])
    print(f"{k[0]} {k[1]:28s} v = {v:+.3f}  (text {txt})")
    assert f"{r3(v):+.3f}" == txt, (k, v)
    in_tex(f"${txt}$")
assert ("CHN", "ICIO2021") not in VEL                                   # China split in the 2021 edition
assert f"{r3(float(VEL[('BOL', 'INE base 1990, 1990-2013')]['v_ols_per_decade'])):+.3f}" == "-0.046"
v17 = VEL[("BOL", "INE base 2017")]
assert f"{r2(float(v17['v_endpoint_per_decade'])):+.2f}" == "-0.12" and f"{r2(float(v17['v_ols_per_decade'])):+.2f}" == "-0.16"
# Korea across the two editions, 1995 -> 2019 (vintage gap at 2011 is small relative to the fall)
assert r3(kor21[2011][0] - kor23[2011][0]) == 0.017 and r3(kor21[1995][0] - kor23[2019][0]) == 0.115

print("\n== SCI with Gamma_ref = 0.50 (gamma_velocity.csv) ==")
ref = 0.50
for c in ("CHL", "ARG", "USA", "KOR", "DEU"):
    r = VEL[(c, "ICIO2021")]
    g0, g1 = float(r["gamma_y0"]), float(r["gamma_y1"])
    sci = (g1 - g0) / (ref - g0) if g0 < ref else None
    print(f"{c}: {g0:.3f} -> {g1:.3f}  SCI = {'undefined' if sci is None else f'{sci:.2f}'}")
    if sci is None:
        assert r["sci_ref050"] == ""
    else:
        assert abs(sci - float(r["sci_ref050"])) < 1e-9
        in_tex(f"${sci:.2f}$")
    in_tex(f"{r3(g0):.3f}") ; in_tex(f"{r3(g1):.3f}")
assert r3(ref - float(VEL[("CHL", "ICIO2021")]["gamma_y0"])) == 0.033
assert VEL[("BOL", "INE base 1990, 1990-2013")]["sci_ref050"] == ""

print("\n== Viz panel: Gamma_t, x = (year - 1988)/4 cm, y = 12.5 (Gamma - 0.30) cm ==")


def coords(series, years):
    return " ".join(f"({(y - 1988) / 4:.2f},{12.5 * (series[y] - 0.30):.2f})" for y in years)


panel = {
    "BOL base 1990, 1988-2009": coords(b90, [y for y in ys if y <= 2009]),
    "BOL base 1990, 2011-2014": coords(b90, [y for y in ys if y >= 2011]),
    "BOL base 2017": coords(b17, sorted(b17)),
    "KOR ICIO 2021": coords({y: v[0] for y, v in kor21.items()}, sorted(kor21)),
    "KOR ICIO 2023": coords({y: v[0] for y, v in kor23.items()}, sorted(kor23)),
}
for k, v in panel.items():
    print(f"{k}: {v}")
    in_tex(v)
print("\nall sourced values match chapters/chapter23.tex")
