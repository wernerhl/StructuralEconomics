"""Chapter 3: reproduces every number in the worked examples, tables, viz panel
and exercise answers.  Convention (Ch1): a_ij = input j per unit of i, L = (I-A')^-1,
W^P = D^-1 A (row-stochastic input shares).  Run: python3 code/ch03.py (numpy only)."""
import numpy as np
np.set_printoptions(precision=4, suppress=True)
def hdr(s): print("\n== " + s)
# ---- Sourced data (data/processed, retrieved 2026-11-05; see data/README.md) ----
import csv, pathlib
DATA = pathlib.Path(__file__).resolve().parent.parent / "data" / "processed"
def rows(name):
    with open(DATA / name, newline="", encoding="utf-8") as f: return list(csv.DictReader(f))
def fact(fid):
    return float(next(r["value"] for r in rows("bolivia_facts.csv") if r["fact_id"] == fid))
def bol_summary(year, agg, norm="input_shares", flows="domestic", base=None):
    base = base or ("2017" if int(year) >= 2017 else "1990")
    return next(r for r in rows("bolivia_io_summary.csv") if (r["base"], r["year"], r["aggregation"], r["flows"], r["normalization"]) == (base, str(year), agg, flows, norm))
def oecd_summary(iso, year=2019, agg="book12_isic", edition="ICIO2023"):
    return next(r for r in rows("oecd_io_summary.csv") if (r["edition"], r["country"], r["year"], r["aggregation"], r["flows"], r["normalization"]) == (edition, iso, str(year), agg, "domestic", "input_shares"))
def bol_pi(sector, year=2019, agg="book12"):
    return float(next(r["pi"] for r in rows("bolivia_io_pi.csv") if (r["base"], r["year"], r["aggregation"], r["flows"], r["normalization"], r["sector"]) == ("2017", str(year), agg, "domestic", "input_shares", sector)))
def sector_fact(sector, col, year=2019):
    return float(next(r[col] for r in rows("bolivia_io_sector_facts.csv") if (r["base"], r["year"], r["book12"]) == ("2017", str(year), sector)))
def check(label, value, text, nd):
    """Assert that the sourced value rounds to the number printed in the text (nd decimals)."""
    print("  %-58s %10.4f  text: %s" % (label, value, text)); assert round(value, nd) == round(text, nd), (label, value, text)

def Linv(A, g=1.0): return np.linalg.inv(np.eye(len(A)) - g * A.T)
def WP(A): return A / A.sum(1, keepdims=True)
def perron(W):
    w, v = np.linalg.eig(W.T); p = np.real(v[:, np.argmin(abs(w - 1))]); return p / p.sum()
def gap(W):
    e = sorted(abs(np.linalg.eigvals(W)), reverse=True); return 1 - e[1]
def cross_mass(A, blocks):
    tot = A.sum(); c = sum(A[i, j] for i in range(len(A)) for j in range(len(A)) if blocks[i] != blocks[j]); return c / tot

hdr("Example two_WP / propagation (rows: Modern, Traditional)")
AA = np.array([[.10, .25], [.20, .10]]); AB = np.array([[.10, .02], [.02, .10]])
for name, A in [("A (connected)", AA), ("B (disconnected)", AB)]:
    L = Linv(A); print(name, "L = (I-A')^-1 =\n", L.round(3), "\n  modern final-demand impulse -> (M, T) =", L[:, 0].round(3),
                       " rho(A) =", round(max(abs(np.linalg.eigvals(A))), 3))
    assert np.allclose(L, sum(np.linalg.matrix_power(A.T, k) for k in range(300)))

hdr("Bolivia stylized example")
Ab = np.array([[.15, .03], [.02, .20]]); Ac = np.array([[.15, .15], [.15, .20]])
for name, A in [("enclave", Ab), ("connected", Ac)]:
    L = Linv(A); col = L[:, 0]; W = WP(A); p = perron(W)
    ev = np.sort(np.linalg.eigvals(A).real)[::-1]
    print(name, "L =\n", L.round(4), "\n  gas impulse ->", col.round(4), "multiplier", col.sum().round(4),
          "share to traditional", round(col[1] / col.sum(), 4), " eig(A) =", ev.round(4), "rho-|l2| (old) =", round(ev[0] - abs(ev[1]), 4))
    print("  W^P = D^-1 A =", W.round(4).tolist(), " Gamma(W^P) =", round(gap(W), 4), " = w12+w21 =", round(W[0, 1] + W[1, 0], 4),
          " pi(W^P) =", p.round(4), " cross-block mass of A =", round(cross_mass(A, [0, 1]), 4))
    assert np.isclose(gap(W), W[0, 1] + W[1, 0])
Lb, Lc = Linv(Ab), Linv(Ac); print("ratio of spillovers connected/enclave =", round(Lc[1, 0] / Lb[1, 0], 2))
m = max(Lb.max(), Lc.max())
for name, L in [("enclave", Lb), ("connected", Lc)]:
    print("heatmap", name, ",".join("%d/%d/%d/%.2f" % (i + 1, j + 1, round(100 * L[i, j] / m), L[i, j]) for i in range(2) for j in range(2)))
A06 = np.diag([0.6, 0.05]); print("diag(0.6,0.05): rho-|l2| =", 0.55, " Gamma of W^P (identity) =", round(gap(WP(A06)), 3))

hdr("Growth-propagation bound")
W = np.array([[.5, .9], [0, .5]]); print("entry 0.9 > rho =", max(abs(np.linalg.eigvals(W))))
C = np.array([[.85, .88, .05], [0, .86, .72], [0, 0, .88]]); rho = max(abs(np.linalg.eigvals(C)))
L = np.linalg.inv(np.eye(3) - C.T); off = (L.sum() - np.trace(L)) / 6
print("3x3 nonnormal: rho =", rho.round(3), "avg off-diag L =", round(off, 2), "old 'bound' ||C||_2/(1-rho) =",
      round(np.linalg.norm(C, 2) / (1 - rho), 2), " ||C||_2 =", round(np.linalg.norm(C, 2), 3), "(new bound needs ||C||<1)")
for name, A in [("Bolivia", Ab), ("connected", Ac), ("A", AA)]:
    n2 = np.linalg.norm(A, 2); L = Linv(A); N = len(A); off = (L.sum() - np.trace(L)) / (N * (N - 1))
    print(name, "avg off-diag =", round(off, 4), "bound g||A||/(1-g||A||) =", round(n2 / (1 - n2), 4)); assert off <= n2 / (1 - n2)
rng = np.random.default_rng(0)
for _ in range(20000):   # random check of the corrected bound
    N = rng.integers(2, 6); A = rng.random((N, N)) * rng.random() * 0.6; n2 = np.linalg.norm(A, 2)
    if n2 >= 0.99: continue
    L = Linv(A); assert (L - np.diag(np.diag(L))).max() <= n2 / (1 - n2) + 1e-12
print("corrected bound holds entrywise on 20000 random draws")

hdr("Exercises")
s = lambda th: th / (1 + th); a = .3
print("3.1 sigma_M, sigma_F =", s(2), s(.5), " sbar =", (s(2) + s(.5)) / 2,
      " K*(sbar) = (sbar(1-a))^(1/(1-a)) ratio vs both theta=2:", round(((s(2) + s(.5)) / 2 / s(2)) ** (1 / (1 - a)), 4))
A2 = np.array([[.10, .20, .05], [.15, .10, .15], [.05, .10, .10]]); L2 = Linv(A2)
print("3.2 L =\n", L2.round(4), "\n  impulse to sector 2 ->", L2[:, 1].round(4), " multipliers", L2.sum(0).round(4),
      " rho(A) =", round(max(abs(np.linalg.eigvals(A2))), 4), " Gamma(W^P) =", round(gap(WP(A2)), 4))
A3 = A2.copy(); A3[2, 0] = A3[2, 1] = A3[0, 2] = A3[1, 2] = .01; L3 = Linv(A3)
print("3.3 L =\n", L3.round(4), "\n  impulse to sector 3 ->", L3[:, 2].round(4), " Gamma(W^P) =", round(gap(WP(A3)), 4),
      " decline factor (spill to 1,2):", (L2[:2, 2] / L3[:2, 2]).round(2))
A4 = np.array([[.15, .10], [.10, .20]]); L4 = Linv(A4); print("3.4 links 0.10: L =", L4.round(4).tolist(), " Gamma(W^P) =", round(gap(WP(A4)), 4))
th = .5
for t in range(1, 41):
    th = .8 * th + .2 * 2 + .1 * (.7 * .1 - .1)
    if t in (1, 5, 10, 20, 40): print("3.7 theta_F(t=%d) = %.4f gap = %.4f" % (t, th, 2 - th))
print("3.7 fixed point with lambda=0.7:", round((.4 + .1 * (.07 - .1)) / .2, 4))

hdr("Proof checks (randomized, seed 3): Prop block_spectral (1)-(3)")
import itertools
rng = np.random.default_rng(3)
for _ in range(300):
    a, b = rng.random(2) * 0.5 + 1e-3; W2 = np.array([[1 - a, a], [b, 1 - b]])          # part (2)
    assert np.isclose(gap(W2), a + b) and np.allclose(perron(W2), np.array([b, a]) / (a + b))
    t = rng.random() * (1 / (a + b)); W2t = np.array([[1 - t * a, t * a], [t * b, 1 - t * b]])
    assert np.isclose(gap(W2t), t * (a + b)) and np.allclose(perron(W2t), perron(W2))
    sizes = rng.integers(1, 4, rng.integers(2, 4)); n = sizes.sum(); Ab = np.zeros((n, n)); o = 0   # part (1)
    for q in sizes:
        Ab[o:o + q, o:o + q] = rng.random((q, q)) + .01; o += q
    ev = np.linalg.eigvals(WP(Ab)); assert np.sum(np.isclose(ev, 1)) == len(sizes) and np.isclose(gap(WP(Ab)), 0)
    n = rng.integers(2, 7); Ksym = rng.random((n, n)) * (rng.random((n, n)) < .7); Ksym = Ksym + Ksym.T + 1e-3 * np.eye(n)  # part (3)
    Wr = Ksym / Ksym.sum(1, keepdims=True); pr = Ksym.sum(1) / Ksym.sum()            # reversible random walk
    assert np.allclose(pr @ Wr, pr) and np.allclose(pr[:, None] * Wr, (pr[:, None] * Wr).T)
    lam = np.sort(np.linalg.eigvals(Wr).real)[::-1]; l2 = lam[1]
    Phi_c = min(sum(pr[i] * Wr[i, j] for i in S for j in range(n) if j not in S) / pr[list(S)].sum()
                for r in range(1, n) for S in itertools.combinations(range(n), r) if pr[list(S)].sum() <= .5 + 1e-12)
    assert Phi_c ** 2 / 2 - 1e-10 <= 1 - l2 <= 2 * Phi_c + 1e-10 and gap(Wr) <= 2 * Phi_c + 1e-10
print("all randomized proof checks passed")
hdr("Sourced Bolivia values quoted in the text (Sec. multi_sector, kuznets, heterogeneity, trickle; summary)")
check("fuel exports % merch. exports 2013 (WDI)", fact("BOL_fuel_exports_share_2013"), 54.5, 1)
check("fuel exports % merch. exports 2019 (WDI)", fact("BOL_fuel_exports_share_2019"), 32.1, 1)
check("fuel exports % merch. exports 2024 (WDI)", fact("BOL_fuel_exports_share_2024"), 18.9, 1)
check("employment share mining & quarrying 2019 (ILOSTAT)", fact("BOL_emp_share_mining_quarrying_incl_hydrocarbons_2019"), 1.3, 1)
check("employment share agriculture 2019 (ILOSTAT)", fact("BOL_emp_share_agriculture_2019"), 27.4, 1)
check("employment share commerce 2019 (ILOSTAT)", fact("BOL_emp_share_commerce_2019"), 17.1, 1)
s19 = bol_summary(2019, "book12"); N = int(s19["N"])
check("pi hydrocarbons 2019, book12 W^P", bol_pi("Hydrocarbons"), 0.039, 3)
check("uniform share 1/N, N = 11", 1 / N, 0.091, 3); assert N == 11 and bol_pi("Hydrocarbons") < 1 / N
check("hydrocarbons % of output exported 2019", 100 * sector_fact("Hydrocarbons", "export_share_of_output"), 70, 0)
check("hydrocarbons domestic inputs per unit output 2019", sector_fact("Hydrocarbons", "dom_int_inputs_per_output"), 0.225, 3)
check("manufacturing domestic inputs per unit output 2019", sector_fact("Manufacturing", "dom_int_inputs_per_output"), 0.558, 3)
check("hydrocarbons compensation % of VA 2019", 100 * sector_fact("Hydrocarbons", "comp_share_of_va"), 6, 0)
check("Gini INE 2005", fact("BOL_INE_gini_2005"), 0.600, 3); check("Gini INE 2019", fact("BOL_INE_gini_2019"), 0.426, 3)
check("Gini PIP 2005", fact("BOL_PIP_gini_2005"), 0.585, 3); check("Gini PIP 2019", fact("BOL_PIP_gini_2019"), 0.416, 3)
check("extreme poverty 2005, old lines (INE)", fact("BOL_INE_extreme_poverty_oldline_2005"), 38.2, 1)
check("extreme poverty 2018, old lines (INE)", fact("BOL_INE_extreme_poverty_oldline_2018"), 15.2, 1)
check("Gamma 2006, book12, base 1990", float(bol_summary(2006, "book12")["Gamma"]), 0.518, 3)
check("Gamma 2013, book12, base 1990", float(bol_summary(2013, "book12")["Gamma"]), 0.484, 3)
check("Gamma 2019, book12", float(s19["Gamma"]), 0.467, 3); check("N_eff 2019, book12", float(s19["Neff"]), 6.48, 2)
check("Gamma 2019, book12_isic", float(bol_summary(2019, "book12_isic")["Gamma"]), 0.475, 3)
for iso, g in [("USA", .463), ("KOR", .511), ("DEU", .544), ("CHL", .557), ("CHN", .582), ("MEX", .649)]:
    check("Gamma 2019 %s, ICIO 2023 book12_isic" % iso, float(oecd_summary(iso)["Gamma"]), g, 3)
gs = [float(bol_summary(2019, agg, nm)["Gamma"]) for agg in ("book12", "detail39")
      for nm in ("input_shares", "input_shares_zero_diag", "sales_shares", "sales_shares_zero_diag")]
check("Gamma 2019 range over 8 variants: min", min(gs), 0.33, 2); check("Gamma 2019 range over 8 variants: max", max(gs), 0.69, 2)
print("\nall assertions passed")
