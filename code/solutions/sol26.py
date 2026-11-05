"""Solutions, Chapter 26 (The Neo-Structuralist Synthesis).
Prints every number used in solutions/sol26.tex.  Run from the repo root:
    python3 code/solutions/sol26.py      (numpy + pandas; reads data/processed)
Sourced values: data/processed (retrieved 2026-11-05; see data/README.md).
Simulation parameters, plan costs and plan targets are illustrative."""
import pathlib
import numpy as np
import pandas as pd

DATA = pathlib.Path(__file__).resolve().parents[2] / "data" / "processed"
eci = pd.read_csv(DATA / "eci_facts.csv")
E = eci.set_index(["country", "year"])
wdi = pd.read_csv(DATA / "wdi_long.csv")
facts = pd.read_csv(DATA / "bolivia_facts.csv").set_index("fact_id").value
gc = pd.read_csv(DATA / "gamma_comparison.csv")
g19 = gc[(gc.year == 2019) & (gc.aggregation == "book12_isic") & (gc.normalization == "input_shares")].copy()
g19["iso"] = g19.country.str[:3]
G = g19.set_index("iso").Gamma


def wdi_val(c, ind, y):
    s = wdi[(wdi.country_iso3 == c) & (wdi.indicator == ind) & (wdi.year == y)].value
    return float(s.iloc[0]) if len(s) and not np.isnan(s.iloc[0]) else np.nan


# ---------------------------------------------------------------------------
print("=== 26.1 Conditions (1)-(2) with sourced 2019 Gamma and ECI; Gamma* = 0.30 (placeholder) ===")
both = sorted(set(G.index) & set(eci[eci.year == 2019].country))
inside, outside = [], []
for c in both:
    e = E.loc[(c, 2019), "eci_hs92"]
    ok1, ok2 = G[c] > 0.30, e > 0
    (inside if ok1 and ok2 else outside).append(c)
    print(f"  {c}: Gamma {G[c]:.3f} ({'>' if ok1 else '<='} 0.30), ECI {e:+.3f} ({'>' if ok2 else '<='} 0)")
print(f"  satisfy (1) and (2): {len(inside)} {inside}; violate (2): {outside}")
assert all(G[c] > 0.30 for c in G.index)
print(f"  all {len(G)} economies: Gamma from {G.min():.3f} ({G.idxmin()}) to {G.max():.3f} ({G.idxmax()})")
assert outside == ["BOL", "CHL", "PER"]
gonly = sorted(set(G.index) - set(both)); eonly = sorted(set(eci[eci.year == 2019].country) - set(G.index))
print(f"  Gamma but no ECI in repo: {gonly}; ECI but no Gamma: {eonly}")
x = np.array([E.loc[(c, 2019), "eci_hs92"] for c in both]); y = np.array([G[c] for c in both])
print(f"  corr(Gamma, ECI) = {np.corrcoef(x, y)[0,1]:.3f} over {len(both)} economies")

# ---------------------------------------------------------------------------
print("\n=== 26.2 Pathology classification inputs ===")
print(f"  Brazil: Gamma {G['BRA']:.3f}, ECI {E.loc[('BRA', 2019), 'eci_hs92']:+.3f}, Gini 2019 (WDI) {wdi_val('BRA', 'SI.POV.GINI', 2019):.1f}")
print(f"  Korea : Gamma {G['KOR']:.3f}, ECI {E.loc[('KOR', 2019), 'eci_hs92']:+.3f}, Gini 2018 (WDI) {wdi_val('KOR', 'SI.POV.GINI', 2018):.1f}")
print(f"  Venezuela: ECI 2019 {E.loc[('VEN', 2019), 'eci_hs92']:+.3f}; Gamma n/a (not in ICIO); "
      f"last WDI Gini year {int(wdi[(wdi.country_iso3=='VEN')&(wdi.indicator=='SI.POV.GINI')].dropna(subset=['value']).year.max())}")
print(f"  Costa Rica: Gamma 2019 {G['CRI']:.3f}; ECI not in repo extract")
print(f"  Bolivia: Gamma {G['BOL']:.3f}, ECI {E.loc[('BOL', 2019), 'eci_hs92']:+.3f}; INE Gini {facts['BOL_INE_gini_2005']:.3f} -> "
      f"{facts['BOL_INE_gini_2019']:.3f}; growth 2006-19 {facts['BOL_gdp_growth_mean_2006_2019']:.2f}%")
assert round(G["CRI"], 3) > 0.30

# ---------------------------------------------------------------------------
print("\n=== 26.3 Complementarity: time to reach the target zone (all z_k >= 1), u = 1, alpha = 1 ===")


def T_analytic(gam, eta):
    seq = sum((np.log(1 + eta / (1 + gam * k)) / eta) if eta > 0 else 1 / (1 + gam * k) for k in range(3))
    a = 2 * gam + eta
    sim = 3 * np.log(1 + a) / a if a > 0 else 3.0
    return seq, sim


def T_numeric(policy, gam, eta, phi=1.0, z0=(0, 0, 0), dt=1e-4):
    z = np.array(z0, float); t = 0.0
    while z.min() < 1:
        need = z < 1
        u = np.zeros(3)
        if policy == "sim":
            u[need] = 1 / need.sum()
        else:
            u[np.argmax(need)] = 1.0
        rate = (1 + gam * (z.sum() - z) + eta * z) * u ** phi
        z = np.minimum(z + dt * rate * need, 1.0)
        t += dt
    return t


rows = []
for gam, eta in ((0.0, 0.0), (1.0, 0.0), (3.0, 0.0), (0.0, 3.0), (1.0, 3.0)):
    seq, sim = T_analytic(gam, eta)
    ns, nm = T_numeric("seq", gam, eta), T_numeric("sim", gam, eta)
    rows.append((gam, eta, seq, sim))
    print(f"  gamma = {gam}, eta = {eta}: sequential {seq:.3f} (numeric {ns:.3f}), simultaneous {sim:.3f} (numeric {nm:.3f})")
    assert abs(ns - seq) < 2e-3 and abs(nm - sim) < 2e-3
print("  T_seq - T_sim on a grid (rows eta, columns gamma):")
grid_g = (0.25, 0.5, 1.0, 3.0, 10.0)
for eta in (0.0, 0.5, 1.0, 3.0, 10.0):
    diffs = [T_analytic(gm, eta)[0] - T_analytic(gm, eta)[1] for gm in grid_g]
    print(f"    eta = {eta:5.2f}: " + "  ".join(f"{d:+.3f}" for d in diffs))
    for gm, d in zip(grid_g, diffs):
        if abs(gm - eta) < 1e-12:
            assert abs(d) < 1e-12
        elif gm > eta:
            assert d > 0
        else:
            assert d < 0
eta_ = 2.0
print(f"  at gamma = eta = {eta_}: T_seq = {T_analytic(eta_, eta_)[0]:.6f}, T_sim = {T_analytic(eta_, eta_)[1]:.6f}, "
      f"ln(1+3 eta)/eta = {np.log(1+3*eta_)/eta_:.6f}")
for phi in (1.5,):
    for gam in (0.0, 1.0, 3.0, 10.0):
        ns, nm = T_numeric("seq", gam, 0.0, phi, dt=1e-3), T_numeric("sim", gam, 0.0, phi, dt=1e-3)
        print(f"  phi = {phi} (increasing returns to intensity), gamma = {gam}: sequential {ns:.3f}, simultaneous {nm:.3f} "
              f"-> {'sequential' if ns < nm else 'simultaneous'} faster")
# Bolivia-like start: condition (1) already met
for gam in (1.0,):
    ns, nm = T_numeric("seq", gam, 0.0, z0=(1, 0, 0)), T_numeric("sim", gam, 0.0, z0=(1, 0, 0))
    print(f"  start z = (1, 0, 0), gamma = {gam}: sequential {ns:.3f}, simultaneous {nm:.3f}")

# ---------------------------------------------------------------------------
print("\n=== 26.4 (a) Gamma-GDP and ECI-growth, cross-country (sourced) ===")
k = wdi[wdi.indicator == "NY.GDP.PCAP.PP.KD"].pivot(index="year", columns="country_iso3", values="value")
kc = wdi[wdi.indicator == "NY.GDP.PCAP.PP.CD"].pivot(index="year", columns="country_iso3", values="value")
cg = [c for c in G.index if c in kc.columns and not np.isnan(kc.loc[2019, c])]
lg = np.log(kc.loc[2019, cg].values.astype(float)); gg = G[cg].values
print(f"  Gamma vs log GDP pc PPP 2019: n = {len(cg)} {cg}, corr = {np.corrcoef(gg, lg)[0,1]:.3f}")
reg = {}
for y0, y1 in ((1995, 2019), (2009, 2019)):
    cs = [c for c in k.columns if (c, y0) in E.index and not np.isnan(k.loc[y0, c])]
    g = np.array([100 * ((k.loc[y1, c] / k.loc[y0, c]) ** (1 / (y1 - y0)) - 1) for c in cs])
    X = np.column_stack([np.ones(len(cs)), [E.loc[(c, y0), "eci_hs92"] for c in cs], np.log([k.loc[y0, c] for c in cs])])
    b = np.linalg.lstsq(X, g, rcond=None)[0]
    r = g - X @ b
    V = (r @ r / (len(g) - 3)) * np.linalg.inv(X.T @ X)
    reg[(y0, y1)] = b[1]
    print(f"  growth {y0}-{y1} (% p.a., GDP pc PPP const.) = {b[0]:.2f} + {b[1]:.3f} ECI_{y0} ({np.sqrt(V[1,1]):.3f}) "
          f"{b[2]:+.3f} ln y_{y0} ({np.sqrt(V[2,2]):.3f}); n = {len(cs)}")
    cs_nc = [c for c in cs if c != "CHN"]
    idx = [cs.index(c) for c in cs_nc]
    bn = np.linalg.lstsq(X[idx], g[idx], rcond=None)[0]
    print(f"    excluding China: ECI coefficient {bn[1]:.3f}")
print(f"  implied growth premium of +1.30 ECI: {1.30*reg[(1995,2019)]:.2f} pp (1995-2019 coef), {1.30*reg[(2009,2019)]:.2f} pp (2009-2019 coef)")

print("\n=== 26.4 (b)-(c) NPV of the plan (costs: Example ch26:ex:plan; ECI path: plan targets) ===")
years = np.arange(2025, 2056)
cost = np.zeros(len(years))
items = [(1.5, 2025, 2034), (2.0, 2030, 2044), (5.0, 2030, 2044), (3.0, 2030, 2044), (1.0, 2030, 2044), (4.0, 2040, 2054)]
for amt, a, b_ in items:
    m = (years >= a) & (years <= b_)
    cost[m] += amt / m.sum()
print(f"  costed items total = {cost.sum():.1f} bn US$ (book rounds to 17)")
assert abs(cost.sum() - 16.5) < 1e-9


def eci_path(t):
    if t <= 2035:
        return -1.30
    if t <= 2045:
        return -1.30 + (t - 2035) * 1.00 / 10
    if t <= 2055:
        return -0.30 + (t - 2045) * 0.30 / 10
    return 0.0


Y0, g0 = facts["BOL_gdp_usd_2024"], facts["BOL_gdp_growth_geo_2006_2024"] / 100


def npv(bcoef, r, T_end=2055):
    yrs = np.arange(2025, T_end + 1)
    base = plan = Y0
    out = 0.0
    for t in yrs:
        base *= 1 + g0
        plan *= 1 + g0 + bcoef / 100 * (eci_path(t - 1) + 1.30)
        c = cost[t - 2025] if t <= 2055 else 0.0
        out += (plan - base - c) / (1 + r) ** (t - 2024)
    return out


def irr(bcoef, T_end=2055):
    lo, hi = 1e-4, 1.0
    if npv(bcoef, lo, T_end) < 0:
        return np.nan
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if npv(bcoef, mid, T_end) > 0:
            lo = mid
        else:
            hi = mid
    return mid


print(f"  baseline: GDP 2024 = {Y0:.1f} bn US$, growth {100*g0:.2f}% (geometric mean 2006-2024)")
for bcoef in (0.0, 0.5, reg[(2009, 2019)], reg[(1995, 2019)]):
    print(f"  b = {bcoef:.3f} pp/ECI: NPV(5%, to 2055) = {npv(bcoef, 0.05):+.1f} bn; IRR = "
          f"{100*irr(bcoef):.1f}%" if not np.isnan(irr(bcoef)) else f"  b = {bcoef:.3f}: NPV(5%) = {npv(bcoef, 0.05):+.1f} bn; IRR n/a")
lo, hi = 0.0, 2.0
for _ in range(60):
    mid = 0.5 * (lo + hi)
    lo, hi = (lo, mid) if npv(mid, 0.05) > 0 else (mid, hi)
print(f"  break-even b at 5% = {mid:.4f} pp per ECI unit")
pv_cost = sum(cost[t - 2025] / 1.05 ** (t - 2024) for t in years)
print(f"  PV of costs at 5% = {pv_cost:.2f} bn")
assert npv(0.0, 0.05) < 0 and npv(reg[(1995, 2019)], 0.05) > 0

# ---------------------------------------------------------------------------
print("\n=== 26.5 Neoliberal counterfactual: extrapolate 1988/1995-2005 velocities ===")
bs = pd.read_csv(DATA / "bolivia_io_summary.csv")
gb = bs[(bs.base == 1990) & (bs.aggregation == "book12") & (bs.flows == "domestic") & (bs.normalization == "input_shares")].set_index("year").Gamma
vG = (gb[2005] - gb[1988]) / 1.7
print(f"  Gamma (INE base 1990): 1988 {gb[1988]:.3f}, 2005 {gb[2005]:.3f}; v = {vG:+.4f} per decade")
for yy in (2013, 2019, 2020):
    print(f"    projected {yy}: {gb[2005] + vG*(yy-2005)/10:.3f}")
print(f"  actual Gamma: 2013 {gb[2013]:.3f} (base 1990); 2019 0.467 (base 2017, not comparable)")
e95, e05 = E.loc[("BOL", 1995), "eci_hs92"], E.loc[("BOL", 2005), "eci_hs92"]
vE = (e05 - e95)
print(f"  ECI 1995 {e95:+.3f}, 2005 {e05:+.3f}; v = {vE:+.3f} per decade; projected 2019 {e05 + vE*1.4:+.3f}, 2020 {e05 + vE*1.5:+.3f}; "
      f"actual 2019 {E.loc[('BOL', 2019), 'eci_hs92']:+.3f}")
pip = pd.read_csv(DATA / "pip_bolivia.csv")
pip = pip[(pip.poverty_line == 3.0) & (pip.reporting_level == "national")].set_index("reporting_year").gini
pre = pip.loc[1997:2005]
print("  PIP Gini 1997-2005 (non-comparable spells): " + ", ".join(f"{y} {v:.3f}" for y, v in pre.items())
      + f"; mean {pre.mean():.3f}; actual 2019 {pip[2019]:.3f}")
assert e05 + vE * 1.4 < E.loc[("BOL", 2019), "eci_hs92"]

# ---------------------------------------------------------------------------
print("\n=== 26.6 Populist counterfactual: size of gas rents, and the fastest observed ECI gains ===")
rent = wdi[(wdi.country_iso3 == "BOL") & (wdi.indicator == "NY.GDP.NGAS.RT.ZS")].set_index("year").value
gdp = wdi[(wdi.country_iso3 == "BOL") & (wdi.indicator == "NY.GDP.MKTP.CD")].set_index("year").value
tot = sum(rent[y] / 100 * gdp[y] for y in range(2006, 2020)) / 1e9
tot14 = sum(rent[y] / 100 * gdp[y] for y in range(2006, 2015)) / 1e9
print(f"  natural gas rents (WDI NY.GDP.NGAS.RT.ZS x NY.GDP.MKTP.CD), 2006-2019: {tot:.2f} bn current US$ (2006-2014: {tot14:.2f})")
assert 10 < tot < 14
best = (-9, None)
for c in sorted(set(eci.country)):
    for y in range(1995, 2011):
        if (c, y) in E.index and (c, y + 14) in E.index:
            d = E.loc[(c, y + 14), "eci_hs92"] - E.loc[(c, y), "eci_hs92"]
            if d > best[0]:
                best = (d, (c, y, y + 14))
print(f"  largest 14-year ECI gain among the {eci.country.nunique()} economies in the repo: {best[0]:+.3f} ({best[1]})")
print(f"  Korea 1995-2019: {E.loc[('KOR', 2019), 'eci_hs92'] - E.loc[('KOR', 1995), 'eci_hs92']:+.3f}")
print(f"  Bolivia 2005 + best gain = {e05 + best[0]:+.3f}")
print(f"  linear scaling of the plan's ECI effect (1.30 per 17 bn): $12bn -> {1.30*tot/17:+.2f}; $50bn -> {1.30*50/17:+.2f}")
print("all chapter-26 solution checks passed")
