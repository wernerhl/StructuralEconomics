"""Solutions to the Chapter 11 exercises (The Psi Ratio: Measuring Formation).
Every number quoted in solutions/sol11.tex is printed here.
Run from the repository root:  python3 code/solutions/sol11.py
numpy, plus pandas to read data/processed/ilostat_bolivia_emp_status.csv (Exercise 11.4)."""
import os
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def pvar(v):
    v = np.asarray(v, float)
    return np.mean((v - v.mean()) ** 2)


# ===========================================================================
print("== Exercise 11.1: cohort-year decomposition by hand")
X = np.array([[0.38, 0.42, 0.44],
              [0.50, 0.55, 0.58],
              [0.60, 0.66, 0.72]])     # rows: born 1940s, 1950s, 1960s; columns: 1980, 1990, 2000
gm = X.mean()
ym = X.mean(0)
cm = X.mean(1)
yr = ym - gm
co = cm - gm
print(f"grand mean = {gm:.4f}")
print("(a) year means (average over cohorts) =", ym, " centred year effects =", yr)
print("(b) cohort effects (cohort mean - grand mean) =", co)
fit = gm + yr[None, :] + co[:, None]
resid = X - fit
print("    two-way fitted values:\n", fit, "\n    residuals:\n", resid)
print(f"    Var(year) = {pvar(yr):.6f}, Var(cohort) = {pvar(co):.6f}  (population variances)")
psi0 = pvar(co) / pvar(yr)
print(f"(c) Psi_var = {psi0:.3f}")
assert abs(psi0 - 8.045) < 1e-3
years = np.array([1980, 1990, 2000])
mid = np.array([1945, 1955, 1965])
print("(e) ages (cohort mid-year 1945/55/65):\n", years[None, :] - mid[:, None])
cidx = np.arange(3.0)
tidx = np.arange(3.0)
for k in [-0.10, -0.05, 0.0, 0.05]:
    co_k, yr_k = co + k * cidx, yr - k * tidx
    fit_k = gm + yr_k[None, :] + co_k[:, None] + k * (tidx[None, :] - cidx[:, None])
    assert np.allclose(fit_k, fit)
    print(f"    k = {k:+.2f}: cohort effects {np.round(co_k - co_k.mean(), 4)}, year effects "
          f"{np.round(yr_k - yr_k.mean(), 4)}, age effects A_a = k*a_index, Psi_var = {pvar(co_k)/pvar(yr_k):.3f}")
ratio = lambda k: pvar(co + k * cidx) / pvar(yr - k * tidx)
kd = -np.cov(co, cidx, bias=True)[0, 1] / pvar(cidx)
ka = np.cov(yr, tidx, bias=True)[0, 1] / pvar(tidx)
print(f"    k_delta = {kd:.4f} (ratio {ratio(kd):.5f});  k_alpha = {ka:.4f} (ratio {ratio(ka):.1f})")
# values of k at which Psi_var = 1
f = lambda k: pvar(co + k * cidx) - pvar(yr - k * tidx)
grid = np.linspace(-1, 1, 200001)
vals = np.array([f(k) for k in grid])
roots = grid[np.where(np.diff(np.sign(vals)))[0]]
print("    Psi_var = 1 at k =", np.round(roots, 4))
cdc = np.cov(co, cidx, bias=True)[0, 1]
cat_ = np.cov(yr, tidx, bias=True)[0, 1]
k1 = (pvar(yr) - pvar(co)) / (2 * (cdc + cat_))
print(f"    exact: Cov(delta,c) = {cdc:.5f}, Cov(alpha,t) = {cat_:.5f}, k* = (Var a - Var d)/(2(Cov+Cov)) = {k1:.5f}")
assert abs(ratio(k1) - 1) < 1e-9
print(f"    limit |k| -> inf: Var(c)/Var(t) = {pvar(cidx)/pvar(tidx):.1f}")

# ===========================================================================
print("\n== Exercise 11.2: epidemiological design")
o = np.linspace(0.20, 0.75, 10)
vo = pvar(o)
vo_s = np.var(o, ddof=1)
b = 0.35
expl = b ** 2 * vo
share = expl / (expl + 0.04)
expl_s = b ** 2 * vo_s
print(f"origin rates {np.round(o, 4)}; spacing {o[1]-o[0]:.4f}")
print(f"(c) Var(origin) = {vo:.6f} (population), explained = beta^2 Var = {expl:.6f}, share = {share:.4f}")
print(f"    with the n-1 variance: Var = {vo_s:.6f}, share = {expl_s/(expl_s+0.04):.4f}")
print(f"    a 10 pp higher origin rate -> {100*b*0.10:.1f} pp; range A..J (0.55) -> {b*0.55:.4f}")
assert abs(share - 0.0862) < 1e-4

# ===========================================================================
print("\n== Exercise 11.3: adoption designs (illustrative gradients)")
ga, gb = 0.10, 0.30
print(f"(b) environment share of biological gradient = {ga/gb:.4f}; genes = {gb-ga:.2f}")
for R in [1.0, 0.5, 0.25]:
    print(f"(c) R_theta = {R}: rho = 0.10 / R_theta = {ga/R:.2f}")
# (d) power calculation for a random-roommate design
z975, z80 = 1.959964, 0.841621
for bpeer, sig in [(0.10, 0.8), (0.10, 1.0), (0.05, 0.8)]:
    n = ((z975 + z80) * sig / bpeer) ** 2
    print(f"(d) peer coefficient {bpeer} SD, residual SD {sig}: n = {n:.1f} -> {int(np.ceil(n))} students")

# ===========================================================================
print("\n== Exercise 11.4: Bolivia self-employment (national series only; EH microdata not fetched)")
try:
    import pandas as pd
    path = os.path.join("data", "processed", "ilostat_bolivia_emp_status.csv")
    d = pd.read_csv(path)
    t = d[(d.sex == "SEX_T") & d.status.isin(["STE_AGGREGATE_SLF", "STE_AGGREGATE_TOTAL"])]
    t = t.pivot_table(index=["year", "source"], columns="status", values="value_thousands").reset_index()
    t["share"] = t.STE_AGGREGATE_SLF / t.STE_AGGREGATE_TOTAL
    w = t[(t.year >= 2000) & (t.year <= 2019)]
    for _, row in w.iterrows():
        print(f"   {int(row.year)} {row.source[:40]:40s} self-employed share = {100*row.share:.1f}%")
    print(f"   2000-2019: {len(w)} survey years; mean {100*w.share.mean():.1f}%, min {100*w.share.min():.1f}% "
          f"({int(w.year[w.share.idxmin()])}), max {100*w.share.max():.1f}% ({int(w.year[w.share.idxmax()])}), "
          f"SD of year shares {100*np.sqrt(pvar(w.share)):.2f} pp")
    lfs = w[w.source.str.startswith("LFS")]
    hies = w[w.source.str.startswith("HIES")]
    print(f"   by survey type: LFS mean {100*lfs.share.mean():.1f}% ({len(lfs)} yrs), HIES mean {100*hies.share.mean():.1f}% ({len(hies)} yrs)")
    assert abs(100 * float(w[w.year == 2019].share.iloc[0]) - 66.6) < 0.05
except FileNotFoundError:
    print("   data/processed/ilostat_bolivia_emp_status.csv not found")
import csv
with open(os.path.join("data", "processed", "bolivia_facts.csv"), newline="", encoding="utf-8") as fh:
    facts = {r["fact_id"]: float(r["value"]) for r in csv.DictReader(fh)}
print("   WDI SL.EMP.SELF.ZS 2019 (ILO modelled):", facts["BOL_self_employed_2019"])

# ===========================================================================
print("\n== Exercise 11.5: WWII (illustrative figures 28 -> 36 -> 33)")
l40, l45, l50 = 28.0, 36.0, 33.0
gain, rev, pers = l45 - l40, l45 - l50, l50 - l40
print(f"wartime gain {gain:.0f} pp; reversal 1945-50 {rev:.0f} pp = {100*rev/gain:.1f}% of the gain; "
      f"persistent {pers:.0f} pp = {100*pers/gain:.1f}%")
for tau in [0.0, 0.2, 0.3, 0.5]:
    print(f"   pre-war trend {tau} pp/yr: formation+persistence residual in 1950 = {pers - 10*tau:.1f} pp")

# ===========================================================================
print("\n== Exercise 11.6: sensitivity of Psi(inf) = k/(1-rho-k), k = mu*beta*eta")
be = 0.5
psi = lambda r, m: m * be / (1 - r - m * be)
for r in [0.50, 0.62, 0.75, 0.90]:
    print(f"(a) rho = {r:.2f}: Psi(inf) = {psi(r, 0.18):.3f}")
for m in [0.10, 0.18, 0.30, 0.50]:
    print(f"(b) mu = {m:.2f}: lambda = rho + mu*be = {0.62 + m*be:.3f}, Psi(inf) = {psi(0.62, m):.3f}")
mu_c = (1 - 0.62) / be
print(f"    unstable for mu >= {mu_c:.2f}")
assert abs(psi(0.90, 0.18) - 9.0) < 1e-9 and abs(psi(0.62, 0.50) - 1.923) < 1e-3
k = 0.18 * be
den = 1 - 0.62 - k
dr = k / den ** 2
dm = be * (1 - 0.62) / den ** 2
print(f"(c) dPsi/drho = {dr:.3f}, dPsi/dmu = {dm:.3f}; elasticities: rho {0.62/den:.3f}, mu {(1-0.62)/den:.3f}")
print("    curve (rho, Psi):", " ".join(f"({r:.2f},{psi(r,0.18):.3f})" for r in np.arange(0.30, 0.8801, 0.02)))
print(f"(d) Psi diverges at rho = 1 - mu*beta*eta = {1-k:.2f}")
for r in [0.50, 0.62, 0.90]:
    lam = r + k
    print(f"    rho = {r}: Psi(T) for T=1,2,5,10 = " +
          ", ".join(f"{k*(1-lam**T)/(1-lam):.3f}" for T in [1, 2, 5, 10]) + f"; Psi(inf) = {psi(r,0.18):.3f}")
    for T in range(1, 30):
        assert k * (1 - lam ** T) / (1 - lam) <= k * T + 1e-12

# ===========================================================================
print("\n== Exercise 11.7: designing an aspiration intervention")
P_ = 2.1
print(f"(a) total long-run effect = (1+Psi) E = {1+P_:.1f} E; E needed for +25 pp = {25/(1+P_):.2f} pp; "
      f"share seen by a fixed-preference evaluation = {1/(1+P_):.3f}")
lam = 0.62
kap = P_ * (1 - lam)
print(f"    one-loop consistency: Psi = kappa/(1-lambda) = 2.1 with lambda = 0.62 needs kappa = {kap:.3f}, rho = {lam-kap:.3f}")
# (c) illustrative concave technologies e(b) = aE sqrt(b), f(b) = aF sqrt(b)
for aE, aF in [(1.0, 1.0), (1.0, 3.0), (1.0, 4.0)]:
    sE = ((1 + P_) * aE) ** 2 / (((1 + P_) * aE) ** 2 + aF ** 2)
    print(f"(c) aE={aE}, aF={aF}: optimal budget share to enforcement = {sE:.3f}")
Tn = np.log(0.1) / np.log(lam)
print(f"(d) generations to 90%: ln(0.1)/ln(0.62) = {Tn:.3f} -> {int(np.ceil(Tn))}; 0.62^5 = {lam**5:.4f}; 0.62^4 = {lam**4:.4f}")
assert int(np.ceil(Tn)) == 5
