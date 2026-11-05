"""Chapter 11: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch11.py
numpy only, plus the csv module for the two sourced national figures (WDI, via
data/processed/bolivia_facts.csv).  All other empirical magnitudes in Chapter 11
are illustrative."""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def pvar(v):
    return np.mean((np.asarray(v) - np.mean(v)) ** 2)


# ---------------------------------------------------------------------------
# Exercise 11.1 and the APC sensitivity (Proposition ch11:prop:apc)
X = np.array([[0.38, 0.42, 0.44],
              [0.50, 0.55, 0.58],
              [0.60, 0.66, 0.72]])          # rows: cohorts 1940s,50s,60s; cols: 1980,1990,2000
years = np.array([1980, 1990, 2000])
cohorts = np.array([1945, 1955, 1965])
gm = X.mean()
yr = X.mean(0) - gm
co = X.mean(1) - gm
print("== Exercise 11.1")
print("year means", X.mean(0), " cohort effects", co)
print("Psi_var = Var(cohort)/Var(year) =", pvar(co) / pvar(yr))
print("ages at observation:\n", years[None, :] - cohorts[:, None])
# observationally equivalent APC decompositions: delta_c + k c, alpha_t - k t, A_a + k a
cidx = np.arange(3.0)
tidx = np.arange(3.0)
print("APC sensitivity (k per decade): ")
for k in [-0.10, -0.05, 0.0, 0.05]:
    co_k = co + k * cidx
    yr_k = yr - k * tidx
    # fitted values identical once age effects A_a = k*(t - c) are included
    fit0 = gm + yr[None, :] + co[:, None]
    fitk = gm + yr_k[None, :] + co_k[:, None] + k * (tidx[None, :] - cidx[:, None])
    assert np.allclose(fit0, fitk)
    print(f"  k={k:+.2f}: Psi_var = {pvar(co_k)/pvar(yr_k):.3f}")

# ---------------------------------------------------------------------------
# Example ch11:ex:psi_lfp (illustrative variance components)
print("\n== Example psi_lfp (illustrative)")
vt, vy, vc = 0.042, 0.013, 0.027
print("Psi_var =", vc / vy, " sqrt =", np.sqrt(vc / vy), " residual share", (vt - vy - vc) / vt)

# Proposition ch11:prop:psivar_map: one-index model -> Psi_var = Psi^2
rng = np.random.default_rng(1)
z = np.cumsum(rng.normal(size=400))
e, f = 0.7, 1.3
a_t = e * z[100:300]
d_c = f * z[100:300]            # same index over a window of the same length & variance
print("one-index check: Psi_var =", pvar(d_c) / pvar(a_t), " (f/e)^2 =", (f / e) ** 2)

# ---------------------------------------------------------------------------
# Example ch11:ex:epi (illustrative generational coefficients)
b2, b3 = 0.40, 0.25
rho_epi = b3 / b2
lam_rel = b2 / rho_epi
print("\n== Example epi: rho = b3/b2 =", rho_epi, " reliability lambda = b2/rho =", lam_rel)
# simulate the derivation: LFP_A = b th_A + e_A, daughters th_d = rho th_A + ..., regression coef = rho*lambda
n = 200000
thA = rng.normal(size=n)
eA = rng.normal(scale=np.sqrt(1 / lam_rel - 1), size=n)    # reliability = Var(th)/(Var(th)+Var(e))
lfpA = thA + eA
thd = rho_epi * thA + rng.normal(scale=0.3, size=n)
lfpd = thd + rng.normal(scale=0.3, size=n)
coef = np.cov(lfpd, lfpA)[0, 1] / np.var(lfpA)
print("   simulated regression coefficient", round(coef, 3), " vs rho*lambda =", rho_epi * lam_rel)
assert abs(coef - rho_epi * lam_rel) < 0.01

# ---------------------------------------------------------------------------
# Example ch11:ex:rho_mu (illustrative)
c1, c2 = 0.35, 0.217
rho_hat = c2 / c1
R = c1 / rho_hat
print("\n== Example rho_mu: rho =", rho_hat, " R_theta =", R, " half-life", np.log(2) / np.log(1 / 0.62))
mu_hat = 0.18
print("   M =", mu_hat / (1 - 0.62), " per-decade 0.62 -> per generation", 0.62 ** 2.5)
print("   cohort gradient (0.65-0.30)/30 =", (0.65 - 0.30) / 30)

# ---------------------------------------------------------------------------
# Network identification (Section ch11:sec:network): BDF rank condition
blk = np.repeat([0, 1, 2], 3)
Wpure = np.where(blk[:, None] == blk[None, :], 1 / 3, 0.0)
Wst = np.where(blk[:, None] == blk[None, :], 0.85 / 3, 0.025)
for nm, W in [("pure block", Wpure), ("Ch10 W^I", Wst)]:
    M = np.stack([np.eye(9).ravel(), W.ravel(), (W @ W).ravel()])
    print(f"{nm}: rank[I, W, W^2] = {np.linalg.matrix_rank(M)}, eigenvalues {np.round(np.unique(np.round(np.linalg.eigvals(W).real, 6)), 4)}")
assert np.allclose(Wpure @ Wpure, Wpure)

# ---------------------------------------------------------------------------
# Proposition ch11:prop:bounds (illustrative sensitivity)
P0 = 2.08
cs_lo, cs_hi = P0 * 0.9, P0 * 1.1
sel = P0 * 0.8
gen = P0 * 0.6
comb = P0 * 0.9 * 0.8 * 0.6
print(f"\n== Bounds: correlated shocks [{cs_lo:.2f},{cs_hi:.2f}], selection {sel:.2f}, gene {gen:.2f}, combined [{comb:.2f},{cs_hi:.2f}]")

# Remark policy shares
for p in [0.5, 1.0]:
    print("share missed when ignoring formation, Psi =", p, ":", p / (1 + p))

# ---------------------------------------------------------------------------
print("\n== Exercises")
# 11.2: ten evenly spaced origin rates 0.20..0.75, beta = 0.35, residual var 0.04
o = np.linspace(0.20, 0.75, 10)
expl = 0.35 ** 2 * pvar(o)
print("11.2 Var(origin) =", pvar(o), " explained =", expl, " share =", expl / (expl + 0.04))
# 11.3 (illustrative): adoptee vs biological-child gradients
print("11.3 ratio adoptee/biological gradient =", 0.10 / 0.30)
# 11.6: one-loop Psi(inf) = k/(1-rho-k), k = mu * beta*eta, beta*eta = 0.5
for r in [0.50, 0.62, 0.75, 0.90]:
    k = 0.18 * 0.5
    print(f"11.6 rho={r}: Psi(inf) = {k/(1-r-k):.3f}")
for m in [0.10, 0.18, 0.30, 0.50]:
    k = m * 0.5
    val = k / (1 - 0.62 - k) if 1 - 0.62 - k > 0 else float('nan')
    print(f"11.6 mu={m}: Psi(inf) = {val:.3f}" + ("" if 1 - 0.62 - k > 0 else "  (unstable: rho+k>=1)"))
print("11.6 divergence at rho = 1 - 0.09 =", 0.91)
# 11.7(d)
print("11.7(d) generations to 90% with lambda=0.62:", np.log(0.1) / np.log(0.62))
print("11.6(b) all four mu values stable; instability when mu*0.5 >= 1-0.62, i.e. mu >=", (1 - 0.62) / 0.5)
assert all(0.62 + m * 0.5 < 1 for m in [0.10, 0.18, 0.30, 0.50])

# ---------------------------------------------------------------------------
# Proof checks (proof pass)
print("\n== Proof checks")
rng = np.random.default_rng(11)
# Prop ch11:prop:challenge: OVB formula
n = 400000
th_ = rng.normal(size=n)
p_ = 0.6 * th_ + rng.normal(size=n)
x_ = 1.5 * p_ + 2.0 * th_ + rng.normal(size=n)
bhat = np.cov(x_, p_)[0, 1] / np.var(p_, ddof=1)
assert abs(bhat - (1.5 + 2.0 * 0.6 / 1.36)) < 0.01   # population: Cov(theta,p)=0.6, Var(p)=1.36
print("OVB formula verified: plim beta_hat =", round(bhat, 3))
# Prop ch11:prop:apc: ratio near 0 at k_delta when delta nearly linear; large at k_alpha when alpha nearly linear
cc = np.arange(8.0)
tt = np.arange(8.0)
d_lin = 0.3 * cc + 0.01 * rng.normal(size=8)
a_wig = rng.normal(size=8)
k_d = -np.cov(d_lin, cc, bias=True)[0, 1] / pvar(cc)
ratio = lambda d, a, k: pvar(d + k * cc) / pvar(a - k * tt)
assert ratio(d_lin, a_wig, k_d) < 1e-3
a_lin = -0.2 * tt + 0.01 * rng.normal(size=8)
d_wig = rng.normal(size=8)
k_a = np.cov(a_lin, tt, bias=True)[0, 1] / pvar(tt)
assert ratio(d_wig, a_lin, k_a) > 1e3
assert abs(ratio(d_wig, a_wig, 1e6) - pvar(cc) / pvar(tt)) < 1e-6
print("APC ratio extremes verified:", ratio(d_lin, a_wig, k_d), ratio(d_wig, a_lin, k_a))
# Prop ch11:prop:id_params: Corr(x_g, x_{g-k}) = rho^k R_theta for AR(1) theta + iid noise
G, rr, bb = 400000, 0.62, 1.0
th = np.zeros(G)
eps = rng.normal(size=G)
for g in range(1, G):
    th[g] = rr * th[g - 1] + 0.3 * eps[g]
xg = 0.2 + bb * th + rng.normal(scale=0.4, size=G)
xg, th = xg[1000:], th[1000:]
Rth = bb ** 2 * np.var(th) / np.var(xg)
for k in (1, 2):
    c = np.corrcoef(xg[k:], xg[:-k])[0, 1]
    assert abs(c - rr ** k * Rth) < 0.01
print("id_params correlations verified: R_theta =", round(Rth, 3))

# ---------------------------------------------------------------------------
# Sourced national figures quoted in Example ch11:prop:bolivia and Exercise 11.7
import csv
import os

_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "processed", "bolivia_facts.csv")
with open(_path, newline="", encoding="utf-8") as _f:
    facts = {r["fact_id"]: float(r["value"]) for r in csv.DictReader(_f)}
for k, text in [("BOL_self_employed_2019", 68.1), ("BOL_female_lfp_2019", 66.2)]:
    print(f"{k}: {facts[k]} -> text {text}")
    assert abs(facts[k] - text) <= 0.05 + 1e-9
assert facts["BOL_self_employed_2019"] > 65      # above the top of the hypothetical departmental range
