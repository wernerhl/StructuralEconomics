"""Chapter 11: reproduces every number used in the worked examples, tables,
viz panels and exercise answers.  Run from the repo root:  python3 code/ch11.py
numpy only.  All empirical magnitudes in Chapter 11 are illustrative."""
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
