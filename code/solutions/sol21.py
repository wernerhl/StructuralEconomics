"""Solutions manual, Chapter 21 (Measuring the Information Network).  Computes and
prints every number used in solutions/sol21.tex, with asserts.

The audience table and population shares are the chapter's illustrative inputs
(Example ch21:ex:bolivia_media).  Sourced values (internet use) come from
data/processed/wdi_long.csv (WDI, retrieved 2026-11-05).
Run from the repo root:  python3 code/solutions/sol21.py   (numpy + pandas).
"""
import itertools
import math
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=150)
PROC = os.path.join("data", "processed")


def perron(W):
    w, v = np.linalg.eig(W.T)
    q = np.real(v[:, np.argmin(abs(w - 1))])
    return q / q.sum()


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def neff(q):
    return 1 / np.sum(q ** 2)


groups = ['Urban NI men', 'Urban NI women', 'Urban indigenous', 'Rural NI', 'Rural indig. men', 'Rural indig. women']
outlets = ['Private TV', 'State', 'Community radio', 'Social media', 'Word of mouth']
S = np.array([[.45, .15, .02, .30, .08],
              [.40, .18, .03, .28, .11],
              [.20, .30, .15, .20, .15],
              [.25, .20, .10, .15, .30],
              [.08, .22, .30, .10, .30],
              [.05, .15, .25, .08, .47]])
p = np.array([0.17, 0.18, 0.25, 0.10, 0.15, 0.15])
assert np.allclose(S.sum(1), 1) and np.isclose(p.sum(), 1)

# =====================================================================
print("=" * 70); print("Exercise 21.1  Media overlap"); print("=" * 70)
O = S @ S.T
print("O = S S':\n", np.round(O, 4))
m = p @ S
Wg = S @ np.diag(1 / m) @ S.T @ np.diag(p)
Wk = np.diag(1 / m) @ S.T @ np.diag(p) @ S
print("W^I_media (population-weighted):\n", np.round(Wg, 4))
assert np.allclose(Wg.sum(1), 1) and np.allclose(Wk.sum(1), 1)
evg, evk = eig_sorted(Wg), eig_sorted(Wk)
print("eigenvalues, group matrix: ", np.round(evg.real, 4), " max |imag| =", f"{abs(evg.imag).max():.1e}")
print("eigenvalues, outlet matrix:", np.round(evk.real, 4))
assert np.allclose(np.sort(evg.real)[1:], np.sort(evk.real), atol=1e-10)   # same nonzero spectrum (+ one extra 0)
Gg, Gk = gap(Wg), gap(Wk)
pg, pk = perron(Wg), perron(Wk)
print(f"Gamma (groups) = {Gg:.4f}; Gamma (outlets) = {Gk:.4f}")
assert np.isclose(Gg, Gk) and np.allclose(pg, p) and np.allclose(pk, m)
print(f"pi (groups) = {np.round(pg,4)} = p; Neff = {neff(pg):.3f} of 6")
print(f"pi (outlets) = reach m = {np.round(m,4)}; Neff = {neff(m):.3f} of 5")
# reversibility (detailed balance) of the group matrix
DB = np.diag(p) @ Wg
assert np.allclose(DB, DB.T)
print("detailed balance p_g W_gg' = p_g' W_g'g holds (reversible)")
Wu = O / O.sum(1, keepdims=True)
pu = perron(Wu)
assert np.allclose(pu, O.sum(1) / O.sum())
print(f"unweighted row-normalisation of O: Gamma = {gap(Wu):.4f}, pi = {np.round(pu,4)} (prop. to row sums of O), Neff = {neff(pu):.3f}")
GAM_U0, NEFF_U0 = gap(Wu), neff(pu)

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.2  Reference-group module: sample size"); print("=" * 70)
z = 1.959964; zb = 0.841621
for p0, h, deff in ((0.30, 0.03, 1.5),):
    n0 = z ** 2 * p0 * (1 - p0) / h ** 2
    print(f"precision: estimate a reference share p = {p0} within +/-{h}: SRS n = {n0:.0f}; deff {deff}: {n0*deff:.0f} per group")
    N_PREC = (n0, n0 * deff)
p1, p2 = 0.20, 0.30
n_test = (z + zb) ** 2 * (p1 * (1 - p1) + p2 * (1 - p2)) / (p2 - p1) ** 2
print(f"test: detect cross-group naming rate {p1} vs {p2} (alpha 0.05 two-sided, power 0.8): n = {n_test:.1f} per group; "
      f"with deff 1.5: {n_test*1.5:.0f}")
assert math.ceil(n_test) == 291
G6 = 6
print(f"six groups x {math.ceil(n_test*1.5)} = {G6*math.ceil(n_test*1.5)} respondents (test design); "
      f"six groups x {N_PREC[1]:.0f} = {G6*N_PREC[1]:.0f} (precision design)")
# construction check on toy responses: alter counts -> row-stochastic W_ref
rng = np.random.default_rng(21)
grp = rng.integers(0, 3, 300); alters = rng.integers(0, 4, (300, 3))       # up to 3 named alters, group codes 0..2, 3 = none
w = rng.uniform(0.5, 2.0, 300)
C = np.zeros((3, 3))
for i in range(300):
    for a in alters[i]:
        if a < 3:
            C[grp[i], a] += w[i]
Wref = C / C.sum(1, keepdims=True)
assert np.allclose(Wref.sum(1), 1)
print("construction check: weighted alter-count matrix row-normalises to a row-stochastic W_ref")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.3  Digital W^I: what N_eff depends on"); print("=" * 70)
# with the population-weighted construction, pi = population base whatever S is
for trial in range(200):
    S_ = rng.dirichlet(np.ones(5), size=6)
    m_ = p @ S_
    W_ = S_ @ np.diag(1 / m_) @ S_.T @ np.diag(p)
    assert np.allclose(perron(W_), p)
print("verified on 200 random audience matrices: group-level pi = population base p, independent of S")
wdi = pd.read_csv(os.path.join(PROC, "wdi_long.csv"))
def wv(c, ind, yr):
    v = wdi[(wdi.country_iso3 == c) & (wdi.indicator == ind) & (wdi.year == yr)]["value"]
    return float(v.iloc[0]) if len(v) else float("nan")
for yr in (2019, 2024):
    print(f"internet users, Bolivia {yr}: {wv('BOL','IT.NET.USER.ZS',yr):.1f}%")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.4  QAP mechanics with six groups"); print("=" * 70)
n_perm = math.factorial(6); n_dyads = 6 * 5 // 2
print(f"label permutations = {n_perm}; dyads = {n_dyads}; smallest attainable one-sided permutation p-value = 1/{n_perm} = {1/n_perm:.4f}")
offd = O[np.triu_indices(6, 1)]
# the number of distinct permuted overlap vectors (automorphisms of O would reduce it)
distinct = len({tuple(np.round(O[np.ix_(pm, pm)][np.triu_indices(6, 1)], 12)) for pm in itertools.permutations(range(6))})
print(f"distinct permuted dyad vectors of the illustrative O: {distinct}")
assert distinct == n_perm

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.5  Assortative mixing parameters"); print("=" * 70)
NI, IND = [0, 1, 3], [2, 4, 5]
URB, RUR = [0, 1, 2], [3, 4, 5]


def lam_blocks(A, B, diag=True):
    cross = O[np.ix_(A, B)].mean()
    def within(X):
        sub = O[np.ix_(X, X)]
        return sub.mean() if diag else sub[~np.eye(len(X), dtype=bool)].mean()
    return cross / np.mean([within(A), within(B)])


LAM = {"gender, urban NI": O[0, 1] / np.mean([O[0, 0], O[1, 1]]),
       "gender, rural indigenous": O[4, 5] / np.mean([O[4, 4], O[5, 5]]),
       "ethnicity": lam_blocks(NI, IND), "urban/rural": lam_blocks(URB, RUR),
       "ethnicity, off-diagonal within": lam_blocks(NI, IND, diag=False),
       "urban/rural, off-diagonal within": lam_blocks(URB, RUR, diag=False)}
for k, v in LAM.items():
    print(f"  lambda {k:34s} = {v:.3f}")
assert round(LAM["ethnicity"], 2) == 0.72 and round(LAM["urban/rural"], 2) == 0.68
# conditional (partial) lambdas: log O_gg' = c + d_eth log l_r + d_urb log l_u + d_gen log l_g on all 21 cells
eth = np.array([0, 0, 1, 0, 1, 1]); urb = np.array([0, 0, 0, 1, 1, 1]); sex = np.array([0, 1, -1, -1, 0, 1])  # -1 = mixed group
rows_, ys_ = [], []
for i in range(6):
    for j in range(i, 6):
        dg = 1.0 if (sex[i] >= 0 and sex[j] >= 0 and sex[i] != sex[j] and eth[i] == eth[j] and urb[i] == urb[j]) else 0.0
        rows_.append([1, float(eth[i] != eth[j]), float(urb[i] != urb[j]), dg]); ys_.append(np.log(O[i, j]))
X = np.array(rows_); yv = np.array(ys_)
beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
res = yv - X @ beta
r2 = 1 - res.var() / yv.var()
PART = {"ethnicity": np.exp(beta[1]), "urban/rural": np.exp(beta[2]), "gender": np.exp(beta[3])}
print(f"partial lambdas from a log-linear fit on the 21 cells (R^2 = {r2:.3f}):", {k: round(v, 3) for k, v in PART.items()})
both = [(i, j) for i in range(6) for j in range(6) if eth[i] != eth[j] and urb[i] != urb[j]]
print(f"pairs differing in both ethnicity and geography: {len(both)//2}; implied joint factor {PART['ethnicity']*PART['urban/rural']:.3f}")

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.6  The operator from data"); print("=" * 70)
iu = np.triu_indices(6, 1)
k_raw = int(np.argmin(O[iu])); i0, j0 = iu[0][k_raw], iu[1][k_raw]
rel = O / np.sqrt(np.outer(np.diag(O), np.diag(O)))
k_rel = int(np.argmin(rel[iu])); i1, j1 = iu[0][k_rel], iu[1][k_rel]
print(f"darkest raw overlap: O[{groups[i0]}, {groups[j0]}] = {O[i0,j0]:.4f}; "
      f"relative to the mean of self-overlaps: {O[i0,j0]/np.mean([O[i0,i0],O[j0,j0]]):.3f}")
print(f"darkest normalised overlap O_gg'/sqrt(O_gg O_g'g'): {groups[i1]} - {groups[j1]} = {rel[i1,j1]:.3f}")
assert (i0, j0) == (0, 5) and (i1, j1) == (0, 5)
O2 = O.copy(); O2[0, 5] += 0.10; O2[5, 0] += 0.10
Wu2 = O2 / O2.sum(1, keepdims=True)
pu2 = perron(Wu2)
G0, G1 = gap(Wu), gap(Wu2)
print(f"Gamma: {G0:.4f} -> {G1:.4f} ({G1-G0:+.4f}); Neff: {neff(pu):.4f} -> {neff(pu2):.4f} ({neff(pu2)-neff(pu):+.4f})")
print(f"pi before {np.round(pu,4)}; after {np.round(pu2,4)}")
assert round(G1, 3) == 0.791 and abs(neff(pu2) - neff(pu)) < 0.02
# first-order prediction of d lambda_2 (perturbation lemma, eigenvalue part)
wr, vr = np.linalg.eig(Wu); wl, vl = np.linalg.eig(Wu.T)
k2 = np.argsort(-abs(wr))[1]; l2 = wr[k2].real
v2 = vr[:, k2].real; u2 = vl[:, np.argmin(abs(wl - wr[k2]))].real
u2 = u2 / (u2 @ v2)
dW = Wu2 - Wu
d1 = u2 @ dW @ v2
lam2_after = eig_sorted(Wu2)[1].real
print(f"lambda_2 = {l2:.4f}; first-order change u'Ev = {d1:+.4f}; exact {lam2_after - l2:+.4f}")
print("lambda_2 right eigenvector (sign-normalised):", np.round(v2 / np.sign(v2[0]) / abs(v2).max(), 3))
# (c) formation: two-block (NI vs indigenous) quotient, pi-weighted, before and after
def quotient(Wx, px, blocks):
    return np.array([[px[A] @ Wx[np.ix_(A, B)].sum(1) / px[A].sum() for B in blocks] for A in blocks])
QB = quotient(Wu, pu, [NI, IND]); QA = quotient(Wu2, pu2, [NI, IND])
l2B, l2A = QB.trace() - 1, QA.trace() - 1
print(f"ethnic quotient lambda_2 (aspiration-gap factor): {l2B:.4f} -> {l2A:.4f}; ratio {l2A/l2B:.4f} "
      f"(gap falls by {100*(1-l2A/l2B):.2f}%)")
lamr_after = O2[np.ix_(NI, IND)].mean() / np.mean([O2[np.ix_(NI, NI)].mean(), O2[np.ix_(IND, IND)].mean()])
fB, fA = (1 - LAM["ethnicity"]) / (1 + LAM["ethnicity"]), (1 - lamr_after) / (1 + lamr_after)
print(f"lambda_ethnicity from O: {LAM['ethnicity']:.4f} -> {lamr_after:.4f}; (1-l)/(1+l): {fB:.4f} -> {fA:.4f} "
      f"(gap falls by {100*(1-fA/fB):.2f}%)")
Psi = 2.0
share_form = Psi / (1 + Psi)
print(f"with Psi = 2: formation share of the outcome gap = {share_form:.4f}; outcome gap falls by "
      f"{100*share_form*(1-l2A/l2B):.2f}% (quotient) or {100*share_form*(1-fA/fB):.2f}% (lambda from O)")
QUOT = (l2B, l2A, fB, fA, lamr_after)

# =====================================================================
print("\n" + "=" * 70); print("Exercise 21.7  Cross-country: sourced inputs available"); print("=" * 70)
for c in ("BOL", "CHL"):
    for yr in (2019, 2024):
        print(f"internet users {c} {yr}: {wv(c,'IT.NET.USER.ZS',yr):.1f}%")
assert "SWE" not in set(wdi.country_iso3)
print("Sweden: not in the WDI extract (wdi_long.csv)")
INET = {(c, yr): wv(c, "IT.NET.USER.ZS", yr) for c in ("BOL", "CHL") for yr in (2019, 2024)}
print(f"Neff/n bounds: 1/n <= Neff/n <= 1; Bolivia illustrative outlets Neff/n = {neff(m)/5:.3f}")
print(f"rank agreement test with 3 countries: P(perfect agreement by chance) = 1/3! = {1/math.factorial(3):.4f}")
print("\nall Chapter 21 solution numbers computed; asserts passed")
