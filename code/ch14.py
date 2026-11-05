"""Chapter 14: reproduces every number in the chapter's examples, figure and
revised exercises, and checks the identities used in the text.
Run from the repo root:  python3 code/ch14.py   (numpy only)."""
import numpy as np
import itertools

np.set_printoptions(precision=3, suppress=True)


def perron(W):
    ev, vl = np.linalg.eig(W.T)
    k = np.argmin(np.abs(ev - 1))
    p = np.real(vl[:, k])
    return p / p.sum()


print("== Example: four alliances (illustrative Delta pi vectors) ==")
pi0 = np.array([0.35, 0.30, 0.20, 0.15])
alliances = {
    "1 Oligarchic-export {F,B}": np.array([+0.05, +0.03, -0.04, -0.04]),
    "2 National-populist {B,M,W}": np.array([-0.10, +0.02, +0.05, +0.03]),
    "3 Techno-bureaucratic {F,B,M}": np.array([+0.03, +0.04, +0.02, -0.09]),
    "4 Left-populist {M,W}": np.array([-0.10, -0.02, +0.05, +0.07]),
}
for k, d in alliances.items():
    assert abs(d.sum()) < 1e-12          # power is constant-sum: sum Delta pi = 0
    assert np.all(pi0 + d > 0)
    print(f"{k}: dpi = {d}, sum = {d.sum():+.2f}, new pi = {pi0 + d}")
# B's ranking at pi_B = 0.30 (mature ISI): alliance 3 > alliance 1 > alliance 2
dB = {k: v[1] for k, v in alliances.items()}
print("Delta pi_B by alliance:", dB)
assert dB["3 Techno-bureaucratic {F,B,M}"] > dB["1 Oligarchic-export {F,B}"] > dB["2 National-populist {B,M,W}"]

print("\n== Exercise 14.1 (Omega = 4, bargaining surplus Phi = 2) ==")
Omega, Phi, N = 4.0, 2.0, 4
x = Omega / N + (pi0 - 1 / N) * Phi
print("x* =", x, " sum =", x.sum())
kappa = 0.5
assert np.all(kappa * x < 1)             # below satiation for U = x - kappa/2 x^2
pi_new = np.full(4, 0.25)
v = Phi * (pi_new - pi0)                  # Delta Omega = 0, Psi = 0
print("v_i for move to uniform pi:", v, "-> v_3, v_4 =", v[2], v[3])
print("agent 2: Delta pi under {1,2} (alliance 1) =", alliances["1 Oligarchic-export {F,B}"][1],
      "; under {2,3,4} (alliance 2) =", alliances["2 National-populist {B,M,W}"][1])

print("\n== Empty core when the surplus term is zero ==")
# TU game v(S) = max over feasible Delta W of sum_{i in S} Phi*Delta pi_i, with Delta Omega = 0.
# Toy feasibility: coalition S can implement any alliance vector whose gainers are exactly S's members.
def v_of(S):
    best = 0.0                              # Delta W = 0 is always feasible
    for d in alliances.values():
        if all((d[i] > 0) <= (i in S) for i in range(4)):   # only members gain
            best = max(best, Phi * sum(d[i] for i in S))
    return best
vN = v_of(set(range(4)))
print(f"v(N) = {abs(vN):.3f}")
positive = [S for r in range(1, 4) for S in itertools.combinations(range(4), r) if v_of(set(S)) > 0]
print("coalitions with v(S) > 0:", positive)
assert abs(vN) < 1e-12 and len(positive) > 0   # => x(N) = 0 but x(S) >= v(S) > 0, x(N\S) >= 0: core empty

print("\n== Appendix 14B: perturbation lemma check ==")
rng = np.random.default_rng(14)
W = rng.uniform(0.1, 1, (4, 4)); W /= W.sum(1)[:, None]
E = rng.normal(size=(4, 4)); E -= E.mean(1)[:, None]          # E 1 = 0
p = perron(W)
Z = np.linalg.inv(np.eye(4) - W + np.outer(np.ones(4), p))
fd = (perron(W + 1e-7 * E) - p) / 1e-7
lam, V = np.linalg.eig(W); o = np.argsort(-np.abs(lam)); lam = lam[o]; V = V[:, o]; U = np.linalg.inv(V)
old = np.real(sum((U[k] @ E @ np.ones(4)) / (1 - lam[k]) * U[k] for k in range(1, 4)))
print("finite difference  :", fd)
print("pi' E Z            :", p @ E @ Z)
print("old formula (=0)   :", old)
assert np.allclose(fd, p @ E @ Z, atol=1e-5) and np.allclose(old, 0, atol=1e-12)

print("\n== Figure: alliance regions in the (pi_F, pi_W) plane ==")
def region(pF, pW):
    if pF >= 0.40:
        return "{F,B}"
    if pW >= 0.15:
        return "{M,W}"
    if pF >= 0.15:
        return "{F,B,M}"
    return "{B,M,W}"
for name, (pF, pW) in {"baseline": (0.35, 0.15), "Bolivia 1985": (0.30, 0.12), "Bolivia 2006": (0.15, 0.20)}.items():
    print(f"{name} (pi_F={pF}, pi_W={pW}) -> {region(pF, pW)}")
assert region(0.30, 0.12) == "{F,B,M}" and region(0.15, 0.20) == "{M,W}"

# ---------------------------------------------------------------------------
print("\n== Sourced: Bolivia's domestic Gamma(W^P) across the nationalization phase ==")
# INE input-output tables, base 1990, book12, domestic flows, input shares (data/processed, csv module)
import csv
import os

_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "processed", "bolivia_io_summary.csv")
with open(_p, newline="", encoding="utf-8") as _f:
    _g = {r["year"]: float(r["Gamma"]) for r in csv.DictReader(_f)
          if r["base"] == "1990" and r["aggregation"] == "book12" and r["flows"] == "domestic"
          and r["normalization"] == "input_shares"}
for y, text in [("2006", 0.518), ("2013", 0.484)]:
    print(f"Gamma {y}: {_g[y]:.4f} -> text {text}")
    assert abs(_g[y] - text) <= 0.0005 + 1e-9
assert _g["2013"] < _g["2006"]      # Phase 3 did not raise domestic Gamma in Bolivia
