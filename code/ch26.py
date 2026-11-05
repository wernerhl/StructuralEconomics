"""Chapter 26: counterexamples to the old proof, viz-panel classification,
Korea check, and plan arithmetic. Run: python3 code/ch26.py (numpy only)."""
import numpy as np


def stats(W):
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi /= pi.sum()
    N = len(W)
    var = np.mean((pi - 1 / N) ** 2)
    return 1 - abs(ev[o[1]]), pi, var, 1 / np.sum(pi ** 2)


print("== Gamma does not bound Var(pi) (Ch6 Prop gap-power-independent) ==")
N = 4
pi = np.array([0.97, 0.01, 0.01, 0.01])
W1 = np.outer(np.ones(N), pi)
G, p, v, ne = stats(W1)
print(f"rank-one 1pi': Gamma = {G:.3f}, Var(pi) = {v:.4f} (max (N-1)/N^2 = {(N-1)/N**2:.4f}), Neff = {ne:.2f}")
assert abs(G - 1) < 1e-9 and v <= (N - 1) / N ** 2
e = 0.01
B = np.full((2, 2), 0.5)
W2 = np.block([[(1 - e) * B, e * B], [e * B, (1 - e) * B]])
G, p, v, ne = stats(W2)
print(f"doubly stochastic weak blocks: Gamma = {G:.3f}, Var(pi) = {v:.2e}, Neff = {ne:.2f}")
assert abs(G - 0.02) < 1e-9 and v < 1e-20

print("\n== Viz panel (illustrative): x = 10*Gamma, y = ECI; zone Gamma>0.30 & ECI>0 ==")
pts = {"BOL": (0.13, -0.8), "BRA": (0.25, 0.8), "COL": (0.20, -0.3), "MEX": (0.28, 1.0),
       "KOR": (0.45, 2.1), "JPN": (0.56, 2.3), "DEU": (0.52, 1.8), "CZE": (0.40, 1.5),
       "CHL": (0.28, 0.3), "ARG": (0.15, 0.5), "VEN": (0.10, -0.5), "MYS": (0.42, 0.8),
       "POL": (0.38, 1.2), "CHN": (0.35, 1.0)}
inside = sorted(k for k, (g, c) in pts.items() if g > 0.30 and c > 0)
print("inside the zone:", inside)
assert "CHL" not in inside and "CHN" in inside
for k, (g, c) in pts.items():
    print(f"  {k}: ({10*g:.1f}, {c:+.1f})")

print("\n== Korea on Chapter 23's illustrative path, Gamma(t) = 0.15 + 0.0075 (t - 1960) ==")
for t in (1965, 1970, 1975, 1980, 1985):
    g = 0.15 + 0.0075 * (t - 1960)
    print(f"  {t}: Gamma = {g:.3f}  {'< Gamma* = 0.30' if g < 0.30 else '>= Gamma*'}")

print("\n== Plan arithmetic ==")
costs = {"broadband": 1.5, "polytechnics": 2, "corridors": 5, "lithium": 3, "agri": 1, "industrial": 4}
tot = sum(costs.values())
print(f"costed items = ${tot}B (rounded to $17B); per year over 30 = ${17/30*1000:.0f}M; share of $40B GDP = {17/30/40*100:.2f}%")
print("Neff after Phase 1: 3.4 + 2.0 (broadband) =", 3.4 + 2.0, "-> text ~5.5 with community journalism")
