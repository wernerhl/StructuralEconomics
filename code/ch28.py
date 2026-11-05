"""Chapter 28: checks behind the Framing Principle and the Bolivia table.
Run from the repo root: python3 code/ch28.py  (numpy only)."""
import numpy as np


def stats(W):
    ev, vl = np.linalg.eig(W.T)
    o = np.argsort(-np.abs(ev))
    pi = np.real(vl[:, o[0]]); pi /= pi.sum()
    return 1 - abs(ev[o[1]]), pi, 1 / np.sum(pi ** 2)


N = 5
s = np.array([0.4, 0.25, 0.15, 0.12, 0.08])          # wealth shares
G, pi, ne = stats(np.outer(np.ones(N), s))
print(f"W_mkt = 1 s': Gamma = {G:.3f}, pi = s: {np.allclose(pi, s)}, Neff = {ne:.2f}")
G, pi, ne = stats(np.outer(np.ones(N), np.eye(N)[0]))
print(f"dictatorship 1 e1': Gamma = {G:.3f}, Neff = {ne:.2f}  -> Gamma = 1 does not identify the benchmark")
# Gamma and pi independent: W = g 1 pi' + (1-g) I
p = np.array([0.6, 0.1, 0.1, 0.1, 0.1])
for g in (0.1, 0.5, 1.0):
    G, pi, ne = stats(g * np.outer(np.ones(N), p) + (1 - g) * np.eye(N))
    print(f"gamma = {g}: Gamma = {G:.3f}, pi = {np.round(pi, 3)}")
    assert abs(G - g) < 1e-9 and np.allclose(pi, p)
# doubly stochastic but not W_CE: uniform pi, Delta != 0
C = np.array([[0.6, 0.3, 0.1], [0.1, 0.6, 0.3], [0.3, 0.1, 0.6]])
G, pi, ne = stats(C)
print(f"circulant: Gamma = {G:.3f}, pi uniform = {np.allclose(pi, 1/3)}, ||Delta||_F = {np.linalg.norm(C - 1/3):.3f}")

print("\nBolivia table (illustrative, from Chapters 23 and 25):")
print("  Delta Gamma 1985-2019 =", round(0.13 - 0.11, 2), "; SCI 1990-2019 =", round((0.13 - 0.10) / (0.50 - 0.10), 3))
f = lambda x: (1 - x) / x
shares = [(f(0.13) - f(g0)) / f(0.13) for g0 in (0.25, 0.35, 0.52)]
print(f"  truncation share of integration deficit: {shares[0]:.2f} - {shares[2]:.2f} (at 0.35: {shares[1]:.2f})")
