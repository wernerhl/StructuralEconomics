"""Chapter 19: reproduces every number used in chapters/chapter19.tex and checks
the lumpability, stochastic-complementation and nested-block results.  numpy only.
Run from the repo root:  python3 code/ch19.py
"""
import contextlib
import io
import os
import runpy

import numpy as np

rng = np.random.default_rng(19)


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def gap(W):
    m = np.sort(abs(np.linalg.eigvals(W)))[::-1]
    return 1 - m[1]


def rownorm(K):
    return K / K.sum(1, keepdims=True)


print("=== Prop. lumpability: Gamma_coarse >= Gamma_fine, pi_bar = block sums ===")
blocks = [0, 0, 0, 1, 1, 2, 2, 2, 2]          # 9 fine nodes, 3 coarse blocks
k, n = 3, len(blocks)
V = np.zeros((n, k)); V[np.arange(n), blocks] = 1
for trial in range(200):
    Wbar = rownorm(rng.random((k, k)))
    W = np.zeros((n, n))
    for i in range(n):                         # each row: split block mass arbitrarily
        for B in range(k):
            idx = [j for j in range(n) if blocks[j] == B]
            w = rng.random(len(idx)); W[i, idx] = Wbar[blocks[i], B] * w / w.sum()
    assert np.allclose(W @ V, V @ Wbar)        # exact lumpability
    assert gap(Wbar) >= gap(W) - 1e-12
    assert np.allclose(perron(W) @ V, perron(Wbar))
print("holds in 200 random exactly-lumpable matrices")
W_ex = W; print(f"last example: Gamma_fine = {gap(W):.3f}, Gamma_coarse = {gap(Wbar):.3f}")

print("\n=== Prop. stochastic complementation (exact nesting identity) ===")
def stoch_complement(W, idx):
    rest = [j for j in range(W.shape[0]) if j not in idx]
    A = W[np.ix_(idx, idx)]; B = W[np.ix_(idx, rest)]
    C = W[np.ix_(rest, idx)]; D = W[np.ix_(rest, rest)]
    return A + B @ np.linalg.solve(np.eye(len(rest)) - D, C)

for trial in range(100):
    n = 10; part = [list(range(0, 4)), list(range(4, 7)), list(range(7, 10))]
    W = rownorm(rng.random((n, n)) ** 3)
    pi = perron(W)
    s = [perron(stoch_complement(W, b)) for b in part]
    Cc = np.array([[s[a] @ W[np.ix_(part[a], part[b])].sum(1) for b in range(3)] for a in range(3)])
    xi = perron(Cc)
    recon = np.concatenate([xi[a] * s[a] for a in range(3)])
    assert np.allclose(recon, pi)
print("pi_i = xi_B * s^(B)_i holds exactly in 100 random matrices (Meyer 1989)")
# Simon-Ando: weak coupling -> s^(B) ~ Perron vector of the diagonal block
n = 10; part = [list(range(0, 4)), list(range(4, 7)), list(range(7, 10))]
W0 = np.zeros((n, n))
for b in part:
    W0[np.ix_(b, b)] = rownorm(rng.random((len(b), len(b))))
for eps in [0.1, 0.01, 0.001]:
    W = (1 - eps) * W0 + eps * rownorm(rng.random((n, n)))
    approx = np.concatenate([perron(W0[np.ix_(b, b)]) for b in part])
    s_true = np.concatenate([perron(stoch_complement(W, b)) for b in part])
    print(f"eps={eps}: max |s^(B) - pi(W_BB)| = {abs(s_true-approx).max():.2e}")

print("\n=== Example: per-capita nesting (illustrative inputs) ===")
piG_bol, piG_deu = 0.0005, 0.08                # Chapter 16 illustrative values
L_bol, L_deu = 6e6, 45e6                       # illustrative employed workforces
# Chapter 17 illustrative sector values (Example ch17:ex:bolivia_pi), imported from
# code/ch17.py so the two chapters cannot drift apart.  The Potosi vendor works in
# "personal services" (traditional/informal block); household production is
# "domestic services".
with contextlib.redirect_stdout(io.StringIO()):
    _c17 = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch17.py"))
_sec = list(_c17["sectors"])
pi17, emp17 = _c17["p"], _c17["emp"]
i_gas, i_inf = _sec.index("Hydrocarbons"), _sec.index("Personal services")
xi_gas, e_gas = pi17[i_gas], emp17[i_gas]     # 0.195 at 1% of employment
xi_inf, e_inf = pi17[i_inf], emp17[i_inf]     # 0.012 at 9% of employment
assert round(xi_gas, 3) == 0.195 and e_gas == 0.01 and round(xi_inf, 3) == 0.012 and e_inf == 0.09
print(f"Ch17 inputs: hydrocarbons pi = {xi_gas:.3f} (emp {e_gas:.2f}); "
      f"personal services pi = {xi_inf:.4f} (emp {e_inf:.2f})")
w_vendor = 0.57 * 0.58                         # Chapter 18: per-capita factor vs non-indigenous man
sector_gas, sector_inf = xi_gas / e_gas, xi_inf / e_inf
print(f"sector per-capita factors: gas {sector_gas:.2f}, personal services {sector_inf:.3f}, "
      f"ratio {sector_gas/sector_inf:.1f}")
print(f"within factor: indigenous woman / non-indigenous man = {w_vendor:.3f} (rounded 0.33)")
r_exec = sector_gas / sector_inf / 0.33
print(f"executive / vendor = {sector_gas/sector_inf:.1f} x {1/0.33:.2f} = {r_exec:.0f}")
country = (piG_deu / L_deu) / (piG_bol / L_bol)
print(f"country per-capita factor Germany/Bolivia = 160 x {L_bol/L_deu:.3f} = {country:.1f}")
r_deu = country * (1.0 / sector_inf) / 0.33
print(f"German engineer (sector factor 1, within 1) / vendor = {country:.1f} x {1/sector_inf:.2f} x {1/0.33:.2f} = {r_deu:.0f}")
print(f"  without dividing by workforce: 160 x {1/sector_inf:.2f} x {1/0.33:.2f} = {160/sector_inf/0.33:.0f}")
print(f"vendor / executive per-capita power = {1/r_exec:.4f} ({100/r_exec:.2f}%)")
print("old text: 0.0112/1.15e-5 =", round(0.0112 / 1.15e-5, 1), "; 160*2.8/0.46 =", round(160 * 2.8 / 0.46, 1),
      "(the 986 vs 974 gap was rounding 1/0.46 = 2.17 to 2.2)")
i_max, i_min = int(np.argmax(pi17)), int(np.argmin(pi17))
print(f"national max/min pi on Ch17 values: {_sec[i_max]}/{_sec[i_min]} = "
      f"{pi17[i_max]:.3f}/{pi17[i_min]:.3f} = {pi17[i_max]/pi17[i_min]:.0f}")
print(f"domestic services (household production) pi = {pi17[_sec.index('Domestic services')]:.3f}")
print(f"Gamma_N/Gamma_G = 0.13/0.20 = {0.13/0.20:.2f} (35% lower)")

print("\n=== Appendix 19B: hierarchical block model ===")
k, m, nn = 3, 3, 4; w2, w1, w0 = 1.0, 0.3, 0.05
N = k * m * nn; K = np.full((N, N), w0)
for b in range(k):
    K[b*m*nn:(b+1)*m*nn, b*m*nn:(b+1)*m*nn] = w1
    for s_ in range(m):
        o = b*m*nn + s_*nn; K[o:o+nn, o:o+nn] = w2
d = K.sum(1)[0]; W = K / d; wbar = d / N
ev = np.sort(np.real(np.linalg.eigvals(W)))[::-1]
vals, counts = np.unique(np.round(ev, 6), return_counts=True)
print("distinct eigenvalues (multiplicity):", {float(v): int(c) for v, c in zip(vals[::-1], counts[::-1])})
lam_block, lam_sub = 1 - w0 / wbar, nn * (w2 - w1) / d
print(f"formula: between blocks 1 - w0/wbar = {lam_block:.3f} (x{k-1}); between sub-blocks "
      f"n(w2-w1)/d = {lam_sub:.3f} (x{k*(m-1)}); within sub-blocks 0 (x{N-k*m})")
assert np.isclose(ev[1], lam_block) and np.isclose(sorted(set(np.round(ev, 8)))[-3], lam_sub)
print(f"old appendix values: {1-(w2-w1)/w2:.3f}, {1-(w1-w0)/w1:.3f};  Gamma = {gap(W):.3f}")

print("\n=== Appendix 19A: Kronecker case ===")
A, B, C = (rownorm(rng.random((q, q))) for q in (3, 4, 2))
Wk = np.kron(np.kron(A, B), C)
assert np.allclose(perron(Wk), np.kron(np.kron(perron(A), perron(B)), perron(C)))
assert np.isclose(gap(Wk), min(gap(A), gap(B), gap(C)))
print("pi and Gamma identities for exact Kronecker products verified")

print("\n=== Proof pass: two-level nesting pi_I = pi_G * pi_N|c * pi_within (Prop. nesting) ===")
def coupling(W, part):
    s = [perron(stoch_complement(W, b)) for b in part]
    Cc = np.array([[s[a] @ W[np.ix_(part[a], part[b])].sum(1) for b in range(len(part))]
                   for a in range(len(part))])
    assert np.allclose(Cc.sum(1), 1)
    return perron(Cc), s
for trial in range(100):
    # 3 countries, each with 2 or 3 sectors of 1-3 individuals
    layout, idx = [], 0
    for c in range(3):
        secs = []
        for s_ in range(rng.integers(2, 4)):
            m_ = int(rng.integers(1, 4)); secs.append(list(range(idx, idx + m_))); idx += m_
        layout.append(secs)
    n = idx
    W = rownorm(rng.random((n, n)) ** 4)
    pi = perron(W)
    countries = [sum(secs, []) for secs in layout]
    for b in countries:                                   # I - W_** invertible, S_B stochastic
        rest = [j for j in range(n) if j not in b]
        assert max(abs(np.linalg.eigvals(W[np.ix_(rest, rest)]))) < 1
        assert np.allclose(stoch_complement(W, b).sum(1), 1)
    xiG, sC = coupling(W, countries)
    assert np.allclose(xiG, [pi[b].sum() for b in countries])
    # pi-weighted aggregate has Perron vector = block sums
    Wt = np.array([[pi[a] @ W[np.ix_(a, b)].sum(1) / pi[a].sum() for b in countries] for a in countries])
    assert np.allclose(perron(Wt), xiG)
    recon = np.zeros(n)
    for c, secs in enumerate(layout):
        Sc = stoch_complement(W, countries[c])
        loc = [[countries[c].index(i) for i in sec] for sec in secs]
        xiN, sW = coupling(Sc, loc)
        for k_, sec in enumerate(secs):
            recon[sec] = xiG[c] * xiN[k_] * sW[k_]
    assert np.allclose(recon, pi)
print("two-level factorization exact in 100 random three-country economies")
