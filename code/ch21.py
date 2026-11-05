"""Chapter 21: media-overlap information matrices.  Reproduces every number in
chapters/chapter21.tex.  numpy only.  Run from the repo root:  python3 code/ch21.py
"""
import numpy as np

np.set_printoptions(precision=3, suppress=True, linewidth=150)
groups = ['Urban NI men', 'Urban NI women', 'Urban indigenous', 'Rural NI',
          'Rural indig. men', 'Rural indig. women']
outlets = ['Private TV', 'State', 'Community radio', 'Social media', 'Word of mouth']
S = np.array([[.45, .15, .02, .30, .08],
              [.40, .18, .03, .28, .11],
              [.20, .30, .15, .20, .15],
              [.25, .20, .10, .15, .30],
              [.08, .22, .30, .10, .30],
              [.05, .15, .25, .08, .47]])
assert np.allclose(S.sum(1), 1)
p = np.array([0.17, 0.18, 0.25, 0.10, 0.15, 0.15])   # illustrative population shares
assert np.isclose(p.sum(), 1)


def perron(W):
    w, v = np.linalg.eig(W.T)
    q = np.real(v[:, np.argmin(abs(w - 1))])
    return q / q.sum()


def gap(W):
    m = np.sort(abs(np.linalg.eigvals(W)))[::-1]
    return 1 - m[1], m


def neff(q):
    return 1 / (q ** 2).sum()


print("=== Raw overlap O = S S' ===")
O = S @ S.T
print(O)
print(f"O11 = {O[0,0]:.3f}, O12 = {O[0,1]:.3f}, O16 = {O[0,5]:.3f}")

print("\n=== Unweighted row-normalisation (as in the earlier draft) ===")
Wu = O / O.sum(1, keepdims=True)
Gu, mu = gap(Wu)
pu = perron(Wu)
assert np.allclose(pu, O.sum(1) / O.sum())          # reversible: pi ~ row sums
print(f"Gamma = {Gu:.2f} (|lambda_2| = {mu[1]:.2f}); pi = {pu.round(3)}; Neff = {neff(pu):.2f} of 6")

print("\n=== Population-weighted group matrix W_g = S diag(m)^-1 S' diag(p) ===")
m = p @ S                                             # outlet reach (audience shares)
Wg = S @ np.diag(1 / m) @ S.T @ np.diag(p)
assert np.allclose(Wg.sum(1), 1)
Gg, mg = gap(Wg)
pg = perron(Wg)
assert np.allclose(pg, p)                             # pi = population shares
print(f"Gamma = {Gg:.2f} (|lambda_2| = {mg[1]:.2f}); pi = population shares; Neff = {neff(p):.2f} of 6")

print("\n=== Outlet projection W_k = diag(m)^-1 S' diag(p) S ===")
Wk = np.diag(1 / m) @ S.T @ np.diag(p) @ S
assert np.allclose(Wk.sum(1), 1)
pk = perron(Wk)
assert np.allclose(pk, m)                             # pi = reach
Gk, _ = gap(Wk)
print("reach / pi:", dict(zip(outlets, m.round(3))), f"Neff = {neff(m):.2f} of 5, Gamma = {Gk:.2f}")
print(f"private TV + state = {m[0]+m[1]:.3f}")
old = np.array([.38, .28, .18, .11, .05])
print(f"(earlier draft's outlet pi: Neff = {neff(old):.2f}, not 3.4)")

print("\n=== Mixing parameters from raw overlaps (cross / within means) ===")
def lam(A, B):
    cross = O[np.ix_(A, B)].mean()
    within = np.mean([O[np.ix_(A, A)].mean(), O[np.ix_(B, B)].mean()])
    return cross / within
print(f"gender, urban NI : {O[0,1]/np.mean([O[0,0],O[1,1]]):.2f}")
print(f"gender, rural ind: {O[4,5]/np.mean([O[4,4],O[5,5]]):.2f}")
print(f"ethnicity        : {lam([0,1,3],[2,4,5]):.2f}")
print(f"urban-rural      : {lam([0,1,2],[3,4,5]):.2f}")
print(f"darkest cross entry: O16 / mean(O11,O66) = {O[0,5]/np.mean([O[0,0],O[5,5]]):.2f}")

print("\n=== Exercise 21.6: add 0.10 to the darkest raw entry (symmetrically) ===")
O2 = O.copy(); O2[0, 5] += 0.10; O2[5, 0] += 0.10
W2 = O2 / O2.sum(1, keepdims=True)
print(f"Gamma: {Gu:.3f} -> {gap(W2)[0]:.3f}; Neff: {neff(pu):.3f} -> {neff(perron(W2)):.3f}")
