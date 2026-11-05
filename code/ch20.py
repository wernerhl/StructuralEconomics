"""Chapter 20: the illustrative 10x10 Bolivia SAM.  Reproduces every number in
chapters/chapter20.tex (balance, eigenvalues, Perron vector, Pyatt-Round and
Stone decompositions, counterfactual, heatmap cells).  numpy only.
Run from the repo root:  python3 code/ch20.py
"""
import numpy as np

np.set_printoptions(precision=3, suppress=True, linewidth=150)
lab = ['Gas', 'Agr', 'Mfg', 'Svc', 'Lab', 'Cap', 'HH-P', 'HH-NP', 'Gov', 'RoW']
# T[i, j] = payment from account j (column) to account i (row); billions of Bs
T = np.array([
    #Gas  Agr  Mfg  Svc  Lab  Cap  HHP  HHNP Gov  RoW
    [2.1, 0.1, 0.8, 0.3, 0.0, 0.0, 0.0, 0.0, 0.2, 4.5],   # Gas
    [0.0, 1.8, 1.5, 0.5, 0.0, 0.0, 0.8, 0.2, 0.1, 0.8],   # Agr
    [0.3, 0.6, 3.2, 1.8, 0.0, 0.0, 1.5, 1.5, 0.4, 2.2],   # Mfg
    [0.2, 0.3, 1.1, 4.5, 0.0, 0.0, 2.8, 3.2, 1.5, 0.5],   # Svc
    [0.8, 2.1, 2.5, 3.8, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0],   # Lab
    [3.6, 0.5, 1.0, 1.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2],   # Cap
    [0.0, 0.0, 0.0, 0.0, 4.2, 0.8, 0.0, 0.0, 1.5, 0.3],   # HH-P
    [0.0, 0.0, 0.0, 0.0, 6.2, 6.3, 0.0, 0.0, 0.5, 0.8],   # HH-NP
    [1.0, 0.3, 0.5, 0.8, 0.0, 0.0, 0.6, 2.4, 0.0, 0.3],   # Gov
    [0.0, 0.0, 0.9, 0.6, 0.0, 0.0, 1.1, 6.5, 0.5, 0.0]])  # RoW
y = T.sum(1)
assert np.allclose(T.sum(1), T.sum(0)), "SAM must balance"
print("account totals y:", dict(zip(lab, y.round(2))))
print("household consumption added: HH-P", T[[1, 2, 3], 6], " HH-NP", T[[1, 2, 3], 7])


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def gap_and_eigs(W):
    ev = np.linalg.eigvals(W)
    ev = ev[np.argsort(-abs(ev))]
    return 1 - abs(ev[1]), ev


def gini(x):
    x = np.sort(x); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


print("\n=== Row-normalised SAM ===")
W = T / y[:, None]
G, ev = gap_and_eigs(W)
print("eigenvalues:", np.round(ev, 3))
print("moduli     :", np.round(abs(ev), 3))
print(f"Gamma_SAM = {G:.3f}")
pi = perron(W)
assert np.allclose(pi, y / y.sum())       # Perron vector = account size
print("pi = y/sum(y):", dict(zip(lab, pi.round(3))), f" Gini(pi) = {gini(pi):.2f}")

# production block: receipt shares (row-normalised T_PP) and input-cost shares W^P = D^-1 A
TPP = T[:4, :4]
G_rows, _ = gap_and_eigs(TPP / TPP.sum(1, keepdims=True))
A = (TPP / y[:4][None, :]).T             # Ch.1 orientation: a_ij = input j per unit of i
WP = A / A.sum(1, keepdims=True)
G_WP, _ = gap_and_eigs(WP)
print(f"I-O block: Gamma(row-normalised T_PP) = {G_rows:.3f}; Gamma(W^P = D^-1 A) = {G_WP:.3f}")
print(f"Gamma_SAM > Gamma_IO in this SAM: {G > G_rows and G > G_WP}")

# pure 3-cycle illustration: production -> factors -> households -> production
C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
print("pure 3-cycle eigenvalue moduli:", np.round(abs(np.linalg.eigvals(C3)), 3), "=> Gamma = 0")

print("\n=== SAM multipliers and decompositions (endogenous: production, factors, households) ===")
e = list(range(8))
An = T[np.ix_(e, e)] / y[e][None, :]
I = np.eye(8)
M = np.linalg.inv(I - An)
At = np.zeros_like(An); At[:4, :4] = An[:4, :4]; At[6:, 6:] = An[6:, 6:]
# block pattern of An - At: factors <- production, households <- factors, production <- households
Off = An - At
assert np.allclose(Off[:4, :4], 0) and np.allclose(Off[4:6, 4:], 0) and np.allclose(Off[6:, :4], 0) \
    and np.allclose(Off[:4, 4:6], 0)
M1 = np.linalg.inv(I - At)
As = M1 @ Off
M2 = I + As + As @ As
M3 = np.linalg.inv(I - As @ As @ As)
assert np.allclose(M, M3 @ M2 @ M1)
# Stone additive decomposition
own = M1 - I                   # transfer (own-block) effects
openl = (M2 - I) @ M1          # open-loop (cross-block) effects
closed = (M3 - I) @ M2 @ M1    # closed-loop (circular) effects
assert np.allclose(M, I + own + openl + closed)
cols = [0, 1, 2, 3, 6, 7]
print("injection into :", [lab[c] for c in cols])
for name, X in [("initial I", I), ("own-block", own), ("open-loop", openl), ("closed-loop", closed), ("total M", M)]:
    print(f"{name:12s}", np.round(X.sum(0)[cols], 2))
print("product of column sums of M1,M2,M3 (NOT a valid decomposition):",
      np.round((M1.sum(0) * M2.sum(0) * M3.sum(0))[cols], 2), " vs true", np.round(M.sum(0)[cols], 2))
tot = M.sum(0)
print(f"gas vs mfg: {(1-tot[0]/tot[2])*100:.0f}% below; gas vs agr: {(1-tot[0]/tot[1])*100:.0f}% below; "
      f"agr/gas = {tot[1]/tot[0]:.2f}")
print(f"poor vs non-poor HH total multiplier: {tot[6]:.2f} vs {tot[7]:.2f}")
print(f"household import (+saving abroad) shares of spending: HH-P {T[9,6]/y[6]:.2f}, HH-NP {T[9,7]/y[7]:.2f}")
print(f"open-loop ratio gas/mfg = {openl.sum(0)[0]/openl.sum(0)[2]:.2f}")
print(f"gas column: capital share {T[5,0]/y[0]:.2f}, labour share {T[4,0]/y[0]:.2f}, "
      f"direct leakage to Gov+RoW {(T[8,0]+T[9,0])/y[0]:.3f}")
print(f"capital income to poor HH: {T[6,5]/y[5]:.2f}; labour income to poor HH: {T[6,4]/y[4]:.3f}")

print("\n=== Counterfactual: shift 1.0 of gas value added from capital to labour ===")
x_exo = T[np.ix_(e, [8, 9])].sum(1)          # injections from Gov and RoW
assert np.allclose(M @ x_exo, y[e])
An2 = An.copy(); An2[4, 0] += 1.0 / y[0]; An2[5, 0] -= 1.0 / y[0]
M_2 = np.linalg.inv(I - An2)
y2 = M_2 @ x_exo
# rebuild the SAM: endogenous columns scale with new activity; exogenous columns unchanged
lk = T[np.ix_([8, 9], e)] / y[e][None, :]     # leakage coefficients to Gov, RoW
T2 = T.copy()
T2[np.ix_(e, e)] = An2 * y2[None, :]
T2[np.ix_([8, 9], e)] = lk * y2[None, :]
W2 = T2 / T2.sum(1, keepdims=True)
G2, _ = gap_and_eigs(W2)
print(f"gas share of VA moved: 1.0/{T[4,0]+T[5,0]:.1f} = {1/(T[4,0]+T[5,0]):.2f} (incl. taxes {1/(T[4,0]+T[5,0]+T[8,0]):.2f})")
print(f"labour->poor share {T[6,4]/y[4]:.2f} => poor HH get {T[6,4]/y[4]:.2f} per unit directly")
print(f"Gamma_SAM: {G:.4f} -> {G2:.4f}  (change {G2-G:+.4f})")
print(f"gas total multiplier: {tot[0]:.2f} -> {M_2.sum(0)[0]:.2f} ({(M_2.sum(0)[0]/tot[0]-1)*100:+.1f}%)")
sh, sh2 = y[6] / (y[6] + y[7]), y2[6] / (y2[6] + y2[7])
print(f"poor-household income share: {sh:.3f} -> {sh2:.3f}")

print("\n=== Heatmap cells: blue!XX with XX = round(100 * T_ij / max T) ===")
H = np.rint(100 * T / T.max()).astype(int)
print(H)
print("tikz rows:", "/".join(",".join(str(v) for v in row) for row in H))
print(f"agriculture labour share of output: {T[4,1]/y[1]:.2f}")

print("\n=== Proof of Prop. 'What the SAM adds', item 3: opposite ordering ===")
# accounts: P1, P2 (production), F (factor), H (household); small intra-production flows
eps = 0.05
# balance by construction: make each column total equal its row total
Tx = np.array([
    [eps, eps, 0.0, 1 - 2 * eps],
    [eps, eps, 0.0, 1 - 2 * eps],
    [1 - 2 * eps, 1 - 2 * eps, 0.0, 0.0],
    [0.0, 0.0, 2 * (1 - 2 * eps), 0.0]])
assert np.allclose(Tx.sum(1), Tx.sum(0))
Gx, evx = gap_and_eigs(Tx / Tx.sum(1, keepdims=True))
Gio, _ = gap_and_eigs(Tx[:2, :2] / Tx[:2, :2].sum(1, keepdims=True))
print(f"Gamma_IO = {Gio:.2f}, Gamma_SAM = {Gx:.2f}  (eigenvalue moduli {np.round(abs(evx),3)})")
assert Gx < Gio
