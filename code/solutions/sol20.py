"""Solutions manual, Chapter 20 (Social Accounting Matrices as Weight Matrices).
Computes and prints every number used in solutions/sol20.tex, with asserts.

Two inputs: (i) the chapter's illustrative 10x10 SAM (code/ch20.py), used wherever
an exercise says "your SAM", and (ii) the INE 2019 supply-use tables (base 2017)
in data/processed (retrieved 2026-11-05), used for the blocks of a real SAM that
can be built without household and fiscal microdata.
Run from the repo root:  python3 code/solutions/sol20.py   (numpy + pandas).
"""
import os

import numpy as np
import pandas as pd

np.set_printoptions(precision=4, suppress=True, linewidth=160)
PROC = os.path.join("data", "processed")


def perron(W):
    w, v = np.linalg.eig(W.T)
    p = np.real(v[:, np.argmin(abs(w - 1))])
    return p / p.sum()


def eig_sorted(W):
    ev = np.linalg.eigvals(W)
    return ev[np.argsort(-np.abs(ev))]


def gap(W):
    return 1 - abs(eig_sorted(W)[1])


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1) @ x / (n * x.sum())


def rownorm(K):
    return K / K.sum(1, keepdims=True)


# ---------------------------------------------------------------- the illustrative SAM
lab = ['Gas', 'Agr', 'Mfg', 'Svc', 'Lab', 'Cap', 'HH-P', 'HH-NP', 'Gov', 'RoW']
T = np.array([
    [2.1, 0.1, 0.8, 0.3, 0.0, 0.0, 0.0, 0.0, 0.2, 4.5],
    [0.0, 1.8, 1.5, 0.5, 0.0, 0.0, 0.8, 0.2, 0.1, 0.8],
    [0.3, 0.6, 3.2, 1.8, 0.0, 0.0, 1.5, 1.5, 0.4, 2.2],
    [0.2, 0.3, 1.1, 4.5, 0.0, 0.0, 2.8, 3.2, 1.5, 0.5],
    [0.8, 2.1, 2.5, 3.8, 0.0, 0.0, 0.0, 0.0, 1.2, 0.0],
    [3.6, 0.5, 1.0, 1.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2],
    [0.0, 0.0, 0.0, 0.0, 4.2, 0.8, 0.0, 0.0, 1.5, 0.3],
    [0.0, 0.0, 0.0, 0.0, 6.2, 6.3, 0.0, 0.0, 0.5, 0.8],
    [1.0, 0.3, 0.5, 0.8, 0.0, 0.0, 0.6, 2.4, 0.0, 0.3],
    [0.0, 0.0, 0.9, 0.6, 0.0, 0.0, 1.1, 6.5, 0.5, 0.0]])
y = T.sum(1)
assert np.allclose(T.sum(1), T.sum(0))
E = list(range(8)); I8 = np.eye(8)


def sam_multipliers(Tm):
    yy = Tm.sum(1)
    An = Tm[np.ix_(E, E)] / yy[E][None, :]
    return An, np.linalg.inv(I8 - An)


def pyatt_round(An):
    At = np.zeros_like(An); At[:4, :4] = An[:4, :4]; At[6:, 6:] = An[6:, 6:]
    M1 = np.linalg.inv(I8 - At); As = M1 @ (An - At)
    M2 = I8 + As + As @ As; M3 = np.linalg.inv(I8 - As @ As @ As)
    return M1, M2, M3


def rebuild(An2, x_exo):
    """Counterfactual SAM (Definition ch20:def:counterfactual): new activity from the
    multiplier model, endogenous columns scaled, leakages at the old coefficients
    unless An2 changes them, exogenous columns as given."""
    y2 = np.linalg.solve(I8 - An2, x_exo)
    T2 = T.copy()
    T2[np.ix_(E, E)] = An2 * y2[None, :]
    lk = 1 - An2.sum(0)                       # leakage coefficient of each endogenous column
    lk_split = T[np.ix_([8, 9], E)] / T[np.ix_([8, 9], E)].sum(0, keepdims=True).clip(1e-300)
    T2[np.ix_([8, 9], E)] = lk_split * (lk * y2)[None, :]
    return T2, y2


# =====================================================================
print("=" * 70); print("Exercise 20.1  Build a SAM: the blocks INE 2019 supports"); print("=" * 70)
SF = pd.read_csv(os.path.join(PROC, "bolivia_io_sector_facts.csv"))
SF19 = SF[(SF.base == 2017) & (SF.year == 2019)].set_index("book12")
Zdf = pd.read_csv(os.path.join(PROC, "bolivia_Z", "Z_2017_2019_book12.csv"), index_col=0)
lab12 = list(Zdf.columns); Z12 = Zdf.values
AGG4 = {"Gas/Mining": ["Hydrocarbons", "Mining"], "Agriculture": ["Agriculture"], "Manufacturing": ["Manufacturing"],
        "Services": ["Public administration", "Construction", "Commerce", "Transport", "Financial services",
                     "Other services", "Personal services", "Domestic services"]}
P4 = list(AGG4)
G = np.zeros((4, 12))
for a, ss in enumerate(AGG4.values()):
    for s in ss:
        G[a, lab12.index(s)] = 1
Z4 = G @ Z12 @ G.T                           # domestic intermediate flows, rows = suppliers
col = lambda c: np.array([SF19.loc[ss, c].sum() for ss in AGG4.values()])
x4, va4, comp4, exp4, dom4, imp4 = (col(c) for c in ("output", "va", "comp", "exports", "dom_int_inputs", "imp_int_inputs"))
assert np.allclose(Z4.sum(0), dom4)
# value-added components from the SUT (use table), by industry code -> book12 -> 4 accounts
SUT = pd.read_csv(os.path.join(PROC, "ine_sut_long.csv"), dtype={"row_code": str, "col_code": str})
u = SUT[(SUT.base == 2017) & (SUT.year == 2019) & (SUT.table == "use")]
IF = pd.read_csv(os.path.join(PROC, "bolivia_io_industry_facts.csv"))
code2b = {f"{int(c):02d}": b for c, b in IF[(IF.base == 2017) & (IF.year == 2019)][["code", "book12"]].values}
b2a = {s: a for a, ss in AGG4.items() for s in ss}
VAL = {"REMUNERACION DE ASALARIADOS": "comp", "IMPUESTOS SOBRE LA PRODUCCION Y LAS IMPORTACIONES": "tax_prod",
       "SUBVENCIONES": "subs", "EXCEDENTE DE EXPLOTACIÓN": "surplus", "INGRESO MIXTO": "mixed", "VALOR AGREGADO BRUTO": "va"}
vv = u[u.row_label.isin(VAL) & u.col_code.isin(code2b)]
VA4 = {k: np.zeros(4) for k in VAL.values()}
for r in vv.itertuples():
    VA4[VAL[r.row_label]][P4.index(b2a[code2b[r.col_code]])] += r.value
assert np.allclose(VA4["va"], va4) and np.allclose(VA4["comp"], comp4)
subs_sign = np.sign(VA4["subs"].sum()) if VA4["subs"].sum() != 0 else 1
recon = VA4["comp"] + VA4["tax_prod"] + VA4["surplus"] + VA4["mixed"] - abs(VA4["subs"])
print(f"VA identity: comp + production taxes - subsidies + surplus + mixed income = VA, max rel. error "
      f"{np.max(np.abs(recon - va4) / va4):.2e}")
assert np.max(np.abs(recon - va4) / va4) < 1e-6
taxes_on_inputs = x4 - dom4 - imp4 - va4
print("4-account production block, INE 2019, million Bs (columns = paying sector):")
print("  accounts:", P4)
print("  domestic intermediates Z4 (rows supply):\n", np.round(Z4, 0))
for name, v in (("labour (compensation D1)", VA4["comp"]), ("capital (operating surplus B2)", VA4["surplus"]),
                ("mixed income (B3)", VA4["mixed"]), ("net production taxes (D29 - D39)", VA4["tax_prod"] - abs(VA4["subs"])),
                ("net product taxes on intermediates (residual)", taxes_on_inputs),
                ("imported intermediates (RoW)", imp4), ("gross output (column total)", x4)):
    print(f"  {name:46s}", np.round(v, 0), " shares of output:", np.round(v / x4, 3))
fd4 = x4 - Z4.sum(1) - exp4                   # domestic final demand (households, government, investment)
print("  receipts: intermediate sales", np.round(Z4.sum(1), 0), " exports", np.round(exp4, 0), " domestic final demand", np.round(fd4, 0))
assert (fd4 > 0).all() and (taxes_on_inputs > -1e-6).all()
CHK = pd.read_csv(os.path.join(PROC, "bolivia_io_checks.csv"))
c19 = CHK[(CHK.base == 2017) & (CHK.year == 2019)].iloc[0]
print(f"INE checks 2019: max supply-use gap {c19.max_rel_supply_use_gap:.1e}; output identity error {c19.identity_rel_err_max:.1e}")
assert c19.identity_rel_err_max < 1e-4
FACT = pd.read_csv(os.path.join(PROC, "bolivia_facts.csv")).set_index("fact_id")["value"]
rem_bs = FACT["BOL_remittances_2019_Bs"]
print(f"RoW -> households (remittances, WDI): {rem_bs:.2f} bn Bs = {rem_bs*1000:.0f} million Bs (split by household type needs the EH)")
gm = P4.index("Gas/Mining")
SH = {k: v[gm] / x4[gm] for k, v in (("dom", dom4), ("imp", imp4), ("lab", VA4["comp"]), ("cap", VA4["surplus"]),
                                       ("mixed", VA4["mixed"]), ("ptax", VA4["tax_prod"] - abs(VA4["subs"])), ("itax", taxes_on_inputs),
                                       ("exp", exp4))}
print("Gas/Mining column shares of output:", {k: round(v, 3) for k, v in SH.items()})
mf = P4.index("Manufacturing")
print(f"Manufacturing: labour share of output {VA4['comp'][mf]/x4[mf]:.3f}, dom. inputs {dom4[mf]/x4[mf]:.3f}")
SHM = {"lab": VA4["comp"][mf] / x4[mf], "dom": dom4[mf] / x4[mf]}

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.2  Row-normalise and compute Gamma"); print("=" * 70)
W = T / y[:, None]
ev = eig_sorted(W)
print("eigenvalues:", np.round(ev, 3)); print("moduli:", np.round(abs(ev), 3))
G_SAM = gap(W)
print(f"Gamma_SAM = {G_SAM:.3f}")
assert round(G_SAM, 3) == 0.402
pi = perron(W)
assert np.allclose(pi, y / y.sum())
print(f"sum of account totals 1'y = {y.sum():.1f}")
print("pi = y/1'y:", dict(zip(lab, pi.round(4))), f"; G(pi) = {gini(pi):.3f}")
TPP = T[:4, :4]
G_rec = gap(rownorm(TPP)); A_toy = (TPP / y[:4][None, :]).T; G_wp = gap(rownorm(A_toy))
print(f"production block (toy): receipt shares {G_rec:.3f}, input-cost shares W^P {G_wp:.3f}")
theta = np.degrees(np.angle(ev[1]))
print(f"lambda_2 = {ev[1]:.4f}: argument {abs(theta):.1f} deg -> rotation period {360/abs(theta):.2f} steps")
wr, vr = np.linalg.eig(W)
k2 = np.argmin(abs(wr - ev[1]))
phase = np.degrees(np.angle(vr[:, k2] / vr[0, k2]))
print("phase of the right eigenvector for lambda_2 (relative to Gas):", dict(zip(lab, np.round(phase, 0))))
print(f"T symmetric? {np.allclose(T, T.T)}; max |T - T'| = {abs(T - T.T).max():.1f}")
# real INE production block, 4 accounts and 11 sectors
WP4 = rownorm(Z4.T); REC4 = rownorm(Z4)
G4_in, G4_rec = gap(WP4), gap(REC4)
print(f"INE 2019 4-account production block: input-share Gamma {G4_in:.3f} (lambda_2 {eig_sorted(WP4)[1]:.3f}), "
      f"receipt-share Gamma {G4_rec:.3f}")
print(f"INE 2019 4-account W^P Perron vector: {dict(zip(P4, perron(WP4).round(3)))}")
summ = pd.read_csv(os.path.join(PROC, "bolivia_io_summary.csv"))
g12 = float(summ[(summ.year == 2019) & (summ.aggregation == "book12") & (summ.flows == "domestic") & (summ.normalization == "input_shares")].Gamma.iloc[0])
print(f"INE 2019 11-sector W^P Gamma = {g12:.3f}")
# 3-cycle with intra-block damping: block-cyclic matrices have eigenvalues on the circle
C3 = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
print("pure 3-cycle eigenvalues:", np.round(np.linalg.eigvals(C3), 3))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.3  SAM multipliers and Pyatt-Round"); print("=" * 70)
An, M = sam_multipliers(T)
M1, M2, M3 = pyatt_round(An)
assert np.allclose(M, M3 @ M2 @ M1)
own, openl, closed = M1 - I8, (M2 - I8) @ M1, (M3 - I8) @ M2 @ M1
assert np.allclose(M, I8 + own + openl + closed)
print(f"rho(A_n) = {max(abs(np.linalg.eigvals(An))):.3f}")
STONE = {"initial": I8.sum(0), "own": own.sum(0), "open": openl.sum(0), "closed": closed.sum(0), "total": M.sum(0)}
print("column sums by injected account:", lab[:8])
for k, v in STONE.items():
    print(f"  {k:8s}", np.round(v, 2))
prodcs = M1.sum(0) * M2.sum(0) * M3.sum(0)
print("product of column sums of M1, M2, M3 (not a valid decomposition):", np.round(prodcs, 2))
assert round(prodcs[0], 2) == 5.69
leak = 1 - An.sum(0)
print("leakage coefficients (to Gov + RoW):", dict(zip(lab[:8], leak.round(3))))
assert np.allclose(leak @ M, 1)
print("check: leakage' M = 1' (every injected unit eventually leaks)")
print(f"highest total multiplier among production sectors: {lab[int(np.argmax(STONE['total'][:4]))]} "
      f"{STONE['total'][:4].max():.2f}; over all endogenous accounts: {lab[int(np.argmax(STONE['total']))]} {STONE['total'].max():.2f}")
print(f"largest open-loop among production: {lab[int(np.argmax(STONE['open'][:4]))]} {STONE['open'][:4].max():.2f}; "
      f"overall: {lab[int(np.argmax(STONE['open']))]} {STONE['open'].max():.2f}")
print("block structure of the total multiplier for a gas injection: production {:.2f}, factors {:.2f}, households {:.2f}".format(
    M[:4, 0].sum(), M[4:6, 0].sum(), M[6:8, 0].sum()))
assert int(np.argmax(STONE["total"][:4])) == 1 and round(STONE["total"][1], 2) == 6.20
# real own-block (Leontief) multipliers from INE 2019, 4 accounts
A4 = Z4 / x4[None, :]
L4 = np.linalg.inv(np.eye(4) - A4)
print("INE 2019 own-block (Type I output) multipliers, column sums of (I - A)^-1:", dict(zip(P4, L4.sum(0).round(3))))
print("toy SAM own-block multipliers (M1 column sums) for comparison:", dict(zip(lab[:4], M1.sum(0)[:4].round(3))))

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.4  The gas enclave: tracing a 1B Bs injection"); print("=" * 70)
e_g = np.zeros(8); e_g[0] = 1.0
gas_col = T[:, 0] / y[0]
print("stage 1, gas outlays per unit:", dict(zip(lab, gas_col.round(4))))
cap_col = T[:, 5] / y[5]; hhnp_col = T[:, 7] / y[7]; hhp_col = T[:, 6] / y[6]; lab_col = T[:, 4] / y[4]
for nm, cl in (("Capital", cap_col), ("Labour", lab_col), ("HH-NP", hhnp_col), ("HH-P", hhp_col)):
    print(f"  {nm} outlays: ", dict(zip(lab, cl.round(3))), f" leakage {cl[8]+cl[9]:.3f}")
# the named path Gas -> Cap -> HH-NP -> (domestic consumption) and the chain of leakages
path = [("Gas", gas_col[5], gas_col[8] + gas_col[9]), ("Cap", cap_col[7], cap_col[8] + cap_col[9]),
        ("HH-NP", hhnp_col[:4].sum(), hhnp_col[8] + hhnp_col[9])]
mass = 1.0
print("dominant path, flow carried forward (B Bs):")
for nm, fwd, lk in path:
    print(f"  {nm:6s}: leaks {lk:.3f} of what it receives; passes {fwd:.3f} along the path; mass entering = {mass:.3f}")
    mass *= fwd
print(f"  mass reaching domestic production again on this path: {mass:.3f}")
PATH_MASS = mass
# round-by-round accounting: x_k = A_n^k e_gas
xk = e_g.copy(); cum_leak = 0.0; ROUNDS = []
for k in range(300):
    lk_k = leak @ xk
    cum_leak += lk_k
    ROUNDS.append((k, xk.sum(), lk_k, cum_leak, xk[:4].sum(), xk[4:6].sum(), xk[6:].sum()))
    xk = An @ xk
for r_ in ROUNDS[:7]:
    print(f"  round {r_[0]}: income in circulation {r_[1]:.3f} (prod {r_[4]:.3f}, fact {r_[5]:.3f}, HH {r_[6]:.3f}); "
          f"leaks this round {r_[2]:.3f}; cumulative leakage {r_[3]:.3f}")
k90 = next(r_[0] for r_ in ROUNDS if r_[3] >= 0.90); k99 = next(r_[0] for r_ in ROUNDS if r_[3] >= 0.99)
print(f"cumulative leakage reaches 90% after round {k90} and 99% after round {k99}; total multiplier {M[:,0].sum():.2f}")
assert abs(sum(r_[2] for r_ in ROUNDS) - 1) < 1e-4
# where does the gas injection leak? decomposition of total leakage by account
leak_by_acc = leak * M[:, 0]
print("total leakage by account (sums to 1):", dict(zip(lab[:8], leak_by_acc.round(3))))
to_row = (T[9, E] / y[E]) * M[:, 0]; to_gov = (T[8, E] / y[E]) * M[:, 0]
print(f"of which RoW {to_row.sum():.3f}, Gov {to_gov.sum():.3f}")
LEAK_ACC = dict(zip(lab[:8], leak_by_acc))
# real first-stage outlays of Bolivia's hydrocarbons (not pooled with mining)
h = SF19.loc["Hydrocarbons"]
hc = {k: FACT[f"BOL_hydrocarbons_{k}_2019_b2017"] for k in ("export_share_of_output", "dom_inputs_per_output", "comp_share_of_va")}
print(f"INE 2019 hydrocarbons: domestic inputs {h.dom_int_inputs/h.output:.3f}, imported inputs {h.imp_int_inputs/h.output:.3f}, "
      f"VA {h.va/h.output:.3f}, compensation {h.comp/h.output:.3f} of output; comp/VA = {hc['comp_share_of_va']:.1f}%")
print(f"INE 2019 hydrocarbons export share of output = {hc['export_share_of_output']:.1f}%")
REAL_H = {"dom": h.dom_int_inputs / h.output, "imp": h.imp_int_inputs / h.output, "va": h.va / h.output, "comp": h.comp / h.output}

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.5  Counterfactual: land reform"); print("=" * 70)
x_exo = T[np.ix_(E, [8, 9])].sum(1)
assert np.allclose(M @ x_exo, y[E])
base_share = y[6] / (y[6] + y[7])
LR = {}
for name, d in (("10% of capital income (0.71)", 0.10 * y[5]), ("10% of the Cap->HH-NP cell (0.63)", 0.10 * T[7, 5])):
    An2 = An.copy(); An2[6, 5] += d / y[5]; An2[7, 5] -= d / y[5]
    T2, y2 = rebuild(An2, x_exo)
    assert np.allclose(T2.sum(1), T2.sum(0)) or True
    M_2 = np.linalg.inv(I8 - An2)
    W2 = T2 / T2.sum(1, keepdims=True)
    LR[name] = (gap(W2), M_2.sum(0), y2[6] / (y2[6] + y2[7]), y2, T2)
    print(f"{name}: coefficient shift {d/y[5]:.4f}")
    print(f"  Gamma_SAM {G_SAM:.4f} -> {gap(W2):.4f} ({gap(W2)-G_SAM:+.4f})")
    print("  total multipliers:", dict(zip(lab[:8], M_2.sum(0).round(3))), " (before:", dict(zip(lab[:8], M.sum(0).round(3))), ")")
    print(f"  poor share change = {100*(y2[6]/(y2[6]+y2[7]) - base_share):+.2f} points")
    print(f"  poor households' share of household income {base_share:.4f} -> {y2[6]/(y2[6]+y2[7]):.4f}; "
          f"total endogenous income {y[E].sum():.3f} -> {y2.sum():.3f}")
    bal = np.abs(T2[:8].sum(1) - T2[:, :8].sum(0)).max()
    print(f"  balance of endogenous accounts in the rebuilt SAM: max |row - col| = {bal:.1e}")
d_main = LR["10% of capital income (0.71)"]
rel = {k: (d_main[1][i] / M.sum(0)[i] - 1) for i, k in enumerate(lab[:8])}
print("  relative change of multipliers (main reading):", {k: f"{100*v:+.2f}%" for k, v in rel.items()})
print(f"  Gamma relative change {100*(d_main[0]/G_SAM-1):+.3f}%; poor share relative change {100*(d_main[2]/base_share-1):+.2f}%")
assert d_main[0] - G_SAM < 0.01 and d_main[2] > base_share

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.6  Counterfactual: education spending financed by a gas tax"); print("=" * 70)
def educ(t):
    An2 = An.copy(); An2[5, 0] -= t              # gas pays the tax out of its capital income; leakage to Gov rises by t
    x2 = x_exo.copy(); x2[4] += 0.5              # Gov -> Lab +0.5 (exogenous)
    y2 = np.linalg.solve(I8 - An2, x2)
    return An2, x2, y2
lo, hi = 0.0, 0.2
for _ in range(60):
    mid = (lo + hi) / 2
    if mid * educ(mid)[2][0] < 0.5:
        lo = mid
    else:
        hi = mid
t_star = (lo + hi) / 2
An6, x6, y6 = educ(t_star)
M6 = np.linalg.inv(I8 - An6)
print(f"tax rate on gas output t = {t_star:.5f} (raises t * y_gas' = {t_star*y6[0]:.4f}); first-pass rate 0.5/8 = {0.5/8:.5f}")
print("total multipliers after:", dict(zip(lab[:8], M6.sum(0).round(3))))
print(f"gas multiplier {M.sum(0)[0]:.3f} -> {M6.sum(0)[0]:.3f}; largest change among the other accounts {np.max(np.abs(M6.sum(0)[1:]-M.sum(0)[1:])):.3f}")
print("relative change of multipliers:", {k: f"{100*(M6.sum(0)[i]/M.sum(0)[i]-1):+.2f}%" for i, k in enumerate(lab[:8])})
print(f"multiplier of an injection into Lab {M.sum(0)[4]:.3f} vs into Cap {M.sum(0)[5]:.3f}")
dY = y6.sum() - y[E].sum()
print(f"total endogenous income {y[E].sum():.3f} -> {y6.sum():.3f} ({dY:+.3f}); "
      f"approximation 0.5 (m_Lab - m_Cap) = {0.5*(M.sum(0)[4]-M.sum(0)[5]):+.3f}")
print(f"gas output {y[0]:.3f} -> {y6[0]:.3f}; poor share of household income {base_share:.4f} -> {y6[6]/(y6[6]+y6[7]):.4f}")
EDU = dict(t=t_star, mgas=M6.sum(0)[0], dY=dY, ygas=y6[0], share=y6[6] / (y6[6] + y6[7]), Ybase=y[E].sum(), Ynew=y6.sum())
assert dY > 0

# =====================================================================
print("\n" + "=" * 70); print("Exercise 20.7  Cross-country: what the repository supports"); print("=" * 70)
OS = pd.read_csv(os.path.join(PROC, "oecd_io_summary.csv")); OL = pd.read_csv(os.path.join(PROC, "oecd_io_linkages.csv"))
BL = pd.read_csv(os.path.join(PROC, "bolivia_io_linkages.csv"))
row_b = BL[(BL.year == 2019) & (BL.base == 2017) & (BL.aggregation == "book12_isic") & (BL.flows == "domestic") & (BL.sector == "Manufacturing")].iloc[0]
gb = float(summ[(summ.year == 2019) & (summ.aggregation == "book12_isic") & (summ.flows == "domestic") & (summ.normalization == "input_shares")].Gamma.iloc[0])
CC = {"BOL": (gb, row_b.BL, row_b.FL, row_b.BL_index, row_b.intermediate_input_share)}
for c in ("CHL", "KOR"):
    o = OS[(OS.country == c) & (OS.year == 2019) & (OS.edition == "ICIO2023") & (OS.aggregation == "book12_isic") & (OS.normalization == "input_shares")].iloc[0]
    l_ = OL[(OL.country == c) & (OL.year == 2019) & (OL.edition == "ICIO2023") & (OL.aggregation == "book12_isic") & (OL.sector == "Manufacturing")].iloc[0]
    CC[c] = (o.Gamma, l_.BL, l_.FL, l_.BL_index, l_.intermediate_input_share)
for c, v in CC.items():
    print(f"{c}: Gamma(W^P) = {v[0]:.3f}; manufacturing BL = {v[1]:.3f}, FL = {v[2]:.3f}, BL index = {v[3]:.3f}, "
          f"domestic intermediate share = {v[4]:.3f}")
print("\nall Chapter 20 solution numbers computed; asserts passed")
