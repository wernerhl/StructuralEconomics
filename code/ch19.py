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

print("\n=== Examples ch19:ex:bolivia_three, ch19:ex:potosi, ch19:ex:vicious (sourced + illustrative) ===")
# Chapter 17's values (INE 2019 W^P, ILOSTAT 2019 employment) are imported from code/ch17.py, which
# reads data/processed, so the two chapters cannot drift apart.  The Potosi vendor works in
# "personal services" (ISIC I, R, S); the gas executive in hydrocarbons, whose employment ILOSTAT
# reports only together with mining (ISIC B).
with contextlib.redirect_stdout(io.StringIO()):
    _c17 = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch17.py"))
chk, rd, FACTS = _c17["check"], _c17["rd"], _c17["FACTS"]
ratio17, pi_emp17, emp17 = _c17["ratio"], _c17["pi_emp"], _c17["emp"]
GWP = {row["year"]: row for row in rd("global_wp_summary.csv")}["2019"]
GPI = {row["country"]: float(row["pi"]) for row in rd("global_wp_pi.csv") if row["year"] == "2019"}

# --- three scales (global trade W^P, 230 economies; national INE W^P, 11 sectors) ---
Gam_G, Gam_N = float(GWP["Gamma"]), _c17["gap"](_c17["W"])
chk("Gamma_G", Gam_G, 0.452, 3); chk("Gamma_N", Gam_N, 0.467, 3)
chk("pi_USA/pi_LAO", float(GWP["pi_USA_over_LAO"]), 193, 0)
assert abs(GPI["USA"] / GPI["LAO"] - float(GWP["pi_USA_over_LAO"])) < 1e-9
p17, sec17 = _c17["p"], _c17["sectors"]
i_max, i_min = int(np.argmax(p17)), int(np.argmin(p17))
print(f"Gamma_G = {Gam_G:.3f} (N = {GWP['N']}), Gamma_N = {Gam_N:.3f} (N = {len(sec17)})")
print(f"global pi_USA/pi_LAO = {float(GWP['pi_USA_over_LAO']):.0f}; national max/min pi = "
      f"{sec17[i_max]}/{sec17[i_min]} = {p17[i_max]:.3f}/{p17[i_min]:.4f} = {p17[i_max]/p17[i_min]:.0f}")
assert sec17[i_max] == "Other services" and sec17[i_min] == "Public administration"
chk("national max/min", p17[i_max] / p17[i_min], 53, 0)
chk("|lambda_2| national", abs(_c17["eig_sorted"](_c17["W"])[1]), 0.53, 2)
chk("Korea Gamma_N", float(_c17["okor"]["Gamma"]), 0.511, 3)

# --- per-capita nesting (Example ch19:ex:potosi) ---
piG_bol, piG_deu = GPI["BOL"], GPI["DEU"]
L_bol, L_deu = FACTS["BOL_labour_force_2019_mn"], FACTS["DEU_labour_force_2019_mn"]
chk("pi_G BOL", piG_bol, 0.00053, 5); chk("pi_G DEU", piG_deu, 0.082, 3)
chk("pi_DEU/pi_BOL", piG_deu / piG_bol, 157, 0)
chk("L_BOL", L_bol, 5.95, 2); chk("L_DEU", L_deu, 43.98, 2)
xi_gas, e_gas = pi_emp17["Hydrocarbons and mining"], emp17["Hydrocarbons and mining"]
xi_inf, e_inf = pi_emp17["Personal services"], emp17["Personal services"]
chk("pi gas+mining", xi_gas, 0.077, 3); chk("emp gas+mining", e_gas, 0.013, 3)
chk("pi personal", xi_inf, 0.021, 3); chk("emp personal", e_inf, 0.107, 3)
sector_gas, sector_inf = ratio17["Hydrocarbons and mining"], ratio17["Personal services"]
chk("sector gas", sector_gas, 5.7, 1); chk("sector inf", sector_inf, 0.20, 2)
sec_ratio = sector_gas / sector_inf
print(f"sector per-capita factors: hydrocarbons+mining {sector_gas:.2f}, personal services {sector_inf:.3f}, "
      f"ratio {sec_ratio:.1f}")
chk("sector ratio", sec_ratio, 29, 0)
w_vendor = 0.57 * 0.58                         # Chapter 18 illustrative: indigenous woman / non-indigenous man
print(f"within factor: indigenous woman / non-indigenous man = {w_vendor:.3f} (rounded 0.33); 1/0.33 = {1/0.33:.2f}")
r_exec = sec_ratio / 0.33
print(f"executive / vendor = {sec_ratio:.1f} x {1/0.33:.2f} = {r_exec:.1f}")
chk("executive/vendor", r_exec, 88, 0)
print(f"vendor / executive per-capita power = {1/r_exec:.4f} ({100/r_exec:.1f}%)")
chk("vendor/executive %", 100 / r_exec, 1.1, 1)
country = (piG_deu / L_deu) / (piG_bol / L_bol)
print(f"country per-capita factor Germany/Bolivia = {piG_deu/piG_bol:.1f} x {L_bol/L_deu:.3f} = {country:.1f}")
chk("country factor", country, 21.2, 1)
chk("1/sector_inf", 1 / sector_inf, 5.05, 2)
r_deu = country * (1.0 / sector_inf) / 0.33
print(f"German engineer (sector factor 1, within 1) / vendor = {country:.1f} x {1/sector_inf:.2f} x {1/0.33:.2f} = {r_deu:.0f}")
chk("engineer/vendor", r_deu, 325, 0)
r_noL = piG_deu / piG_bol / sector_inf / 0.33
print(f"  without dividing by workforce: {piG_deu/piG_bol:.0f} x {1/sector_inf:.2f} x {1/0.33:.2f} = {r_noL:.0f}")
chk("without workforce (hundreds)", r_noL, 2400, -2)

# --- vicious circle (Example ch19:ex:vicious) and Prebisch box ---
eci = {(row["country"], row["year"]): float(row["eci_hs92"]) for row in rd("eci_facts.csv") if row["eci_hs92"]}
chk("ECI 2019", eci[("BOL", "2019")], -1.30, 2); chk("ECI 2006", eci[("BOL", "2006")], -1.21, 2)
chk("ECI 2018", eci[("BOL", "2018")], -1.46, 2)
chk("Gamma 2006", float(_c17["summ"](2006)["Gamma"]), 0.518, 3)
hc_exp = FACTS["BOL_hydrocarbons_export_share_of_output_2019_b2017"]
chk("hydro exported", hc_exp, 70, 0); chk("pi hydro", p17[sec17.index("Hydrocarbons")], 0.039, 3)
chk("pi_URY", GPI["URY"], 0.00055, 5)
print(f"ECI {eci[('BOL','2006')]:+.3f} (2006) -> {eci[('BOL','2018')]:+.3f} (2018) -> {eci[('BOL','2019')]:+.3f} (2019); "
      f"pi_URY = {GPI['URY']:.5f}")
print(f"Gamma_N / Gamma_G = {Gam_N/Gam_G:.2f}")
chk("Gamma_N/Gamma_G", Gam_N / Gam_G, 1.03, 2)
chk("Gamma 1990", float(_c17["summ"](1990)["Gamma"]), 0.567, 3); chk("Gamma 2013", float(_c17["summ"](2013)["Gamma"]), 0.484, 3)
chk("vendor share 1/88", r_exec, 88, 0)
for _k, _v in {"dom_inputs_per_output": 0.225, "comp_share_of_va": 6.0}.items():
    assert FACTS[f"BOL_hydrocarbons_{_k}_2019_b2017"] == _v
chk("USA Gamma", float(_c17["oecd"]("USA", 2019)["Gamma"]), 0.463, 3)

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
