"""Diagnostics of docs/conventions.md applied to real input-output data.

CONVENTIONS (read this before using the outputs)
------------------------------------------------
Data sources publish intermediate flows Z with ROWS = SUPPLIERS (seller s,
buyer s'): Z[s, s'] = purchases by s' of s's output.  Total output x.

* Statistical (data) coefficients:  A_data[s, s'] = Z[s, s'] / x[s']
  (column s' = recipe of s'),  Leontief inverse L = (I - A_data)^-1.
* Book convention (Chapter 1):  a_ij = input of j per unit of i, row i is
  i's recipe, so  A_book = A_data'  and the book's L = (I - A_book')^-1
  = (I - A_data)^-1  -- the SAME matrix L.  Linkages (Rasmussen, raw sums):
      BL_j = sum_i L[i, j]   (column sums;  code/ch17.py: L.sum(0))
      FL_i = sum_j L[i, j]   (row sums;     code/ch17.py: L.sum(1))
  They are computed from A only, never from W^P.  mean(BL) = mean(FL).
  Normalised Rasmussen indices BL_j / mean(BL), FL_i / mean(FL) are also given.
* W^P = D^-1 A_book, D = diag(A_book 1):  W^P[s, s'] = Z[s', s] / sum_s'' Z[s'', s]
  = input-cost share of s' in s's intermediate purchases (row-stochastic).
* Sales-share alternative B[s, s'] = Z[s, s'] / sum_s'' Z[s, s''].
* pi = left Perron vector of W (pi' W = pi', sum 1).  Gamma = 1 - |lambda_2|.
  N_eff = 1 / sum pi^2.  Power Gini as in code/ch06.py.
  d_pi = 0.5 * ||pi - 1/N||_1.  Var(pi) = mean((pi - 1/N)^2).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def row_normalize(K: np.ndarray) -> np.ndarray:
    s = K.sum(1, keepdims=True)
    if np.any(s <= 0):
        raise ValueError("row with zero sum: drop or aggregate that sector first")
    return K / s


def wp_input_shares(Z: np.ndarray, zero_diag: bool = False) -> np.ndarray:
    """W^P (book): row s = shares of s's intermediate inputs by supplier."""
    K = Z.T.copy().astype(float)
    if zero_diag:
        np.fill_diagonal(K, 0.0)
    return row_normalize(K)


def sales_shares(Z: np.ndarray, zero_diag: bool = False) -> np.ndarray:
    K = Z.copy().astype(float)
    if zero_diag:
        np.fill_diagonal(K, 0.0)
    return row_normalize(K)


def spectrum(W: np.ndarray):
    ev = np.linalg.eigvals(W)
    order = np.argsort(-np.abs(ev))
    return ev[order]


def gap(W: np.ndarray) -> float:
    ev = spectrum(W)
    return float(1 - np.abs(ev[1]))


def perron(W: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(W.T)
    k = np.argmax(vals.real)
    p = np.abs(vecs[:, k].real)
    return p / p.sum()


def neff(p):
    return float(1 / np.sum(np.asarray(p) ** 2))


def gini(p):
    p = np.sort(np.asarray(p, float))
    N = len(p)
    return float(np.sum((2 * np.arange(1, N + 1) - N - 1) * p) / (N * p.sum()))


def diagnostics(W: np.ndarray, labels=None) -> dict:
    ev = spectrum(W)
    p = perron(W)
    N = len(p)
    l2 = ev[1]
    absl2 = float(abs(l2))
    return {
        "N": N,
        "lambda2_re": float(l2.real), "lambda2_im": float(l2.imag), "abs_lambda2": absl2,
        "abs_lambda3": float(abs(ev[2])) if N > 2 else np.nan,
        "Gamma": 1 - absl2,
        "half_life": float(np.log(2) / -np.log(absl2)) if 0 < absl2 < 1 else np.nan,
        "Neff": neff(p), "Neff_over_N": neff(p) / N,
        "power_gini": gini(p),
        "d_pi": float(0.5 * np.abs(p - 1 / N).sum()),
        "var_pi": float(np.mean((p - 1 / N) ** 2)),
        "pi_max": float(p.max()),
        "pi_argmax": (labels[int(p.argmax())] if labels is not None else int(p.argmax())),
        "pi": p,
    }


def leontief_linkages(Z: np.ndarray, x: np.ndarray, labels) -> pd.DataFrame:
    """Rasmussen linkages from A (data convention; identical L to the book's)."""
    x = np.asarray(x, float)
    A = Z / x[None, :]
    rho = max(abs(np.linalg.eigvals(A)))
    if rho >= 1:
        raise ValueError(f"spectral radius of A = {rho} >= 1 (Hawkins-Simon fails)")
    L = np.linalg.inv(np.eye(len(x)) - A)
    BL, FL = L.sum(0), L.sum(1)
    # Ghosh (supply-driven) forward linkage, for reference only
    B = Z / x[:, None]
    G = np.linalg.inv(np.eye(len(x)) - B)
    df = pd.DataFrame({
        "sector": labels,
        "output": x,
        "output_share": x / x.sum(),
        "BL": BL, "FL": FL,
        "BL_index": BL / BL.mean(), "FL_index": FL / FL.mean(),
        "FL_ghosh": G.sum(1), "FL_ghosh_index": G.sum(1) / G.sum(1).mean(),
        "intermediate_input_share": Z.sum(0) / x,
    })
    df["type"] = np.select(
        [(df.BL_index >= 1) & (df.FL_index >= 1), (df.BL_index >= 1), (df.FL_index >= 1)],
        ["key sector", "backward", "forward"], default="weak")
    df.attrs["rho_A"] = float(rho)
    return df


def aggregate(Z: np.ndarray, x: np.ndarray, labels, mapping: dict, order=None,
              extra: dict | None = None):
    """Sum flows into groups. mapping: label -> group. extra: name -> vector to sum too."""
    groups = order or sorted(set(mapping[l] for l in labels))
    G = np.zeros((len(groups), len(labels)))
    for k, l in enumerate(labels):
        G[groups.index(mapping[l]), k] = 1.0
    Zg = G @ Z @ G.T
    xg = G @ x
    ex = {k: G @ np.asarray(v, float) for k, v in (extra or {}).items()}
    return Zg, xg, groups, ex


def full_report(Z, x, labels, tag: dict) -> tuple[list[dict], pd.DataFrame, pd.DataFrame]:
    """Run every variant on one (Z, x). Returns summary rows, linkages, pi table."""
    rows, pis = [], []
    variants = [
        ("input_shares", False, wp_input_shares),
        ("input_shares_zero_diag", True, wp_input_shares),
        ("sales_shares", False, sales_shares),
        ("sales_shares_zero_diag", True, sales_shares),
    ]
    for name, zd, f in variants:
        # drop (iteratively) sectors with no intermediate purchases or no
        # intermediate sales on this variant's matrix: their W rows are undefined
        keep = np.ones(len(labels), bool)
        Zv = Z.astype(float).copy()
        if zd:
            np.fill_diagonal(Zv, 0.0)
        while True:
            Zs = Zv[np.ix_(keep, keep)]
            bad = (Zs.sum(0) <= 0) | (Zs.sum(1) <= 0)
            if not bad.any():
                break
            idx = np.where(keep)[0]
            keep[idx[bad]] = False
        Zk, lab_k = Z[np.ix_(keep, keep)], [l for l, k in zip(labels, keep) if k]
        W = f(Zk, zero_diag=zd)
        d = diagnostics(W, lab_k)
        p = d.pop("pi")
        rows.append({**tag, "normalization": name, "dropped_sectors": ";".join(
            l for l, k in zip(labels, keep) if not k), **d})
        pis.append(pd.DataFrame({**tag, "normalization": name, "sector": lab_k, "pi": p}))
    link = leontief_linkages(Z, x, labels)
    for k, v in tag.items():
        link.insert(0, k, v)
    return rows, link, pd.concat(pis, ignore_index=True)
