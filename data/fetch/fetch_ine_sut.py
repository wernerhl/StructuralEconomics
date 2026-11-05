"""INE Bolivia supply-use tables (Cuadros de Oferta y Utilizacion, COU).

Two vintages, both at CURRENT prices:

1. Base 2017 national accounts ("Cambio de Ano de Referencia"),
   https://www.ine.gob.bo/referencia2017/cuadros_oferta_utilizacion.html
   Cuadros 4.2.1-4.2.16: supply (odd) and use (even) tables 2017-2024,
   42 products x 39 industries ("2 digitos"), millions of Bs.  2023-2024 are
   marked (P) preliminary in the files.
2. Base 1990 "Matrices de Insumo Producto" (INE web, Estadisticas economicas >
   PIB y cuentas nacionales > Matrices > Matrices de insumo producto,
   current prices, thousands of Bs), 1988-2009 and 2011-2014: one sheet with
   the production (make/supply) matrix 35 products x 35 activities and the use
   matrix at purchasers' prices, plus FISIM ("imputacion bancaria") dummy.
   https://www.ine.gob.bo/index.php/estadisticas-economicas/pib-y-cuentas-nacionales/matrices/matrices-de-insumo-producto/

Output (tidy, long): data/processed/ine_sut_long.csv with columns
  base, year, table (supply|use), row_code, row_label, col_code, col_label,
  value, unit
All raw xlsx in data/raw/ine_cou/ and data/raw/ine_mip90/ (each < 0.2 MB).

Run: python3 data/fetch/fetch_ine_sut.py
"""
from __future__ import annotations

import re
import warnings

import openpyxl
import pandas as pd

from common import PROC, RAW, download

warnings.filterwarnings("ignore", category=UserWarning)

B17 = "https://www.ine.gob.bo/referencia2017/CUADROS/pagina_web/4.2.{}.xlsx"
# odd -> supply, even -> use; pairs (2k-1, 2k) for year 2016+k
B17_YEARS = {2017 + k: (2 * k + 1, 2 * k + 2) for k in range(8)}

B90 = ("https://www.ine.gob.bo/index.php/descarga/406/matrices-de-insumo-producto-corriente/"
       "{id}/matriz-de-insumo-producto-en-miles-de-bolivianos{y}.xlsx")
B90_IDS = {y: 44576 + i for i, y in enumerate(list(range(1988, 2010)) + [2011, 2012, 2013, 2014])}


def _s(c):
    return re.sub(r"\s+", " ", str(c)).strip() if c is not None else ""


def parse_b17(path, year, table):
    ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
    R = [list(r) for r in ws.iter_rows(values_only=True)]
    title = _s(R[7][1])
    assert str(year) in title, (title, year)
    unit = _s(R[8][1])
    codes, l12, l13 = R[10], R[11], R[12]
    cols = {}
    group = ""
    for j in range(3, len(codes)):
        lab = _s(l13[j]) or _s(l12[j])
        if _s(l12[j]) and not _s(l13[j]):
            group = _s(l12[j])
        if not lab:
            continue
        code = _s(codes[j]) if _s(codes[j]) and re.match(r"^\d+$", _s(codes[j])) else ""
        cols[j] = (code, lab)
    out = []
    for r in R[13:]:
        rc, rl = _s(r[1]), _s(r[2])
        if not rl or rl.startswith("Fuente") or rc.startswith("Fuente") or rc.startswith("("):
            continue
        for j, (cc, cl) in cols.items():
            v = r[j] if j < len(r) else None
            if isinstance(v, (int, float)):
                out.append(dict(base=2017, year=year, table=table, row_code=rc, row_label=rl,
                                col_code=cc, col_label=cl, value=float(v), unit=unit))
    return out


PROD1 = re.compile(r"^\s*0?1\.\s")


def _parse_b90_sheet(ws, year):
    """Generic parser for the base-1990 layouts.

    A sheet holds one or more blocks; each block has a label column whose
    cells contain product labels ' 1. PRODUCTOS ...', '35. SERVICIOS ...'
    followed by total/VA rows.  Value columns of a block lie to the right of its
    label column up to the next label column; in the 'MIP' layout the columns
    to the LEFT of the first label column (VBP, imports, taxes, margins, total
    supply by product) also belong to that label column.
    Blocks are classified as:
      make   = production matrix (products x activities, basic prices) [MP sheets
               or the left block of the 2011-2014 combined sheet]
      supply = product-level supply aggregates (VBP, imports, taxes, margins)
      use    = use matrix at purchasers' prices + final demand + VA rows
    """
    R = [list(r) for r in ws.iter_rows(values_only=True)]
    width = max(len(r) for r in R)
    R = [r + [None] * (width - len(r)) for r in R]
    unit = _s(R[1][0])
    i0 = next(i for i, r in enumerate(R)
              if any(isinstance(c, str) and PROD1.match(_s(c)) for c in r)
              and sum(isinstance(c, (int, float)) for c in r) >= 5)
    labcols = [j for j, c in enumerate(R[i0]) if isinstance(c, str) and PROD1.match(_s(c))]
    hdr = R[2:i0]

    def clab(j):
        return " / ".join(_s(h[j]) for h in hdr if _s(h[j]))

    out = []
    for b, lc in enumerate(labcols):
        right_end = labcols[b + 1] if b + 1 < len(labcols) else width
        blocks = []
        if b == 0 and lc > 0:
            blocks.append(("left", range(0, lc)))
        blocks.append(("right", range(lc + 1, right_end)))
        for side, cols in blocks:
            labels = " ".join(clab(j) for j in cols).upper()
            if "CONSUMO INTERMEDIO" in labels or "DEMANDA" in labels:
                table = "use"
            elif side == "left" or "OFERTA TOTAL" in labels or "IMPORTACI" in labels:
                table = "supply"
            else:
                table = "make"
            for i in range(i0, len(R)):
                rl = _s(R[i][lc])
                if not rl or rl.upper().startswith(("FUENTE", "(P)", "PIB", "P I B", "NOTA")):
                    continue
                rm = re.match(r"^(\d+)\.", rl)
                for j in cols:
                    v = R[i][j]
                    cl = clab(j)
                    if isinstance(v, (int, float)) and cl and not cl.upper().startswith("PRODUCTOS - ACTIV"):
                        cm = re.match(r"^(\d+)\.?(\s|$)", cl)
                        out.append(dict(base=1990, year=year, table=table,
                                        row_code=rm.group(1) if rm else "", row_label=rl,
                                        col_code=cm.group(1) if cm else "", col_label=cl,
                                        value=float(v), unit=unit))
    return out


def parse_b90(path, year):
    wb = openpyxl.load_workbook(path, data_only=True)
    names = [w.title for w in wb.worksheets]
    # current-price sheets only: 'MIPyyyy' and 'MPyyyy' (the 'MIPyy90'/'MPyy90'
    # sheets in the 2002/2003/2006 files are at constant 1990 prices)
    keep = [n for n in names if n in (f"MIP{year}", f"MP{year}")]
    out = []
    for n in keep:
        out += _parse_b90_sheet(wb[n], year)
    if not keep:
        raise ValueError(f"no current-price sheet in {names}")
    return out


def main():
    rows = []
    d17 = RAW / "ine_cou"
    for y, (s, u) in B17_YEARS.items():
        ps = download(B17.format(s), d17 / f"4.2.{s}.xlsx")
        pu = download(B17.format(u), d17 / f"4.2.{u}.xlsx")
        rows += parse_b17(ps, y, "supply") + parse_b17(pu, y, "use")
        print("base2017", y, "ok")
    d90 = RAW / "ine_mip90"
    for y, i in B90_IDS.items():
        p = download(B90.format(id=i, y=y), d90 / f"mip_corr_{y}.xlsx")
        try:
            rows += parse_b90(p, y)
            print("base1990", y, "ok")
        except Exception as e:  # noqa: BLE001
            print("base1990", y, "PARSE FAILED", e)
    df = pd.DataFrame(rows)
    df.to_csv(PROC / "ine_sut_long.csv", index=False)
    print(df.groupby(["base", "table"]).size())


if __name__ == "__main__":
    main()
