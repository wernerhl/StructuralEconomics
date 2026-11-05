"""INE Bolivia: official monetary poverty, extreme poverty and Gini (Encuesta
de Hogares), from the page "Pobreza y Desigualdad - Encuestas de Hogares":
https://www.ine.gob.bo/index.php/estadisticas-economicas/encuestas-de-hogares/

Two poverty-line vintages exist and must NOT be spliced:
  * 2005-2018 tables: old basic baskets (UDAPE/INE methodology);
  * 2016-2025 tables: lines updated with the Encuesta de Presupuestos
    Familiares 2015-2016 (Nota tecnica LP01/2020); "-rev25" files use the
    2025 revision of the population projections (EPP) as expansion factors.
No survey was carried out in 2010.

Outputs:
  data/raw/ine_pobreza/*.xlsx
  data/processed/ine_poverty_gini_long.csv  (table, series_vintage, area,
      indicator, unit, year, value)
  data/processed/ine_poverty_notes.csv      (footnotes of every table, verbatim)

Run: python3 data/fetch/fetch_ine_poverty.py
"""
from __future__ import annotations

import re

import openpyxl
import pandas as pd

from common import PROC, RAW, download

NUBE = "https://nube.ine.gob.bo/index.php/s/{}/download"
FILES = {
    # name: (share id, INE table code, vintage)
    "pobreza_extrema_area_2005_2018": ("6ffM4Tt50ILU23N", "3.06.01.02", "old_lines_2005_2018"),
    "pobreza_area_2005_2018": ("RkqFdKgwXMq1Jh9", "3.06.01.01", "old_lines_2005_2018"),
    "pobreza_extrema_area_2016_2025": ("S9q946tPUVX3P8i", "3.06.01.02", "EPF2015_lines_2016_2025"),
    "pobreza_extrema_area_rev25": ("ybB63ivM10ViiJH", "3.06.01.02", "EPF2015_lines_2016_2025_rev25"),
    "pobreza_area_2016_2025": ("MwlP5BcDCs1yIg1", "3.06.01.01", "EPF2015_lines_2016_2025"),
    "pobreza_area_rev25": ("pmC0tIulettaTEG", "3.06.01.01", "EPF2015_lines_2016_2025_rev25"),
    "gini_area_2005_2025": ("HUlsrhvWNgm5hvC", "3.06.01.03", "gini_2005_2025"),
    "gini_area_rev25": ("0ap0uccAMh9zMX5", "3.06.01.03", "gini_2005_2025_rev25"),
    "percentil_2005_2025": ("WwHJTCz1GFtDAvJ", "3.06.01.04", "percentile_ratio_2005_2025"),
    "lineas_pobreza_rev25": ("bNBX2DDzNkLrxJd", "3.06.01.05", "poverty_lines_rev25"),
}


def parse(path, table, vintage):
    ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    title = " ".join(str(c) for r in rows[:3] for c in r if c)[:300]
    hdr_i = next(i for i, r in enumerate(rows) if r and isinstance(r[0], str)
                 and r[0].strip().upper().startswith(("INDICADOR", "ÁREA", "AREA", "DESCRIPCI")))
    hdr = rows[hdr_i]
    years = {}
    for j, c in enumerate(hdr):
        m = re.match(r"\s*(\d{4})", str(c)) if c is not None else None
        if m:
            years[j] = int(m.group(1))
    has_unit = hdr[1] is not None and "UNIDAD" in str(hdr[1]).upper()
    out, notes, area = [], [], "Bolivia"
    for r in rows[hdr_i + 1:]:
        if not r or r[0] is None:
            continue
        lab = str(r[0]).strip()
        vals = {years[j]: r[j] for j in years if j < len(r) and isinstance(r[j], (int, float))}
        if not vals:
            if re.match(r"(Fuente|\(|Nota|Metodolog|La Canast|El c|No se|A partir|Serie|Elab)", lab):
                notes.append(lab)
            else:
                area = lab  # section header (BOLIVIA / Área Urbana / Área Rural / dept)
            continue
        indicator = lab if has_unit else "Gini" if "gini" in vintage else lab
        sub_area = area if has_unit else lab
        for y, v in vals.items():
            out.append({"table": table, "series_vintage": vintage, "title": title,
                        "area": sub_area, "indicator": indicator,
                        "unit": (r[1] if has_unit else ""), "year": y, "value": v})
    return out, [{"table": table, "series_vintage": vintage, "note": n} for n in notes]


def main():
    d = RAW / "ine_pobreza"
    allrows, allnotes = [], []
    for name, (sid, table, vintage) in FILES.items():
        p = download(NUBE.format(sid), d / f"{name}.xlsx")
        try:
            r, n = parse(p, table, vintage)
        except Exception as e:  # noqa: BLE001
            print("PARSE FAILED", name, e)
            continue
        print(name, len(r), "values")
        allrows += r
        allnotes += n
    pd.DataFrame(allrows).to_csv(PROC / "ine_poverty_gini_long.csv", index=False)
    pd.DataFrame(allnotes).to_csv(PROC / "ine_poverty_notes.csv", index=False)


if __name__ == "__main__":
    main()
