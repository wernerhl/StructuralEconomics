"""OECD Inter-Country Input-Output (ICIO) tables, 2023 edition, "SML" version
(76 economies + rest of world, 45 ISIC Rev.4 industries, 1995-2020, USD million,
basic prices).  Cite: OECD (2023), OECD Inter-Country Input-Output Database,
http://oe.cd/icio .

The www.oecd.org landing page is behind a Cloudflare browser challenge from this
machine (HTTP 403), but the file server stats.oecd.org/wbos is reachable.  The
five year-block file IDs below are those listed by the iotr R package
(https://rdrr.io/github/okrebs/iotr/man/io_dl_icio_v2.html); each was verified
by its Content-Disposition file name (e.g. "2016-2020_SML.zip").

Bolivia is NOT a separate economy in ICIO (it is inside ROW).

Raw zips (~85 MB each) -> data/raw/oecd_icio/  (git-ignored; re-downloaded)
Processed: data/processed/oecd_icio_domestic/<ISO3>.csv.gz : domestic block
  Z_cc (rows = supplying industry, cols = using industry) and gross output, long
  format (year, row, col, value) for the comparator economies below; plus
  data/processed/oecd_icio_output.csv.gz (year, country, industry, output).

Run: python3 data/fetch/fetch_oecd_icio.py [--blocks 2016-2020 ...]
"""
from __future__ import annotations

import argparse
import io
import zipfile

import numpy as np
import pandas as pd

from common import PROC, RAW, download

BASE = "https://stats.oecd.org/wbos/fileview2.aspx?IDFile={}"
BLOCKS = {
    "1995-2000": "d26ad811-5b58-4f0c-a4e3-06a1469e475c",
    "2001-2005": "7cb93dae-e491-4cfd-ac67-889eb7016a4a",
    "2006-2010": "ea165bfb-3a85-4e0a-afee-6ba8e6c16052",
    "2011-2015": "1f791bc6-befb-45c5-8b34-668d08a1702a",
    "2016-2020": "d1ab2315-298c-4e93-9a81-c6f2273139fe",
}
# 2021 edition (1995-2018, 66 economies + ROW, 45 industries with numeric ISIC
# codes, e.g. "01T02"; China and Mexico split into CN1/CN2, MX1/MX2), used only
# to extend trajectories before 2011.  IDs from pymrio 0.7.0
# (pymrio/tools/iodownloader.py); verified by Content-Disposition
# ("ICIO_1995-1999.zip").  Do not mix editions in one series without a break.
BLOCKS_2021 = {
    "1995-1999": "91d8e84b-7406-46b9-af5f-ec096242755c",
    "2000-2004": "8adf89dd-18b4-40fe-bc7f-c822052eb961",
    "2005-2009": "fe218690-0a3b-44aa-a82c-b3e3da6d24db",
    "2010-2014": "2c2f499f-5703-4034-9457-2f7518e8f2fc",
    "2015-2018": "59a3d7f2-3f23-40d5-95ca-48da84c0f861",
}
SPLITS = {"CHN": ["CN1", "CN2"], "MEX": ["MX1", "MX2"]}
COMPARATORS = ["KOR", "DEU", "USA", "CHL", "CHN", "MEX", "BRA", "ARG", "COL", "PER",
               "JPN", "NOR", "CRI", "IND", "IDN", "ZAF", "ESP", "VNM", "THA", "MYS", "POL", "CZE"]


def read_year(zf: zipfile.ZipFile, name: str) -> pd.DataFrame:
    with zf.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f), index_col=0)
    return df


def domestic_blocks(df: pd.DataFrame, countries):
    rows = df.index.astype(str)
    cols = df.columns.astype(str)
    out_key = next((k for k in ("OUT", "OUTPUT") if k in df.index), None)
    out_row = df.loc[out_key] if out_key else None
    res = {}
    for c in countries:
        colset = set(cols)
        ri = [r for r in rows if r.startswith(c + "_") and r in colset]  # drops e.g. ARG_TAXSUB
        if not ri and c in SPLITS:  # 2021 edition: sum the split blocks
            parts = SPLITS[c]
            r0 = [r for r in rows if r.startswith(parts[0] + "_") and r in colset]
            if not r0:
                continue
            inds = [r.split("_", 1)[1] for r in r0]
            Z = sum(df.loc[[f"{a}_{i}" for i in inds], [f"{b}_{i}" for i in inds]].values.astype(float)
                    for a in parts for b in parts)
            x = sum(out_row[[f"{a}_{i}" for i in inds]].values.astype(float) for a in parts)
            res[c] = (inds, Z, x)
            continue
        if not ri:
            continue
        inds = [r.split("_", 1)[1] for r in ri]
        ci = [f"{c}_{i}" for i in inds]  # industry columns only (drop final demand)
        Z = df.loc[ri, ci].values.astype(float)
        x = out_row[ci].values.astype(float) if out_row is not None else None
        res[c] = (inds, Z, x)
    return res


def all_country_codes(df):
    cs = {r.split("_", 1)[0] for r in df.index.astype(str) if "_" in r}
    cs |= {c for c, parts in SPLITS.items() if parts[0] in cs}
    return sorted(cs - {p for ps in SPLITS.values() for p in ps})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edition", default="2023", choices=["2023", "2021"])
    ap.add_argument("--blocks", nargs="*", default=None)
    a = ap.parse_args()
    blocks = BLOCKS if a.edition == "2023" else BLOCKS_2021
    a.blocks = a.blocks or list(blocks)
    sfx = "" if a.edition == "2023" else "2021"
    raw = RAW / "oecd_icio"
    outdir = PROC / f"oecd_icio{sfx}_domestic"
    outdir.mkdir(parents=True, exist_ok=True)
    long = {c: [] for c in COMPARATORS}
    outputs = []
    failed = []
    for blk in a.blocks:
        try:
            fname = f"{blk}_SML.zip" if a.edition == "2023" else f"ICIO2021_{blk}.zip"
            p = download(BASE.format(blocks[blk]), raw / fname, timeout=1200)
        except Exception as e:  # noqa: BLE001
            print("DOWNLOAD FAILED", blk, e)
            failed.append(blk)
            continue
        with zipfile.ZipFile(p) as zf:
            for name in sorted(n for n in zf.namelist() if n.endswith(".csv")):
                year = int(next(t for t in name.replace(".", "_").split("_") if t.isdigit() and len(t) == 4))
                df = read_year(zf, name)
                dblocks = domestic_blocks(df, all_country_codes(df))
                for c, (inds, Z, x) in dblocks.items():
                    outputs += [dict(year=year, country=c, industry=i, output=v) for i, v in zip(inds, x)]
                    if c in long:
                        r, k = np.nonzero(Z)
                        long[c].append(pd.DataFrame({"year": year, "row": np.array(inds)[r],
                                                     "col": np.array(inds)[k], "value": Z[r, k]}))
                print(blk, name, "countries:", len(dblocks))
    for c, parts in long.items():
        if parts:
            f = outdir / f"{c}.csv.gz"
            new = pd.concat(parts)
            if f.exists():  # merge with blocks fetched earlier
                old = pd.read_csv(f)
                new = pd.concat([old[~old.year.isin(new.year.unique())], new])
            new.sort_values(["year", "row", "col"]).to_csv(f, index=False)
    print("failed blocks:", failed)
    fo = PROC / f"oecd_icio{sfx}_output.csv.gz"
    new = pd.DataFrame(outputs)
    if fo.exists():
        old = pd.read_csv(fo)
        new = pd.concat([old[~old.year.isin(new.year.unique())], new])
    new.to_csv(fo, index=False)


if __name__ == "__main__":
    main()
