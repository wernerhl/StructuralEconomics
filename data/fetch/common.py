"""Shared helpers for the data/fetch scripts.

All downloads go to data/raw/<source>/, tidy outputs to data/processed/.
Every download is logged (URL, UTC time, HTTP status, bytes, sha256) in
data/raw/download_log.csv so that the README can cite exact retrieval details.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import os
import time
from pathlib import Path

import requests

DATA = Path(__file__).resolve().parents[1]
RAW = DATA / "raw"
PROC = DATA / "processed"
RAW.mkdir(exist_ok=True)
PROC.mkdir(exist_ok=True)
LOG = RAW / "download_log.csv"

# The session proxy re-terminates TLS; requests honours REQUESTS_CA_BUNDLE /
# SSL_CERT_FILE if set.  Fall back to the proxy bundle when present.
_CA = os.environ.get("REQUESTS_CA_BUNDLE") or os.environ.get("SSL_CERT_FILE")
if not _CA and Path("/root/.ccr/ca-bundle.crt").exists():
    _CA = "/root/.ccr/ca-bundle.crt"

UA = {"User-Agent": "StructuralEconomics-data-pass/1.0 (research; python-requests)"}


def _log(url: str, status, path: Path | None, nbytes: int, note: str = "") -> None:
    new = not LOG.exists()
    sha = ""
    if path is not None and path.exists():
        sha = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    with LOG.open("a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(["utc_time", "url", "status", "bytes", "sha256_16", "local_path", "note"])
        w.writerow([
            dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            url, status, nbytes, sha,
            str(path.relative_to(DATA)) if path is not None else "", note,
        ])


def download(url: str, dest: Path, *, timeout: int = 180, retries: int = 3,
             force: bool = False, params: dict | None = None) -> Path:
    """Download url -> dest (skips if dest exists and force is False)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not force:
        return dest
    last = None
    for k in range(retries):
        try:
            with requests.get(url, params=params, headers=UA, timeout=timeout,
                              stream=True, verify=_CA or True) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                n = 0
                with tmp.open("wb") as f:
                    for chunk in r.iter_content(1 << 16):
                        f.write(chunk)
                        n += len(chunk)
                tmp.replace(dest)
                _log(r.url, r.status_code, dest, n)
                return dest
        except Exception as e:  # noqa: BLE001
            last = e
            _log(url, f"ERR:{type(e).__name__}", None, 0, str(e)[:200])
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"download failed: {url}: {last}")


def get_json(url: str, *, params: dict | None = None, timeout: int = 120, retries: int = 3):
    last = None
    for k in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout, verify=_CA or True)
            r.raise_for_status()
            _log(r.url, r.status_code, None, len(r.content), "json")
            return r.json()
        except Exception as e:  # noqa: BLE001
            last = e
            _log(url, f"ERR:{type(e).__name__}", None, 0, str(e)[:200])
            time.sleep(3 * (k + 1))
    raise RuntimeError(f"json fetch failed: {url}: {last}")
