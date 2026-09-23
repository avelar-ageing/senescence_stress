"""Compare two rerun_outputs directories, file by file, for every output in pipeline.tsv.

Usage (from bulk_RNA_seq/):
    .venv/bin/python env/compare_outputs.py <reference_rerun_dir> <new_rerun_dir> [--rtol 1e-6] [--atol 1e-9]

CSV: same columns and rows in the same order; numeric columns equal within
|a - b| <= atol + rtol * |b|; other columns identical (NA == NA).
Other files (.png, .rds): md5 only, reported but never counted as a failure, because
PNG rendering depends on fonts and RDS files embed R-version metadata.
Exit code 1 if any CSV differs or is missing.
"""
import argparse
import hashlib
import os
import sys

import numpy as np
import pandas as pd


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def outputs(tsv):
    for ln in open(tsv):
        if ln.startswith("#") or ln.startswith("stage\t") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        for o in f[3].split():
            yield f[1], o


def compare_csv(a, b, rtol, atol):
    x, y = pd.read_csv(a, low_memory=False), pd.read_csv(b, low_memory=False)
    if list(x.columns) != list(y.columns):
        return f"columns differ: {sorted(set(x.columns) ^ set(y.columns))}"
    if len(x) != len(y):
        return f"rows {len(x)} vs {len(y)}"
    bad = []
    for c in x.columns:
        u, v = x[c], y[c]
        if pd.api.types.is_numeric_dtype(u) and pd.api.types.is_numeric_dtype(v):
            u, v = u.to_numpy(float), v.to_numpy(float)
            same = (np.isnan(u) & np.isnan(v)) | (np.abs(u - v) <= atol + rtol * np.abs(v))
            if not same.all():
                i = np.flatnonzero(~same)
                bad.append(f"{c}: {len(i)} values, max |diff| {np.nanmax(np.abs(u[i] - v[i])):.3g}")
        else:
            same = (u.isna() & v.isna()) | (u.astype(str) == v.astype(str))
            if not same.all():
                bad.append(f"{c}: {int((~same).sum())} values")
    return "; ".join(bad)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    ap.add_argument("--pipeline", default="pipeline.tsv")
    a = ap.parse_args()
    fail = 0
    for script, o in outputs(a.pipeline):
        r, n = os.path.join(a.ref, o), os.path.join(a.new, o)
        if not os.path.exists(n) or not os.path.exists(r):
            print(f"MISSING  {o}  ({'new' if not os.path.exists(n) else 'reference'})  [{script}]")
            fail += o.endswith(".csv")
            continue
        if o.endswith(".csv"):
            d = compare_csv(r, n, a.rtol, a.atol)
            print(f"{'DIFFERS' if d else 'same   '}  {o}" + (f"  {d}  [{script}]" if d else ""))
            fail += bool(d)
        else:
            s = md5(r) == md5(n)
            print(f"{'same   ' if s else 'md5 ≠  '}  {o}" + ("" if s else "  (not a failure: binary/figure)"))
    print(f"== {fail} CSV file(s) differ or are missing")
    sys.exit(1 if fail else 0)


if __name__ == "__main__":
    main()
