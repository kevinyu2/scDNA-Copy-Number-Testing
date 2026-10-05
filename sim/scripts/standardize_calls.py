#!/usr/bin/env python3
"""
Convert a caller's output to the common calls table that eval_calls.py reads.

Output columns (TSV, one row per cell x bin):
  cell     cell name without the _sim suffix (batchN_cellM)
  chrom    chr-prefixed
  start    bin start (the caller's own coordinates)
  end      bin end
  CN_A     copy number of haplotype A
  CN_B     copy number of haplotype B
  plus caller-specific extras (kept for later metrics, ignored by eval_calls.py
  unless named there):
    hiscanner: gamma
    chisel:    cluster, rdr, baf, a_count, b_count, clone (from clones/mapping.tsv)

Usage
  standardize_calls.py hiscanner --calls-dir <hiscanner>/output/final_calls --out calls.tsv
  standardize_calls.py chisel --calls <chisel>/calls/calls.tsv \
      --barcodes <prep>/barcodedcells.info.tsv [--clones <chisel>/clones/mapping.tsv] --out calls.tsv
"""
import argparse
import glob
import os
import sys

import pandas as pd

COLS = ["cell", "chrom", "start", "end", "CN_A", "CN_B"]


def norm_cell(name):
    name = str(name)
    return name[:-4] if name.endswith("_sim") else name


def norm_chrom(c):
    c = str(c)
    return c if c.startswith("chr") else f"chr{c}"


def from_hiscanner(a):
    skip = ("_seg.txt", "_seg_merged.txt", "_seg_initial_anno.txt")
    files = sorted(f for f in glob.glob(os.path.join(a.calls_dir, "*.txt")) if not f.endswith(skip))
    if not files:
        sys.exit(f"ERROR: no HiScanner call files in {a.calls_dir}")
    parts, empty = [], []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        if not {"CHROM", "START", "END", "CN_A", "CN_B"}.issubset(df.columns):
            empty.append(os.path.basename(f)); continue
        df = df.dropna(subset=["CN_A", "CN_B"])
        if df.empty:
            empty.append(os.path.basename(f)); continue
        out = pd.DataFrame({
            "cell": norm_cell(os.path.basename(f)[:-4]),
            "chrom": df.CHROM.map(norm_chrom),
            "start": df.START.astype("int64"), "end": df.END.astype("int64"),
            "CN_A": df.CN_A.astype(int), "CN_B": df.CN_B.astype(int),
            "gamma": df["gamma"] if "gamma" in df.columns else float("nan"),
        })
        parts.append(out)
    if empty:
        print(f"[standardize_calls] {len(empty)} HiScanner files without calls, e.g. {empty[:3]}", file=sys.stderr)
    if not parts:
        sys.exit("ERROR: no usable HiScanner calls")
    return pd.concat(parts, ignore_index=True)


def from_chisel(a):
    calls = pd.read_csv(a.calls, sep="\t")
    calls.columns = [c.lstrip("#") for c in calls.columns]
    need = {"CHR", "START", "END", "CELL", "CN_STATE"}
    if not need.issubset(calls.columns):
        sys.exit(f"ERROR: {a.calls} lacks columns {sorted(need - set(calls.columns))}")

    bc = pd.read_csv(a.barcodes, sep="\t")
    bc.columns = [c.lstrip("#") for c in bc.columns]
    if bc.BARCODE.duplicated().any():
        sys.exit(f"ERROR: duplicated barcodes in {a.barcodes}")
    to_cell = dict(zip(bc.BARCODE.astype(str), bc.CELL.map(norm_cell)))

    calls["CELL"] = calls.CELL.astype(str)
    unknown = sorted(set(calls.CELL) - set(to_cell))
    if unknown:
        sys.exit(f"ERROR: {len(unknown)} CHISEL barcodes not in {a.barcodes}, e.g. {unknown[:3]}")

    ab = calls.CN_STATE.astype(str).str.split("|", expand=True)
    if ab.shape[1] != 2:
        sys.exit(f"ERROR: unexpected CN_STATE values in {a.calls}, e.g. {calls.CN_STATE.iloc[0]}")

    out = pd.DataFrame({
        "cell": calls.CELL.map(to_cell),
        "chrom": calls.CHR.map(norm_chrom),
        "start": calls.START.astype("int64"), "end": calls.END.astype("int64"),
        "CN_A": ab[0].astype(int), "CN_B": ab[1].astype(int),
    })
    for src, dst in (("CLUSTER", "cluster"), ("RDR", "rdr"), ("BAF", "baf"), ("A_COUNT", "a_count"), ("B_COUNT", "b_count")):
        if src in calls.columns:
            out[dst] = calls[src].values

    if a.clones and os.path.exists(a.clones) and os.path.getsize(a.clones) > 0:
        cl = pd.read_csv(a.clones, sep="\t")
        cl.columns = [c.lstrip("#") for c in cl.columns]
        if {"CELL", "CLONE"}.issubset(cl.columns):
            out["clone"] = out.cell.map(dict(zip(cl.CELL.astype(str).map(to_cell), cl.CLONE)))

    n_bc = len(bc)
    print(f"[standardize_calls] CHISEL: {out.cell.nunique()}/{n_bc} prepped cells have calls "
          f"(the rest failed CHISEL's minreads filter)", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="caller", required=True)
    h = sub.add_parser("hiscanner")
    h.add_argument("--calls-dir", required=True)
    h.add_argument("--out", required=True)
    c = sub.add_parser("chisel")
    c.add_argument("--calls", required=True)
    c.add_argument("--barcodes", required=True)
    c.add_argument("--clones")
    c.add_argument("--out", required=True)
    a = ap.parse_args()

    df = from_hiscanner(a) if a.caller == "hiscanner" else from_chisel(a)
    df = df.sort_values(["cell", "chrom", "start"]).reset_index(drop=True)
    df = df[COLS + [c for c in df.columns if c not in COLS]]
    tmp = a.out + ".tmp"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, a.out)
    print(f"[standardize_calls] wrote {a.out}: {len(df)} rows, {df.cell.nunique()} cells", file=sys.stderr)


if __name__ == "__main__":
    main()
