#!/usr/bin/env python3
"""event_separability.py -- can two groups of cells be told apart on given events at all? (diagnosis)

    python3 test/event_separability.py --run RUN --barcodes BARCODES --truth TRUTH \\
        --seg mat:chr1:5056000-5551000 --seg mat:chr1:153022000-153523000 \\
        --group-a "mat:chr1:153022000-153523000=0" --group-b "pat:chr1:10008200-12008300=2" [--exclude-b ...]

For every --seg (haplotype : tree segment as in truth_hg38.tsv's seg column), the segment's lifted hg38
span is taken from the truth, and every bin of the run lying at least --min-frac inside it is used.
Per cell and segment: mean RDR over those bins, and the pooled BAF (sum B / sum (A + B)) as distance
from 0.5. Group A = cells whose CN in the given haplotype segment equals the value (all --group-a
conditions); group B likewise with --group-b, minus cells matching any --exclude-b. Reports per
segment the bins used, both groups' quantiles and the AUC of each feature (1 or 0 = perfectly
separable per cell, 0.5 = not at all), and the AUC of a combined score (the features standardized,
summed with the sign that separates). This uses the truth, so it measures what is possible, not
what a method finds.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def norm_cell(n):
    n = str(n)
    return n[:-4] if n.endswith("_sim") else n


def auc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    r = pd.Series(np.r_[pos, neg]).rank().values
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def missing(t, hap, seg):
    """Error listing the truth segments of this haplotype that overlap the requested coordinates."""
    try:
        lo, hi = (int(x) for x in seg.split(":")[1].split("-"))
        segs = sorted(set(t[t.hap == hap].seg))
        near = [x for x in segs if ":" in x and "-" in x and int(x.split(":")[1].split("-")[0]) < hi and int(x.split(":")[1].split("-")[1]) > lo]
    except ValueError:
        near = []
    sys.exit(f"ERROR: no truth rows for {hap} {seg}" + (f"; overlapping {hap} segments: {', '.join(near[:8])}" if near else ""))


def cells_with(t, conds):
    keep = None
    for c in conds:
        (hs, val) = c.rsplit("=", 1)
        hap, seg = hs.split(":", 1)
        g = t[(t.hap == hap) & (t.seg == seg)].groupby("cell").cn.max()
        if g.empty:
            missing(t, hap, seg)
        s = set(g[g == int(val)].index)
        keep = s if keep is None else keep & s
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--seg", action="append", required=True, help="hap:seg, e.g. mat:chr1:5056000-5551000 (repeatable)")
    ap.add_argument("--group-a", action="append", required=True, help="hap:seg=CN condition for group A (repeatable, all must hold)")
    ap.add_argument("--group-b", action="append", required=True, help="hap:seg=CN condition for group B (repeatable)")
    ap.add_argument("--exclude-b", action="append", default=[], help="hap:seg=CN: cells matching this are removed from group B (repeatable)")
    ap.add_argument("--min-frac", type=float, default=0.5, help="bins at least this share inside the segment (default 0.5)")
    a = ap.parse_args()

    t = pd.read_csv(a.truth, sep="\t", dtype={"cell": str})
    t["cell"] = t.cell.map(norm_cell)
    A = cells_with(t, a.group_a)
    B = cells_with(t, a.group_b) - A
    for ex in a.exclude_b:
        B -= cells_with(t, [ex])
    print(f"group A: {len(A)} cells; group B: {len(B)} cells")

    bc = pd.read_csv(a.barcodes, sep="\t")
    bc.columns = [c.lstrip("#") for c in bc.columns]
    to_cell = dict(zip(bc.BARCODE.astype(str), bc.CELL.map(norm_cell)))
    c = pd.read_csv(f"{a.run}/combo/combo.tsv", sep="\t", header=None, comment="#", usecols=[0, 1, 2, 3, 6, 7, 8],
                    names=["chrom", "s", "e", "bc", "rdr", "A", "B"], dtype={"chrom": str, "bc": str})
    c["cell"] = c.bc.map(to_cell)
    c = c[c.cell.isin(A | B)]

    feats = {}
    for spec in a.seg:
        hap, seg = spec.split(":", 1)
        rows = t[(t.hap == hap) & (t.seg == seg)]
        if rows.empty:
            missing(t, hap, seg)
        lo, hi = int(rows.start.min()), int(rows.end.max())
        ov = (np.minimum(c.e, hi) - np.maximum(c.s, lo)).clip(lower=0) / (c.e - c.s)
        sub = c[ov >= a.min_frac]
        nb = sub[["s", "e"]].drop_duplicates()
        print(f"\n== {spec}: lifted to {lo:,}-{hi:,} ({(hi - lo) / 1e6:.2f} Mb); bins >= {a.min_frac:g} inside: {len(nb)} "
              f"({(nb.e - nb.s).sum() / 1e6:.2f} Mb)")
        if nb.empty:
            continue
        g = sub.groupby("cell").agg(rdr=("rdr", "mean"), A=("A", "sum"), B=("B", "sum"))
        g["baf_dist"] = (g.B / (g.A + g.B).replace(0, np.nan) - 0.5).abs()
        g["snp_reads"] = g.A + g.B
        for f in ("rdr", "baf_dist"):
            ga, gb = g.loc[g.index.isin(A), f].dropna(), g.loc[g.index.isin(B), f].dropna()
            print(f"  {f:8s}  A quantiles {np.round(ga.quantile([.1, .5, .9]).values, 3)}  B quantiles "
                  f"{np.round(gb.quantile([.1, .5, .9]).values, 3)}  AUC {auc(ga, gb):.3f}")
            feats[f"{spec} {f}"] = g[f]
        print(f"  SNP reads per cell over these bins: median {g.snp_reads.median():.0f}")
    if len(feats) > 1:
        F = pd.DataFrame(feats).dropna()
        Fz = (F - F.mean()) / F.std()
        sign = np.sign([auc(Fz.loc[Fz.index.isin(A), col], Fz.loc[Fz.index.isin(B), col]) - 0.5 for col in Fz])
        score = (Fz * sign).sum(axis=1)
        print(f"\ncombined score of {F.shape[1]} features: AUC {auc(score[score.index.isin(A)], score[score.index.isin(B)]):.3f} "
              f"({len(F)} cells with every feature)")


if __name__ == "__main__":
    main()
