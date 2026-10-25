#!/usr/bin/env python3
"""oracle_clusters.py -- the bin clusters a perfect clustering method would find (for diagnosis).

    python3 scripts/oracle_clusters.py --truth eval/truth_hg38.tsv --run chisel_pbins8/perfect [--out DIR]

Two stretches of genome belong to the same oracle cluster when every cell has the same copy numbers
on both haplotypes in both: the profile of a stretch is (mat CN, pat CN) of every cell. Regions
without any CNA (every cell 1|1) are therefore one cluster, and two CNAs carried by the same cells
with the same copy numbers share a cluster.

The truth (lift_truth.py's truth_hg38.tsv) cuts every tree segment into the pieces that align to
hg38, separately for the maternal and paternal haplotype (their breakpoints land in different hg38
places). Per haplotype the pieces are made to tile the chromosome: unaligned gaps are split at their
middle between the neighbouring pieces and the first/last pieces are extended to the chromosome
ends; overlapping pieces (duplicated mappings) keep the earlier piece. Overlaying both haplotypes
gives stretches with one maternal and one paternal segment, hence one profile.

Every bin of the run (the bins of its combo/combo.tsv, or of --bins) goes to the cluster covering
most of its base pairs.

Outputs (DIR, default RUN/oracle):
  clusters.tsv  #CHR START END CLUSTER FRAC ALIGNED NORMAL   one line per bin (Caller --clusters
                format plus): FRAC = share of the bin's bp in its cluster, ALIGNED = share of the
                bin's bp that aligned on both haplotypes (the rest was gap-filled), NORMAL = 1 for
                the cluster without CNAs
  summary.tsv   one line per cluster: bins, Mb, the tree segments in it, cells whose maternal /
                paternal CN is not 1, and the distinct copy-number states
"""
import argparse
import os
import sys
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd

INF = 1 << 62


def chrom_key(c):
    c = str(c)
    c = c[3:] if c.lower().startswith("chr") else c
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def norm_chrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def read_bins(args):
    if args.bins:
        path = args.bins
    else:
        path = os.path.join(args.run, "combo", "combo.tsv")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: no bins file {path}")
    bins = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
                       names=["chrom", "start", "end"], dtype={"chrom": str}).drop_duplicates()
    bins["_order"] = bins.chrom.map(chrom_key)
    bins = bins.sort_values(["_order", "start"]).drop(columns="_order")
    return bins.reset_index(drop=True), path


def tile(pieces, notes):
    """pieces: DataFrame start, end, seg of one haplotype and chromosome -> list of
    (start, end, seg, aligned) tiling [0, INF): gaps split at the middle, ends extended."""
    p = pieces.sort_values(["start", "end"]).to_records(index=False)
    clean = []
    for s, e, seg in zip(p["start"], p["end"], p["seg"]):
        if clean and s < clean[-1][1]:                       # overlap (duplicated mapping): keep earlier
            notes["overlap_bp"] += min(e, clean[-1][1]) - s
            s = clean[-1][1]
        if e > s:
            clean.append([int(s), int(e), seg])
    out = []
    for i, (s, e, seg) in enumerate(clean):
        left = 0 if i == 0 else (clean[i - 1][1] + s) // 2
        right = INF if i == len(clean) - 1 else (e + clean[i + 1][0]) // 2
        if left < s:
            out.append((left, s, seg, False))
        out.append((s, e, seg, True))
        if e < right:
            out.append((e, right, seg, False))
    return out


def overlay(a, b):
    """Two tilings of [0, INF) -> list of (start, end, seg_a, seg_b, aligned in both)."""
    out, i, j, pos = [], 0, 0, 0
    while i < len(a) and j < len(b):
        end = min(a[i][1], b[j][1])
        if end > pos:
            out.append((pos, end, a[i][2], b[j][2], a[i][3] and b[j][3]))
            pos = end
        if a[i][1] == end:
            i += 1
        if b[j][1] == end:
            j += 1
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--truth", required=True, help="eval/truth_hg38.tsv (lift_truth.py)")
    ap.add_argument("--run", help="CHISEL run folder whose bins to use (its combo/combo.tsv), e.g. chisel_pbins8/perfect")
    ap.add_argument("--bins", help="Instead of --run: a file whose first 3 columns are CHR START END (e.g. a combo.tsv)")
    ap.add_argument("--out", help="Output folder (default: RUN/oracle)")
    a = ap.parse_args()
    if not a.run and not a.bins:
        ap.error("give --run or --bins")
    out_dir = a.out or os.path.join(a.run, "oracle")
    os.makedirs(out_dir, exist_ok=True)

    t = pd.read_csv(a.truth, sep="\t", dtype={"chrom": str})
    need = {"cell", "hap", "chrom", "start", "end", "cn", "seg"}
    if not need <= set(t.columns):
        sys.exit(f"ERROR: {a.truth} lacks columns {sorted(need - set(t.columns))}")
    cells = sorted(t.cell.unique())
    haps = sorted(t.hap.unique())
    if haps != ["mat", "pat"]:
        sys.exit(f"ERROR: expected haplotypes mat and pat in {a.truth}, found {haps}")

    # copy-number profile of every (haplotype, segment): CN of every cell
    cn = t.groupby(["hap", "seg", "cell"]).cn.first().unstack("cell").reindex(columns=cells)
    if cn.isna().any().any():
        bad = cn[cn.isna().any(axis=1)].index[:3].tolist()
        sys.exit(f"ERROR: some cells have no copy number for segments {bad}")
    profile = {k: tuple(int(v) for v in row) for k, row in zip(cn.index, cn.values)}

    # the pieces are the same in every cell: take them once per haplotype
    pieces = t[["hap", "chrom", "start", "end", "seg"]].drop_duplicates()
    per_cell = t.groupby("cell").size()
    if per_cell.nunique() != 1:
        print(f"[oracle] WARNING: cells have different numbers of truth pieces ({per_cell.min()}-{per_cell.max()})", file=sys.stderr)

    bins, bins_path = read_bins(a)
    notes = defaultdict(int)
    stretches = {}                       # normalized chrom -> list of (start, end, mat_seg, pat_seg, aligned)
    for chrom in pieces.chrom.unique():
        til = {h: tile(pieces[(pieces.chrom == chrom) & (pieces.hap == h)][["start", "end", "seg"]], notes) for h in ("mat", "pat")}
        if not til["mat"] or not til["pat"]:
            print(f"[oracle] WARNING: {chrom} lacks a haplotype in the truth; skipped", file=sys.stderr)
            continue
        stretches[norm_chrom(chrom)] = overlay(til["mat"], til["pat"])

    # oracle cluster of every stretch, and majority cluster of every bin
    key_of = {}                          # (mat profile, pat profile) -> cluster key
    rows, segs_in = [], defaultdict(lambda: (set(), set()))
    for b in bins.itertuples(index=False):
        st = stretches.get(norm_chrom(b.chrom))
        if st is None:
            sys.exit(f"ERROR: bin chromosome {b.chrom} is not in the truth")
        bp, aligned = defaultdict(int), 0
        for s, e, ms, ps, al in st:
            if e <= b.start:
                continue
            if s >= b.end:
                break
            ov = min(e, b.end) - max(s, b.start)
            k = (profile[("mat", ms)], profile[("pat", ps)])
            bp[k] += ov
            aligned += ov if al else 0
            segs_in[k][0].add(ms)
            segs_in[k][1].add(ps)
        k, kbp = max(bp.items(), key=lambda x: x[1])
        size = b.end - b.start
        if k not in key_of:
            key_of[k] = len(key_of)          # clusters numbered in genome order
        is_normal = int(all(v == 1 for v in k[0]) and all(v == 1 for v in k[1]))
        rows.append((b.chrom, b.start, b.end, key_of[k], kbp / size, aligned / size, is_normal))
    res = pd.DataFrame(rows, columns=["#CHR", "START", "END", "CLUSTER", "FRAC", "ALIGNED", "NORMAL"])
    res.to_csv(os.path.join(out_dir, "clusters.tsv"), sep="\t", index=False, float_format="%.4f")

    # summary
    summ = []
    for k, cid in key_of.items():
        sel = res.CLUSTER == cid
        mat, pat = np.array(k[0]), np.array(k[1])
        states = pd.Series([f"{m}|{p}" for m, p in zip(mat, pat)]).value_counts()
        normal = (mat == 1).all() and (pat == 1).all()
        def seg_list(s):
            s = sorted(s, key=lambda x: int(x.split(":")[1].split("-")[0]))
            return ",".join(s) if len(s) <= 4 else f"{','.join(s[:3])},... ({len(s)} segments)"
        summ.append(OrderedDict([
            ("CLUSTER", cid), ("BINS", int(sel.sum())), ("MB", round((res.END[sel] - res.START[sel]).sum() / 1e6, 3)),
            ("MEAN_FRAC", round(res.FRAC[sel].mean(), 3)),
            ("CELLS_MAT_NOT1", int((mat != 1).sum())), ("CELLS_PAT_NOT1", int((pat != 1).sum())),
            ("STATES", "normal (all 1|1)" if normal else ", ".join(f"{s} x{n}" for s, n in states.items())),
            ("MAT_SEGMENTS", seg_list(segs_in[k][0])), ("PAT_SEGMENTS", seg_list(segs_in[k][1]))]))
    summ = pd.DataFrame(summ).sort_values("CLUSTER")
    summ.to_csv(os.path.join(out_dir, "summary.tsv"), sep="\t", index=False)

    print(f"[oracle] bins: {len(res)} from {bins_path}", file=sys.stderr)
    print(f"[oracle] {len(key_of)} oracle clusters (cells: {len(cells)}); bins whose majority cluster has < 75% of their bp: "
          f"{int((res.FRAC < 0.75).sum())}; gap-filled bp: {100 * (1 - (res.ALIGNED * (res.END - res.START)).sum() / (res.END - res.START).sum()):.1f}%"
          + (f"; overlapping truth pieces trimmed: {notes['overlap_bp']:,} bp" if notes["overlap_bp"] else ""), file=sys.stderr)
    with pd.option_context("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 60):
        print(summ.drop(columns=["MAT_SEGMENTS", "PAT_SEGMENTS"]).to_string(index=False), file=sys.stderr)
    print(f"[oracle] wrote {out_dir}/clusters.tsv and summary.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
