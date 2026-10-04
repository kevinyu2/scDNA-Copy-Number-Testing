#!/usr/bin/env python3
"""
Lift the simulator's true copy numbers (cn_mat/*_mat.tsv, *_pat.tsv) from
HG002 haplotype-assembly coordinates onto hg38, the coordinates HiScanner
reports in.

Each haplotype has its own coordinate system, so maternal segments are lifted
through a maternal->hg38 alignment and paternal through paternal->hg38
(minimap2 PAF with --cs, produced by eval_run.sh).

Output (long format, 0-based half-open hg38 coordinates):
    cell  hap  chrom  start  end  cn  orig_start  orig_end  start_gap  end_gap
start_gap / end_gap: distance (bp) from the boundary to the nearest aligned
base; >0 means the boundary fell in an unaligned region (e.g. centromere) and
was interpolated.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

CS_RE = re.compile(r":(\d+)|\*([a-z]{2})|\+([a-z]+)|-([a-z]+)|~[a-z]{2}(\d+)[a-z]{2}")


class LiftMap:
    """Piecewise-linear query->target map built from gap-free aligned blocks."""

    def __init__(self, paf, query_chrom, target_chrom):
        qs, qe, ts = [], [], []
        n_aln = 0
        with open(paf) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12 or p[0] != query_chrom or p[5] != target_chrom or p[4] != "+":
                    continue
                tags = {t[:2]: t[5:] for t in p[12:]}
                if tags.get("tp") != "P" or "cs" not in tags:
                    continue
                n_aln += 1
                q, t = int(p[2]), int(p[7])
                bq = bt = None                      # open block start
                for m in CS_RE.finditer(tags["cs"]):
                    match, sub, ins, dele, intron = m.groups()
                    if match is not None or sub is not None:
                        n = int(match) if match is not None else 1
                        if bq is None:
                            bq, bt = q, t
                        q += n
                        t += n
                    else:
                        if bq is not None:          # close block
                            qs.append(bq); qe.append(q); ts.append(bt)
                            bq = None
                        if ins is not None:
                            q += len(ins)           # bases only in query
                        elif dele is not None:
                            t += len(dele)          # bases only in target
                        elif intron is not None:
                            t += int(intron)
                if bq is not None:
                    qs.append(bq); qe.append(q); ts.append(bt)

        if not qs:
            sys.exit(f"ERROR: no primary '+' alignments of {query_chrom} -> {target_chrom} with cs tags in {paf}")

        order = np.argsort(qs, kind="stable")
        self.qs = np.asarray(qs, dtype=np.int64)[order]
        self.qe = np.asarray(qe, dtype=np.int64)[order]
        self.ts = np.asarray(ts, dtype=np.int64)[order]
        self.n_aln = n_aln

    def lift(self, q):
        """Map a 0-based query position. Returns (target_pos, gap_to_nearest_aligned_base)."""
        i = np.searchsorted(self.qs, q, side="right") - 1
        if i >= 0 and q < self.qe[i]:
            return int(self.ts[i] + (q - self.qs[i])), 0
        # Unaligned: interpolate between the block ending before q and the next block
        if i < 0:
            j = 0
            return int(max(0, self.ts[j] - (self.qs[j] - q))), int(self.qs[j] - q)
        lq = self.qe[i]
        lt = self.ts[i] + (self.qe[i] - self.qs[i])
        if i + 1 >= len(self.qs):
            return int(lt + (q - lq)), int(q - lq)
        rq, rt = self.qs[i + 1], self.ts[i + 1]
        frac = (q - lq) / max(1, rq - lq)
        return int(round(lt + frac * (rt - lt))), int(min(q - lq, rq - q))


def read_truth(cn_mat_dir, hap):
    files = sorted(glob.glob(os.path.join(cn_mat_dir, f"*{hap}.tsv")))
    if not files:
        sys.exit(f"ERROR: no *{hap}.tsv files in {cn_mat_dir}")
    rows = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        if "cell" not in df.columns:
            sys.exit(f"ERROR: {f} has no 'cell' column")
        seg_cols = [c for c in df.columns if re.fullmatch(r"[^:]+:\d+-\d+", str(c))]
        long = df.melt(id_vars="cell", value_vars=seg_cols, var_name="seg", value_name="cn").dropna(subset=["cn"])
        parts = long["seg"].str.extract(r"^([^:]+):(\d+)-(\d+)$")
        long["chrom"] = parts[0]
        # Simulator segments are 1-based inclusive -> 0-based half-open
        long["orig_start"] = parts[1].astype(np.int64) - 1
        long["orig_end"] = parts[2].astype(np.int64)
        long["hap"] = hap
        rows.append(long[["cell", "hap", "chrom", "orig_start", "orig_end", "cn"]])
    out = pd.concat(rows, ignore_index=True)
    print(f"[lift_truth] {hap}: {len(files)} file(s), {out.cell.nunique()} cells, "
          f"{out[['orig_start', 'orig_end']].drop_duplicates().shape[0]} distinct segments", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cn-mat-dir", required=True)
    ap.add_argument("--mat-paf", required=True)
    ap.add_argument("--pat-paf", required=True)
    ap.add_argument("--query-chrom", default="chr1", help="contig name in the HG002 haplotype FASTAs")
    ap.add_argument("--target-chrom", default="chr1", help="contig name in hg38")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = []
    for hap, paf in (("mat", a.mat_paf), ("pat", a.pat_paf)):
        truth = read_truth(a.cn_mat_dir, hap)
        truth = truth[truth.chrom == a.query_chrom]
        lm = LiftMap(paf, a.query_chrom, a.target_chrom)
        print(f"[lift_truth] {hap}: {lm.n_aln} alignment(s), {len(lm.qs)} aligned blocks", file=sys.stderr)

        # Lift each distinct boundary once
        bounds = np.unique(np.concatenate([truth.orig_start.values, truth.orig_end.values]))
        lifted = {int(b): lm.lift(int(b)) for b in bounds}
        truth["start"] = truth.orig_start.map(lambda b: lifted[b][0])
        truth["end"] = truth.orig_end.map(lambda b: lifted[b][0])
        truth["start_gap"] = truth.orig_start.map(lambda b: lifted[b][1])
        truth["end_gap"] = truth.orig_end.map(lambda b: lifted[b][1])
        truth["chrom"] = a.target_chrom

        bad = truth.end <= truth.start
        if bad.any():
            print(f"[lift_truth] {hap}: {bad.sum()} segment rows collapsed or inverted after liftover "
                  f"(very short segments or rearranged regions); dropped", file=sys.stderr)
            truth = truth[~bad]
        far = (truth[["start_gap", "end_gap"]].max(axis=1) > 50_000)
        if far.any():
            n = truth.loc[far, ["orig_start", "orig_end"]].drop_duplicates().shape[0]
            print(f"[lift_truth] {hap}: {n} segment(s) have a boundary >50 kb from aligned sequence "
                  f"(interpolated; check start_gap/end_gap)", file=sys.stderr)
        out.append(truth)

    res = pd.concat(out, ignore_index=True)
    res = res[["cell", "hap", "chrom", "start", "end", "cn", "orig_start", "orig_end", "start_gap", "end_gap"]]
    res.sort_values(["cell", "hap", "start"]).to_csv(a.out, sep="\t", index=False)
    print(f"[lift_truth] wrote {a.out}: {res.cell.nunique()} cells", file=sys.stderr)


if __name__ == "__main__":
    main()