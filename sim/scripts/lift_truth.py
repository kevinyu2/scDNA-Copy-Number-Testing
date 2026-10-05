#!/usr/bin/env python3
"""
Lift the simulator's true copy numbers (cn_mat/*_mat.tsv, *_pat.tsv) from
HG002 haplotype-assembly coordinates onto hg38, the coordinates HiScanner
reports in. Nothing is interpolated.

Each truth segment is cut into the pieces that actually align to hg38
(minimap2 PAF with --cs, one per haplotype). Sequence with no alignment
(centromere, large insertions) is simply absent from the output, so it
shows up as a blank in plots. If two parts of a haplotype align onto the
same hg38 region, both pieces are kept.

Aligned blocks separated by small indels (<= --merge-gap bp in both
genomes) are joined into one piece, so a line is not broken at every
small indel.

Output (0-based half-open hg38 coordinates), one row per aligned piece:
    cell  hap  chrom  start  end  cn  seg  orig_start  orig_end
seg is the cn_mat column the piece came from; orig_start/orig_end are the
piece's coordinates in the haplotype assembly.
"""
import argparse
import bisect
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

CS_RE = re.compile(r":(\d+)|\*([a-z]{2})|\+([a-z]+)|-([a-z]+)|~[a-z]{2}(\d+)[a-z]{2}")


def cs_blocks(cs, q, t):
    """Gap-free aligned blocks (q_start, q_end, t_start) from a cs string."""
    out = []
    bq = bt = None
    for m in CS_RE.finditer(cs):
        match, sub, ins, dele, intron = m.groups()
        if match is not None or sub is not None:
            n = int(match) if match is not None else 1
            if bq is None:
                bq, bt = q, t
            q += n
            t += n
        else:
            if bq is not None:
                out.append((bq, q, bt))
                bq = None
            if ins is not None:
                q += len(ins)
            elif dele is not None:
                t += len(dele)
            elif intron is not None:
                t += int(intron)
    if bq is not None:
        out.append((bq, q, bt))
    return out


class AlignedPieces:
    def __init__(self, paf, query_chrom, target_chrom, merge_gap):
        alns = []
        with open(paf) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12 or p[0] != query_chrom or p[5] != target_chrom or p[4] != "+":
                    continue
                tags = {x[:2]: x[5:] for x in p[12:]}
                if tags.get("tp") != "P" or "cs" not in tags:
                    continue
                alns.append((int(p[3]) - int(p[2]), int(p[2]), int(p[7]), tags["cs"]))
        if not alns:
            sys.exit(f"ERROR: no primary '+' alignments of {query_chrom} -> {target_chrom} with cs tags in {paf}")

        # Longest alignments first; later ones may not reuse query bases already
        # placed, so each haplotype base maps to at most one hg38 position.
        alns.sort(reverse=True)
        taken = []                      # sorted, non-overlapping (q_start, q_end)
        blocks = []
        for _, q0, t0, cs in alns:
            for bq, eq, bt in cs_blocks(cs, q0, t0):
                for s, e, tt in self._free_parts(taken, bq, eq, bt):
                    blocks.append((s, e, tt))
                    bisect.insort(taken, (s, e))
        blocks.sort()
        self.qs = np.array([b[0] for b in blocks], dtype=np.int64)
        self.qe = np.array([b[1] for b in blocks], dtype=np.int64)
        self.ts = np.array([b[2] for b in blocks], dtype=np.int64)

        # Join blocks into pieces across small indels
        regions = []                    # [q_start, q_end, first_block, last_block]
        for i in range(len(blocks)):
            if regions:
                r = regions[-1]
                j = r[3]
                dq = self.qs[i] - self.qe[j]
                dt = self.ts[i] - (self.ts[j] + self.qe[j] - self.qs[j])
                if 0 <= dq <= merge_gap and 0 <= dt <= merge_gap:
                    r[1] = self.qe[i]; r[3] = i
                    continue
            regions.append([self.qs[i], self.qe[i], i, i])
        self.regions = regions
        self.rq0 = np.array([r[0] for r in regions], dtype=np.int64)
        self.n_aln = len(alns)
        self.aligned_bp = int((self.qe - self.qs).sum())

    @staticmethod
    def _free_parts(taken, s, e, t):
        """Parts of [s,e) not already in `taken`, with their target starts."""
        parts, cur = [], s
        k = max(0, bisect.bisect_left(taken, (s, s)) - 1)
        while k < len(taken) and taken[k][0] < e:
            ts_, te_ = taken[k]
            if te_ > cur:
                if ts_ > cur:
                    parts.append((cur, min(ts_, e), t + (cur - s)))
                cur = max(cur, te_)
            k += 1
        if cur < e:
            parts.append((cur, e, t + (cur - s)))
        return [p for p in parts if p[1] > p[0]]

    def _map(self, q, lo, hi, side):
        """Map q inside region blocks lo..hi; q in a small internal gap snaps to a block edge."""
        i = lo + np.searchsorted(self.qs[lo:hi + 1], q, side="right") - 1
        i = max(i, lo)
        if self.qs[i] <= q < self.qe[i]:
            return int(self.ts[i] + (q - self.qs[i]))
        if side == "end" or i == hi:
            return int(self.ts[i] + (min(q, self.qe[i]) - self.qs[i]))
        return int(self.ts[i + 1])

    def pieces(self, a, b):
        """hg38 intervals for the aligned parts of haplotype interval [a, b)."""
        out = []
        k = max(0, np.searchsorted(self.rq0, a, side="right") - 1)
        while k < len(self.regions) and self.regions[k][0] < b:
            r0, r1, lo, hi = self.regions[k]
            qa, qb = max(a, r0), min(b, r1)
            if qb > qa:
                ta, tb = self._map(qa, lo, hi, "start"), self._map(qb, lo, hi, "end")
                if tb > ta:
                    out.append((ta, tb, qa, qb))
            k += 1
        return out


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
        long["seg_start"] = parts[1].astype(np.int64) - 1      # 1-based inclusive -> 0-based half-open
        long["seg_end"] = parts[2].astype(np.int64)
        rows.append(long[["cell", "chrom", "seg", "seg_start", "seg_end", "cn"]])
    out = pd.concat(rows, ignore_index=True)
    print(f"[lift_truth] {hap}: {len(files)} file(s), {out.cell.nunique()} cells", file=sys.stderr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cn-mat-dir", required=True)
    ap.add_argument("--mat-paf", required=True)
    ap.add_argument("--pat-paf", required=True)
    ap.add_argument("--query-chrom", default="chr1", help="contig name in the HG002 haplotype FASTAs")
    ap.add_argument("--target-chrom", default="chr1", help="contig name in hg38")
    ap.add_argument("--merge-gap", type=int, default=1000,
                    help="join aligned blocks separated by indels up to this size (bp)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    out = []
    for hap, paf in (("mat", a.mat_paf), ("pat", a.pat_paf)):
        truth = read_truth(a.cn_mat_dir, hap)
        truth = truth[truth.chrom == a.query_chrom]
        ap_ = AlignedPieces(paf, a.query_chrom, a.target_chrom, a.merge_gap)
        span = int(truth.seg_end.max())
        print(f"[lift_truth] {hap}: {ap_.n_aln} alignment(s); {ap_.aligned_bp / span:.1%} of the haplotype aligns; "
              f"{len(ap_.regions)} aligned pieces after joining indels <= {a.merge_gap} bp", file=sys.stderr)

        # Lift each distinct segment once, then attach every cell's CN
        segs = truth[["seg_start", "seg_end"]].drop_duplicates()
        lifted = []
        for s, e in segs.itertuples(index=False):
            for ta, tb, qa, qb in ap_.pieces(int(s), int(e)):
                lifted.append((s, e, ta, tb, qa, qb))
        lf = pd.DataFrame(lifted, columns=["seg_start", "seg_end", "start", "end", "orig_start", "orig_end"])
        res = truth.merge(lf, on=["seg_start", "seg_end"])
        res["hap"] = hap
        res["chrom"] = a.target_chrom
        out.append(res[["cell", "hap", "chrom", "start", "end", "cn", "seg", "orig_start", "orig_end"]])

    res = pd.concat(out, ignore_index=True).sort_values(["cell", "hap", "start"])
    res.to_csv(a.out, sep="\t", index=False)
    print(f"[lift_truth] wrote {a.out}: {res.cell.nunique()} cells, {len(res)} rows", file=sys.stderr)


if __name__ == "__main__":
    main()