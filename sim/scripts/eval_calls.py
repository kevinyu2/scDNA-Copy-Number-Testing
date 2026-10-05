#!/usr/bin/env python3
"""
Compare a caller's haplotype-specific copy numbers with simulated truth.

Inputs
  --truth        truth_hg38.tsv from lift_truth.py (aligned pieces only; no interpolation)
  --calls        common calls table from standardize_calls.py
                 (cell, chrom, start, end, CN_A, CN_B [, gamma, ...])
  --caller-name  label used in plots and the summary (e.g. "CHISEL (ugp)")

Truth for a bin (used only for scoring and the heatmaps; the track plots draw
the truth at its own breakpoints)
  Of the truth that aligns inside the bin, the copy number covering the most
  bases is the bin's true value.
    purity_*   fraction of the bin's aligned truth with that value
               (low = the bin straddles a real breakpoint)
    covered_*  fraction of the bin that has any aligned truth
  A bin is scored when purity >= --purity and covered >= --min-covered on
  both haplotypes. Unaligned gaps do not count against purity.
  Bins whose length-weighted mean truth differs from the majority truth by
  more than --max-mean-dev are also skipped: small high-CN pieces (ecDNA)
  barely move the majority but dominate the read depth callers measure.

Haplotype orientation (which called haplotype, A or B, is maternal)
  actual  one assignment for all cells (--orientation global), or one per
          cell (--orientation cell). Phase switch errors count as errors.
  ideal   "ideal phasing interpretation": for every cell and bin separately,
          whichever assignment gives the smaller |mat diff| + |pat diff|.
          Upper bound on what the same calls could score if every phase
          switch were corrected.

Outputs (in --out-dir)
  per_cell_tracks.pdf       maternal/paternal predicted vs true (actual orientation);
                            shows all data, scoring is not marked
  per_cell_tracks_ideal.pdf same, ideal phasing interpretation (purple ticks = A/B swapped)
  confusion_all_cells.pdf true vs predicted per bin, all cells pooled (actual + ideal)
  bins.tsv, cells.tsv, summary.txt
  sections.tsv            per cn_mat column (and haplotype), pooled over cells:
                          % of compared bases called correctly, mean (signed and
                          absolute) CN difference, % of bases that are gain/loss false
                          positives / false negatives, and number of cells gained /
                          lost / normal (true and called); actual and ideal phasing
  sections_by_cell.tsv    the same per cell x column
  bin_event_counts.tsv    per bin: number of cells with a gain / loss / normal state,
                          true and predicted, per haplotype and for total CN
  event_counts.pdf        those counts along the chromosome
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm


# ---------------------------------------------------------------- helpers

def norm_cell(name):
    return name[:-4] if name.endswith("_sim") else name


def norm_chrom(c):
    c = str(c)
    return c if c.startswith("chr") else f"chr{c}"


def merge_equal_neighbours(g, max_gap):
    """Join touching truth pieces with the same CN (the simulator splits every
    cell at every CNV boundary in the tree). Gaps > max_gap and overlaps are
    left alone, so blanks and duplicated mappings survive."""
    g = g.sort_values("start")
    out = []
    for s, e, c in zip(g.start.values, g.end.values, g.cn.values):
        if out and out[-1][2] == c and 0 <= s - out[-1][1] <= max_gap:
            out[-1][1] = e
        else:
            out.append([int(s), int(e), c])
    return pd.DataFrame(out, columns=["start", "end", "cn"])


def bin_truth(bins, segs):
    """Per bin: majority CN, purity among aligned truth, covered fraction, mean CN."""
    n = len(bins)
    if segs is None or segs.empty:
        nan = np.full(n, np.nan)
        return nan, nan, np.zeros(n), nan
    bs, be = bins[:, 0][:, None], bins[:, 1][:, None]
    ov = np.clip(np.minimum(be, segs.end.values[None, :]) - np.maximum(bs, segs.start.values[None, :]), 0, None).astype(float)
    cn = segs.cn.values.astype(float)
    blen = np.maximum((bins[:, 1] - bins[:, 0]).astype(float), 1)
    covered = ov.sum(axis=1)
    values = np.unique(cn)
    per_val = np.stack([ov[:, cn == v].sum(axis=1) for v in values], axis=1)   # n x k
    best = per_val.argmax(axis=1)
    has = covered > 0
    majority = np.where(has, values[best], np.nan)
    purity = np.where(has, per_val[np.arange(n), best] / np.maximum(covered, 1), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(has, (ov * cn).sum(axis=1) / covered, np.nan)
    return majority, purity, np.minimum(covered / blen, 1.0), mean


def step_xy(starts, ends, vals, max_gap=1000):
    """Step-line coordinates, broken (NaN) at gaps > max_gap and at overlaps,
    so blanks stay blank and overlapping pieces are drawn separately."""
    xs, ys = [], []
    prev_end = None
    for s, e, v in zip(starts, ends, vals):
        if prev_end is not None and (s - prev_end > max_gap or s < prev_end):
            xs.append(np.nan); ys.append(np.nan)
        xs += [s, e]; ys += [v, v]
        prev_end = e
    return np.asarray(xs, dtype=float) / 1e6, np.asarray(ys, dtype=float)


class YScale:
    """Linear up to `linear_max`; each higher value gets its own evenly spaced
    row above a break mark (not to scale)."""

    def __init__(self, values, linear_max):
        v = np.unique(np.asarray([x for x in values if np.isfinite(x)], dtype=float))
        low = v[v <= linear_max]
        self.top = max(4.0, float(low.max()) if low.size else 0.0)
        self.high = [float(x) for x in v[v > linear_max]]
        self.pos = {h: self.top + 1.5 + k for k, h in enumerate(self.high)}

    def __call__(self, arr):
        arr = np.asarray(arr, dtype=float)
        out = arr.copy()
        for h, p in self.pos.items():
            out[arr == h] = p
        return out

    def decorate(self, ax):
        step = 1 if self.top <= 8 else 2
        ticks = list(np.arange(0, self.top + 0.01, step))
        labels = [f"{int(t)}" for t in ticks]
        ticks += [self.pos[h] for h in self.high]
        labels += [f"{int(h)}" for h in self.high]
        ax.set_yticks(ticks); ax.set_yticklabels(labels)
        ymax = (self.pos[self.high[-1]] if self.high else self.top) + 0.5
        ax.set_ylim(-0.4, ymax)
        if self.high:
            yb = self.top + 0.75
            ax.axhline(yb, color="grey", lw=0.6, ls=(0, (2, 3)))
            kw = dict(transform=ax.get_yaxis_transform(), color="k", clip_on=False, lw=0.9)
            for dx in (-0.006, 0.006):
                ax.plot([dx - 0.006, dx + 0.006], [yb - 0.25, yb + 0.25], **kw)
            ax.text(1.002, yb, " axis break:\n values above\n not to scale", transform=ax.get_yaxis_transform(),
                    fontsize=6, va="center", ha="left", color="grey")


def event_stats(true, pred, normal):
    """Gain = CN > normal, loss = CN < normal. Returns TP/FP/FN counts per event."""
    ok = true.notna() & pred.notna()
    t, p = true[ok], pred[ok]
    out = {}
    for ev, ft in (("gain", lambda x: x > normal), ("loss", lambda x: x < normal)):
        tt, pp = ft(t), ft(p)
        out[f"{ev}_tp"] = int((tt & pp).sum())
        out[f"{ev}_fp"] = int((~tt & pp).sum())
        out[f"{ev}_fn"] = int((tt & ~pp).sum())
    return out



# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls", required=True, help="common calls table from standardize_calls.py")
    ap.add_argument("--caller-name", default="caller", help="label for plots and summary")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--purity", type=float, default=0.9)
    ap.add_argument("--min-covered", type=float, default=0.5, help="min fraction of a bin with aligned truth to score it")
    ap.add_argument("--orientation", choices=["global", "cell"], default="global")
    ap.add_argument("--cells-per-page", type=int, default=4)
    ap.add_argument("--max-cells-plot", type=int, default=0)
    ap.add_argument("--max-cn", type=int, default=8, help="heatmap axis cap; larger values lumped")
    ap.add_argument("--y-linear-max", type=float, default=8, help="track plots: CN above this goes above an axis break")
    ap.add_argument("--merge-gap", type=int, default=1000)
    ap.add_argument("--max-mean-dev", type=float, default=0.5,
                    help="skip bins where length-weighted mean truth differs from majority truth by more than this "
                         "(small high-CN pieces such as ecDNA dominate read depth there)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    LABEL = a.caller_name

    truth = pd.read_csv(a.truth, sep="\t")
    truth["cell"] = truth["cell"].map(norm_cell)
    truth["chrom"] = truth["chrom"].map(norm_chrom)
    tgroups = {k: merge_equal_neighbours(g, a.merge_gap) for k, g in truth.groupby(["cell", "hap", "chrom"])}
    truth_cells = set(truth.cell)

    allcalls = pd.read_csv(a.calls, sep="\t")
    missing = {"cell", "chrom", "start", "end", "CN_A", "CN_B"} - set(allcalls.columns)
    if missing:
        sys.exit(f"ERROR: {a.calls} lacks columns {sorted(missing)} (make it with standardize_calls.py)")
    allcalls = allcalls.dropna(subset=["CN_A", "CN_B"])
    allcalls["cell"] = allcalls["cell"].map(norm_cell)
    allcalls["chrom"] = allcalls["chrom"].map(norm_chrom)
    called = set(allcalls.cell)
    no_calls = sorted(truth_cells - called)        # e.g. dropped by the caller's read filters
    no_truth = sorted(called - truth_cells)
    rows = []
    for cell, calls in allcalls.groupby("cell"):
        if cell not in truth_cells:
            continue
        for chrom, cb in calls.groupby("chrom"):
            cb = cb.sort_values("start")
            b = cb[["start", "end"]].values.astype(np.int64)
            m = bin_truth(b, tgroups.get((cell, "mat", chrom)))
            p = bin_truth(b, tgroups.get((cell, "pat", chrom)))
            rows.append(pd.DataFrame({
                "cell": cell, "chrom": chrom, "start": b[:, 0], "end": b[:, 1],
                "CN_A": cb.CN_A.values.astype(int), "CN_B": cb.CN_B.values.astype(int),
                "true_mat": m[0], "purity_mat": m[1], "covered_mat": m[2], "truemean_mat": m[3],
                "true_pat": p[0], "purity_pat": p[1], "covered_pat": p[2], "truemean_pat": p[3],
                "gamma": cb["gamma"].values if "gamma" in cb.columns else np.nan,
            }))
    if not rows:
        sys.exit(f"ERROR: no cells with both {LABEL} calls and truth (check cell / chromosome names)")
    bins = pd.concat(rows, ignore_index=True)

    has_truth = bins.true_mat.notna() & bins.true_pat.notna()
    covered_ok = (bins.covered_mat >= a.min_covered) & (bins.covered_pat >= a.min_covered)
    straddles = (bins.purity_mat < a.purity) | (bins.purity_pat < a.purity)
    mean_dev = ((bins.truemean_mat - bins.true_mat).abs() > a.max_mean_dev) | \
               ((bins.truemean_pat - bins.true_pat).abs() > a.max_mean_dev)
    bins["straddles"] = has_truth & straddles
    bins["depth_mismatch"] = has_truth & ~straddles & mean_dev
    bins["scored"] = has_truth & covered_ok & ~straddles & ~mean_dev

    # ---------------- actual orientation ----------------
    def agree(df, swap):
        pm, pp = (df.CN_B, df.CN_A) if swap else (df.CN_A, df.CN_B)
        return ((pm == df.true_mat) & (pp == df.true_pat)).mean()

    sc = bins[bins.scored]
    if a.orientation == "global":
        s0, s1 = agree(sc, False), agree(sc, True)
        bins["swap"] = bool(s1 > s0)
        orient_note = (f"actual: one orientation for all cells, {'A=pat, B=mat' if s1 > s0 else 'A=mat, B=pat'} "
                       f"(exact match {max(s0, s1):.3f} vs {min(s0, s1):.3f})")
    else:
        sw = {c: agree(g, True) > agree(g, False) for c, g in sc.groupby("cell")}
        bins["swap"] = bins.cell.map(sw).fillna(False).astype(bool)
        orient_note = f"actual: one orientation per cell ({sum(sw.values())}/{len(sw)} cells A=pat)"
    bins["pred_mat"] = np.where(bins.swap, bins.CN_B, bins.CN_A)
    bins["pred_pat"] = np.where(bins.swap, bins.CN_A, bins.CN_B)

    # ---------------- ideal phasing interpretation ----------------
    # A haplotype with no aligned truth in the bin contributes nothing, so the
    # other haplotype still decides; ties (and no truth at all) keep the actual orientation.
    d_keep = (bins.pred_mat - bins.true_mat).abs().fillna(0) + (bins.pred_pat - bins.true_pat).abs().fillna(0)
    d_flip = (bins.pred_pat - bins.true_mat).abs().fillna(0) + (bins.pred_mat - bins.true_pat).abs().fillna(0)
    flip = d_flip < d_keep
    bins["ideal_flipped"] = flip
    bins["pred_mat_ideal"] = np.where(flip, bins.pred_pat, bins.pred_mat)
    bins["pred_pat_ideal"] = np.where(flip, bins.pred_mat, bins.pred_pat)

    for k in ("true", "pred"):
        bins[f"{k}_total"] = bins[f"{k}_mat"] + bins[f"{k}_pat"]
    bins["truemean_total"] = bins.truemean_mat + bins.truemean_pat
    bins["true_minor"] = bins[["true_mat", "true_pat"]].min(axis=1)
    bins["pred_minor"] = bins[["pred_mat", "pred_pat"]].min(axis=1)
    bins.to_csv(os.path.join(a.out_dir, "bins.tsv"), sep="\t", index=False)

    # ---------------- per-section (cn_mat column) scores ----------------
    # Base-pair overlap between each lifted truth piece and the caller's bins;
    # no bins are filtered here. A base counts once it has both truth and a call.
    if "seg" in truth.columns:
        recs = []
        for (cell, chrom), cb in bins.groupby(["cell", "chrom"]):
            bs, be = cb.start.values, cb.end.values
            for h in ("mat", "pat"):
                tg = truth[(truth.cell == cell) & (truth.hap == h) & (truth.chrom == chrom)]
                if tg.empty:
                    continue
                ts, te, cn = tg.start.values, tg.end.values, tg.cn.values.astype(float)
                ov = np.clip(np.minimum(te[:, None], be[None, :]) - np.maximum(ts[:, None], bs[None, :]), 0, None).astype(float)
                rec = pd.DataFrame({"cell": cell, "hap": h, "seg": tg.seg.values,
                                    "lifted_bp": (te - ts).astype(float), "true_cn": cn,
                                    "compared_bp": ov.sum(axis=1)})
                for suffix in ("", "_ideal"):
                    pred = cb[f"pred_{h}{suffix}"].values.astype(float)
                    diff = pred[None, :] - cn[:, None]
                    rec[f"correct_bp{suffix}"] = (ov * (diff == 0)).sum(axis=1)
                    rec[f"diff_bp{suffix}"] = (ov * diff).sum(axis=1)
                    rec[f"absdiff_bp{suffix}"] = (ov * np.abs(diff)).sum(axis=1)
                    # Gains/losses per haplotype: gain = CN > 1, loss = CN 0
                    for ev, t_ev, p_ev in (("gain", cn > 1, pred > 1), ("loss", cn < 1, pred < 1)):
                        rec[f"{ev}_fp_bp{suffix}"] = (ov * (~t_ev[:, None] & p_ev[None, :])).sum(axis=1)
                        rec[f"{ev}_fn_bp{suffix}"] = (ov * (t_ev[:, None] & ~p_ev[None, :])).sum(axis=1)
                        rec[f"pred{ev}_bp{suffix}"] = (ov * p_ev[None, :]).sum(axis=1)
                    rec[f"prednormal_bp{suffix}"] = (ov * (pred == 1)[None, :]).sum(axis=1)
                recs.append(rec)
        if recs:
            sec = pd.concat(recs, ignore_index=True)
            sums = ["lifted_bp", "compared_bp"] + [
                f"{k}{x}" for x in ("", "_ideal")
                for k in ("correct_bp", "diff_bp", "absdiff_bp", "gain_fp_bp", "gain_fn_bp", "loss_fp_bp", "loss_fn_bp",
                          "predgain_bp", "predloss_bp", "prednormal_bp")]
            # per cell x section
            byc = sec.groupby(["cell", "hap", "seg"], sort=False).agg({**{c: "sum" for c in sums}, "true_cn": "first"}).reset_index()

            def rates(df):
                cmp_ = df.compared_bp.where(df.compared_bp > 0)
                out = pd.DataFrame(index=df.index)
                for x in ("", "_ideal"):
                    out[f"pct_correct{x}"] = 100 * df[f"correct_bp{x}"] / cmp_
                    out[f"mean_diff{x}"] = df[f"diff_bp{x}"] / cmp_
                    out[f"mean_abs_diff{x}"] = df[f"absdiff_bp{x}"] / cmp_
                    for k in ("gain_fp", "gain_fn", "loss_fp", "loss_fn"):
                        out[f"pct_{k}{x}"] = 100 * df[f"{k}_bp{x}"] / cmp_
                return out

            def state(cn):
                return np.select([cn > 1, cn < 1], ["gain", "loss"], default="normal")

            def pred_state(df, x):
                st = np.array(["gain", "loss", "normal"])[
                    np.argmax(df[[f"predgain_bp{x}", f"predloss_bp{x}", f"prednormal_bp{x}"]].values, axis=1)]
                return np.where(df.compared_bp > 0, st, "na")

            byc["true_state"] = state(byc.true_cn.values)
            for x in ("", "_ideal"):
                byc[f"pred_state{x}"] = pred_state(byc, x)   # majority call over the section's compared bases
            byc = pd.concat([byc[["cell", "hap", "seg", "true_cn", "true_state", "pred_state", "pred_state_ideal",
                                  "lifted_bp", "compared_bp"]], rates(byc)], axis=1)
            byc.to_csv(os.path.join(a.out_dir, "sections_by_cell.tsv"), sep="\t", index=False)

            # per section, pooled over cells (base-pair weighted)
            agg = sec.groupby(["hap", "seg"], sort=False)
            col = agg[sums].sum()
            col["n_cells"] = agg.cell.nunique()
            col["mean_true_cn"] = byc.groupby(["hap", "seg"], sort=False).true_cn.mean()
            col["true_cn_values"] = byc.groupby(["hap", "seg"], sort=False).true_cn.apply(
                lambda v: ",".join(str(int(x)) for x in sorted(v.unique())))
            # number of cells whose section is gained / lost / normal (true, and HiScanner's majority call)
            gk = byc.groupby(["hap", "seg"], sort=False)
            for src, lab in (("true_state", "true"), ("pred_state", "pred"), ("pred_state_ideal", "pred_ideal")):
                for ev in ("gain", "loss", "normal"):
                    col[f"n_cells_{lab}_{ev}"] = gk[src].apply(lambda v, ev=ev: int((v == ev).sum()))
            col = col.reset_index()
            se = col.seg.str.extract(r":(\d+)-(\d+)$").astype(np.int64)
            col["seg_len"] = se[1] - se[0] + 1
            col["lifted_bp_per_cell"] = col.lifted_bp / col.n_cells
            col["frac_compared"] = col.compared_bp / col.lifted_bp
            ncols = [f"n_cells_{lab}_{ev}" for lab in ("true", "pred", "pred_ideal") for ev in ("gain", "loss", "normal")]
            col = pd.concat([col[["hap", "seg", "seg_len", "lifted_bp_per_cell", "frac_compared", "n_cells",
                                  "mean_true_cn", "true_cn_values"] + ncols], rates(col)], axis=1)
            col["seg_start"] = se[0]
            col = col.sort_values(["hap", "seg_start"]).drop(columns="seg_start")
            col.to_csv(os.path.join(a.out_dir, "sections.tsv"), sep="\t", index=False)
    else:
        print("[eval] truth file has no 'seg' column (made by an older lift_truth.py); skipping per-section scores",
              file=sys.stderr)

    # ---------------- metrics ----------------
    def hap_events(p, suffix):
        """Gains/losses per haplotype (gain = CN > 1, loss = CN 0), mat and pat pooled."""
        tr = pd.concat([p.true_mat, p.true_pat], ignore_index=True)
        pr = pd.concat([p[f"pred_mat{suffix}"], p[f"pred_pat{suffix}"]], ignore_index=True)
        return {f"hap_{k}{suffix}": v for k, v in event_stats(tr, pr, 1).items()}

    def metrics(p):
        if len(p) == 0:
            return {}
        return {
            "acc_mat": np.mean(p.pred_mat == p.true_mat),
            "acc_pat": np.mean(p.pred_pat == p.true_pat),
            "acc_both": np.mean((p.pred_mat == p.true_mat) & (p.pred_pat == p.true_pat)),
            "acc_mat_ideal": np.mean(p.pred_mat_ideal == p.true_mat),
            "acc_pat_ideal": np.mean(p.pred_pat_ideal == p.true_pat),
            "acc_both_ideal": np.mean((p.pred_mat_ideal == p.true_mat) & (p.pred_pat_ideal == p.true_pat)),
            "acc_total": np.mean(p.pred_total == p.true_total),
            "mae_hap": np.mean((p.pred_mat - p.true_mat).abs() + (p.pred_pat - p.true_pat).abs()) / 2,
            "mae_hap_ideal": np.mean((p.pred_mat_ideal - p.true_mat).abs() + (p.pred_pat_ideal - p.true_pat).abs()) / 2,
            "mae_total": np.mean((p.pred_total - p.true_total).abs()),
            "frac_bins_flipped_ideal": np.mean(p.ideal_flipped),
            **hap_events(p, ""),
            **hap_events(p, "_ideal"),
            **{f"total_{k}": v for k, v in event_stats(p.true_total, p.pred_total, 2).items()},
        }

    cells = pd.DataFrame([dict(cell=c, n_bins=len(g), n_scored=int(g.scored.sum()),
                               **metrics(g[g.scored]), gamma=g.gamma.iloc[0], swapped=bool(g.swap.iloc[0]))
                          for c, g in bins.groupby("cell")])
    cells.to_csv(os.path.join(a.out_dir, "cells.tsv"), sep="\t", index=False)

    mt = metrics(bins[bins.scored])
    n_bins, n_sc = len(bins), int(bins.scored.sum())
    summ = [
        f"cells evaluated: {bins.cell.nunique()}",
        f"caller: {LABEL}",
        f"cells without {LABEL} calls (skipped): {len(no_calls)}" + (f"  e.g. {no_calls[:5]}" if no_calls else ""),
        f"cells with calls but no truth (skipped): {len(no_truth)}" + (f"  e.g. {no_truth[:5]}" if no_truth else ""),
        f"bins: {n_bins}   scored: {n_sc} ({n_sc / n_bins:.1%})",
        f"  not scored, straddle a true breakpoint: {int(bins.straddles.sum())}",
        f"  not scored, small high-CN truth pieces (e.g. ecDNA; mean vs majority > {a.max_mean_dev}): {int(bins.depth_mismatch.sum())}",
        f"  not scored, < {a.min_covered:.0%} of bin has aligned truth: {int((has_truth & ~covered_ok & ~straddles & ~mean_dev).sum() + (~has_truth).sum())}",
        orient_note,
        "",
        "accuracy over scored bins (exact match)        actual    ideal phasing",
        f"  maternal                                      {mt['acc_mat']:.3f}     {mt['acc_mat_ideal']:.3f}",
        f"  paternal                                      {mt['acc_pat']:.3f}     {mt['acc_pat_ideal']:.3f}",
        f"  both haplotypes                               {mt['acc_both']:.3f}     {mt['acc_both_ideal']:.3f}",
        f"  mean |error| per haplotype                    {mt['mae_hap']:.3f}     {mt['mae_hap_ideal']:.3f}",
        f"  total CN (orientation-free)                   {mt['acc_total']:.3f}",
        f"  total CN mean |error|                         {mt['mae_total']:.3f}",
        f"  bins the ideal interpretation flips           {mt['frac_bins_flipped_ideal']:.3f}",
        "",
        "gains / losses over scored bins (haplotype: gain = CN > 1, loss = CN 0; total: gain > 2, loss < 2)",
        f"  {'':34s}{'TP':>9s}{'FP':>9s}{'FN':>9s}{'precision':>11s}{'recall':>8s}",
    ]
    scored_bins = bins[bins.scored]

    def ev_line(label, d, ev):
        tp, fp, fn = d[f"{ev}_tp"], d[f"{ev}_fp"], d[f"{ev}_fn"]
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        return f"  {label:34s}{tp:9d}{fp:9d}{fn:9d}{prec:11.3f}{rec:8.3f}"

    for lab, t_, p_, nrm in (("maternal", "true_mat", "pred_mat", 1), ("paternal", "true_pat", "pred_pat", 1),
                             ("maternal (ideal phasing)", "true_mat", "pred_mat_ideal", 1),
                             ("paternal (ideal phasing)", "true_pat", "pred_pat_ideal", 1),
                             ("total CN", "true_total", "pred_total", 2)):
        d = event_stats(scored_bins[t_], scored_bins[p_], nrm)
        summ.append(ev_line(f"{lab} gain", d, "gain"))
        summ.append(ev_line(f"{lab} loss", d, "loss"))
    summ += [
        "  (TP/FP/FN count bins; FP = called but not true, FN = true but not called)",
        "",
        f"per-cell both-haplotype accuracy: actual median {cells.acc_both.median():.3f}, "
        f"ideal median {cells.acc_both_ideal.median():.3f}",
    ]
    if bins.gamma.notna().any():
        summ.append(f"cells {LABEL} called as WGD-scaled (gamma > 3): {(cells.gamma > 3).sum()}/{len(cells)}")
    text = "\n".join(summ)
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)

    # ---------------- cells gained / lost / normal at each bin ----------------
    # Population view over all bins (not only scored ones); truth = bin majority.
    def classify(x, normal):
        return np.select([x > normal, x < normal, x == normal], ["gain", "loss", "normal"], default="na")

    cnt_cols = []
    for lab, col, nrm in (("mat_true", "true_mat", 1), ("mat_pred", "pred_mat", 1), ("mat_pred_ideal", "pred_mat_ideal", 1),
                          ("pat_true", "true_pat", 1), ("pat_pred", "pred_pat", 1), ("pat_pred_ideal", "pred_pat_ideal", 1),
                          ("total_true", "true_total", 2), ("total_pred", "pred_total", 2)):
        bins[f"_c_{lab}"] = classify(bins[col].values.astype(float), nrm)
        cnt_cols.append(lab)
    g = bins.groupby(["chrom", "start", "end"], sort=True)
    bc = g.size().rename("n_cells").to_frame()
    for lab in cnt_cols:
        for ev in ("gain", "loss", "normal"):
            bc[f"{lab}_{ev}"] = g[f"_c_{lab}"].apply(lambda v, ev=ev: int((v == ev).sum()))
    bc = bc.reset_index()
    bins.drop(columns=[c for c in bins.columns if c.startswith("_c_")], inplace=True)
    bc.to_csv(os.path.join(a.out_dir, "bin_event_counts.tsv"), sep="\t", index=False)

    with PdfPages(os.path.join(a.out_dir, "event_counts.pdf")) as pdf:
        for suffix, heading in (("", "actual orientation"), ("_ideal", "ideal phasing interpretation")):
            rows_ = [("mat", "maternal (gain = CN > 1, loss = CN 0)"), ("pat", "paternal (gain = CN > 1, loss = CN 0)")]
            if suffix == "":
                rows_.append(("total", "total CN (gain > 2, loss < 2)"))
            for chrom, cc in bc.groupby("chrom"):
                fig, axs = plt.subplots(len(rows_), 1, figsize=(12.5, 3.2 * len(rows_)), squeeze=False)
                for ax, (h, title) in zip(axs[:, 0], rows_):
                    for ev, colr in (("gain", "tab:red"), ("loss", "tab:blue")):
                        tx, ty = step_xy(cc.start.values, cc.end.values, cc[f"{h}_true_{ev}"].values.astype(float), 1000)
                        px, py = step_xy(cc.start.values, cc.end.values, cc[f"{h}_pred{suffix}_{ev}"].values.astype(float), 1000)
                        ax.plot(tx, ty, color=colr, lw=3, alpha=0.35, label=f"true {ev}")
                        ax.plot(px, py, color=colr, lw=1.2, label=f"{LABEL} {ev}")
                    ax.set_ylabel("cells", fontsize=8)
                    top = max(1, *(cc[f"{h}_{k}_{ev}"].max() for k in ("true", f"pred{suffix}") for ev in ("gain", "loss")))
                    ax.set_ylim(-0.02 * top, top * 1.15)        # headroom so lines at the max stay visible
                    ax.set_title(title, fontsize=9, loc="left")
                    ax.grid(axis="y", alpha=0.3)
                    ax.legend(loc="upper right", fontsize=7, ncol=4)
                axs[-1, 0].set_xlabel(f"{chrom} position (Mb, hg38)", fontsize=8)
                fig.suptitle(f"Cells with a gain or loss at each bin - {heading} "
                             f"(normal = cells called minus gains and losses)", fontsize=10)
                fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)

    # ---------------- heatmaps ----------------
    cap = a.max_cn

    def confusion(ax, t, pr, title, vmax):
        ok = t.notna() & pr.notna()
        exact = float(np.mean(t[ok].astype(int).values == pr[ok].astype(int).values)) if ok.any() else np.nan
        t = np.clip(t[ok].astype(int), 0, cap); pr = np.clip(pr[ok].astype(int), 0, cap)
        m = np.zeros((cap + 1, cap + 1), dtype=np.int64)
        np.add.at(m, (pr.values, t.values), 1)
        im = ax.imshow(np.where(m > 0, m, np.nan), origin="lower", cmap="viridis", norm=LogNorm(vmin=1, vmax=vmax))
        for i in range(cap + 1):
            for j in range(cap + 1):
                if m[i, j]:
                    ax.text(j, i, f"{m[i, j]}", ha="center", va="center", fontsize=6,
                            color="white" if m[i, j] < vmax ** 0.5 else "black")
        lab = [str(i) for i in range(cap)] + [f"{cap}+"]
        ax.set_xticks(range(cap + 1)); ax.set_xticklabels(lab)
        ax.set_yticks(range(cap + 1)); ax.set_yticklabels(lab)
        ax.plot([-0.5, cap + 0.5], [-0.5, cap + 0.5], color="red", lw=0.8, ls="--")
        # exact match uses real values, so e.g. 38 vs 40 is wrong even though both sit in the "cap+" box
        ax.set_title(f"{title}\nexact match {exact:.3f} (n={m.sum()})", fontsize=9)
        ax.set_xlabel("true copy number"); ax.set_ylabel(f"{LABEL} copy number")
        return im

    def heat_page(pdf, df, panels, heading):
        fig, axs = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2))
        vmax = max(1, len(df))
        for ax, (t, pr, title) in zip(np.atleast_1d(axs), panels):
            im = confusion(ax, df[t], df[pr], title, vmax)
        fig.colorbar(im, ax=axs, shrink=0.8, label="bins (log)")
        fig.suptitle(heading, fontsize=10)
        pdf.savefig(fig); plt.close(fig)

    scored = bins[bins.scored]
    with PdfPages(os.path.join(a.out_dir, "confusion_all_cells.pdf")) as pdf:
        heat_page(pdf, scored,
                  [("true_mat", "pred_mat", "maternal"), ("true_pat", "pred_pat", "paternal"),
                   ("true_total", "pred_total", "total"), ("true_minor", "pred_minor", "minor allele (orientation-free)")],
                  f"All cells, scored bins - actual orientation\n{orient_note}")
        heat_page(pdf, scored,
                  [("true_mat", "pred_mat_ideal", "maternal (ideal phasing)"),
                   ("true_pat", "pred_pat_ideal", "paternal (ideal phasing)"),
                   ("true_total", "pred_total", "total (unchanged)")],
                  "All cells, scored bins - ideal phasing interpretation\n"
                  "each cell/bin uses whichever haplotype assignment fits the truth better")
        heat_page(pdf, bins[has_truth],
                  [("true_mat", "pred_mat", "maternal"), ("true_pat", "pred_pat", "paternal"),
                   ("true_total", "pred_total", "total")],
                  "All bins with any aligned truth (includes breakpoint and partly unaligned bins) - actual orientation")
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        for ax, h in zip(axs, ("mat", "pat", "total")):
            x, y = bins[f"truemean_{h}"], bins[f"pred_{h}"]
            ok = x.notna()
            ax.hist2d(x[ok].clip(0, cap), y[ok].clip(0, cap),
                      bins=[np.arange(0, cap + 0.26, 0.25), np.arange(-0.5, cap + 1)], norm=LogNorm(), cmap="viridis")
            ax.plot([0, cap], [0, cap], color="red", lw=0.8, ls="--")
            ax.set_xlabel("true CN (length-weighted mean over the bin's aligned truth)")
            ax.set_ylabel(f"{LABEL} CN"); ax.set_title(h)
        fig.suptitle("All bins - fractional truth for bins that straddle breakpoints", fontsize=10)
        pdf.savefig(fig); plt.close(fig)

    # ---------------- per-cell tracks ----------------
    plot_cells = sorted(cells.cell)
    if a.max_cells_plot > 0:
        plot_cells = plot_cells[: a.max_cells_plot]
    cst = cells.set_index("cell")

    def plot_tracks(path, suffix, heading):
        cpp = max(1, a.cells_per_page)
        with PdfPages(path) as pdf:
            for k in range(0, len(plot_cells), cpp):
                chunk = plot_cells[k:k + cpp]
                fig, axs = plt.subplots(2 * len(chunk), 1, figsize=(12.5, 4.0 * len(chunk)), squeeze=False)
                for i, cell in enumerate(chunk):
                    cb = bins[bins.cell == cell].sort_values("start")
                    chrom = cb.chrom.iloc[0]
                    st = cst.loc[cell]
                    for j, (h, name, col) in enumerate((("mat", "maternal", "tab:red"), ("pat", "paternal", "tab:blue"))):
                        ax = axs[2 * i + j, 0]
                        pcol = f"pred_{h}{suffix}"
                        tg = tgroups.get((cell, h, chrom), pd.DataFrame(columns=["start", "end", "cn"]))
                        ys = YScale(list(tg.cn.values) + list(cb[pcol].values), a.y_linear_max)
                        tx, ty = step_xy(tg.start.values, tg.end.values, tg.cn.values.astype(float), a.merge_gap)
                        px, py = step_xy(cb.start.values, cb.end.values, cb[pcol].values.astype(float), 1000)
                        ax.plot(tx, ys(ty), color="black", lw=3, alpha=0.45, solid_capstyle="butt", zorder=2,
                                label="true")
                        ax.plot(px, ys(py), color=col, lw=1.3, zorder=3, label=LABEL)
                        ys.decorate(ax)
                        if suffix == "_ideal":
                            fl = cb[cb.ideal_flipped]
                            if len(fl):
                                ax.scatter((fl.start + fl.end) / 2e6, np.full(len(fl), -0.25), marker="|", s=40,
                                           color="tab:purple", zorder=4, label="A/B swapped here")
                        ax.set_ylabel(f"{name}\nCN", fontsize=8)
                        ax.grid(axis="y", alpha=0.3)
                        ax.legend(loc="upper left", fontsize=6, ncol=3)
                        if j == 0:
                            ax.set_title(f"{cell}   both-hap acc {st.get('acc_both' + suffix, np.nan):.2f}   "
                                         f"total acc {st.get('acc_total', np.nan):.2f}" +
                                         (f"   gamma {st.gamma:.2f}" if pd.notna(st.gamma) else "") + f"   [{heading}]",
                                         fontsize=9, loc="left")
                        else:
                            ax.set_xlabel(f"{chrom} position (Mb, hg38)", fontsize=8)
                fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)

    plot_tracks(os.path.join(a.out_dir, "per_cell_tracks.pdf"), "", "actual orientation")
    plot_tracks(os.path.join(a.out_dir, "per_cell_tracks_ideal.pdf"), "_ideal", "ideal phasing interpretation")

    print(f"\nwrote {a.out_dir}/per_cell_tracks.pdf, per_cell_tracks_ideal.pdf, confusion_all_cells.pdf, event_counts.pdf, bins.tsv, cells.tsv, summary.txt")


if __name__ == "__main__":
    main()