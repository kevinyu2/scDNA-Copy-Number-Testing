#!/usr/bin/env python3
"""
Compare a caller's haplotype-specific copy numbers with simulated truth.

Inputs
  --truth        truth_hg38.tsv from lift_truth.py (aligned pieces only; no interpolation)
  --calls        common calls table from standardize_calls.py
                 (cell, chrom, start, end, CN_A, CN_B [, gamma, ...])
  --caller-name  label used in plots and the summary (e.g. "CHISEL (ugp)")

Truth for a bin (used only for scoring; the track plots draw
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
  bin     "ideal phasing": one assignment per genomic bin, shared by all cells,
          whichever gives the smaller |mat diff| + |pat diff| summed over cells.
          Germline phase comes from the bulk and is read the same way in every
          cell, so a phase switch error flips a bin in all cells at once; this is
          what the same calls would score with every switch corrected.
          (Column suffix _bin.)
  ideal   "allele-specific": for every cell and bin separately, whichever
          assignment fits better. With two haplotypes this equals comparing the
          sorted (major, minor) pairs, i.e. allele-specific CN scoring. It also
          forgives events put on the wrong haplotype in some cells, so it is an
          upper bound, not a phasing correction. (Column suffix _ideal, kept for
          compatibility.)

Outputs (in --out-dir)
  summary.txt               accuracy over scored bins and gain / normal / loss
                            TP / TN / FP / FN tables, for all three orientation levels
  confusion_segments_cn.pdf     true CN vs called CN (CN 0-4 exact, then groups of 5)
  confusion_segments_events.pdf gain / normal / loss (3x3; gain = CN > 1, normal = 1, loss = 0)
                            Both: all cn_mat segments pooled, then one page per segment;
                            maternal and paternal x actual / ideal phasing / allele-specific.
                            --conf-unit bp (default, base-pair weighted, in cell-Mb) or cells
                            (each cell x segment once, majority call)
  section_confusion.tsv     the numbers behind it, exact CNs, both units
                            (hap, seg, level, true_cn, pred_cn, cells, bp)
  per_cell_tracks.pdf       maternal/paternal predicted vs true (actual orientation);
                            shows all data, scoring is not marked
  per_cell_tracks_bin.pdf   same, ideal phasing (per bin; purple ticks = bin flipped vs actual)
  per_cell_tracks_ideal.pdf same, allele-specific (purple ticks = A/B swapped)

--levels allele_specific (used for HiScanner) keeps only the allele-specific level in
every output: such callers report CN_A|CN_B as major|minor (BAF mirrored per bin), so
there is no haplotype orientation to score. per_cell_tracks.pdf and _bin.pdf are not made.
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


EVENTS = ("gain", "normal", "loss")
CONF = ("tp", "tn", "fp", "fn")


def event_masks(x, normal):
    """gain = CN > normal, normal = CN == normal, loss = CN < normal."""
    return {"gain": x > normal, "normal": x == normal, "loss": x < normal}


def event_stats(true, pred, normal):
    """One-vs-rest TP/TN/FP/FN counts for each of gain / normal / loss."""
    ok = true.notna() & pred.notna()
    mt, mp = event_masks(true[ok], normal), event_masks(pred[ok], normal)
    out = {}
    for ev in EVENTS:
        tt, pp = mt[ev], mp[ev]
        out[f"{ev}_tp"] = int((tt & pp).sum())
        out[f"{ev}_tn"] = int((~tt & ~pp).sum())
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
    ap.add_argument("--conf-unit", choices=["cells", "bp"], default="bp",
                    help="confusion_segments_*.pdf: base pairs (default) or cell x segments (majority call)")
    ap.add_argument("--y-linear-max", type=float, default=8, help="track plots: CN above this goes above an axis break")
    ap.add_argument("--merge-gap", type=int, default=1000)
    ap.add_argument("--max-mean-dev", type=float, default=0.5,
                    help="skip bins where length-weighted mean truth differs from majority truth by more than this "
                         "(small high-CN pieces such as ecDNA dominate read depth there)")
    ap.add_argument("--levels", choices=["all", "allele_specific"], default="all",
                    help="all: actual / ideal phasing / allele-specific; allele_specific: only that level "
                         "(for callers whose A|B is major|minor, e.g. HiScanner)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    PHASED = a.levels == "all"
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

    # ---------------- ideal phasing: one orientation per genomic bin ----------------
    # Chosen from the raw A/B calls (independent of --orientation), summing the error
    # over every cell; ties keep each cell's actual orientation.
    dk = (bins.CN_A - bins.true_mat).abs().fillna(0) + (bins.CN_B - bins.true_pat).abs().fillna(0)
    df_ = (bins.CN_B - bins.true_mat).abs().fillna(0) + (bins.CN_A - bins.true_pat).abs().fillna(0)
    bkey = [bins.chrom, bins.start, bins.end]
    sk, sf = dk.groupby(bkey).transform("sum"), df_.groupby(bkey).transform("sum")
    swap_bin = np.where(sf < sk, True, np.where(sf > sk, False, bins.swap.values))
    bins["bin_flipped"] = swap_bin != bins.swap.values          # differs from the actual orientation
    bins["pred_mat_bin"] = np.where(swap_bin, bins.CN_B, bins.CN_A)
    bins["pred_pat_bin"] = np.where(swap_bin, bins.CN_A, bins.CN_B)

    # ---------------- allele-specific: one orientation per cell x bin ----------------
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

    # ---------------- copy-number confusion per cn_mat segment ----------------
    # Each lifted truth piece is compared with the caller's bins by base-pair overlap
    # (no bins filtered); true CN vs called CN, per haplotype.
    # --conf-unit cells: each cell x segment counts once, called as the CN that covers
    #                    most of its compared bases (ties: lowest CN).
    # --conf-unit bp:    every compared base counts (summed over cells, in cell-Mb).
    # Axes: CN 0-4 exactly, then groups of 5 (5-9, 10-14, ...). "exact match" in the
    # panel titles always uses exact values, not the groups.
    LEVELS = [("", "actual", "actual"), ("_bin", "ideal_phasing", "ideal phasing (per bin)"),
              ("_ideal", "allele_specific", "allele-specific")]
    if not PHASED:
        LEVELS = [lv for lv in LEVELS if lv[1] == "allele_specific"]
    if "seg" in truth.columns:
        parts = []
        # index the truth once (a per-cell filter of the whole table is quadratic in cells)
        truth_by = {k: g for k, g in truth.groupby(["cell", "hap", "chrom"], sort=False)}
        for (cell, chrom), cb in bins.groupby(["cell", "chrom"]):
            bs, be = cb.start.values, cb.end.values
            for h in ("mat", "pat"):
                tg = truth_by.get((cell, h, chrom))
                if tg is None or tg.empty:
                    continue
                ts, te = tg.start.values, tg.end.values
                ov = np.clip(np.minimum(te[:, None], be[None, :]) - np.maximum(ts[:, None], bs[None, :]), 0, None)
                i, j = np.nonzero(ov)
                if len(i) == 0:
                    continue
                for x, lvl, _ in LEVELS:
                    pv = cb[f"pred_{h}{x}"].values[j]
                    ok = ~np.isnan(pv.astype(float))
                    d = pd.DataFrame({"seg": tg.seg.values[i][ok], "true_cn": tg.cn.values[i][ok].astype(int),
                                      "pred_cn": pv[ok].astype(int), "bp": ov[i, j][ok].astype(float)})
                    d = d.groupby(["seg", "true_cn", "pred_cn"], as_index=False).bp.sum()
                    d.insert(0, "level", lvl); d.insert(0, "hap", h); d.insert(0, "cell", cell)
                    parts.append(d)
        long = pd.concat(parts, ignore_index=True)      # cell, hap, level, seg, true_cn, pred_cn, bp

        # cells: majority called CN per cell x haplotype x segment x level
        maj = (long.sort_values(["bp", "pred_cn"], ascending=[False, True])
               .drop_duplicates(["cell", "hap", "level", "seg", "true_cn"]))
        key = ["hap", "seg", "level", "true_cn", "pred_cn"]
        conf = (long.groupby(key).bp.sum().to_frame()
                .join(maj.groupby(key).size().rename("cells"), how="outer").fillna(0).reset_index())
        conf["cells"] = conf.cells.astype(int)
        conf[["hap", "seg", "level", "true_cn", "pred_cn", "cells", "bp"]].to_csv(
            os.path.join(a.out_dir, "section_confusion.tsv"), sep="\t", index=False)

        unit = a.conf_unit
        val = "cells" if unit == "cells" else "bp"

        # Two views of the same comparison:
        #   cn     true CN vs called CN (0-4 exact, then groups of 5)
        #   events gain / normal / loss per haplotype (gain = CN > 1, normal = 1, loss = 0)
        STATES = ["loss", "normal", "gain"]

        def group(v, kind):
            v = int(v)
            if kind == "events":
                return "loss" if v < 1 else "normal" if v == 1 else "gain"
            return v if v < 5 else 5 * (v // 5)          # 0..4 exact, then 5-9 -> 5, 10-14 -> 10, ...

        def glabel(g):
            return g if isinstance(g, str) else (f"{g}" if g < 5 else f"{g}-{g + 4}")

        def fmt(v):
            return f"{int(v):,}" if unit == "cells" else f"{v / 1e6:,.1f}"

        def draw(ax, df, cats, title, kind):
            idx = {c: k for k, c in enumerate(cats)}
            m = np.zeros((len(cats), len(cats)))
            for t, p_, v in zip(df.true_cn, df.pred_cn, df[val]):
                m[idx[group(t, kind)], idx[group(p_, kind)]] += v
            tot = m.sum(axis=1, keepdims=True)
            with np.errstate(invalid="ignore", divide="ignore"):
                frac = np.where(tot > 0, m / tot, np.nan)
            ax.imshow(np.where(tot > 0, frac, np.nan), cmap="Blues", vmin=0, vmax=1)
            n = len(cats)
            fs = 7.5 if n <= 6 else 6 if n <= 10 else 4.5 if n <= 16 else 3.2
            for r_ in range(n):
                for c_ in range(n):
                    if m[r_, c_] > 0 or (kind == "events" and tot[r_, 0] > 0):
                        ax.text(c_, r_, f"{100 * frac[r_, c_]:.1f}%\n{fmt(m[r_, c_])}", ha="center", va="center",
                                fontsize=fs, color="white" if frac[r_, c_] > 0.6 else "black")
                    elif kind == "events":
                        ax.text(c_, r_, "-", ha="center", va="center", fontsize=8, color="0.6")
            labs = [glabel(c) for c in cats]
            lfs = 8 if n <= 10 else 6 if n <= 16 else 4.5
            ax.set_xticks(range(n)); ax.set_xticklabels(labs, fontsize=lfs, rotation=90 if n > 5 else 0)
            ax.set_yticks(range(n)); ax.set_yticklabels(labs, fontsize=lfs)
            what = "CN" if kind == "cn" else "state"
            ax.set_xlabel(f"{LABEL} {what}", fontsize=8)
            ax.set_ylabel(f"true {what}", fontsize=8)
            tw = df[val].sum()
            if kind == "cn":                              # exact values, not the groups
                score = df.loc[df.true_cn == df.pred_cn, val].sum() / tw if tw else np.nan
                ax.set_title(f"{title}\nexact match {score:.3f}   (n = {fmt(tw)})", fontsize=8.5)
            else:
                score = np.trace(m) / m.sum() if m.sum() else np.nan
                ax.set_title(f"{title}\nagreement {score:.3f}   (n = {fmt(tw)})", fontsize=8.5)

        how = ("each box: % of the true row, then cells (each cell x segment once, called by the CN "
               "covering most of its compared bases)" if unit == "cells" else
               "each box: % of the true row, then cell-Mb (Mb summed over cells; e.g. 10 cells x 2 Mb = 20)")

        def page(pdf, df, heading, kind):
            haps = [(h, n_) for h, n_ in (("mat", "maternal"), ("pat", "paternal")) if (df.hap == h).any()]
            if not haps:
                return
            nl = len(LEVELS)
            fig, axs = plt.subplots(len(haps), nl, figsize=(13 if nl == 3 else 7.5, 4.9 * len(haps) + 1.2),
                                    squeeze=False)
            for r_, (h, hname) in enumerate(haps):
                dh = df[df.hap == h]                      # one set of axis groups per haplotype row
                if kind == "events":
                    cats = STATES
                else:
                    cats = sorted({group(v, kind) for v in np.concatenate([dh.true_cn.values, dh.pred_cn.values])})
                for c_, (_, lvl, lname) in enumerate(LEVELS):
                    draw(axs[r_, c_], df[(df.hap == h) & (df.level == lvl)], cats, f"{hname} - {lname}", kind)
            sub = ("CN 0-4 exact, then groups of 5" if kind == "cn" else
                   "per haplotype: gain = CN > 1, normal = 1, loss = 0")
            if nl < 3:                                    # narrow figure: wrap the long explanation
                import textwrap
                how_txt = "\n".join(textwrap.wrap(how, 95))
            else:
                how_txt = how
            fig.suptitle(f"{heading}\n{how_txt}\n{sub}", fontsize=9.5 if nl == 3 else 8.5)
            fig.tight_layout(rect=[0, 0, 1, 1 - (0.75 if nl == 3 else 1.0) / fig.get_figheight()], h_pad=2.5)
            pdf.savefig(fig); plt.close(fig)

        def cn_range(v):
            lo, hi = int(v.min()), int(v.max())
            return f"{lo}" if lo == hi else f"{lo}-{hi}"

        segs = sorted(conf.seg.unique(), key=lambda v: int(v.rsplit(":", 1)[1].split("-")[0]))
        for kind, fname in (("cn", "confusion_segments_cn.pdf"), ("events", "confusion_segments_events.pdf")):
            with PdfPages(os.path.join(a.out_dir, fname)) as pdf:
                page(pdf, conf, "All cn_mat segments pooled, all cells", kind)
                for sg in segs:
                    g = maj[(maj.seg == sg) & (maj.level == LEVELS[0][1])]
                    rng = ", ".join(f"{h} CN {cn_range(g[g.hap == h].true_cn)}" for h in ("mat", "pat") if (g.hap == h).any())
                    page(pdf, conf[conf.seg == sg], f"Segment {sg}  ({rng}; {g.cell.nunique()} cells)", kind)
    else:
        print("[eval] truth file has no 'seg' column (made by an older lift_truth.py); skipping confusion_segments_*.pdf",
              file=sys.stderr)

    # ---------------- metrics ----------------
    def hap_events(p, suffix):
        """Gain / normal / loss TP/TN/FP/FN per haplotype (gain = CN > 1, normal = 1, loss = 0), mat and pat pooled."""
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
            "acc_mat_bin": np.mean(p.pred_mat_bin == p.true_mat),
            "acc_pat_bin": np.mean(p.pred_pat_bin == p.true_pat),
            "acc_both_bin": np.mean((p.pred_mat_bin == p.true_mat) & (p.pred_pat_bin == p.true_pat)),
            "acc_mat_ideal": np.mean(p.pred_mat_ideal == p.true_mat),
            "acc_pat_ideal": np.mean(p.pred_pat_ideal == p.true_pat),
            "acc_both_ideal": np.mean((p.pred_mat_ideal == p.true_mat) & (p.pred_pat_ideal == p.true_pat)),
            "acc_total": np.mean(p.pred_total == p.true_total),
            "mae_hap": np.mean((p.pred_mat - p.true_mat).abs() + (p.pred_pat - p.true_pat).abs()) / 2,
            "mae_hap_bin": np.mean((p.pred_mat_bin - p.true_mat).abs() + (p.pred_pat_bin - p.true_pat).abs()) / 2,
            "mae_hap_ideal": np.mean((p.pred_mat_ideal - p.true_mat).abs() + (p.pred_pat_ideal - p.true_pat).abs()) / 2,
            "mae_total": np.mean((p.pred_total - p.true_total).abs()),
            "frac_bins_flipped_bin": np.mean(p.bin_flipped),
            "frac_bins_flipped_ideal": np.mean(p.ideal_flipped),
            **hap_events(p, ""),
            **hap_events(p, "_bin"),
            **hap_events(p, "_ideal"),
            **{f"total_{k}": v for k, v in event_stats(p.true_total, p.pred_total, 2).items()},
        }

    cells = pd.DataFrame([dict(cell=c, n_bins=len(g), n_scored=int(g.scored.sum()),
                               **metrics(g[g.scored]), gamma=g.gamma.iloc[0], swapped=bool(g.swap.iloc[0]))
                          for c, g in bins.groupby("cell")])

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
        orient_note if PHASED else
        "levels: allele-specific only (the caller reports CN_A|CN_B as major|minor, so A/B carries no phase)",
        "",
    ] + ([
        "accuracy over scored bins (exact match)        actual    ideal phasing   allele-specific",
        "                                                          (per bin)       (per cell x bin)",
        f"  maternal                                      {mt['acc_mat']:.3f}     {mt['acc_mat_bin']:.3f}           {mt['acc_mat_ideal']:.3f}",
        f"  paternal                                      {mt['acc_pat']:.3f}     {mt['acc_pat_bin']:.3f}           {mt['acc_pat_ideal']:.3f}",
        f"  both haplotypes                               {mt['acc_both']:.3f}     {mt['acc_both_bin']:.3f}           {mt['acc_both_ideal']:.3f}",
        f"  mean |error| per haplotype                    {mt['mae_hap']:.3f}     {mt['mae_hap_bin']:.3f}           {mt['mae_hap_ideal']:.3f}",
        f"  scored cell-bins re-oriented vs actual                  {mt['frac_bins_flipped_bin']:.3f}           {mt['frac_bins_flipped_ideal']:.3f}",
        f"  total CN (orientation-free)                   {mt['acc_total']:.3f}",
        f"  total CN mean |error|                         {mt['mae_total']:.3f}",
        f"  genomic bins whose per-bin orientation differs from actual: "
        f"{bins[bins.bin_flipped].groupby(['chrom', 'start', 'end']).ngroups}/{bins.groupby(['chrom', 'start', 'end']).ngroups}",
        "  (actual -> ideal phasing = cost of phase switch errors; ideal phasing -> allele-specific = events",
        "   put on different haplotypes in different cells, which phasing cannot fix)",
    ] if PHASED else [
        "accuracy over scored bins (exact match)        allele-specific (sorted major/minor, per cell x bin)",
        f"  maternal                                      {mt['acc_mat_ideal']:.3f}",
        f"  paternal                                      {mt['acc_pat_ideal']:.3f}",
        f"  both haplotypes                               {mt['acc_both_ideal']:.3f}",
        f"  mean |error| per haplotype                    {mt['mae_hap_ideal']:.3f}",
        f"  total CN (orientation-free)                   {mt['acc_total']:.3f}",
        f"  total CN mean |error|                         {mt['mae_total']:.3f}",
        "  (maternal/paternal = the true haplotype each called allele was matched to, whichever fits better)",
    ]) + [
        "",
        "gain / normal / loss over scored bins, one-vs-rest (haplotype: gain = CN > 1, normal = 1, loss = 0;",
        "total: gain > 2, normal = 2, loss < 2). Counts are cell x bin pairs.",
    ]
    scored_bins = bins[bins.scored]

    def ev_line(label, d, ev):
        tp, tn, fp, fn = (d[f"{ev}_{c}"] for c in CONF)
        prec = tp / (tp + fp) if tp + fp else float("nan")
        rec = tp / (tp + fn) if tp + fn else float("nan")
        spec = tn / (tn + fp) if tn + fp else float("nan")
        return f"    {label:18s}{tp:10d}{tn:10d}{fp:10d}{fn:10d}{prec:11.3f}{rec:8.3f}{spec:13.3f}"

    head = f"    {'':18s}{'TP':>10s}{'TN':>10s}{'FP':>10s}{'FN':>10s}{'precision':>11s}{'recall':>8s}{'specificity':>13s}"
    blocks = [(f"{lvl} - {hname}", t_, f"pred_{h}{sfx}", 1)
              for lvl, sfx in ((("actual", ""), ("ideal phasing (per bin)", "_bin"),
                                ("allele-specific (per cell x bin)", "_ideal")) if PHASED else
                               (("allele-specific (per cell x bin)", "_ideal"),))
              for h, hname, t_ in (("mat", "maternal", "true_mat"), ("pat", "paternal", "true_pat"))]
    blocks.append(("total CN (orientation-free)", "true_total", "pred_total", 2))
    for title, t_, p_, nrm in blocks:
        d = event_stats(scored_bins[t_], scored_bins[p_], nrm)
        summ += ["", f"  {title}", head] + [ev_line(ev, d, ev) for ev in EVENTS]
    summ += [
        "",
        "  (one-vs-rest per class: FP = called this class but truth is not; FN = truth is this class but not called)",
        "",
        (f"per-cell both-haplotype accuracy: actual median {cells.acc_both.median():.3f}, "
         f"ideal phasing median {cells.acc_both_bin.median():.3f}, "
         f"allele-specific median {cells.acc_both_ideal.median():.3f}") if PHASED else
        f"per-cell both-haplotype accuracy (allele-specific): median {cells.acc_both_ideal.median():.3f}",
    ]
    if bins.gamma.notna().any():
        summ.append(f"cells {LABEL} called as WGD-scaled (gamma > 3): {(cells.gamma > 3).sum()}/{len(cells)}")
    text = "\n".join(summ)
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)

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
                        if suffix in ("_bin", "_ideal"):
                            fl = cb[cb.bin_flipped if suffix == "_bin" else cb.ideal_flipped]
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

    track_sets = [("per_cell_tracks.pdf", "", "actual orientation"),
                  ("per_cell_tracks_bin.pdf", "_bin", "ideal phasing (per bin)"),
                  ("per_cell_tracks_ideal.pdf", "_ideal", "allele-specific (per cell x bin)")]
    made = []
    for fname, sfx, heading in track_sets:
        path = os.path.join(a.out_dir, fname)
        if PHASED or sfx == "_ideal":
            plot_tracks(path, sfx, heading)
            made.append(fname)
        elif os.path.exists(path):          # from an earlier run with all levels
            os.remove(path)

    print(f"\nwrote {a.out_dir}/summary.txt, confusion_segments_cn.pdf, confusion_segments_events.pdf, section_confusion.tsv, "
          + ", ".join(made))


if __name__ == "__main__":
    main()