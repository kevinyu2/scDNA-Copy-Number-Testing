#!/usr/bin/env python3
"""
Compare HiScanner haplotype-specific copy numbers with simulated truth.

Inputs
  --truth      truth_hg38.tsv from lift_truth.py (aligned pieces only; no interpolation)
  --calls-dir  HiScanner output/final_calls (per-cell <cell>.txt with CN_A, CN_B)

Truth for a bin (used only for scoring and the heatmaps; the track plots draw
the truth at its own breakpoints)
  Of the truth that aligns inside the bin, the copy number covering the most
  bases is the bin's true value.
    purity_*   fraction of the bin's aligned truth with that value
               (low = the bin straddles a real breakpoint)
    covered_*  fraction of the bin that has any aligned truth
  A bin is scored when purity >= --purity and covered >= --min-covered on
  both haplotypes. Unaligned gaps do not count against purity.

Haplotype orientation (which HiScanner haplotype is maternal)
  actual  one assignment for all cells (--orientation global), or one per
          cell (--orientation cell). Phase switch errors count as errors.
  ideal   "ideal phasing interpretation": for every cell and bin separately,
          whichever assignment gives the smaller |mat diff| + |pat diff|.
          Upper bound on what the same calls could score if every phase
          switch were corrected.

Outputs (in --out-dir)
  per_cell_tracks.pdf     maternal/paternal predicted vs true along the chromosome
  confusion_all_cells.pdf true vs predicted per bin, all cells pooled (actual + ideal)
  bins.tsv, cells.tsv, summary.txt
"""
import argparse
import glob
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


def load_calls(path):
    df = pd.read_csv(path, sep="\t")
    if not {"CHROM", "START", "END", "CN_A", "CN_B"}.issubset(df.columns):
        return None
    df = df.dropna(subset=["CN_A", "CN_B"]).copy()
    df["CHROM"] = df["CHROM"].map(norm_chrom)
    return df


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--purity", type=float, default=0.9)
    ap.add_argument("--min-covered", type=float, default=0.5, help="min fraction of a bin with aligned truth to score it")
    ap.add_argument("--orientation", choices=["global", "cell"], default="global")
    ap.add_argument("--cells-per-page", type=int, default=4)
    ap.add_argument("--max-cells-plot", type=int, default=0)
    ap.add_argument("--max-cn", type=int, default=8, help="heatmap axis cap; larger values lumped")
    ap.add_argument("--y-linear-max", type=float, default=8, help="track plots: CN above this goes above an axis break")
    ap.add_argument("--merge-gap", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    truth = pd.read_csv(a.truth, sep="\t")
    truth["cell"] = truth["cell"].map(norm_cell)
    truth["chrom"] = truth["chrom"].map(norm_chrom)
    tgroups = {k: merge_equal_neighbours(g, a.merge_gap) for k, g in truth.groupby(["cell", "hap", "chrom"])}
    truth_cells = set(truth.cell)

    files = sorted(f for f in glob.glob(os.path.join(a.calls_dir, "*.txt"))
                   if not os.path.basename(f).endswith(("_seg.txt", "_seg_merged.txt", "_seg_initial_anno.txt")))
    rows, no_calls, no_truth = [], [], []
    for f in files:
        cell = norm_cell(os.path.basename(f)[:-4])
        calls = load_calls(f)
        if calls is None or calls.empty:
            no_calls.append(cell); continue
        if cell not in truth_cells:
            no_truth.append(cell); continue
        for chrom, cb in calls.groupby("CHROM"):
            b = cb[["START", "END"]].values.astype(np.int64)
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
        sys.exit("ERROR: no cells with both HiScanner calls and truth (check cell / chromosome names)")
    bins = pd.concat(rows, ignore_index=True)

    has_truth = bins.true_mat.notna() & bins.true_pat.notna()
    covered_ok = (bins.covered_mat >= a.min_covered) & (bins.covered_pat >= a.min_covered)
    straddles = (bins.purity_mat < a.purity) | (bins.purity_pat < a.purity)
    bins["scored"] = has_truth & covered_ok & ~straddles
    bins["straddles"] = has_truth & straddles

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
    d_keep = (bins.pred_mat - bins.true_mat).abs() + (bins.pred_pat - bins.true_pat).abs()
    d_flip = (bins.pred_pat - bins.true_mat).abs() + (bins.pred_mat - bins.true_pat).abs()
    flip = d_flip < d_keep                      # ties keep the actual orientation
    bins["ideal_flipped"] = flip
    bins["pred_mat_ideal"] = np.where(flip, bins.pred_pat, bins.pred_mat)
    bins["pred_pat_ideal"] = np.where(flip, bins.pred_mat, bins.pred_pat)

    for k in ("true", "pred"):
        bins[f"{k}_total"] = bins[f"{k}_mat"] + bins[f"{k}_pat"]
    bins["truemean_total"] = bins.truemean_mat + bins.truemean_pat
    bins["true_minor"] = bins[["true_mat", "true_pat"]].min(axis=1)
    bins["pred_minor"] = bins[["pred_mat", "pred_pat"]].min(axis=1)
    bins.to_csv(os.path.join(a.out_dir, "bins.tsv"), sep="\t", index=False)

    # ---------------- metrics ----------------
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
        }

    cells = pd.DataFrame([dict(cell=c, n_bins=len(g), n_scored=int(g.scored.sum()),
                               **metrics(g[g.scored]), gamma=g.gamma.iloc[0], swapped=bool(g.swap.iloc[0]))
                          for c, g in bins.groupby("cell")])
    cells.to_csv(os.path.join(a.out_dir, "cells.tsv"), sep="\t", index=False)

    mt = metrics(bins[bins.scored])
    n_bins, n_sc = len(bins), int(bins.scored.sum())
    summ = [
        f"cells evaluated: {bins.cell.nunique()}",
        f"cells without HiScanner calls (skipped): {len(no_calls)}" + (f"  e.g. {no_calls[:5]}" if no_calls else ""),
        f"cells with calls but no truth (skipped): {len(no_truth)}" + (f"  e.g. {no_truth[:5]}" if no_truth else ""),
        f"bins: {n_bins}   scored: {n_sc} ({n_sc / n_bins:.1%})",
        f"  not scored, straddle a true breakpoint: {int(bins.straddles.sum())}",
        f"  not scored, < {a.min_covered:.0%} of bin has aligned truth: {int((has_truth & ~covered_ok & ~straddles).sum() + (~has_truth).sum())}",
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
        f"per-cell both-haplotype accuracy: actual median {cells.acc_both.median():.3f}, "
        f"ideal median {cells.acc_both_ideal.median():.3f}",
    ]
    if bins.gamma.notna().any():
        summ.append(f"cells HiScanner called as WGD-scaled (gamma > 3): {(cells.gamma > 3).sum()}/{len(cells)}")
    text = "\n".join(summ)
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)

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
        ax.set_xlabel("true copy number"); ax.set_ylabel("HiScanner copy number")
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
            ax.set_ylabel("HiScanner CN"); ax.set_title(h)
        fig.suptitle("All bins - fractional truth for bins that straddle breakpoints", fontsize=10)
        pdf.savefig(fig); plt.close(fig)

    # ---------------- per-cell tracks ----------------
    plot_cells = sorted(cells.cell)
    if a.max_cells_plot > 0:
        plot_cells = plot_cells[: a.max_cells_plot]
    cst = cells.set_index("cell")
    cpp = max(1, a.cells_per_page)
    with PdfPages(os.path.join(a.out_dir, "per_cell_tracks.pdf")) as pdf:
        for k in range(0, len(plot_cells), cpp):
            chunk = plot_cells[k:k + cpp]
            fig, axs = plt.subplots(2 * len(chunk), 1, figsize=(12.5, 4.0 * len(chunk)), squeeze=False)
            for i, cell in enumerate(chunk):
                cb = bins[bins.cell == cell].sort_values("start")
                chrom = cb.chrom.iloc[0]
                st = cst.loc[cell]
                for j, (h, name, col) in enumerate((("mat", "maternal", "tab:red"), ("pat", "paternal", "tab:blue"))):
                    ax = axs[2 * i + j, 0]
                    tg = tgroups.get((cell, h, chrom), pd.DataFrame(columns=["start", "end", "cn"]))
                    ys = YScale(list(tg.cn.values) + list(cb[f"pred_{h}"].values), a.y_linear_max)
                    tx, ty = step_xy(tg.start.values, tg.end.values, tg.cn.values.astype(float), a.merge_gap)
                    px, py = step_xy(cb.start.values, cb.end.values, cb[f"pred_{h}"].values.astype(float), 1000)
                    ax.plot(tx, ys(ty), color="black", lw=6, alpha=0.25, solid_capstyle="butt", label="true (aligned regions)")
                    ax.plot(px, ys(py), color=col, lw=1.3, label="HiScanner")
                    first = True
                    for _, r in cb[cb.straddles].iterrows():
                        ax.axvspan(r.start / 1e6, r.end / 1e6, color="grey", alpha=0.15, lw=0,
                                   label="straddles true breakpoint (not scored)" if first else None)
                        first = False
                    ys.decorate(ax)
                    ax.set_ylabel(f"{name}\nCN", fontsize=8)
                    ax.grid(axis="y", alpha=0.3)
                    ax.legend(loc="upper left", fontsize=6, ncol=3)
                    if j == 0:
                        ax.set_title(f"{cell}   both-hap acc {st.get('acc_both', np.nan):.2f} "
                                     f"(ideal phasing {st.get('acc_both_ideal', np.nan):.2f})   "
                                     f"total acc {st.get('acc_total', np.nan):.2f}   gamma {st.gamma:.2f}",
                                     fontsize=9, loc="left")
                    else:
                        ax.set_xlabel(f"{chrom} position (Mb, hg38)", fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

    print(f"\nwrote {a.out_dir}/per_cell_tracks.pdf, confusion_all_cells.pdf, bins.tsv, cells.tsv, summary.txt")


if __name__ == "__main__":
    main()