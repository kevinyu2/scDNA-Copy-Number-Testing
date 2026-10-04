#!/usr/bin/env python3
"""
Compare HiScanner haplotype-specific copy numbers with simulated truth.

Inputs
  --truth      truth_hg38.tsv from lift_truth.py
  --calls-dir  HiScanner output/final_calls (per-cell <cell>.txt with CN_A, CN_B)

What counts as "true" for a bin
  Each bin is compared against the true segments it overlaps:
    true_*      majority state: CN of the segment covering most of the bin
    purity_*    fraction of the bin covered by that majority segment
    truemean_*  length-weighted mean CN over the bin
  Bins whose purity is below --purity (they straddle a true breakpoint) are
  kept in the table but left out of the accuracy numbers and the main
  confusion matrix; a second page shows all bins.

Haplotype orientation
  HiScanner's A/B haplotypes come from population phasing, so which one is
  maternal is arbitrary. --orientation global (default) picks one A/B->mat/pat
  assignment for all cells, which is honest because every cell shares the
  same phasing; local phase switch errors then show up as errors.
  --orientation cell picks the best assignment per cell (more forgiving).
  Allele-specific (major/minor) agreement is also reported; it needs no
  orientation at all.

Outputs (in --out-dir)
  per_cell_tracks.pdf     maternal/paternal predicted vs true along the chromosome
  confusion_all_cells.pdf true vs predicted per bin, all cells pooled
  bins.tsv                per-bin comparison table
  cells.tsv               per-cell summary
  summary.txt             overall numbers
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


def norm_cell(name):
    return name[:-4] if name.endswith("_sim") else name


def norm_chrom(c):
    c = str(c)
    return c if c.startswith("chr") else f"chr{c}"


def bin_truth(bins, segs):
    """bins: (n,2) starts/ends; segs: DataFrame start,end,cn -> majority, purity, mean, covered."""
    bs, be = bins[:, 0][:, None], bins[:, 1][:, None]
    ss, se = segs.start.values[None, :], segs.end.values[None, :]
    ov = np.clip(np.minimum(be, se) - np.maximum(bs, ss), 0, None).astype(float)   # n x m
    blen = (bins[:, 1] - bins[:, 0]).astype(float)
    covered = ov.sum(axis=1)
    cn = segs.cn.values.astype(float)
    best = ov.argmax(axis=1)
    majority = np.where(covered > 0, cn[best], np.nan)
    purity = np.where((blen > 0) & (covered > 0), ov[np.arange(len(best)), best] / blen, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(covered > 0, (ov * cn).sum(axis=1) / covered, np.nan)
    return majority, purity, mean, covered / np.maximum(blen, 1)


def load_calls(path):
    df = pd.read_csv(path, sep="\t")
    need = {"CHROM", "START", "END", "CN_A", "CN_B"}
    if not need.issubset(df.columns):
        return None
    df = df.dropna(subset=["CN_A", "CN_B"]).copy()
    df["CHROM"] = df["CHROM"].map(norm_chrom)
    return df


def step_xy(starts, ends, vals):
    x = np.empty(2 * len(starts)); y = np.empty(2 * len(starts))
    x[0::2], x[1::2] = starts, ends
    y[0::2], y[1::2] = vals, vals
    return x / 1e6, y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--calls-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--purity", type=float, default=0.9, help="min fraction of a bin in one true segment to score it")
    ap.add_argument("--orientation", choices=["global", "cell"], default="global")
    ap.add_argument("--cells-per-page", type=int, default=4)
    ap.add_argument("--max-cells-plot", type=int, default=0, help="limit per-cell track plots (0 = all)")
    ap.add_argument("--max-cn", type=int, default=8, help="confusion-matrix axis cap; larger values are lumped")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    truth = pd.read_csv(a.truth, sep="\t")
    truth["cell"] = truth["cell"].map(norm_cell)
    truth["chrom"] = truth["chrom"].map(norm_chrom)
    tgroups = {k: g.sort_values("start") for k, g in truth.groupby(["cell", "hap", "chrom"])}
    truth_cells = set(truth.cell)

    files = sorted(f for f in glob.glob(os.path.join(a.calls_dir, "*.txt"))
                   if not os.path.basename(f).endswith(("_seg.txt", "_seg_merged.txt", "_seg_initial_anno.txt")))
    rows, skipped_nocalls, skipped_notruth = [], [], []
    for f in files:
        cell = norm_cell(os.path.basename(f)[:-4])
        calls = load_calls(f)
        if calls is None or calls.empty:
            skipped_nocalls.append(cell); continue
        if cell not in truth_cells:
            skipped_notruth.append(cell); continue
        for chrom, cb in calls.groupby("CHROM"):
            tm, tp = tgroups.get((cell, "mat", chrom)), tgroups.get((cell, "pat", chrom))
            if tm is None or tp is None:
                continue
            b = cb[["START", "END"]].values.astype(np.int64)
            m_maj, m_pur, m_mean, m_cov = bin_truth(b, tm)
            p_maj, p_pur, p_mean, p_cov = bin_truth(b, tp)
            d = pd.DataFrame({
                "cell": cell, "chrom": chrom, "start": b[:, 0], "end": b[:, 1],
                "CN_A": cb.CN_A.values.astype(int), "CN_B": cb.CN_B.values.astype(int),
                "true_mat": m_maj, "purity_mat": m_pur, "truemean_mat": m_mean,
                "true_pat": p_maj, "purity_pat": p_pur, "truemean_pat": p_mean,
                "gamma": cb["gamma"].values if "gamma" in cb.columns else np.nan,
            })
            rows.append(d)
    if not rows:
        sys.exit("ERROR: no cells with both HiScanner calls and truth (check cell names / chromosome names)")
    bins = pd.concat(rows, ignore_index=True)
    bins["pure"] = (bins.purity_mat >= a.purity) & (bins.purity_pat >= a.purity)
    no_truth = bins.true_mat.isna() | bins.true_pat.isna()

    # ---------------- orientation ----------------
    def agree(df, swap):
        pm, pp = (df.CN_B, df.CN_A) if swap else (df.CN_A, df.CN_B)
        return ((pm == df.true_mat) & (pp == df.true_pat)).mean()

    pure = bins[bins.pure]
    if a.orientation == "global":
        s0, s1 = agree(pure, False), agree(pure, True)
        swap = s1 > s0
        bins["swap"] = swap
        orient_note = (f"global orientation: {'A=pat, B=mat' if swap else 'A=mat, B=pat'} "
                       f"(exact-match rate {max(s0, s1):.3f} vs {min(s0, s1):.3f} for the other way)")
    else:
        sw = {c: agree(g, True) > agree(g, False) for c, g in pure.groupby("cell")}
        bins["swap"] = bins.cell.map(sw).fillna(False).astype(bool)
        orient_note = f"per-cell orientation: {sum(sw.values())}/{len(sw)} cells use A=pat"
    bins["pred_mat"] = np.where(bins.swap, bins.CN_B, bins.CN_A)
    bins["pred_pat"] = np.where(bins.swap, bins.CN_A, bins.CN_B)
    bins["pred_total"] = bins.pred_mat + bins.pred_pat
    bins["true_total"] = bins.true_mat + bins.true_pat
    bins["truemean_total"] = bins.truemean_mat + bins.truemean_pat
    bins["pred_major"] = bins[["pred_mat", "pred_pat"]].max(axis=1)
    bins["pred_minor"] = bins[["pred_mat", "pred_pat"]].min(axis=1)
    bins["true_major"] = bins[["true_mat", "true_pat"]].max(axis=1)
    bins["true_minor"] = bins[["true_mat", "true_pat"]].min(axis=1)
    bins.to_csv(os.path.join(a.out_dir, "bins.tsv"), sep="\t", index=False)

    # ---------------- per-cell summary ----------------
    def cell_stats(g):
        p = g[g.pure]
        n = len(p)
        f = lambda s: s.mean() if n else np.nan
        return pd.Series({
            "n_bins": len(g), "n_scored": n,
            "acc_mat": f(p.pred_mat == p.true_mat),
            "acc_pat": f(p.pred_pat == p.true_pat),
            "acc_both": f((p.pred_mat == p.true_mat) & (p.pred_pat == p.true_pat)),
            "acc_total": f(p.pred_total == p.true_total),
            "acc_allele_specific": f((p.pred_major == p.true_major) & (p.pred_minor == p.true_minor)),
            "mae_total": f((p.pred_total - p.true_total).abs()),
            "gamma": g.gamma.iloc[0],
            "swapped": bool(g.swap.iloc[0]),
        })
    cells = pd.DataFrame([dict(cell=c, **cell_stats(g)) for c, g in bins.groupby("cell")])
    cells.to_csv(os.path.join(a.out_dir, "cells.tsv"), sep="\t", index=False)

    p = bins[bins.pure]
    summ = [
        f"cells evaluated: {bins.cell.nunique()}",
        f"cells without HiScanner calls (skipped): {len(skipped_nocalls)}" + (f"  e.g. {skipped_nocalls[:5]}" if skipped_nocalls else ""),
        f"cells with calls but no truth (skipped): {len(skipped_notruth)}" + (f"  e.g. {skipped_notruth[:5]}" if skipped_notruth else ""),
        f"bins: {len(bins)}  scored (purity >= {a.purity}): {len(p)}  "
        f"mixed (straddle a true breakpoint): {len(bins) - len(p) - no_truth.sum()} "
        f"({(len(bins) - len(p) - no_truth.sum()) / len(bins):.1%})  outside truth span: {no_truth.sum()}",
        orient_note,
        "",
        "accuracy over scored bins (exact match):",
        f"  maternal          {np.mean(p.pred_mat == p.true_mat):.3f}",
        f"  paternal          {np.mean(p.pred_pat == p.true_pat):.3f}",
        f"  both haplotypes   {np.mean((p.pred_mat == p.true_mat) & (p.pred_pat == p.true_pat)):.3f}",
        f"  total CN          {np.mean(p.pred_total == p.true_total):.3f}",
        f"  allele-specific   {np.mean((p.pred_major == p.true_major) & (p.pred_minor == p.true_minor)):.3f}  (major/minor, no orientation needed)",
        f"  total CN MAE      {np.mean((p.pred_total - p.true_total).abs()):.3f}",
        "",
        f"per-cell both-haplotype accuracy: median {cells.acc_both.median():.3f}, "
        f"10th pct {cells.acc_both.quantile(0.1):.3f}, min {cells.acc_both.min():.3f}",
    ]
    if bins.gamma.notna().any():
        hi = (cells.gamma > 3).sum()
        summ.append(f"cells HiScanner called as WGD-scaled (gamma > 3): {hi}/{len(cells)}")
    text = "\n".join(summ)
    with open(os.path.join(a.out_dir, "summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)

    # ---------------- confusion matrices ----------------
    cap = a.max_cn

    def confusion(ax, t, pr, title, vmax):
        ok = t.notna() & pr.notna()          # bins outside the truth's span have no truth
        t = np.clip(t[ok].astype(int), 0, cap); pr = np.clip(pr[ok].astype(int), 0, cap)
        m = np.zeros((cap + 1, cap + 1), dtype=np.int64)
        np.add.at(m, (pr, t), 1)
        im = ax.imshow(np.where(m > 0, m, np.nan), origin="lower", cmap="viridis",
                       norm=LogNorm(vmin=1, vmax=vmax))
        for i in range(cap + 1):
            for j in range(cap + 1):
                if m[i, j]:
                    ax.text(j, i, f"{m[i, j]}", ha="center", va="center", fontsize=6,
                            color="white" if m[i, j] < m.max() ** 0.5 else "black")
        lab = [str(i) for i in range(cap)] + [f"{cap}+"]
        ax.set_xticks(range(cap + 1)); ax.set_xticklabels(lab)
        ax.set_yticks(range(cap + 1)); ax.set_yticklabels(lab)
        ax.plot([-0.5, cap + 0.5], [-0.5, cap + 0.5], color="red", lw=0.8, ls="--")
        acc = np.trace(m) / max(1, m.sum())
        ax.set_title(f"{title}\nexact match {acc:.3f} (n={m.sum()})", fontsize=9)
        ax.set_xlabel("true copy number"); ax.set_ylabel("HiScanner copy number")
        return im

    with PdfPages(os.path.join(a.out_dir, "confusion_all_cells.pdf")) as pdf:
        for subset, label in ((bins[bins.pure], f"bins with purity >= {a.purity}"),
                              (bins, "all bins (majority truth, includes breakpoint bins)")):
            fig, axs = plt.subplots(1, 4, figsize=(20, 5.2))
            vmax = max(1, len(subset))       # shared color scale across panels
            confusion(axs[0], subset.true_mat, subset.pred_mat, "maternal", vmax)
            confusion(axs[1], subset.true_pat, subset.pred_pat, "paternal", vmax)
            confusion(axs[2], subset.true_total, subset.pred_total, "total", vmax)
            im = confusion(axs[3], subset.true_minor, subset.pred_minor, "minor allele (orientation-free)", vmax)
            fig.colorbar(im, ax=axs, shrink=0.8, label="bins (log)")
            fig.suptitle(f"All cells pooled - {label}\n{orient_note}", fontsize=10)
            pdf.savefig(fig); plt.close(fig)

        # Continuous view: length-weighted mean truth vs prediction
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        for ax, h in zip(axs, ("mat", "pat", "total")):
            x, y = bins[f"truemean_{h}"], bins[f"pred_{h}"]
            ok = x.notna()
            ax.hist2d(x[ok].clip(0, cap), y[ok].clip(0, cap), bins=[np.arange(0, cap + 0.26, 0.25), np.arange(-0.5, cap + 1)],
                      norm=LogNorm(), cmap="viridis")
            ax.plot([0, cap], [0, cap], color="red", lw=0.8, ls="--")
            ax.set_xlabel("true CN (length-weighted mean over bin)"); ax.set_ylabel("HiScanner CN")
            ax.set_title(h)
        fig.suptitle("All bins - fractional truth for bins that straddle breakpoints", fontsize=10)
        pdf.savefig(fig); plt.close(fig)

    # ---------------- per-cell tracks ----------------
    plot_cells = cells.sort_values("cell").cell.tolist()
    if a.max_cells_plot > 0:
        plot_cells = plot_cells[: a.max_cells_plot]
    cstats = cells.set_index("cell")
    cpp = max(1, a.cells_per_page)
    with PdfPages(os.path.join(a.out_dir, "per_cell_tracks.pdf")) as pdf:
        for k in range(0, len(plot_cells), cpp):
            chunk = plot_cells[k:k + cpp]
            fig, axs = plt.subplots(2 * len(chunk), 1, figsize=(12, 2.0 * 2 * len(chunk)), squeeze=False)
            for i, cell in enumerate(chunk):
                cb = bins[bins.cell == cell].sort_values("start")
                chrom = cb.chrom.iloc[0]
                st = cstats.loc[cell]
                for j, (h, name, col) in enumerate((("mat", "maternal", "tab:red"), ("pat", "paternal", "tab:blue"))):
                    ax = axs[2 * i + j, 0]
                    tg = tgroups[(cell, h, chrom)]
                    tx, ty = step_xy(tg.start.values, tg.end.values, tg.cn.values)
                    px, py = step_xy(cb.start.values, cb.end.values, cb[f"pred_{h}"].values)
                    # Truth as a wide translucent band, prediction as a thin line on top:
                    # agreement = line runs inside the band
                    ax.plot(tx, ty, color="black", lw=6, alpha=0.25, solid_capstyle="butt", label="true")
                    ax.plot(px, py, color=col, lw=1.3, label="HiScanner")
                    # Bins that straddle a true breakpoint
                    mixed = cb[cb[f"purity_{h}"].notna() & (cb[f"purity_{h}"] < a.purity)]
                    for _, r in mixed.iterrows():
                        ax.axvspan(r.start / 1e6, r.end / 1e6, color="grey", alpha=0.15, lw=0)
                    ymax = max(4, int(np.nanmax(cb[f"pred_{h}"])), int(np.nanmax(cb[f"true_{h}"]))) + 1
                    ax.set_ylim(-0.3, ymax + 0.3)
                    ax.set_yticks(range(0, ymax + 1, 1 if ymax <= 8 else 2))
                    ax.set_ylabel(f"{name}\nCN", fontsize=8)
                    ax.grid(axis="y", alpha=0.3)
                    if j == 0:
                        ax.set_title(f"{cell}   both-hap acc {st.acc_both:.2f}   total acc {st.acc_total:.2f}"
                                     f"   gamma {st.gamma:.2f}", fontsize=9, loc="left")
                        ax.legend(loc="upper right", fontsize=7, ncol=2)
                    if j == 1:
                        ax.set_xlabel(f"{chrom} position (Mb, hg38)", fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig); plt.close(fig)

    print(f"\nwrote {a.out_dir}/per_cell_tracks.pdf, confusion_all_cells.pdf, bins.tsv, cells.tsv, summary.txt")


if __name__ == "__main__":
    main()