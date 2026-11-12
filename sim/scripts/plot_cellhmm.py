#!/usr/bin/env python3
"""Where the states of a per-cell HMM run sit in RDR x BAF, against the data at single-cell and bin resolution.

Works on both per-cell HMM runs:
  --cellhmm K     chisel_hmm's CellHMM: <run>/calls/calls.tsv + <run>/calls/cellhmm.tsv
  --bincellhmm K  RunpengLuo/chisel bin_cell_hmm: <run>/calls/calls.tsv + <run>/calls/bincellhmm.npz
(found automatically from --run). Every point is one cell in one bin: x = RDR, y = BAF = B / (A + B) (phased, not
mirrored), CLUSTER of calls.tsv = the HMM state.

Page 1, overview:
  1. every cell-bin with SNP reads as a 2D histogram; the copy-number grid ((a + b) / ploidy, b / (a + b), grey
     rings), the grid points the states started from (open markers) and where they ended (filled), with an arrow;
  2. segments: runs of consecutive bins of one cell in one state, reads pooled (RDR = sum RDR x NORM_COUNT /
     sum NORM_COUNT, BAF = sum B / sum (A + B)); far less noisy than single cell-bins;
  3. the RDR of the cell-bins decoded to every state, with its fitted mean.
Pages 2+: single cells along the genome, RDR and BAF of every bin coloured by its state, the state's fitted mean
drawn over every run of the state (default: the cells with the most bins outside their most common state, and
some typical ones; or --cells).
Prints per state: its labels, initial grid point, fitted means, pooled RDR / BAF of its cell-bins, share.

  python3 scripts/plot_cellhmm.py --run <PROJECT_DIR>/calls/chisel_pbins8_cellhmm7/ugp
  -> <run>/calls/cellhmm_states.pdf (all pages) and .png (page 1)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

# categorical slots in fixed order; states beyond 8 reuse a hue with another marker shape
HUES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SHAPES = ["o", "s", "^", "D"]
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
INK, MUTED, GRID = "#2b2b2a", "#77766f", "#e4e3de"


def style(k):
    return HUES[k % len(HUES)], SHAPES[(k // len(HUES)) % len(SHAPES)]


def cn_str(v):
    return "{}|{}".format(int(v[0]), int(v[1]))


def read_states(run_calls, states_arg):
    """Per state: CN (label), INIT (initial grid point), MEAN_RDR, MEAN_BAF (NaN where not one shared mean)."""
    path = states_arg
    if path is None:
        for name in ("cellhmm.tsv", "bincellhmm.npz"):
            if os.path.exists(os.path.join(run_calls, name)):
                path = os.path.join(run_calls, name)
                break
    if path is None:
        sys.exit("ERROR: no cellhmm.tsv or bincellhmm.npz in {} (give --states)".format(run_calls))
    if path.endswith(".npz"):
        z = np.load(path, allow_pickle=True)
        st = pd.DataFrame({"CN": [cn_str(v) for v in z["state_cn"]], "INIT": [cn_str(v) for v in z["init_states"]],
                           "MEAN_RDR": z["means"][:, 0], "MEAN_BAF": z["means"][:, 1],
                           "SD_RDR": np.sqrt(z["variances"][:, 0]), "SD_BAF": np.sqrt(z["variances"][:, 1])})
        st.index.name = "STATE"
        return st, path
    st = pd.read_csv(path, sep="\t", comment="#").set_index("STATE")
    st["MEAN_RDR"] = pd.to_numeric(st["MEAN_RDR"], errors="coerce")       # 'gamma*N' with --means grid
    st = st.rename(columns={"INIT_CN": "INIT"})
    return st[["CN", "INIT", "MEAN_RDR", "MEAN_BAF"]], path


def read_calls(path):
    want = {"CHR", "START", "END", "CELL", "NORM_COUNT", "RDR", "A_COUNT", "B_COUNT", "CLUSTER", "CN_STATE"}
    df = pd.read_csv(path, sep="\t", usecols=lambda c: c.lstrip("#") in want,
                     dtype={"#CHR": str, "CHR": str, "CELL": str, "CN_STATE": str})
    df.columns = [c.lstrip("#") for c in df.columns]
    df["N"] = df.A_COUNT + df.B_COUNT
    df["BAF"] = np.where(df.N > 0, df.B_COUNT / df.N.where(df.N > 0, 1), np.nan)
    return df


def chrom_key(c):
    c = str(c)
    c = c[3:] if c.lower().startswith("chr") else c
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24}.get(c.upper(), 99))


def segments(df):
    """Runs of consecutive bins of one cell in one state (within a chromosome), reads pooled."""
    d = df.sort_values(["CELL", "CHR", "START"], kind="mergesort")
    new = (d.CELL.ne(d.CELL.shift()) | d.CHR.ne(d.CHR.shift()) | d.CLUSTER.ne(d.CLUSTER.shift())).to_numpy()
    d = d.assign(seg=np.cumsum(new), rn=d.RDR * d.NORM_COUNT)
    g = d.groupby("seg").agg(state=("CLUSTER", "first"), bins=("CLUSTER", "size"), rn=("rn", "sum"),
                             norm=("NORM_COUNT", "sum"), B=("B_COUNT", "sum"), N=("N", "sum"))
    g = g[(g.norm > 0) & (g.N > 0)]
    return g.assign(rdr=g.rn / g.norm, baf=g.B / g.N)


def pick_cells(df, n, given):
    if given:
        missing = [c for c in given if c not in set(df.CELL)]
        if missing:
            sys.exit("ERROR: cells not in calls.tsv: {}".format(", ".join(missing[:5])))
        return given
    top = df.groupby("CELL").CLUSTER.agg(lambda s: 1.0 - s.value_counts(normalize=True).iloc[0]).sort_values()
    k = max(1, n // 2)                                     # the most aberrant, then evenly spread over the rest
    rest = top.index[:-k]
    spread = list(rest[np.linspace(0, len(rest) - 1, n - k).round().astype(int)]) if len(rest) and n > k else []
    return list(top.index[-k:][::-1]) + [c for c in spread if c not in top.index[-k:]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", help="a --cellhmm / --bincellhmm run folder (with calls/)")
    ap.add_argument("--calls", help="calls.tsv (default: <run>/calls/calls.tsv)")
    ap.add_argument("--states", help="cellhmm.tsv or bincellhmm.npz (default: found in <run>/calls/)")
    ap.add_argument("-o", "--out", help="output prefix (default: <run>/calls/cellhmm_states)")
    ap.add_argument("--ploidy", type=float, default=2.0, help="grid RDR of a|b = (a + b) / ploidy (default: 2)")
    ap.add_argument("--maxcn", type=int, default=6, help="largest a + b on the grid (default: 6)")
    ap.add_argument("--cells", nargs="+", help="cells for the genome pages (default: --ncells picked)")
    ap.add_argument("--ncells", type=int, default=8, help="cells on the genome pages (default: 8, 4 per page)")
    ap.add_argument("--min-seg-bins", type=int, default=3, help="segments with fewer bins are not drawn (default: 3)")
    ap.add_argument("--max-points", type=int, default=100000, help="segments drawn at most (random; default: 100000)")
    ap.add_argument("--rdr-max", type=float, default=None, help="RDR axis limit (default: 99.5th percentile, at least the largest mean)")
    ap.add_argument("--title", default=None)
    a = ap.parse_args()
    if not (a.run or (a.calls and a.states)):
        ap.error("give --run, or --calls and --states")
    calls = a.calls or os.path.join(a.run, "calls", "calls.tsv")
    st, spath = read_states(os.path.dirname(os.path.abspath(calls)), a.states)
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(calls)), "cellhmm_states")
    df = read_calls(calls)
    seg = segments(df)
    K = len(st)
    kind = "bin_cell_hmm (RunpengLuo)" if spath.endswith(".npz") else "CellHMM"

    ok = df.NORM_COUNT > 0
    pooled = pd.DataFrame({
        "rdr": (df.RDR * df.NORM_COUNT)[ok].groupby(df.CLUSTER[ok]).sum() / df.NORM_COUNT[ok].groupby(df.CLUSTER[ok]).sum(),
        "baf": df.B_COUNT.groupby(df.CLUSTER).sum() / df.N.groupby(df.CLUSTER).sum().replace(0, np.nan),
        "share": df.CLUSTER.value_counts(normalize=True)}).reindex(st.index)
    called = df.groupby("CLUSTER").CN_STATE.agg(lambda s: ", ".join("{} {:.0%}".format(k, v) for k, v in
                                                                     s.value_counts(normalize=True).head(3).items()))
    tab = pd.DataFrame({"LABEL": st.CN, "INIT": st.INIT, "MEAN_RDR": st.MEAN_RDR, "MEAN_BAF": st.MEAN_BAF,
                        "POOLED_RDR": pooled.rdr, "POOLED_BAF": pooled.baf, "SHARE": pooled.share,
                        "CALLED_AS": called.reindex(st.index)})
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(tab.to_string(float_format=lambda v: "{:.3f}".format(v)))

    mx = st.MEAN_RDR.where(st.MEAN_RDR.notna(), pooled.rdr)          # where to draw a state
    my = st.MEAN_BAF
    rmax = a.rdr_max or max(float(np.nanpercentile(df.RDR[ok], 99.5)), float(np.nanmax(mx)) * 1.15)
    gpos = {}
    for t in range(1, a.maxcn + 1):
        for x in range(t + 1):
            gpos["{}|{}".format(x, t - x)] = (t / a.ploidy, (t - x) / float(t))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap, LogNorm
    from matplotlib.lines import Line2D
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.titlesize": 10, "axes.titlecolor": INK})
    title = a.title or "{}, K = {}: {}".format(kind, K, os.path.dirname(os.path.dirname(os.path.abspath(calls))))

    def legend_handles():
        h = [Line2D([], [], marker=style(k)[1], color=style(k)[0], linestyle="none", markersize=8, markeredgecolor="white",
                    label="{}: {}  ({:.1%})".format(k, st.CN[k], pooled.share[k] if np.isfinite(pooled.share[k]) else 0))
             for k in st.index]
        h += [Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor=MUTED, linestyle="none", markersize=8,
                     label="grid a|b"),
              Line2D([], [], marker="o", markerfacecolor="none", markeredgecolor=INK, linestyle="none", markersize=8,
                     markeredgewidth=1.5, label="start (k-means)")]
        return h

    # ---------------- page 1
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw={"width_ratios": [1, 1, 1.05]})

    def decorate(axis, ttl):
        axis.set_xlim(0, rmax)
        axis.set_ylim(-0.03, 1.03)
        axis.set_xlabel("RDR")
        axis.set_ylabel("BAF  (B / (A + B))")
        axis.set_title(ttl, loc="left")
        axis.grid(color=GRID, linewidth=0.6)
        axis.set_axisbelow(True)
        for lab, (gx, gy) in gpos.items():
            if gx <= rmax:
                axis.scatter(gx, gy, s=70, facecolors="none", edgecolors=MUTED, linewidths=0.8, zorder=3)
                axis.annotate(lab, (gx, gy), xytext=(6, -10), textcoords="offset points", fontsize=7, color=MUTED, zorder=3)

    def means(axis, arrows=False):
        placed = []                                         # stack labels of states that sit close together
        for k in st.index:
            if not np.isfinite(mx[k]):
                continue
            c, m = style(k)
            if arrows and st.INIT[k] in gpos:
                sx, sy = gpos[st.INIT[k]]
                axis.scatter(sx, sy, s=120, marker=m, facecolors="none", edgecolors=c, linewidths=1.8, zorder=4)
                if np.hypot((mx[k] - sx) / rmax, my[k] - sy) > 0.01:
                    axis.annotate("", xy=(mx[k], my[k]), xytext=(sx, sy), zorder=4,
                                  arrowprops=dict(arrowstyle="-|>", color=c, lw=1.5, shrinkA=6, shrinkB=7))
            axis.scatter(mx[k], my[k], s=150, marker=m, color=c, edgecolors="white", linewidths=2, zorder=5)
            near = sum(1 for (px, py) in placed if abs(px - mx[k]) / rmax < 0.08 and abs(py - my[k]) < 0.06)
            placed.append((mx[k], my[k]))
            axis.annotate("{}: {}".format(k, st.CN[k]), (mx[k], my[k]), xytext=(8, 6 + 11 * near), textcoords="offset points",
                          fontsize=8, color=INK, fontweight="bold", zorder=6)

    sel = ok & (df.N > 0)
    h = ax[0].hist2d(df.RDR[sel].clip(upper=rmax), df.BAF[sel], bins=[160, 101], range=[[0, rmax], [0, 1]],
                     cmap=LinearSegmentedColormap.from_list("blues", BLUES), norm=LogNorm(), zorder=1)
    fig.colorbar(h[3], ax=ax[0], pad=0.01, label="cell-bins")
    decorate(ax[0], "Cell-bins ({:,} with SNP reads): start -> fitted mean".format(int(sel.sum())))
    means(ax[0], arrows=True)

    s = seg[seg.bins >= a.min_seg_bins]
    if len(s) > a.max_points:
        s = s.sample(a.max_points, random_state=0)
    decorate(ax[1], "Segments of one state in one cell, reads pooled ({:,} with >= {} bins)".format(
        int((seg.bins >= a.min_seg_bins).sum()), a.min_seg_bins))
    for k in st.index:
        t = s[s.state == k]
        if len(t):
            c, m = style(k)
            ax[1].scatter(t.rdr.clip(upper=rmax), t.baf, s=4 + 2 * np.sqrt(t.bins), marker=m, color=c, alpha=0.35,
                          linewidths=0, zorder=2)
    means(ax[1])

    edges = np.linspace(0, rmax, 121)
    for k in st.index:
        r = df.RDR[ok & (df.CLUSTER == k)].clip(upper=rmax)
        if not len(r):
            continue
        c, m = style(k)
        cnt, _ = np.histogram(r, bins=edges)
        ax[2].step(edges[:-1], cnt, where="post", color=c, linewidth=2 if k < len(HUES) else 1.2,
                   linestyle="-" if k < len(HUES) else "--")
        if np.isfinite(mx[k]):
            ax[2].axvline(mx[k], color=c, linewidth=1, linestyle=":")
    ax[2].set_yscale("log")
    ax[2].set_xlim(0, rmax)
    ax[2].set_xlabel("RDR")
    ax[2].set_ylabel("cell-bins")
    ax[2].set_title("RDR of the cell-bins of each state (dotted: fitted mean)", loc="left")
    ax[2].grid(color=GRID, linewidth=0.6)
    ax[2].set_axisbelow(True)
    fig.legend(handles=legend_handles(), loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
               title="state: label (share)")
    fig.suptitle(title, x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pages = [fig]

    # ---------------- pages 2+: single cells along the genome
    cells = pick_cells(df, a.ncells, a.cells)
    chroms = sorted(df.CHR.unique(), key=chrom_key)
    span = df.groupby("CHR").END.max()
    off, acc = {}, 0.0
    for c in chroms:
        off[c] = acc
        acc += float(span[c])
    sub = df[df.CELL.isin(cells)].copy()
    sub["x"] = (sub.CHR.map(off) + (sub.START + sub.END) / 2.0) / 1e6
    sub["x0"], sub["x1"] = (sub.CHR.map(off) + sub.START) / 1e6, (sub.CHR.map(off) + sub.END) / 1e6
    for p0 in range(0, len(cells), 4):
        group = cells[p0:p0 + 4]
        fig, axes = plt.subplots(2 * len(group), 1, figsize=(15, 2.3 * 2 * len(group)), sharex=True, squeeze=False)
        for i, cell in enumerate(group):
            d = sub[sub.CELL == cell].sort_values(["x"])
            ar, ab = axes[2 * i, 0], axes[2 * i + 1, 0]
            for k in st.index:
                t = d[d.CLUSTER == k]
                if not len(t):
                    continue
                c, m = style(k)
                ar.scatter(t.x, t.RDR.clip(upper=rmax), s=6, marker=m, color=c, linewidths=0, alpha=0.8, zorder=2)
                ab.scatter(t.x, t.BAF, s=6, marker=m, color=c, linewidths=0, alpha=0.8, zorder=2)
            run = (d.CLUSTER.ne(d.CLUSTER.shift()) | d.CHR.ne(d.CHR.shift())).cumsum()
            for _, r in d.groupby(run):
                k = int(r.CLUSTER.iloc[0])
                c = style(k)[0]
                if np.isfinite(st.MEAN_RDR[k]):
                    ar.hlines(st.MEAN_RDR[k], r.x0.min(), r.x1.max(), color=INK, linewidth=2.2, zorder=3)
                    ar.hlines(st.MEAN_RDR[k], r.x0.min(), r.x1.max(), color=c, linewidth=1.4, zorder=4)
                ab.hlines(st.MEAN_BAF[k], r.x0.min(), r.x1.max(), color=INK, linewidth=2.2, zorder=3)
                ab.hlines(st.MEAN_BAF[k], r.x0.min(), r.x1.max(), color=c, linewidth=1.4, zorder=4)
            states_here = d.CLUSTER.value_counts(normalize=True)
            ar.set_title("{}: states {}".format(cell, ", ".join("{} ({}) {:.0%}".format(k, st.CN[k], v)
                                                             for k, v in states_here.head(5).items())), loc="left")
            ar.set_ylim(0, rmax)
            ab.set_ylim(-0.03, 1.03)
            ar.set_ylabel("RDR")
            ab.set_ylabel("BAF")
            for axis in (ar, ab):
                axis.grid(color=GRID, linewidth=0.6)
                axis.set_axisbelow(True)
                for c in chroms[1:]:
                    axis.axvline(off[c] / 1e6, color=MUTED, linewidth=0.8)
        axes[-1, 0].set_xlabel("genome position (Mb){}".format("" if len(chroms) == 1 else ", chromosomes " + " ".join(chroms)))
        fig.legend(handles=legend_handles()[:K], loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
                   title="state: label (share)")
        fig.suptitle("Single cells: every bin coloured by its state; lines = the state's fitted mean over each run",
                     x=0.01, ha="left", fontsize=11, color=INK)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        pages.append(fig)

    pages[0].savefig(out + ".png", dpi=150, bbox_inches="tight")
    with PdfPages(out + ".pdf") as pdf:
        for f in pages:
            pdf.savefig(f, bbox_inches="tight")
    print("wrote {0}.pdf ({1} pages) and {0}.png (overview)".format(out, len(pages)), file=sys.stderr)


if __name__ == "__main__":
    main()
