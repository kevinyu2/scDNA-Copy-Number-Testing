#!/usr/bin/env python3
"""Detection of the planted events: how often is each event found, by subclone size x CNA size.

Reads section_confusion.tsv (written by eval_calls.py) and the tree's design table. For every event, among the
carrying cells (true CN on the event haplotype = the event's CN) and only bins with >= --seg-min-frac (the
`min75` rows) of their truth inside the segment, it reports the percentage whose called CN shows the event:
AMP: called CN >= true CN of the carriers; DEL: called CN <= true CN of the carriers.
Then a matrix: x = subclone size (cells carrying the event, pooled over batches), y = CNA size (bp), one colour
and one label per cell = mean over the events at that (size, subclone), with the number of events.

  eval_detection.py --eval-dir eval/ugp_chisel_pbins8 --design data/trees/sens_amp_b20.design.tsv
  eval_detection.py --compare eval/A eval/B eval/C --names oracle hmm plain --design ... -o compare   # side by side

--level   allele_specific (default; each cell x bin oriented to the truth, so phasing does not matter),
          actual, ideal_phasing
--metric  bins (default: % of cell x bin pairs), bp (% of bases), cells (% of cells, by majority call)
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

SEG = re.compile(r"^(?:chr)?(\w+):(\d+)-(\d+)$")


def per_event(conf_path, design, level, tag, metric, tol=2):
    conf = pd.read_csv(conf_path, sep="\t")
    conf = conf[(conf.bins == tag) & (conf.level == level)].copy()
    if conf.empty:
        sys.exit("ERROR: no rows for bins={} level={} in {}".format(tag, level, conf_path))
    m = conf.seg.str.extract(SEG)
    conf["c"], conf["s"], conf["e"] = m[0], pd.to_numeric(m[1]), pd.to_numeric(m[2])
    col = {"bins": "cell_bins", "bp": "bp", "cells": "cells"}[metric]
    hap_of = {"M": "mat", "P": "pat"}
    out = []
    for ev in design.itertuples():
        g = conf[(conf.hap == hap_of[ev.hap]) & (conf.c == str(ev.chrom).replace("chr", ""))
                 & (conf.s >= ev.start - tol) & (conf.e <= ev.end + 1 + tol)]
        if g.empty:
            out.append(dict(event_id=ev.event_id, n=0, found=0.0, pct=np.nan))
            continue
        amp = ev.type == "AMP"
        cn = g.true_cn.max() if amp else g.true_cn.min()          # the carriers' CN
        car = g[g.true_cn == cn]
        tot = car[col].sum()
        hit = car[(car.pred_cn >= cn) if amp else (car.pred_cn <= cn)][col].sum()
        out.append(dict(event_id=ev.event_id, n=float(tot), found=float(hit), pct=100.0 * hit / tot if tot > 0 else np.nan))
    return design.merge(pd.DataFrame(out), on="event_id")


def matrix(ev):
    xs = sorted(ev.clade_cells.unique())
    ys = sorted(ev.size_bp.unique())
    z = np.full((len(ys), len(xs)), np.nan)
    cnt = np.zeros((len(ys), len(xs)), dtype=int)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            s = ev[(ev.size_bp == y) & (ev.clade_cells == x)].pct.dropna()
            cnt[i, j] = int(((ev.size_bp == y) & (ev.clade_cells == x)).sum())
            if len(s):
                z[i, j] = s.mean()
    return xs, ys, z, cnt


def lab_size(v):
    return "{:g} Mb".format(v / 1e6) if v >= 1e6 else "{:g} kb".format(v / 1e3)


def draw(ax, ev, title):
    xs, ys, z, cnt = matrix(ev)
    im = ax.imshow(z, origin="lower", aspect="auto", cmap="viridis", vmin=0, vmax=100)
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([str(int(x)) for x in xs])
    ax.set_yticks(range(len(ys)))
    ax.set_yticklabels([lab_size(y) for y in ys])
    ax.set_xlabel("subclone size (cells carrying the event)")
    ax.set_ylabel("CNA size")
    ax.set_title(title, fontsize=9)
    for i in range(len(ys)):
        for j in range(len(xs)):
            if cnt[i, j] == 0:
                ax.add_patch(__import__("matplotlib").patches.Rectangle((j - .5, i - .5), 1, 1, fc="#dddddd", ec="white"))
                continue
            v = z[i, j]
            t = "n/a" if v != v else "{:.0f}%".format(v)
            ax.text(j, i, "{}\n(n={})".format(t, cnt[i, j]), ha="center", va="center", fontsize=7,
                    color="white" if (v != v or v < 55) else "black")
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eval-dir", help="an eval folder with section_confusion.tsv")
    ap.add_argument("--compare", nargs="+", help="several eval folders: one heatmap each, side by side")
    ap.add_argument("--names", nargs="+", help="titles for --compare")
    ap.add_argument("--design", required=True, help="the tree's .design.tsv")
    ap.add_argument("--level", default="allele_specific", choices=["allele_specific", "actual", "ideal_phasing"])
    ap.add_argument("--metric", default="bins", choices=["bins", "bp", "cells"])
    ap.add_argument("--seg-min-frac", type=float, default=0.75)
    ap.add_argument("--title", default=None)
    ap.add_argument("-o", "--out", default=None, help="output prefix (default: <eval-dir>/detection)")
    a = ap.parse_args()
    dirs = a.compare or ([a.eval_dir] if a.eval_dir else None)
    if not dirs:
        ap.error("give --eval-dir or --compare")
    names = a.names or [os.path.basename(os.path.normpath(d)) for d in dirs]
    design = pd.read_csv(a.design, sep="\t")
    tag = "min{:g}".format(100 * a.seg_min_frac)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    evs = []
    for d, nm in zip(dirs, names):
        e = per_event(os.path.join(d, "section_confusion.tsv"), design, a.level, tag, a.metric)
        e["run"] = nm
        evs.append(e)
    allev = pd.concat(evs, ignore_index=True)
    prefix = a.out or os.path.join(dirs[0], "detection")
    allev.to_csv(prefix + ".events.tsv", sep="\t", index=False, float_format="%.4g")
    n = len(evs)
    fig, axes = plt.subplots(1, n, figsize=(4.6 * n + 1.2, 4.6), squeeze=False)
    for ax, e, nm in zip(axes[0], evs, names):
        im = draw(ax, e, nm)
    fig.colorbar(im, ax=axes[0].tolist(), label="% of carrier {} with the event called".format(a.metric), shrink=0.85)
    fig.suptitle(a.title or "Detection ({}, bins >= {:g}% inside the event, {})".format(a.level, 100 * a.seg_min_frac, a.metric), fontsize=10)
    fig.savefig(prefix + ".png", dpi=160, bbox_inches="tight")
    print("wrote {0}.png and {0}.events.tsv".format(prefix))
    print(allev.groupby("run").pct.agg(["mean", "count"]).to_string())


if __name__ == "__main__":
    main()
