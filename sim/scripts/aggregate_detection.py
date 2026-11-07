#!/usr/bin/env python3
"""Pool the detection of several runs (trees / seeds) into one subclone size x CNA size heatmap.

Takes the detection.events.tsv files that eval_detection.py writes (one row per planted event: size_bp,
clade_cells, pct, n, type, ...), pools events with the same (CNA size, subclone size) over all files and draws
the matrix of eval_detection.py with the pooled number of events. Several groups (e.g. oracle / hmm40 / plain)
are drawn side by side.

  aggregate_detection.py -g oracle eval/*/pbins8_oracle_nc/detection.events.tsv ... -o pooled
  aggregate_detection.py --group oracle A/detection.events.tsv B/detection.events.tsv \\
                         --group hmm40  C/detection.events.tsv D/detection.events.tsv  -o pooled
  --type AMP|DEL   only events of that type          --weighted  the table's mean weights events by their n (cell x bin pairs); the heatmap is a plain mean
Also writes <out>.tsv (per box: mean, SD, events, files).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_detection import draw  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-g", "--group", nargs="+", action="append", metavar=("NAME", "FILE"), required=True,
                    help="NAME followed by its detection.events.tsv files (repeat for several groups)")
    ap.add_argument("--type", choices=["AMP", "DEL"], default=None)
    ap.add_argument("--weighted", action="store_true")
    ap.add_argument("--title", default=None)
    ap.add_argument("-o", "--out", default="detection_pooled")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups, rows = [], []
    for g in a.group:
        name, files = g[0], g[1:]
        if not files:
            ap.error("group {} has no files".format(name))
        parts = []
        for f in files:
            d = pd.read_csv(f, sep="\t")
            d["file"] = f
            parts.append(d)
        ev = pd.concat(parts, ignore_index=True)
        if a.type:
            ev = ev[ev.type == a.type]
        groups.append((name, ev))
        for (sz, cc), x in ev.groupby(["size_bp", "clade_cells"]):
            ok = x.dropna(subset=["pct"])
            w = np.average(ok.pct, weights=ok.n) if a.weighted and ok.n.sum() > 0 else (ok.pct.mean() if len(ok) else np.nan)
            rows.append(dict(group=name, size_bp=sz, clade_cells=cc, mean_pct=w, sd=ok.pct.std() if len(ok) > 1 else np.nan,
                             events=len(x), events_scored=len(ok), files=x.file.nunique()))
    pd.DataFrame(rows).to_csv(a.out + ".tsv", sep="\t", index=False, float_format="%.4g")
    fig, axes = plt.subplots(1, len(groups), figsize=(4.8 * len(groups) + 1.2, 4.8), squeeze=False)
    for ax, (name, ev) in zip(axes[0], groups):
        im = draw(ax, ev, "{} ({} events, {} files)".format(name, len(ev), ev.file.nunique()))
    fig.colorbar(im, ax=axes[0].tolist(), label="% of carrier bins with the event called", shrink=0.85)
    fig.suptitle(a.title or "Detection pooled over runs" + (" ({} only)".format(a.type) if a.type else ""), fontsize=10)
    fig.savefig(a.out + ".png", dpi=160, bbox_inches="tight")
    print("wrote {0}.png and {0}.tsv".format(a.out))


if __name__ == "__main__":
    main()
