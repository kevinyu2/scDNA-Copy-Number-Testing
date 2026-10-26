#!/usr/bin/env python3
"""select_bins.py -- which bins would a variability threshold keep? (diagnosis)

    python3 test/select_bins.py PREFIX_bins.tsv --rdr-pct 80 --baf-pct 90 [--metric dispersion] [--combine or] [--sweep]

Reads the per-bin table written by test/bin_variability.py (label CNA / no CNA, RDR and BAF
variance and dispersion). A bin is selected when its RDR metric is at or above the --rdr-pct
percentile of all bins OR (--combine and: AND) its BAF metric is at or above the --baf-pct
percentile. Percentiles are over all bins, without the labels, as a real method would have to
pick them. --rdr-pct 100 or --baf-pct 100 switches that metric off. Prints, for the selected and
the not-selected bins, how many are CNA and no CNA (bins and Mb), and the share of CNA bins kept
(recall) and of selected bins that are CNA (precision). --sweep adds a table over percentiles.
"""
import argparse
import sys

import numpy as np
import pandas as pd


def select(df, rcol, bcol, rpct, bpct, combine):
    def over(col, pct):
        if pct >= 100:
            return pd.Series(False, index=df.index)
        thr = np.nanpercentile(df[col], pct)
        return (df[col] >= thr).fillna(False)
    r, b = over(rcol, rpct), over(bcol, bpct)
    if combine == "and":
        r = r if rpct < 100 else pd.Series(True, index=df.index)
        b = b if bpct < 100 else pd.Series(True, index=df.index)
        return r & b
    return r | b


def table(df, sel):
    mb = (df.end - df.start) / 1e6
    rows = []
    for name, s in (("selected", sel), ("not selected", ~sel)):
        rows.append({"": name, "CNA bins": int((s & (df.label == "CNA")).sum()),
                     "no-CNA bins": int((s & (df.label == "no CNA")).sum()),
                     "CNA Mb": round(mb[s & (df.label == "CNA")].sum(), 1),
                     "no-CNA Mb": round(mb[s & (df.label == "no CNA")].sum(), 1)})
    t = pd.DataFrame(rows)
    cna = df.label == "CNA"
    recall = (sel & cna).sum() / max(cna.sum(), 1)
    precision = (sel & cna).sum() / max(sel.sum(), 1)
    return t, recall, precision


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("bins", help="PREFIX_bins.tsv from test/bin_variability.py")
    ap.add_argument("--rdr-pct", type=float, default=80, help="RDR metric percentile threshold, 0-100 (100 = off; default 80)")
    ap.add_argument("--baf-pct", type=float, default=100, help="BAF metric percentile threshold, 0-100 (100 = off; default 100)")
    ap.add_argument("--metric", choices=["dispersion", "variance"], default="dispersion", help="default: dispersion")
    ap.add_argument("--combine", choices=["or", "and"], default="or", help="select bins over either threshold (or) or both (and); default or")
    ap.add_argument("--sweep", action="store_true", help="also print recall/precision over a grid of percentiles")
    a = ap.parse_args()

    df = pd.read_csv(a.bins, sep="\t")
    rcol, bcol = f"RDR {a.metric}", f"BAF {a.metric}"
    for c in (rcol, bcol, "label", "start", "end"):
        if c not in df.columns:
            sys.exit(f"ERROR: {a.bins} has no column '{c}' (run bin_variability.py with --run so it has the data metrics)")
    n_cna, n_all = int((df.label == "CNA").sum()), len(df)
    print(f"{n_all} bins: {n_cna} CNA ({100 * n_cna / n_all:.1f}%), {n_all - n_cna} no CNA")

    sel = select(df, rcol, bcol, a.rdr_pct, a.baf_pct, a.combine)
    t, rec, prec = table(df, sel)
    rule = []
    if a.rdr_pct < 100:
        rule.append(f"{rcol} >= {a.rdr_pct:g}th pct ({np.nanpercentile(df[rcol], a.rdr_pct):.3g})")
    if a.baf_pct < 100:
        rule.append(f"{bcol} >= {a.baf_pct:g}th pct ({np.nanpercentile(df[bcol], a.baf_pct):.3g})")
    print(f"\nselected if {f' {a.combine.upper()} '.join(rule) if rule else 'nothing (both thresholds off)'}")
    print(t.to_string(index=False))
    print(f"recall (CNA bins kept) {100 * rec:.1f}%   precision (selected bins that are CNA) {100 * prec:.1f}%   "
          f"selected {int(sel.sum())} of {n_all} bins")

    if a.sweep:
        grid = [50, 60, 70, 75, 80, 85, 90, 95, 100]
        rows = []
        for rp in grid:
            for bp in grid:
                s = select(df, rcol, bcol, rp, bp, a.combine)
                if s.sum() == 0:
                    continue
                _, r, p = table(df, s)
                rows.append({"rdr_pct": rp, "baf_pct": bp, "selected bins": int(s.sum()), "recall %": round(100 * r, 1),
                             "precision %": round(100 * p, 1)})
        sw = pd.DataFrame(rows)
        print(f"\n== sweep ({a.combine}; 100 = metric off)")
        with pd.option_context("display.max_rows", 200, "display.width", 150):
            print(sw.to_string(index=False))


if __name__ == "__main__":
    main()
