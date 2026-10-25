#!/usr/bin/env python3
"""cluster_ari.py -- how close CHISEL's bin clusters are to the oracle clusters (for diagnosis).

    python3 scripts/cluster_ari.py ORACLE RUN [RUN ...] [--min-frac 0.75] [--events]

ORACLE: clusters.tsv from oracle_clusters.py (e.g. chisel_pbins8_oracle/perfect/oracle/clusters.tsv).
RUN:    a CHISEL run folder (reads calls/calls.tsv) or a calls.tsv; all runs must use the oracle's bins.

Per run: the number of bin clusters, the adjusted Rand index (ARI) of its bin clusters against the
oracle's, over all bins and over the bins that lie >= --min-frac in their oracle cluster (no
boundary bins), and how many oracle CNA clusters were kept apart from the normal cluster: an
oracle cluster counts as merged into normal when most of its bins sit in the predicted cluster that
holds most of the normal bins. --events lists every oracle cluster per run.
ARI: over all pairs of bins, how often both clusterings agree on "same cluster" vs "different",
corrected for chance: 1 = identical partitions, 0 = chance level. Dominated by the large normal
cluster, hence the per-event columns.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd


def ari(a, b):
    a = pd.factorize(pd.Series(a))[0]
    b = pd.factorize(pd.Series(b))[0]
    table = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(table, (a, b), 1)
    comb = lambda x: x * (x - 1) / 2.0
    sij, sa, sb, n = comb(table).sum(), comb(table.sum(1)).sum(), comb(table.sum(0)).sum(), comb(len(a))
    expected = sa * sb / n if n else 0.0
    mx = (sa + sb) / 2.0
    return 1.0 if mx == expected else (sij - expected) / (mx - expected)


def read_calls(path):
    f = os.path.join(path, "calls", "calls.tsv") if os.path.isdir(path) else path
    if not os.path.isfile(f):
        sys.exit(f"ERROR: no calls file {f}")
    head = open(f).readline().lstrip("#").strip().split("\t")
    col = head.index("CLUSTER") if "CLUSTER" in head else 10
    d = pd.read_csv(f, sep="\t", comment="#", header=None, usecols=[0, 1, 2, col], names=["chrom", "start", "end", "pred"],
                    dtype={"chrom": str})
    d = d.drop_duplicates(["chrom", "start", "end"])
    return d, f


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("oracle")
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--min-frac", type=float, default=0.75, help="bins at least this much inside their oracle cluster, for the second ARI (default 0.75)")
    ap.add_argument("--events", action="store_true", help="list every oracle cluster for every run")
    a = ap.parse_args()

    o = pd.read_csv(a.oracle, sep="\t", dtype={"#CHR": str}).rename(
        columns={"#CHR": "chrom", "START": "start", "END": "end", "CLUSTER": "true"})
    if "NORMAL" in o.columns and (o.NORMAL == 1).any():
        normal_id = o[o.NORMAL == 1].true.iloc[0]          # the oracle cluster without CNAs
    else:
        normal_id = o.groupby("true").size().idxmax()      # older oracle files: the largest cluster
        print(f"[ari] no NORMAL column: taking the largest oracle cluster ({normal_id}) as normal", file=sys.stderr)
    summ = []
    for run in a.runs:
        c, f = read_calls(run)
        m = o.merge(c, on=["chrom", "start", "end"], how="left")
        if m.pred.isna().any():
            sys.exit(f"ERROR: {int(m.pred.isna().sum())} oracle bins are not in {f} (different bins?)")
        pure = m.FRAC >= a.min_frac
        normal_pred = m[m.true == normal_id].pred.value_counts().idxmax()
        ev = []
        for t, g in m.groupby("true"):
            vc = g.pred.value_counts()
            ev.append({"oracle": t, "bins": len(g), "mb": (g.end - g.start).sum() / 1e6, "main_pred": vc.index[0],
                       "share": vc.iloc[0] / len(g), "n_pred": len(vc),
                       "merged_into_normal": t != normal_id and vc.index[0] == normal_pred})
        ev = pd.DataFrame(ev)
        cna = ev[ev.oracle != normal_id]
        summ.append({"run": run, "clusters": m.pred.nunique(), "oracle_clusters": m.true.nunique(),
                     "ARI": ari(m.true, m.pred), f"ARI_frac>={a.min_frac:g}": ari(m.true[pure], m.pred[pure]),
                     "CNA_clusters_kept_apart": f"{int((~cna.merged_into_normal).sum())}/{len(cna)}",
                     "CNA_Mb_merged": round(cna[cna.merged_into_normal].mb.sum(), 2)})
        if a.events:
            print(f"\n== {run}  (oracle normal = cluster {normal_id}, its main predicted cluster = {normal_pred})")
            print(ev.round(3).to_string(index=False))
    s = pd.DataFrame(summ)
    with pd.option_context("display.width", 250, "display.max_columns", 20, "display.max_colwidth", 70):
        print("\n" + s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
