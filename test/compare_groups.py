#!/usr/bin/env python3
"""compare_groups.py -- run the fork's divisive cell grouping on a run and score it against the true clones.

    python3 test/compare_groups.py --run RUN --barcodes BARCODES --truth TRUTH --cc CELLCLUSTERER \\
        [--permq 90 --floormargin 1.2 --minsize 14 --seed 12 ...] [--groups FILE ...]

--run RUN: a CHISEL run folder (reads RUN/combo/combo.tsv); --combo FILE instead to name the file.
--cc: the fork's CellClusterer.py (chisel_hmm src/chisel/CellClusterer.py; it runs under Python 3 too).
The grouping is the one chisel_calling runs before calling (calls/cellgroups.tsv), with the fork's
defaults unless options are given; its split tree is printed.
--groups FILE (repeatable): also score existing groupings, e.g. RUN/calls/cellgroups.tsv (CELL GROUP)
or RUN/clones/mapping.tsv (#CELL CLUSTER CLONE); --no-run scores only these.
--barcodes: <project>/chisel/prep/barcodedcells.info.tsv; --truth: <project>/eval/truth_hg38.tsv. True
groups are cells with identical copy numbers in every tree segment (as test/cluster_cells.py).

Per grouping: true groups x clusters, ARI, and the two checks that matter here:
  mixed clusters  clusters whose cells are < --purity (default 90%) one true group: clones merged
  lost groups     true groups that are not the majority of any cluster: a clone not found
Oversplitting (one true group in several clusters) is listed but not counted as a failure.
Noisy cells set apart (None) are left out of the checks and counted.
"""
import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_cells import ari, norm_cell, true_groups   # noqa: E402


def load_cc(path):
    spec = importlib.util.spec_from_file_location("CellClusterer", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read_groups(path):
    g = pd.read_csv(path, sep="\t", dtype=str)
    g.columns = [c.lstrip("#") for c in g.columns]
    col = "GROUP" if "GROUP" in g.columns else "CLONE" if "CLONE" in g.columns else g.columns[-1]
    return pd.Series(g[col].values, index=g.CELL.astype(str))


def score(name, groups, to_cell, truth, desc, purity):
    """groups: barcode -> cluster (str; 'None' = set apart). Prints the table and checks; returns True if no clone merged or lost."""
    d = pd.DataFrame({"barcode": groups.index, "cluster": groups.astype(str).values})
    d["cell"] = d.barcode.map(to_cell)
    d["group"] = d.cell.map(truth)
    unknown = int(d.group.isna().sum())
    noisy = d.cluster.isin(["None", "nan", "-1"])
    r = d[~noisy & d.group.notna()]
    print(f"\n== {name}: {r.cluster.nunique()} clusters, {int(noisy.sum())} cells set apart (None)"
          + (f", {unknown} cells not in the truth" if unknown else "") + f"; ARI {ari(r.group, r.cluster):.3f}")
    print(pd.crosstab(r.group, r.cluster, margins=True).to_string())
    mixed, majority = [], {}
    for c, s in r.groupby("cluster").group:
        top = s.value_counts()
        majority[c] = top.index[0]
        if top.iloc[0] / len(s) < purity:
            mixed.append(f"cluster {c} ({len(s)} cells): " + ", ".join(f"{g} {n}" for g, n in top.items()))
    lost = sorted(set(r.group) - set(majority.values()), key=lambda x: int(x[1:]))
    split = {g: sorted(c for c, m in majority.items() if m == g) for g in set(majority.values())}
    split = {g: c for g, c in split.items() if len(c) > 1}
    print("  mixed clusters (clones merged): " + ("; ".join(mixed) if mixed else "none"))
    print("  lost true groups (no cluster of their own): " + (", ".join(f"{g} ({desc[g]})" for g in lost) if lost else "none"))
    print("  oversplit true groups (fine): " + (", ".join(f"{g} -> clusters {', '.join(c)}" for g, c in sorted(split.items())) if split else "none"))
    ok = not mixed and not lost
    print(f"  => {'OK' if ok else 'PROBLEM'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", help="CHISEL run folder (reads combo/combo.tsv)")
    ap.add_argument("--combo", help="instead of --run: a combo.tsv (or calls.tsv)")
    ap.add_argument("--cc", help="the fork's CellClusterer.py")
    ap.add_argument("--barcodes", required=True, help="<project>/chisel/prep/barcodedcells.info.tsv")
    ap.add_argument("--truth", required=True, help="<project>/eval/truth_hg38.tsv")
    ap.add_argument("--groups", action="append", default=[], help="existing grouping(s) to score as well (CELL GROUP or Cloner mapping.tsv)")
    ap.add_argument("--no-run", action="store_true", help="only score --groups, do not run the grouping")
    ap.add_argument("--purity", type=float, default=0.9, help="a cluster with less of one true group is 'mixed' (default 0.9)")
    ap.add_argument("--seed", type=int, default=12)
    ap.add_argument("--minsize", type=int, default=14)
    ap.add_argument("--rdrweight", type=float, default=0.5)
    ap.add_argument("--noisescale", action="store_true")
    ap.add_argument("--perm", type=int, default=20)
    ap.add_argument("--permq", type=float, default=None, help="default: the fork's (CellClusterer.PERM_Q)")
    ap.add_argument("--floormargin", type=float, default=None, help="default: the fork's (CellClusterer.FLOOR_MARGIN)")
    ap.add_argument("--noreassign", action="store_true")
    ap.add_argument("--out", help="write the new grouping here (CELL GROUP, as calls/cellgroups.tsv)")
    a = ap.parse_args()

    bc = pd.read_csv(a.barcodes, sep="\t")
    bc.columns = [c.lstrip("#") for c in bc.columns]
    to_cell = dict(zip(bc.BARCODE.astype(str), bc.CELL.map(norm_cell)))
    truth, desc = true_groups(a.truth)
    print("true groups: " + "; ".join(f"{g} {int((truth == g).sum())} cells: {desc[g]}"
                                      for g in sorted(set(truth), key=lambda x: int(x[1:]))))
    results = []
    if not a.no_run:
        if not a.cc or not (a.run or a.combo):
            ap.error("give --cc and --run (or --combo), or --no-run with --groups")
        cc = load_cc(a.cc)
        combo = a.combo or os.path.join(a.run, "combo", "combo.tsv")
        q = cc.PERM_Q if a.permq is None else a.permq
        m = cc.FLOOR_MARGIN if a.floormargin is None else a.floormargin
        groups, log = cc.cell_groups(combo, seed=a.seed, min_size=a.minsize, rdr_weight=a.rdrweight,
                                     noise_scale=a.noisescale, n_perm=a.perm, q=q, margin=m,
                                     reassign=not a.noreassign)
        print(f"\n== split tree ({combo}; permq {q:g}, floormargin {m:g})")
        print("\n".join(log))
        g = pd.Series({str(k): str(v) for k, v in groups.items()})
        if a.out:
            pd.DataFrame({"#CELL": g.index, "GROUP": g.values}).to_csv(a.out, sep="\t", index=False)
        results.append(score(f"new grouping (permq {q:g}, floormargin {m:g})", g, to_cell, truth, desc, a.purity))
    for f in a.groups:
        results.append(score(f, read_groups(f), to_cell, truth, desc, a.purity))
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
