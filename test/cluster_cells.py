#!/usr/bin/env python3
"""cluster_cells.py -- cluster cells on RDR and BAF of selected bins, compare with the true clones.

    python3 test/cluster_cells.py --bins PREFIX_bins.tsv --run RUN --barcodes BARCODES --truth TRUTH --k 3 \\
        [--rdr-pct 70 --baf-pct 70 --combine or] [--rdr-weight 0.5] [--linkage average] [--clones MAPPING]

A proof of concept for replacing CHISEL's Cloner (Hamming distance on called states, all bins):
1. bins: chosen as test/select_bins.py does (RDR or BAF dispersion at or above a percentile of all
   bins; no labels used), from the PREFIX_bins.tsv of test/bin_variability.py;
2. features: every cell's RDR and BAF (B / (A + B), 0.5 without SNP reads) on those bins, from
   RUN/combo/combo.tsv;
3. distance between cells: w x Euclidean distance of RDR vectors + (1 - w) x Euclidean distance of
   BAF vectors (--rdr-weight w, default 0.5), each divided by its mean over all pairs first so the two
   count equally whatever their scales (--no-normalize to skip);
4. hierarchical clustering (--linkage, default average), cut into --k clusters;
5. compared with the true clonal structure: cells with identical copy numbers (both haplotypes, every
   tree segment) in TRUTH (eval/truth_hg38.tsv) form a group. Prints groups x clusters, the adjusted
   Rand index, and the same for CHISEL's clones with --clones RUN/clones/mapping.tsv.
CHISEL barcodes are matched to simulated cell names with BARCODES (<project>/chisel/prep/barcodedcells.info.tsv).
--out PREFIX writes PREFIX_clusters.tsv (cell, barcode, cluster, true group) and PREFIX_mapping.tsv
(CHISEL Cloner's mapping.tsv format: #CELL CLUSTER CLONE, by barcode).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from select_bins import select   # noqa: E402  same rule as test/select_bins.py


def nchrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def norm_cell(name):
    name = str(name)
    return name[:-4] if name.endswith("_sim") else name


def ari(a, b):
    a = pd.factorize(pd.Series(a))[0]
    b = pd.factorize(pd.Series(b))[0]
    t = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(t, (a, b), 1)
    comb = lambda x: x * (x - 1) / 2.0
    sij, sa, sb, n = comb(t).sum(), comb(t.sum(1)).sum(), comb(t.sum(0)).sum(), comb(len(a))
    exp = sa * sb / n if n else 0.0
    mx = (sa + sb) / 2.0
    return 1.0 if mx == exp else (sij - exp) / (mx - exp)


def true_groups(truth):
    t = pd.read_csv(truth, sep="\t", dtype={"cell": str})
    t["cell"] = t.cell.map(norm_cell)
    prof = t.groupby(["cell", "hap", "seg"]).cn.max().unstack(["hap", "seg"])
    key = prof.astype(int).astype(str).agg("/".join, axis=1)
    order = key.value_counts().index
    names = {k: f"G{i}" for i, k in enumerate(order)}
    desc = {}
    for k in order:
        row = prof[key == k].iloc[0]
        alt = []
        for (h, s), v in row.items():
            if v != 1:
                lo, hi = s.split(":")[1].split("-")
                alt.append((h, int(lo), int(hi), int(v)))
        alt.sort(key=lambda x: (x[0], x[1]))
        merged = []                                      # adjacent segments, same haplotype and CN
        for h, lo, hi, v in alt:
            if merged and merged[-1][0] == h and merged[-1][3] == v and lo - merged[-1][2] <= 1:
                merged[-1][2] = hi
            else:
                merged.append([h, lo, hi, v])
        merged.sort(key=lambda x: x[1])
        desc[names[k]] = ", ".join(f"{h} {lo / 1e6:.1f}-{hi / 1e6:.1f} Mb CN{v}" for h, lo, hi, v in merged) or "no CNA"
    return key.map(names), desc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bins", required=True, help="PREFIX_bins.tsv from test/bin_variability.py")
    ap.add_argument("--run", help="CHISEL run folder (reads combo/combo.tsv)")
    ap.add_argument("--combo", help="instead of --run: a combo.tsv")
    ap.add_argument("--barcodes", required=True, help="<project>/chisel/prep/barcodedcells.info.tsv")
    ap.add_argument("--truth", required=True, help="eval/truth_hg38.tsv")
    ap.add_argument("--k", type=int, required=True, help="number of cell clusters")
    ap.add_argument("--rdr-pct", type=float, default=70, help="bin selection: RDR dispersion percentile (100 = off; default 70)")
    ap.add_argument("--baf-pct", type=float, default=70, help="bin selection: BAF dispersion percentile (100 = off; default 70)")
    ap.add_argument("--combine", choices=["or", "and"], default="or")
    ap.add_argument("--metric", choices=["dispersion", "variance"], default="dispersion")
    ap.add_argument("--rdr-weight", type=float, default=0.5, help="weight of the RDR distance; BAF gets 1 - w (default 0.5)")
    ap.add_argument("--no-normalize", action="store_true", help="do not scale the RDR and BAF distances to mean 1 before mixing")
    ap.add_argument("--linkage", default="average", choices=["average", "complete", "weighted", "single", "ward"],
                    help="hierarchical linkage (default average; CHISEL's Cloner uses weighted)")
    ap.add_argument("--clones", help="CHISEL clones/mapping.tsv to compare as well")
    ap.add_argument("--out", help="output prefix for the cluster assignments")
    a = ap.parse_args()
    if not (a.run or a.combo):
        ap.error("give --run or --combo")
    if not 0 <= a.rdr_weight <= 1:
        ap.error("--rdr-weight must be in [0, 1]")

    # 1. bins
    b = pd.read_csv(a.bins, sep="\t", dtype={"chrom": str})
    b["chrom"] = b.chrom.map(nchrom)
    sel = select(b, f"RDR {a.metric}", f"BAF {a.metric}", a.rdr_pct, a.baf_pct, a.combine)
    chosen = b[sel]
    print(f"bins selected: {len(chosen)} of {len(b)} ({int((chosen.label == 'CNA').sum())} CNA, "
          f"{int((chosen.label == 'no CNA').sum())} no CNA; labels only reported, not used)")
    if len(chosen) < 2:
        sys.exit("ERROR: fewer than 2 bins selected")

    # 2. features
    combo = a.combo or os.path.join(a.run, "combo", "combo.tsv")
    cols = ["chrom", "start", "end", "cell", "norm", "count", "rdr", "a", "b", "baf"]
    d = pd.read_csv(combo, sep="\t", header=None, comment="#", names=cols, usecols=[0, 1, 2, 3, 6, 9],
                    dtype={"chrom": str, "cell": str})
    d["chrom"] = d.chrom.map(nchrom)
    d = d.merge(chosen[["chrom", "start", "end"]], on=["chrom", "start", "end"], how="inner")
    rdr = d.pivot_table(index="cell", columns=["chrom", "start"], values="rdr")
    baf = d.pivot_table(index="cell", columns=["chrom", "start"], values="baf").reindex(index=rdr.index)
    if rdr.isna().any().any() or baf.isna().any().any():
        print(f"[cluster_cells] {int(rdr.isna().sum().sum() + baf.isna().sum().sum())} missing cell-bin values set to the bin mean")
        rdr, baf = rdr.fillna(rdr.mean()), baf.fillna(baf.mean())
    print(f"cells: {rdr.shape[0]}, features: {rdr.shape[1]} RDR + {baf.shape[1]} BAF")

    # 3. distance, 4. clustering
    dr, db = pdist(rdr.values), pdist(baf.values)
    if not a.no_normalize:
        dr, db = dr / max(dr.mean(), 1e-12), db / max(db.mean(), 1e-12)
    dist = a.rdr_weight * dr + (1 - a.rdr_weight) * db
    if a.linkage == "ward":
        print("[cluster_cells] note: ward assumes Euclidean distances; the mixed distance is not exactly that")
    Z = linkage(dist, method=a.linkage)
    clus = pd.Series(fcluster(Z, t=a.k, criterion="maxclust"), index=rdr.index, name="cluster")

    # 5. comparison
    bc = pd.read_csv(a.barcodes, sep="\t")
    bc.columns = [c.lstrip("#") for c in bc.columns]
    to_cell = dict(zip(bc.BARCODE.astype(str), bc.CELL.map(norm_cell)))
    unknown = [c for c in clus.index if c not in to_cell]
    if unknown:
        sys.exit(f"ERROR: {len(unknown)} barcodes not in {a.barcodes}, e.g. {unknown[:3]}")
    groups, desc = true_groups(a.truth)
    res = pd.DataFrame({"barcode": clus.index, "cell": [to_cell[c] for c in clus.index], "cluster": clus.values})
    res["group"] = res.cell.map(groups)
    missing = res.group.isna().sum()
    if missing:
        print(f"[cluster_cells] {missing} cells not in the truth, left out of the comparison")
    r = res.dropna(subset=["group"])
    print(f"\n== true groups x clusters (k = {a.k}, {a.linkage} linkage, RDR weight {a.rdr_weight:g}); "
          f"ARI {ari(r.group, r.cluster):.3f}")
    print(pd.crosstab(r.group, r.cluster, margins=True).to_string())
    print("\n".join(f"  {g}: {desc[g]}" for g in sorted(set(r.group), key=lambda x: int(x[1:]))))
    if a.clones:
        m = pd.read_csv(a.clones, sep="\t")
        m.columns = [c.lstrip("#") for c in m.columns]
        m["cell"] = m.CELL.astype(str).map(to_cell)
        m["group"] = m.cell.map(groups)
        m = m.dropna(subset=["group"])
        print(f"\n== for comparison, CHISEL's clones ({a.clones}); ARI {ari(m.group, m.CLONE.astype(str)):.3f} "
              f"(None counted as one cluster)")
        print(pd.crosstab(m.group, m.CLONE, margins=True).to_string())
    if a.out:
        res.to_csv(f"{a.out}_clusters.tsv", sep="\t", index=False)
        with open(f"{a.out}_mapping.tsv", "w") as f:
            f.write("#CELL\tCLUSTER\tCLONE\n")
            for bcode, c in zip(res.barcode, res.cluster):
                f.write(f"{bcode}\t{c}\tClone{c}\n")
        print(f"\n(wrote {a.out}_clusters.tsv and {a.out}_mapping.tsv)")


if __name__ == "__main__":
    main()
