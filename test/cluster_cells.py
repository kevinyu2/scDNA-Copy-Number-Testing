#!/usr/bin/env python3
"""cluster_cells.py -- cluster cells on RDR and BAF of selected bins, compare with the true clones.

    python3 test/cluster_cells.py --bins PREFIX_bins.tsv --run RUN --barcodes BARCODES --truth TRUTH --k 3 \\
        [--rdr-pct 70 --baf-pct 70 --combine or | --all-bins] [--rdr-weight 0.5] [--noise-scale]
        [--pcs 10] [--method kmeans] [--plot PCS.pdf] [--clones MAPPING]

A proof of concept for replacing CHISEL's Cloner (Hamming distance on called states, all bins):
1. bins: chosen as test/select_bins.py does (RDR or BAF dispersion at or above a percentile of all
   bins; no labels used), from the PREFIX_bins.tsv of test/bin_variability.py; --all-bins: every bin;
2. features: every cell's RDR and BAF (B / (A + B)) on those bins, from RUN/combo/combo.tsv, each bin
   centred on its mean over cells (cells without SNP reads in a bin get the bin's mean BAF);
   --noise-scale: each value divided by its expected counting noise (RDR: Poisson on the cell's and
   the normal's reads; BAF: binomial on the cell's SNP reads), so bins are in units of noise;
3. the RDR block and the BAF block are each scaled to total variance 1 and weighted (--rdr-weight w,
   BAF 1 - w), so squared distances between cells mix the two as w : 1 - w; within a block, bins
   keep their own variance (bins where cells differ count more than noise-only bins);
4. --pcs N > 0: PCA (SVD) to the top N components (default 10; 0 = no PCA);
5. --method kmeans (default; k-means++ starts, best of 20) or hierarchical (ward, average, complete,
   weighted) on Euclidean distances, cut into --k clusters. Average/weighted linkage on noisy cells
   tends to split off single outlier cells first; k-means and ward split large groups first;
6. compared with the true clonal structure: cells with identical copy numbers (both haplotypes, every
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


def kmeans(X, K, rng, n_init=20, iters=300):
    """Lloyd's k-means with k-means++ starts, best of n_init; returns labels 0..K-1."""
    n = X.shape[0]
    best, best_sse = None, np.inf
    for _ in range(n_init):
        C = [X[rng.randint(n)]]
        d2 = ((X - C[0]) ** 2).sum(1)
        for _ in range(1, K):
            C.append(X[rng.choice(n, p=d2 / d2.sum()) if d2.sum() > 0 else rng.randint(n)])
            d2 = np.minimum(d2, ((X - C[-1]) ** 2).sum(1))
        C = np.array(C)
        for _ in range(iters):
            D = ((X[:, None, :] - C[None, :, :]) ** 2).sum(2)
            lab = D.argmin(1)
            newC = np.array([X[lab == k].mean(0) if (lab == k).any() else X[D.min(1).argmax()] for k in range(K)])
            if np.allclose(newC, C):
                break
            C = newC
        sse = ((X - C[lab]) ** 2).sum()
        if sse < best_sse:
            best, best_sse = lab, sse
    return best


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
    ap.add_argument("--all-bins", action="store_true", help="use every bin (ignore the percentile selection)")
    ap.add_argument("--rdr-weight", type=float, default=0.5, help="weight of RDR in squared distances; BAF gets 1 - w (default 0.5)")
    ap.add_argument("--noise-scale", action="store_true", help="divide values by their expected counting noise")
    ap.add_argument("--pcs", type=int, default=10, help="principal components to keep (0 = no PCA; default 10)")
    ap.add_argument("--method", default="kmeans", choices=["kmeans", "ward", "average", "complete", "weighted"],
                    help="kmeans (default) or hierarchical linkage (CHISEL's Cloner uses weighted)")
    ap.add_argument("--seed", type=int, default=0, help="k-means random seed (default 0)")
    ap.add_argument("--plot", help="PC1 vs PC2 scatter (coloured by true group, marker by cluster) to this PDF/PNG")
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
    sel = pd.Series(True, index=b.index) if a.all_bins else select(b, f"RDR {a.metric}", f"BAF {a.metric}", a.rdr_pct, a.baf_pct, a.combine)
    chosen = b[sel]
    print(f"bins {'(all)' if a.all_bins else 'selected'}: {len(chosen)} of {len(b)} ({int((chosen.label == 'CNA').sum())} CNA, "
          f"{int((chosen.label == 'no CNA').sum())} no CNA; labels only reported, not used)")
    if len(chosen) < 2:
        sys.exit("ERROR: fewer than 2 bins selected")

    # 2. features
    combo = a.combo or os.path.join(a.run, "combo", "combo.tsv")
    cols = ["chrom", "start", "end", "cell", "norm", "count", "rdr", "a", "b", "baf"]
    d = pd.read_csv(combo, sep="\t", header=None, comment="#", names=cols, dtype={"chrom": str, "cell": str})
    d["chrom"] = d.chrom.map(nchrom)
    d = d.merge(chosen[["chrom", "start", "end"]], on=["chrom", "start", "end"], how="inner")
    d["n"] = d.a + d.b
    d.loc[d.n == 0, "baf"] = np.nan                                   # no SNP reads: no BAF
    piv = (lambda v: d.pivot_table(index="cell", columns=["chrom", "start"], values=v))
    rdr = piv("rdr")
    baf = piv("baf").reindex(index=rdr.index, columns=rdr.columns)
    R = rdr - rdr.mean()
    B = baf - baf.mean()
    if a.noise_scale:
        cnt, nrm, n = piv("count").reindex_like(rdr), piv("norm").reindex_like(rdr), piv("n").reindex_like(rdr)
        rsd = np.sqrt(rdr ** 2 * (1.0 / cnt.where(cnt > 0) + 1.0 / nrm.clip(lower=1)))
        p = (piv("b").reindex_like(rdr).sum() / n.sum().replace(0, np.nan)).clip(0.01, 0.99)
        bsd = np.sqrt(p * (1 - p) / n.where(n > 0))
        R, B = R / rsd, B / bsd
    R, B = R.fillna(0.0).values, B.fillna(0.0).values                 # missing: the bin's mean
    print(f"cells: {R.shape[0]}, features: {R.shape[1]} RDR + {B.shape[1]} BAF"
          + (" (in units of counting noise)" if a.noise_scale else ""))

    # 3. mix the blocks, 4. PCA, 5. clustering
    tr, tb = max((R ** 2).sum() / R.shape[0], 1e-12), max((B ** 2).sum() / B.shape[0], 1e-12)
    X = np.hstack([R * np.sqrt(a.rdr_weight / tr), B * np.sqrt((1 - a.rdr_weight) / tb)])
    pcs = None
    if a.pcs > 0:
        U, S, Vt = np.linalg.svd(X - X.mean(axis=0), full_matrices=False)
        npc = min(a.pcs, len(S))
        explained = S ** 2 / (S ** 2).sum()
        X = pcs = U[:, :npc] * S[:npc]
        print(f"PCA: top {npc} components explain {100 * explained[:npc].sum():.1f}% of the variance "
              f"(PC1-5: {', '.join(f'{100 * e:.1f}%' for e in explained[:5])})")
    if a.method == "kmeans":
        lab = kmeans(X, a.k, np.random.RandomState(a.seed))
    else:
        lab = fcluster(linkage(pdist(X), method=a.method), t=a.k, criterion="maxclust") - 1
    clus = pd.Series(lab + 1, index=rdr.index, name="cluster")

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
    how = f"{a.method}{'' if a.method == 'kmeans' else ' linkage'}, {f'{a.pcs} PCs' if a.pcs > 0 else 'no PCA'}"
    print(f"\n== true groups x clusters (k = {a.k}, {how}, RDR weight {a.rdr_weight:g}); "
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
    if a.plot:
        if pcs is None or pcs.shape[1] < 2:
            print("[cluster_cells] --plot needs --pcs >= 2")
        else:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            grp = res.set_index("barcode").reindex(rdr.index).group.fillna("?")
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
            marks = "o^sDvP*Xh<>p"
            for ax, (i, j) in zip(axes, ((0, 1), (0, 2) if pcs.shape[1] > 2 else (0, 1))):
                for gi, g in enumerate(sorted(grp.unique(), key=lambda x: (x == "?", x))):
                    labelled = False
                    for ci, c in enumerate(sorted(clus.unique())):
                        m = ((grp == g) & (clus == c)).values
                        if m.any():
                            ax.scatter(pcs[m, i], pcs[m, j], s=9, color=f"C{gi % 10}", marker=marks[ci % len(marks)],
                                       alpha=0.6, label=None if labelled else f"{g} ({int((grp == g).sum())} cells)")
                            labelled = True
                ax.set_xlabel(f"PC{i + 1}")
                ax.set_ylabel(f"PC{j + 1}")
            axes[0].legend(fontsize=7, title="true group (colour); marker = cluster", title_fontsize=7)
            fig.suptitle(f"cells in PC space: {how}, k = {a.k}", fontsize=9)
            fig.tight_layout()
            fig.savefig(a.plot)
            print(f"(plot: {a.plot})")
    if a.out:
        res.to_csv(f"{a.out}_clusters.tsv", sep="\t", index=False)
        with open(f"{a.out}_mapping.tsv", "w") as f:
            f.write("#CELL\tCLUSTER\tCLONE\n")
            for bcode, c in zip(res.barcode, res.cluster):
                f.write(f"{bcode}\t{c}\tClone{c}\n")
        print(f"\n(wrote {a.out}_clusters.tsv and {a.out}_mapping.tsv)")


if __name__ == "__main__":
    main()
