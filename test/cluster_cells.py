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
   --k MIN-MAX (e.g. 2-10): every K in the range is clustered and scored by the mean silhouette of
   the cells in the clustering space (the PCs, or the mixed features without PCA; Euclidean); the K
   with the highest silhouette among those whose clusters all have at least --min-size cells is kept
   (default 14, as CHISEL's Cloner's clone size): without that, the silhouette tends to grow with K by
   isolating noisy cells (a cell alone in its cluster scores 0, more than it would inside a group).
   The table also shows each K's ARI against the true groups, to judge the choice (--plot-k draws
   both against K);
   --method divisive: top-down instead, K is not given. Starting from all cells, a group is split in
   two (PCA of the group's own cells, 2-means on its signal PCs) when it still has structure: at
   least one PC whose variance beats the noise floor of a permutation test (every feature's values
   shuffled across the group's cells, which keeps each feature's noise but destroys shared
   structure; floor = --perm-q percentile of the shuffled top PC variance over --perm permutations,
   Horn's parallel analysis; a signal PC must beat the floor by --floor-margin, default 1.5x, since
   every group of the tree is tested and PCs just above the floor are usually one noisy feature), the split is bimodal (cells projected on the line between the two
   halves' centres: two Gaussians fit better than one by BIC, and they are apart, Ashman's
   D = |mu1 - mu2| sqrt(2 / (s1^2 + s2^2)) >= --ashman, default 2; BIC alone also prefers two
   Gaussians for one heavy-tailed hump, e.g. from cells of different depth) and both halves have
   >= --min-size cells. When 2-means cuts off fewer than --min-size cells, those cells (outliers) are set aside
   and the split is retried on the rest (up to --peel times); set-aside cells join the nearer half.
   Groups that do not split are the clusters. Then every cell is reassigned to the nearest cluster centre (repeated until
   stable; --no-reassign to skip), which corrects cells sent the wrong way by an early split. The
   split tree is printed;
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
from scipy.spatial.distance import pdist, squareform

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


def silhouette(Dm, labels):
    """Mean silhouette for a precomputed distance matrix (cells alone in their cluster score 0)."""
    labels = np.asarray(labels)
    uniq, lab = np.unique(labels, return_inverse=True)
    n, k = len(lab), len(uniq)
    if not 2 <= k <= n - 1:
        return float("nan")
    onehot = np.zeros((n, k))
    onehot[np.arange(n), lab] = 1.0
    sizes = onehot.sum(0)
    sums = Dm.dot(onehot)
    own = sizes[lab]
    a_ = sums[np.arange(n), lab] / np.maximum(own - 1, 1)
    other = sums / sizes[None, :]
    other[np.arange(n), lab] = np.inf
    b_ = other.min(1)
    s = (b_ - a_) / np.maximum(np.maximum(a_, b_), 1e-300)
    s[own == 1] = 0.0
    return float(s.mean())


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


def signal_pcs(X, rng, n_perm, q, max_pcs=20, margin=1.0):
    """Horn's parallel analysis: PCs of X (cells x features) whose variance beats the q-th percentile
    of the top PC variance after shuffling every feature independently across cells.
    Returns (number of signal PCs, PC scores (cells x max_pcs), real variances, noise floor)."""
    from scipy.sparse.linalg import svds
    Xc = X - X.mean(axis=0)
    n = Xc.shape[0]
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2 / max(n - 1, 1)
    top = []
    for _ in range(n_perm):
        idx = rng.rand(*Xc.shape).argsort(axis=0)                    # an independent shuffle per feature
        Xp = Xc[idx, np.arange(Xc.shape[1])]
        top.append(svds(Xp, k=1, return_singular_vectors=False)[0] ** 2 / max(n - 1, 1))
    floor = float(np.percentile(top, q)) * margin
    m = int((var > floor).sum())
    k = min(max_pcs, len(S))
    return m, U[:, :k] * S[:k], var, floor


def bimodal(y, min_d=2.0):
    """1-D: do two Gaussians fit better than one by BIC, and are they apart (Ashman's D >= min_d)?
    Returns (decision, BIC1 - BIC2, D)."""
    y = np.asarray(y, float)
    n = len(y)
    m, v = y.mean(), max(y.var(), 1e-12)
    ll1 = -0.5 * n * (np.log(2 * np.pi * v) + 1)
    # EM for two components, started from the split at the median
    lab = y > np.median(y)
    w = np.array([1 - lab.mean(), lab.mean()])
    mu = np.array([y[~lab].mean(), y[lab].mean()])
    va = np.array([max(y[~lab].var(), 1e-6 * v), max(y[lab].var(), 1e-6 * v)])
    for _ in range(200):
        dens = w / np.sqrt(2 * np.pi * va) * np.exp(-0.5 * (y[:, None] - mu) ** 2 / va)
        tot = np.maximum(dens.sum(1, keepdims=True), 1e-300)
        r = dens / tot
        nk = np.maximum(r.sum(0), 1e-9)
        w, mu = nk / n, (r * y[:, None]).sum(0) / nk
        va = np.maximum((r * (y[:, None] - mu) ** 2).sum(0) / nk, 1e-6 * v)
    dens = w / np.sqrt(2 * np.pi * va) * np.exp(-0.5 * (y[:, None] - mu) ** 2 / va)
    ll2 = np.log(np.maximum(dens.sum(1), 1e-300)).sum()
    delta = (-2 * ll1 + 2 * np.log(n)) - (-2 * ll2 + 5 * np.log(n))
    d = abs(mu[0] - mu[1]) * np.sqrt(2.0 / (va[0] + va[1]))
    return bool(delta > 0 and d >= min_d), float(delta), float(d)


def divisive(X, a, rng):
    """Top-down splitting with a permutation test per group; returns labels 0..K-1."""
    leaves, queue, log = [], [(np.arange(X.shape[0]), "root", 0)], []
    while queue:
        idx, name, depth = queue.pop(0)
        pad = "  " * depth
        if len(idx) < 2 * a.min_size or len(leaves) + len(queue) + 1 >= a.max_clusters:
            why = f"< {2 * a.min_size} cells" if len(idx) < 2 * a.min_size else f"--max-clusters {a.max_clusters}"
            log.append(f"{pad}{name}: {len(idx)} cells -> cluster ({why})")
            leaves.append(idx)
            continue
        core, peeled, verdict, split = idx, np.array([], dtype=int), None, None
        for attempt in range(a.peel + 1):
            if len(core) < 2 * a.min_size:
                verdict = f"{len(core)} cells left after setting aside {len(peeled)} outliers -> cluster"
                break
            m, scores, var, floor = signal_pcs(X[core], rng, a.perm, a.perm_q, margin=a.floor_margin)
            if m == 0:
                tops = ", ".join(f"{v:.3g}" for v in var[:3])
                verdict = f"PC variances {tops} vs noise floor {floor:.3g}: no signal PC -> cluster"
                break
            lab = kmeans(scores[:, :m], 2, rng)
            sizes = np.bincount(lab, minlength=2)
            if sizes.min() < a.min_size:                      # outliers: set them aside, retry on the rest
                small = int(sizes.argmin())
                peeled = np.concatenate([peeled, core[lab == small]])
                core = core[lab != small]
                continue
            c0, c1 = scores[lab == 0, :m].mean(0), scores[lab == 1, :m].mean(0)
            y = scores[:, :m].dot((c1 - c0) / max(np.linalg.norm(c1 - c0), 1e-12))
            ok, delta, dd = bimodal(y, a.ashman)
            if not ok:
                verdict = (f"{m} signal PC(s) (PC1 {var[0]:.3g} vs floor {floor:.3g}) but the split {sizes[0]}/{sizes[1]} is not "
                           f"bimodal (BIC 1 - 2 Gaussians {delta:.3g}, Ashman's D {dd:.2f}) -> cluster")
                break
            split = (core, lab, m, var, floor, sizes, (delta, dd), scores)
            break
        else:
            verdict = f"still cutting off < {a.min_size} cells after {a.peel} retries -> cluster"
        peel_note = f" ({len(peeled)} outlier cells set aside)" if len(peeled) else ""
        if split is None:
            log.append(f"{pad}{name}: {len(idx)} cells, {verdict}{peel_note}")
            leaves.append(idx)
            continue
        core, lab, m, var, floor, sizes, delta, scores = split
        halves = [core[lab == 0], core[lab == 1]]
        if len(peeled):                                       # outliers join the nearer half (feature space)
            cen = [X[h].mean(0) for h in halves]
            near = np.argmin([((X[peeled] - c) ** 2).sum(1) for c in cen], axis=0)
            halves = [np.concatenate([halves[k], peeled[near == k]]) for k in (0, 1)]
        log.append(f"{pad}{name}: {len(idx)} cells, {m} signal PC(s) (PC1 {var[0]:.3g} vs floor {floor:.3g}), bimodal "
                   f"(BIC 1 - 2 Gaussians {delta[0]:.3g}, Ashman's D {delta[1]:.2f}) -> split {len(halves[0])} / {len(halves[1])}{peel_note}")
        queue.append((halves[0], name + ".0", depth + 1))
        queue.append((halves[1], name + ".1", depth + 1))
    lab = np.empty(X.shape[0], dtype=int)
    for i, idx in enumerate(leaves):
        lab[idx] = i
    print("\n== divisive splits")
    print("\n".join(log))
    if not a.no_reassign and len(leaves) > 1:
        # nearest leaf centre in a global PCA space, repeated until stable (k-means started at the leaves)
        Xc = X - X.mean(axis=0)
        U, S, _ = np.linalg.svd(Xc, full_matrices=False)
        P = U[:, :min(20, len(S))] * S[:min(20, len(S))]
        moved_total = 0
        for _ in range(100):
            C = np.array([P[lab == k].mean(0) for k in range(len(leaves))])
            new = ((P[:, None, :] - C[None, :, :]) ** 2).sum(2).argmin(1)
            moved = int((new != lab).sum())
            if not moved:
                break
            moved_total += moved
            lab = new
            keep = np.unique(lab)                                         # drop emptied clusters
            lab = np.searchsorted(keep, lab)
            leaves = keep
        print(f"reassignment to the nearest cluster centre (top 20 global PCs): {moved_total} moves, {len(np.unique(lab))} clusters")
    return lab


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
    ap.add_argument("--bins", help="PREFIX_bins.tsv from test/bin_variability.py of the same run (needed for the percentile selection; with --all-bins the run's own bins are used)")
    ap.add_argument("--run", help="CHISEL run folder (reads combo/combo.tsv)")
    ap.add_argument("--combo", help="instead of --run: a combo.tsv")
    ap.add_argument("--barcodes", required=True, help="<project>/chisel/prep/barcodedcells.info.tsv")
    ap.add_argument("--truth", required=True, help="eval/truth_hg38.tsv")
    ap.add_argument("--k", help="number of cell clusters, or MIN-MAX to choose by silhouette (e.g. 2-10); not used by --method divisive")
    ap.add_argument("--plot-k", help="with --k MIN-MAX: silhouette and ARI against K to this PDF/PNG")
    ap.add_argument("--min-size", type=int, default=14, help="with --k MIN-MAX: only Ks whose clusters all have at least this many cells (default 14; 1 = no limit)")
    ap.add_argument("--rdr-pct", type=float, default=70, help="bin selection: RDR dispersion percentile (100 = off; default 70)")
    ap.add_argument("--baf-pct", type=float, default=70, help="bin selection: BAF dispersion percentile (100 = off; default 70)")
    ap.add_argument("--combine", choices=["or", "and"], default="or")
    ap.add_argument("--metric", choices=["dispersion", "variance"], default="dispersion")
    ap.add_argument("--all-bins", action="store_true", help="use every bin (ignore the percentile selection)")
    ap.add_argument("--rdr-weight", type=float, default=0.5, help="weight of RDR in squared distances; BAF gets 1 - w (default 0.5)")
    ap.add_argument("--noise-scale", action="store_true", help="divide values by their expected counting noise")
    ap.add_argument("--pcs", type=int, default=10, help="principal components to keep (0 = no PCA; default 10)")
    ap.add_argument("--perm", type=int, default=20, help="divisive: permutations per split test (default 20)")
    ap.add_argument("--perm-q", type=float, default=95, help="divisive: noise floor = this percentile of the permuted top PC variance (default 95)")
    ap.add_argument("--max-clusters", type=int, default=30, help="divisive: stop splitting at this many clusters (default 30)")
    ap.add_argument("--no-reassign", action="store_true", help="divisive: keep the split tree's assignment (no final reassignment)")
    ap.add_argument("--floor-margin", type=float, default=1.5, help="divisive: a signal PC must exceed the permutation noise floor by this factor (default 1.5)")
    ap.add_argument("--ashman", type=float, default=2.0, help="divisive: minimum Ashman's D of the two halves (default 2)")
    ap.add_argument("--peel", type=int, default=5, help="divisive: times to set aside small outlier groups and retry a split (default 5)")
    ap.add_argument("--method", default="kmeans", choices=["kmeans", "ward", "average", "complete", "weighted", "divisive"],
                    help="kmeans (default) or hierarchical linkage (CHISEL's Cloner uses weighted)")
    ap.add_argument("--seed", type=int, default=0, help="k-means random seed (default 0)")
    ap.add_argument("--plot", help="PC1 vs PC2 scatter (coloured by true group, marker by cluster) to this PDF/PNG")
    ap.add_argument("--clones", help="CHISEL clones/mapping.tsv to compare as well")
    ap.add_argument("--out", help="output prefix for the cluster assignments")
    a = ap.parse_args()
    if not (a.run or a.combo):
        ap.error("give --run or --combo")
    if not a.bins and not a.all_bins:
        ap.error("selecting bins needs --bins (PREFIX_bins.tsv of this run); or use --all-bins")
    if a.method == "divisive":
        a.k = a.k or "1"
    elif not a.k:
        ap.error("--k is needed (or --method divisive)")
    try:
        kmin, kmax = (int(x) for x in a.k.split("-")) if "-" in a.k else (int(a.k), int(a.k))
    except ValueError:
        ap.error("--k needs a number or MIN-MAX, e.g. 3 or 2-10")
    if not 1 <= kmin <= kmax:
        ap.error("--k needs 1 <= MIN <= MAX")
    if kmin < kmax and kmin < 2:
        ap.error("choosing K by silhouette needs MIN >= 2 (the silhouette is undefined for one cluster)")
    if not 0 <= a.rdr_weight <= 1:
        ap.error("--rdr-weight must be in [0, 1]")

    # 1. bins (2. features read first: the run's own bins)
    combo = a.combo or os.path.join(a.run, "combo", "combo.tsv")
    cols = ["chrom", "start", "end", "cell", "norm", "count", "rdr", "a", "b", "baf"]
    d = pd.read_csv(combo, sep="\t", header=None, comment="#", names=cols, dtype={"chrom": str, "cell": str})
    d["chrom"] = d.chrom.map(nchrom)
    run_bins = d[["chrom", "start", "end"]].drop_duplicates()
    if a.bins:
        b = pd.read_csv(a.bins, sep="\t", dtype={"chrom": str})
        b["chrom"] = b.chrom.map(nchrom)
        common = b.merge(run_bins, on=["chrom", "start", "end"], how="inner")
        if len(common) < 0.95 * len(run_bins):
            msg = (f"only {len(common)} of the run's {len(run_bins)} bins are in {a.bins}: it was made for other bins "
                   f"(another run or project?)")
            if not a.all_bins:
                sys.exit(f"ERROR: {msg}; rerun test/bin_variability.py on this run")
            print(f"[cluster_cells] WARNING: {msg}; using the run's own bins")
            b = None
    else:
        b = None
    if b is None:                                       # --all-bins without a matching bins file
        b = run_bins.assign(label="?").reset_index(drop=True)
    sel = pd.Series(True, index=b.index) if a.all_bins else select(b, f"RDR {a.metric}", f"BAF {a.metric}", a.rdr_pct, a.baf_pct, a.combine)
    chosen = b[sel]
    if (chosen.label != "?").any():
        print(f"bins {'(all)' if a.all_bins else 'selected'}: {len(chosen)} of {len(b)} ({int((chosen.label == 'CNA').sum())} CNA, "
              f"{int((chosen.label == 'no CNA').sum())} no CNA; labels only reported, not used)")
    else:
        print(f"bins (all of the run's): {len(chosen)}")
    if len(chosen) < 2:
        sys.exit("ERROR: fewer than 2 bins selected")

    # 2. features
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
    X = Xfeat = np.hstack([R * np.sqrt(a.rdr_weight / tr), B * np.sqrt((1 - a.rdr_weight) / tb)])
    pcs = None
    if a.pcs > 0:
        U, S, Vt = np.linalg.svd(X - X.mean(axis=0), full_matrices=False)
        npc = min(a.pcs, len(S))
        explained = S ** 2 / (S ** 2).sum()
        X = pcs = U[:, :npc] * S[:npc]
        print(f"PCA: top {npc} components explain {100 * explained[:npc].sum():.1f}% of the variance "
              f"(PC1-5: {', '.join(f'{100 * e:.1f}%' for e in explained[:5])})")
    def cluster(k):
        if a.method == "kmeans":
            return kmeans(X, k, np.random.RandomState(a.seed))
        return fcluster(Z, t=k, criterion="maxclust") - 1
    Z = None if a.method in ("kmeans", "divisive") else linkage(pdist(X), method=a.method)

    # true groups (for the K table and the comparison)
    bc = pd.read_csv(a.barcodes, sep="\t")
    bc.columns = [c.lstrip("#") for c in bc.columns]
    to_cell = dict(zip(bc.BARCODE.astype(str), bc.CELL.map(norm_cell)))
    unknown = [c for c in rdr.index if c not in to_cell]
    if unknown:
        sys.exit(f"ERROR: {len(unknown)} barcodes not in {a.barcodes}, e.g. {unknown[:3]}")
    groups, desc = true_groups(a.truth)
    cell_group = pd.Series([groups.get(to_cell[c]) for c in rdr.index], index=rdr.index)
    known = cell_group.notna().values

    if a.method == "divisive":
        lab = divisive(Xfeat, a, np.random.RandomState(a.seed))
        a.k = int(lab.max()) + 1
    elif kmin < kmax:
        Dm = squareform(pdist(X))
        rows = []
        for k in range(kmin, kmax + 1):
            lab = cluster(k)
            rows.append({"K": k, "silhouette": silhouette(Dm, lab), "ARI vs truth": ari(cell_group[known], lab[known]),
                         "smallest cluster": int(np.bincount(lab).min()), "_lab": lab})
        free = max(rows, key=lambda r: r["silhouette"])
        ok = [r for r in rows if r["smallest cluster"] >= a.min_size]
        if ok:
            best = max(ok, key=lambda r: r["silhouette"])
        else:
            best = free
            print(f"[cluster_cells] no K in {kmin}-{kmax} has all clusters >= {a.min_size} cells; using the highest silhouette")
        tab = pd.DataFrame(rows).drop(columns="_lab")
        tab[f"all >= {a.min_size}"] = np.where(tab["smallest cluster"] >= a.min_size, "yes", "")
        tab["chosen"] = np.where(tab.K == best["K"], "*", "")
        tab.loc[(tab.K == free["K"]) & (tab.K != best["K"]), "chosen"] = "(best silhouette)"
        print(f"\n== choosing K by silhouette ({'on ' + str(a.pcs) + ' PCs' if a.pcs > 0 else 'no PCA'}, {a.method})")
        print(tab.round(4).to_string(index=False))
        if a.plot_k:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(5.5, 3.6))
            ax.plot(tab.K, tab.silhouette, "o-", label="silhouette")
            ax.plot(tab.K, tab["ARI vs truth"], "s--", label="ARI vs truth")
            ax.axvline(best["K"], color="0.6", lw=0.8, ls=":")
            ax.set_xlabel("K")
            ax.legend(fontsize=8)
            ax.set_title(f"chosen K = {best['K']} (clusters >= {a.min_size} cells; highest silhouette overall: K = {free['K']})", fontsize=8)
            fig.tight_layout()
            fig.savefig(a.plot_k)
            print(f"(plot: {a.plot_k})")
        a.k, lab = best["K"], best["_lab"]
    else:
        a.k, lab = kmin, cluster(kmin)
    clus = pd.Series(lab + 1, index=rdr.index, name="cluster")

    # 6. comparison
    res = pd.DataFrame({"barcode": clus.index, "cell": [to_cell[c] for c in clus.index], "cluster": clus.values})
    res["group"] = res.cell.map(groups)
    missing = res.group.isna().sum()
    if missing:
        print(f"[cluster_cells] {missing} cells not in the truth, left out of the comparison")
    r = res.dropna(subset=["group"])
    how = ("divisive (permutation-tested splits)" if a.method == "divisive" else
           f"{a.method}{'' if a.method == 'kmeans' else ' linkage'}, {f'{a.pcs} PCs' if a.pcs > 0 else 'no PCA'}")
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
