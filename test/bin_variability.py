#!/usr/bin/env python3
"""bin_variability.py -- do bins with CNAs vary more across cells than bins without? (diagnosis)

    python3 test/bin_variability.py --oracle ORACLE_DIR --run RUN_DIR [--eval EVAL_DIR ...] [--out PREFIX]

Every bin is labelled by the oracle cluster covering most of its bp (oracle_clusters.py's
clusters.tsv): "no CNA" (the cluster where every cell is 1|1) or "CNA". --min-frac F keeps only
bins with at least F of their bp in that cluster (no boundary bins).

PREFIX_data.pdf  from RUN/combo/combo.tsv, per bin across cells:
    RDR variance, RDR dispersion = variance / expected counting noise, mean over cells of
        RDR^2 (1/COUNT + 1/NORM_COUNT)  (Poisson on the cell's and the normal's reads)
    BAF variance (B / (A + B), cells with SNP reads), BAF dispersion = variance / expected binomial
        noise, mean over cells of p (1 - p) / n with p the bin's pooled BAF and n the cell's SNP reads
    A dispersion near 1 = all cells look alike up to counting noise; above 1 = cells differ.
PREFIX_calls.pdf  from each EVAL_DIR/calls.tsv (eval's standardized calls), per bin across cells,
    for the scored (corrected by default) and raw calls: variance of the total copy number, and
    the share of cells whose state (A|B, ordered as CHISEL's Cloner compares them) is not the bin's
    most common one -- what the bin adds to Cloner's Hamming distances.
PREFIX_bins.tsv  every bin with its label and all metrics.

Per metric, printed and in the panel titles: medians, AUC (probability that a CNA bin scores higher
than a no-CNA bin; 0.5 = no separation) and the share of CNA bins above the 95th percentile of
no-CNA bins (what a threshold set on CNA-free bins would keep). Clonal CNAs (every cell) do not
vary between cells, so they score like no-CNA bins on purpose: they do not separate clones either.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd


def nchrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie), from ranks (Mann-Whitney)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[np.isfinite(pos)], neg[np.isfinite(neg)]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    r = pd.Series(np.concatenate([pos, neg])).rank().values
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def bin_labels(oracle, min_frac):
    f = os.path.join(oracle, "clusters.tsv") if os.path.isdir(oracle) else oracle
    o = pd.read_csv(f, sep="\t", dtype={"#CHR": str}).rename(columns={"#CHR": "chrom", "START": "start", "END": "end"})
    if "NORMAL" not in o.columns:
        sys.exit(f"ERROR: {f} has no NORMAL column; rerun oracle_clusters.py")
    o["chrom"] = o.chrom.map(nchrom)
    o["label"] = np.where(o.NORMAL == 1, "no CNA", "CNA")
    o = o[o.FRAC >= min_frac] if min_frac > 0 else o
    return o[["chrom", "start", "end", "label", "FRAC"]]


def data_metrics(combo, lab):
    cols = ["chrom", "start", "end", "cell", "norm", "count", "rdr", "a", "b", "baf"]
    d = pd.read_csv(combo, sep="\t", header=None, comment="#", names=cols, dtype={"chrom": str, "cell": str})
    d["chrom"] = d.chrom.map(nchrom)
    d = d.merge(lab[["chrom", "start", "end"]], on=["chrom", "start", "end"], how="inner")
    if d.empty:
        sys.exit(f"ERROR: no bin of {combo} matches the oracle's bins")
    key = ["chrom", "start", "end"]
    r = d[d["count"] > 0].copy()
    r["noise"] = r.rdr ** 2 * (1.0 / r["count"] + 1.0 / r["norm"].clip(lower=1))
    g = r.groupby(key)
    out = pd.DataFrame({"RDR variance": g.rdr.var(), "RDR dispersion": g.rdr.var() / g.noise.mean()})
    s = d[(d.a + d.b) > 0].copy()
    s["n"] = s.a + s.b
    pooled = s.groupby(key).agg(sb=("b", "sum"), sn=("n", "sum"))
    s = s.merge((pooled.sb / pooled.sn).rename("p").reset_index(), on=key)
    s["noise"] = s.p * (1 - s.p) / s.n
    g = s.groupby(key)
    out["BAF variance"] = g.baf.var()
    out["BAF dispersion"] = g.baf.var() / g.noise.mean().where(lambda x: x > 0)
    return out.reset_index()


def call_metrics(ev_dir, lab):
    f = os.path.join(ev_dir, "calls.tsv")
    if not os.path.isfile(f):
        sys.exit(f"ERROR: no {f}")
    e = pd.read_csv(f, sep="\t", dtype={"cell": str, "chrom": str})
    e["chrom"] = e.chrom.map(nchrom)
    e = e.merge(lab[["chrom", "start", "end"]], on=["chrom", "start", "end"], how="inner")
    key = ["chrom", "start", "end"]
    out = None
    kinds = [("scored", "CN_A", "CN_B")] + ([("raw", "CN_A_raw", "CN_B_raw")] if "CN_A_raw" in e.columns else [])
    for name, ca, cb in kinds:
        e["tot"] = e[ca] + e[cb]
        e["state"] = e[ca].astype(int).astype(str) + "|" + e[cb].astype(int).astype(str)
        g = e.groupby(key)
        offmode = 1.0 - g.state.agg(lambda s: s.value_counts().iloc[0] / len(s))
        m = pd.DataFrame({f"{name}: total CN variance": g.tot.var(), f"{name}: share off the modal state": offmode})
        out = m if out is None else out.join(m)
    return out.reset_index()


def summarize(df, metrics):
    rows = []
    for m in metrics:
        cna, nrm = df.loc[df.label == "CNA", m].dropna(), df.loc[df.label == "no CNA", m].dropna()
        thr = np.nanpercentile(nrm, 95) if len(nrm) else np.nan
        rows.append({"metric": m, "median no CNA": nrm.median(), "median CNA": cna.median(), "AUC": auc(cna, nrm),
                     "CNA bins > 95th pct of no-CNA": (cna > thr).mean() if len(cna) else np.nan})
    return pd.DataFrame(rows)


def panels(df, metrics, summ, title, path, ncol=2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(metrics)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.4 * nrow), squeeze=False)
    for ax, m in zip(axes.flat, metrics):
        v = df[[m, "label"]].dropna()
        pos = v[m] > 0
        logx = ("variance" in m or "dispersion" in m) and pos.mean() > 0.95
        if logx and len(v[m][pos]):                    # log axis only for values spanning a wide range
            q = np.nanpercentile(v[m][pos], [0.5, 99.5])
            logx = q[1] / max(q[0], 1e-300) > 20
        vals = v[m][pos] if logx else v[m]
        lo, hi = (np.nanpercentile(vals, 0.5), np.nanpercentile(vals, 99.5)) if len(vals) else (0, 1)
        if hi <= lo:
            hi = lo + 1
        bins = np.logspace(np.log10(lo), np.log10(hi), 50) if logx else np.linspace(lo, hi, 50)
        for lab_, col in (("no CNA", "0.55"), ("CNA", "crimson")):
            x = v.loc[(v.label == lab_) & (pos if logx else True), m]
            ax.hist(x.clip(lo, hi), bins=bins, density=True, histtype="stepfilled", alpha=0.45, color=col,
                    label=f"{lab_} ({len(x)} bins)")
        if logx:
            ax.set_xscale("log")
        if "dispersion" in m:
            ax.axvline(1.0, color="k", lw=0.8, ls="--")
        s = summ[summ.metric == m].iloc[0]
        ax.set_title(f"{m}\nAUC {s['AUC']:.2f}; CNA bins above 95th pct of no-CNA: {100 * s['CNA bins > 95th pct of no-CNA']:.0f}%", fontsize=8)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--oracle", required=True, help="oracle folder (clusters.tsv from oracle_clusters.py)")
    ap.add_argument("--run", help="CHISEL run folder whose combo/combo.tsv to read (RDR/BAF panels)")
    ap.add_argument("--combo", help="instead of --run: a combo.tsv")
    ap.add_argument("--eval", action="append", default=[], help="eval folder(s) of the run(s) to compare (calls panels); repeatable")
    ap.add_argument("--min-frac", type=float, default=0.0, help="only bins with at least this share of their bp in their oracle cluster")
    ap.add_argument("--out", default="bin_variability", help="output prefix (default: bin_variability)")
    a = ap.parse_args()
    if not (a.run or a.combo or a.eval):
        ap.error("give --run/--combo and/or --eval")

    lab = bin_labels(a.oracle, a.min_frac)
    print(f"bins: {(lab.label == 'CNA').sum()} CNA, {(lab.label == 'no CNA').sum()} no CNA"
          + (f" (only bins >= {a.min_frac:g} in their oracle cluster)" if a.min_frac > 0 else ""))
    tables = [lab]
    pd.set_option("display.width", 200)
    if a.run or a.combo:
        combo = a.combo or os.path.join(a.run, "combo", "combo.tsv")
        dm = data_metrics(combo, lab)
        df = lab.merge(dm, on=["chrom", "start", "end"])
        metrics = ["RDR variance", "RDR dispersion", "BAF variance", "BAF dispersion"]
        summ = summarize(df, metrics)
        print(f"\n== data ({combo})\n" + summ.round(3).to_string(index=False))
        panels(df, metrics, summ, f"per-bin variability across cells: {combo}", f"{a.out}_data.pdf")
        tables.append(dm)
        print(f"(plot: {a.out}_data.pdf)")
    if a.eval:
        allm, seen = [], {}
        for ev in a.eval:
            cm = call_metrics(ev, lab)
            tag = os.path.basename(os.path.normpath(ev))
            seen[tag] = seen.get(tag, 0) + 1
            tag = tag if seen[tag] == 1 else f"{tag} ({seen[tag]})"     # same folder name twice
            cm = cm.rename(columns={c: f"{tag} | {c}" for c in cm.columns if c not in ("chrom", "start", "end")})
            tables.append(cm)
            allm += [c for c in cm.columns if c not in ("chrom", "start", "end")]
        df = lab
        for t in tables[1:]:
            df = df.merge(t, on=["chrom", "start", "end"], how="left")
        summ = summarize(df, allm)
        print("\n== calls\n" + summ.round(3).to_string(index=False))
        ncol = 4 if any("raw:" in m for m in allm) else 2
        panels(df, allm, summ, "per-bin variability of the calls across cells (scored = corrected unless CHISEL_EVAL_CN=raw)",
               f"{a.out}_calls.pdf", ncol=ncol)
        print(f"(plot: {a.out}_calls.pdf)")
    out = tables[0]
    for t in tables[1:]:
        out = out.merge(t, on=["chrom", "start", "end"], how="left")
    out.to_csv(f"{a.out}_bins.tsv", sep="\t", index=False, float_format="%.6g")
    print(f"(per-bin table: {a.out}_bins.tsv)")


if __name__ == "__main__":
    main()
