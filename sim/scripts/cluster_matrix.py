#!/usr/bin/env python3
"""cluster_matrix.py -- which predicted bin clusters the oracle clusters ended up in (for diagnosis).

    python3 scripts/cluster_matrix.py ORACLE RUN [--eval EVAL_DIR] [--plot out.pdf] [--mb]

ORACLE   oracle folder from oracle_clusters.py (clusters.tsv, summary.tsv, profiles.tsv), e.g.
         chisel_pbins8_oracle/perfect/oracle
RUN      a CHISEL run folder (reads calls/calls.tsv) or its calls.tsv, on the same bins

1. Contingency matrix: rows = oracle clusters (labelled with their true states), columns =
   predicted clusters, entries = bins (--mb: Mb). Columns are matched to rows with the Hungarian
   algorithm (the one-to-one matching that puts the most bins on the diagonal); predicted clusters
   left unmatched follow, largest first. Printed, written to RUN-name_matrix.tsv next to --plot or
   in the current folder, and drawn with --plot.

2. With --eval EVAL_DIR (the run's eval folder; its calls.tsv has simulated cell names and the
   scored copy numbers, plus the per-cell raw calls when corrected ones were scored), per oracle
   cluster and cell: the call covering most of the cluster's bp in that cell vs the cell's true
   (mat, pat). Reported: share of cells with the right allele-specific state (a|b in either order)
   and the right total, for the corrected and the raw calls, and the most common wrong calls. Raw
   right but corrected wrong = the clone consensus (CHISEL's Cloner) changed it.

3. With --eval: true cell groups (cells with identical copy numbers in every oracle cluster, i.e.
   the tree's leaves) vs CHISEL's clones.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def nchrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def read_oracle(path):
    d = path if os.path.isdir(path) else os.path.dirname(path)
    f = os.path.join(d, "clusters.tsv") if os.path.isdir(path) else path
    o = pd.read_csv(f, sep="\t", dtype={"#CHR": str}).rename(
        columns={"#CHR": "chrom", "START": "start", "END": "end", "CLUSTER": "true"})
    o["chrom"] = o.chrom.map(nchrom)
    summ = os.path.join(d, "summary.tsv")
    summ = pd.read_csv(summ, sep="\t") if os.path.isfile(summ) else None
    prof = os.path.join(d, "profiles.tsv")
    prof = pd.read_csv(prof, sep="\t", dtype={"CELL": str}) if os.path.isfile(prof) else None
    return o, summ, prof, d


def read_run(path):
    f = os.path.join(path, "calls", "calls.tsv") if os.path.isdir(path) else path
    if not os.path.isfile(f):
        sys.exit(f"ERROR: no calls file {f}")
    head = open(f).readline().lstrip("#").strip().split("\t")
    col = head.index("CLUSTER") if "CLUSTER" in head else 10
    c = pd.read_csv(f, sep="\t", comment="#", header=None, usecols=[0, 1, 2, col],
                    names=["chrom", "start", "end", "pred"], dtype={"chrom": str})
    c = c.drop_duplicates(["chrom", "start", "end"])
    c["chrom"] = c.chrom.map(nchrom)
    return c


def label(summ, t, short=False):
    if summ is None:
        return f"O{t}"
    r = summ[summ.CLUSTER == t]
    if r.empty:
        return f"O{t}"
    r = r.iloc[0]
    st = str(r.STATES)
    st = st if len(st) <= 34 else st[:31] + "..."
    return f"O{t} {st}" if short else f"O{t} ({r.MB:g} Mb) {st}"


def matrix(m, weight):
    rows = sorted(m.true.unique())
    cols = sorted(m.pred.unique())
    tab = m.pivot_table(index="true", columns="pred", values=weight, aggfunc="sum", fill_value=0)
    tab = tab.reindex(index=rows, columns=cols, fill_value=0)
    r, c = linear_sum_assignment(-tab.values)               # maximize the bins on the diagonal
    match = {tab.index[i]: tab.columns[j] for i, j in zip(r, c)}
    order = [match[t] for t in tab.index if t in match]
    rest = sorted([p for p in tab.columns if p not in order], key=lambda p: -tab[p].sum())
    tab = tab[order + rest]
    matched = sum(tab.loc[t, match[t]] for t in match)
    return tab, match, matched


def per_cell_calls(o, prof, ev_dir):
    f = os.path.join(ev_dir, "calls.tsv")
    if not os.path.isfile(f):
        sys.exit(f"ERROR: no {f} (the standardized calls eval writes)")
    e = pd.read_csv(f, sep="\t", dtype={"cell": str, "chrom": str})
    e["chrom"] = e.chrom.map(nchrom)
    e = e.merge(o[["chrom", "start", "end", "true"]], on=["chrom", "start", "end"], how="inner")
    if e.empty:
        sys.exit(f"ERROR: the bins of {f} do not match the oracle's")
    e["bp"] = e.end - e.start
    kinds = [("scored", "CN_A", "CN_B")]
    if "CN_A_raw" in e.columns:
        kinds.append(("raw", "CN_A_raw", "CN_B_raw"))
    out = {}
    for name, ca, cb in kinds:
        e["state"] = e[ca].astype(int).astype(str) + "|" + e[cb].astype(int).astype(str)
        maj = (e.groupby(["true", "cell", "state"]).bp.sum().reset_index()
               .sort_values("bp", ascending=False).drop_duplicates(["true", "cell"]))
        maj = maj.merge(prof.rename(columns={"CLUSTER": "true", "CELL": "cell"}), on=["true", "cell"], how="inner")
        ab = maj.state.str.split("|", expand=True).astype(int)
        maj["call_lo"], maj["call_hi"] = ab.min(axis=1), ab.max(axis=1)
        maj["true_lo"], maj["true_hi"] = maj[["MAT", "PAT"]].min(axis=1), maj[["MAT", "PAT"]].max(axis=1)
        maj["ok_state"] = (maj.call_lo == maj.true_lo) & (maj.call_hi == maj.true_hi)
        maj["ok_total"] = (maj.call_lo + maj.call_hi) == (maj.true_lo + maj.true_hi)
        out[name] = maj
    clones = e.drop_duplicates("cell").set_index("cell").clone if "clone" in e.columns else None
    return out, clones


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("oracle")
    ap.add_argument("run")
    ap.add_argument("--eval", help="the run's eval folder (per-cell call accuracy, clones vs true groups)")
    ap.add_argument("--plot", help="draw the matrix to this PDF/PNG")
    ap.add_argument("--mb", action="store_true", help="matrix in Mb instead of bins")
    a = ap.parse_args()

    o, summ, prof, odir = read_oracle(a.oracle)
    c = read_run(a.run)
    m = o.merge(c, on=["chrom", "start", "end"], how="left")
    if m.pred.isna().any():
        sys.exit(f"ERROR: {int(m.pred.isna().sum())} oracle bins are not in the run's calls (different bins?)")
    m["pred"] = m.pred.astype(int)
    m["bins"] = 1
    m["mb"] = (m.end - m.start) / 1e6
    weight = "mb" if a.mb else "bins"
    tab, match, matched = matrix(m, weight)
    tot = tab.values.sum()

    run_name = os.path.normpath(a.run).replace(os.sep, "_").strip("_")[-60:]
    disp = tab.copy()
    disp.index = [label(summ, t, short=True) for t in tab.index]
    disp.columns = [f"P{p}" + ("" if p in match.values() else "*") for p in tab.columns]
    fmt = (lambda v: f"{v:.2f}" if a.mb else f"{int(v)}")
    with pd.option_context("display.width", 300, "display.max_columns", 60, "display.max_colwidth", 45):
        print(f"\n== {a.run}: oracle clusters (rows) x predicted clusters (columns), {'Mb' if a.mb else 'bins'}; "
              f"columns Hungarian-matched to rows, * = unmatched")
        print(disp.apply(lambda col: col.map(fmt)).to_string())
    print(f"\nbins on the matched diagonal: {matched:.4g} of {tot:.4g} ({100 * matched / tot:.1f}%)")
    # one-to-one matching can pair an oracle cluster with a predicted cluster that is mostly another
    # oracle cluster's: show the matched cluster's purity, and flag clusters merged into normal
    normal = o[o.NORMAL == 1].true.iloc[0] if "NORMAL" in o.columns and (o.NORMAL == 1).any() else None
    normal_main = tab.loc[normal].idxmax() if normal is not None else None
    per = []
    for t in tab.index:
        row = tab.loc[t]
        main = row.idxmax()
        merged = normal is not None and t != normal and main == normal_main
        per.append((label(summ, t), int(m[m.true == t].bins.sum()), match.get(t, "-"),
                    f"{100 * row[match[t]] / row.sum():.1f}%" if t in match else "-",
                    f"{100 * tab.loc[t, match[t]] / tab[match[t]].sum():.1f}%" if t in match else "-",
                    int((row > 0).sum()), "MERGED INTO NORMAL" if merged else ""))
    print(pd.DataFrame(per, columns=["oracle cluster", "bins", "matched P", "share of row in P", "purity of P",
                                     "P clusters used", ""]).to_string(index=False))
    print("(purity of P = share of the matched predicted cluster's bins that belong to this oracle cluster)")
    out_tsv = (os.path.splitext(a.plot)[0] if a.plot else run_name) + "_matrix.tsv"
    full = pd.concat({"bins": matrix(m, "bins")[0], "mb": matrix(m, "mb")[0]}, axis=1)
    full.to_csv(out_tsv, sep="\t")
    print(f"(matrix in bins and Mb: {out_tsv})")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        frac = tab.div(tab.sum(axis=1), axis=0).values
        fig, ax = plt.subplots(figsize=(1.2 + 0.55 * tab.shape[1] + 3.5, 1.0 + 0.45 * tab.shape[0]))
        ax.imshow(frac, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                v = tab.values[i, j]
                if v:
                    ax.text(j, i, fmt(v), ha="center", va="center", fontsize=7,
                            color="white" if frac[i, j] > 0.6 else "black")
            p = match.get(tab.index[i])
            if p is not None:
                j = list(tab.columns).index(p)
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec="crimson", lw=1.2))
        ax.set_xticks(range(tab.shape[1]))
        ax.set_xticklabels(disp.columns, fontsize=7)
        ax.set_yticks(range(tab.shape[0]))
        ax.set_yticklabels([label(summ, t) for t in tab.index], fontsize=7)
        ax.set_xlabel("predicted cluster (Hungarian-matched to the rows; * unmatched)")
        ax.set_title(f"{a.run}\n{'Mb' if a.mb else 'bins'} per oracle x predicted cluster; colour = share of the row; "
                     f"{100 * matched / tot:.1f}% on the matched diagonal", fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot)
        print(f"(plot: {a.plot})")

    if a.eval:
        if prof is None:
            sys.exit(f"ERROR: {odir}/profiles.tsv is missing: rerun oracle_clusters.py (it writes it now) with the same --out")
        res, clones = per_cell_calls(o, prof, a.eval)
        rows = []
        for t in sorted(prof.CLUSTER.unique()):
            r = {"oracle cluster": label(summ, t)}
            for name, maj in res.items():
                g = maj[maj.true == t]
                if g.empty:
                    continue
                changed = g[(g.MAT != 1) | (g.PAT != 1)]
                r[f"{name}: state ok"] = f"{100 * g.ok_state.mean():.1f}%"
                r[f"{name}: total ok"] = f"{100 * g.ok_total.mean():.1f}%"
                r[f"{name}: CNA cells ok"] = f"{100 * changed.ok_state.mean():.1f}% of {len(changed)}" if len(changed) else "-"
                wrong = g[~g.ok_state].state.value_counts()
                r[f"{name}: top wrong calls"] = ", ".join(f"{s} x{n}" for s, n in wrong.head(3).items()) or "-"
            rows.append(r)
        print(f"\n== per oracle cluster: each cell's call (covering most of the cluster's bp) vs its true state ({a.eval})")
        print("   scored = what eval scored (corrected = clone consensus by default); raw = the cell's own call before cloning")
        with pd.option_context("display.width", 300, "display.max_columns", 30, "display.max_colwidth", 50):
            print(pd.DataFrame(rows).to_string(index=False))
        if clones is not None:
            groups = prof.assign(st=prof.MAT.astype(str) + "|" + prof.PAT.astype(str)).pivot(index="CELL", columns="CLUSTER", values="st")
            key = groups.apply(lambda r: "/".join(r.values), axis=1)
            names = {k: f"G{i}" for i, k in enumerate(key.value_counts().index)}
            grp = key.map(names)
            both = pd.DataFrame({"group": grp}).join(clones.rename("clone"), how="inner")
            print("\n== true cell groups (identical copy numbers in every oracle cluster) vs CHISEL clones (cells)")
            print(pd.crosstab(both.group, both.clone, margins=True).to_string())
            desc = []
            for k, gname in names.items():
                parts = [f"O{t}:{s}" for t, s in zip(groups.columns, k.split("/")) if s != "1|1"]
                desc.append(f"  {gname}: {', '.join(parts) if parts else 'no CNA'}")
            print("\n".join(desc))


if __name__ == "__main__":
    main()
