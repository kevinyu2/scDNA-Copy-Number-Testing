#!/usr/bin/env python3
"""hmm_noise_plot.py -- does the HMM's noise model match the data, at every bin length?

    python3 test/hmm_noise_plot.py RUN/calls/hmm_noise.tsv [more hmm_noise.tsv ...] [--out noise.pdf]
        [--state all|largest|K] [--nbins 12] [--label NAME ...]

hmm_noise.tsv is written by the chisel_hmm fork's HMM next to calls/hmm_selection.tsv: one line per bin
and feature (each sample's BAF and RDR) with the observation X, its state's MEAN, the pure counting
variance COUNT_VAR (Poisson reads, binomial SNP reads), the variance the model assumes MODEL_VAR
(counting noise: tau2 + phi * the sample's counting variance + the normal's; legacy: the state's own
variance) and RESID2 = (X - MEAN)^2.

Page 1, one panel per feature: x = bin length, y = variance (both log). Faint dots: RESID2 of every bin.
Blue: the mean RESID2 in bins of length (quantiles), i.e. the observed variance. Orange: the mean
MODEL_VAR there. Dashed grey: the mean COUNT_VAR (phi = 1, tau2 = 0). If the model is right, blue and
orange agree at every length; if blue rises above orange for short bins, short bins are noisier than the
model allows (what made the legacy HMM split normal by length).
Page 2: observed / model per length bin (1 = right), every feature.
Printed: per feature, phi and tau2 (counting noise) and observed / model for the shortest and longest
third of the bins, with the rough uncertainty of such a ratio, sqrt(2 / bins in a third): with few bins
(a small state, --state K) a ratio of 0.7 can be chance.
--state: all bins (default), only the state with most bins ('largest', usually normal; the cleanest
check, as a misassigned bin inflates residuals) or one state number. Several reports (e.g. the legacy and
the counting run of the same data) give one set of pages each, in one PDF.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

OBS, MODEL, COUNT, DOTS = "#2a78d6", "#eb6834", "#52514e", "#9ec5f4"   # reference palette: blue, orange, secondary ink


def load(path, state):
    d = pd.read_csv(path, sep="\t")
    d.columns = [c.lstrip("#") for c in d.columns]
    need = {"LENGTH", "STATE", "SAMPLE", "FEATURE", "MODEL_VAR", "RESID2"}
    if not need <= set(d.columns):
        sys.exit(f"ERROR: {path} lacks columns {sorted(need - set(d.columns))} (an hmm_noise.tsv from the chisel_hmm HMM?)")
    if "VALID" in d.columns:
        d = d[d.VALID == 1]
    d["COUNT_VAR"] = pd.to_numeric(d.get("COUNT_VAR", np.nan), errors="coerce")
    if state == "largest":
        k = d.drop_duplicates(["CHR", "START"]).STATE.value_counts().idxmax()
        d = d[d.STATE == k]
        note = f"state {k} only (the one with most bins)"
    elif state != "all":
        d = d[d.STATE == int(state)]
        note = f"state {state} only"
    else:
        note = "all bins"
    if d.empty:
        sys.exit(f"ERROR: no rows left in {path} ({note})")
    d["feat"] = d.SAMPLE.astype(str) + " " + d.FEATURE.astype(str)
    return d, note


def by_length(g, nbins):
    q = np.unique(np.quantile(g.LENGTH, np.linspace(0, 1, nbins + 1)))
    if len(q) < 3:
        q = np.array([g.LENGTH.min(), g.LENGTH.max() + 1])
    cut = pd.cut(g.LENGTH, q, include_lowest=True)
    t = g.groupby(cut, observed=True).agg(length=("LENGTH", "median"), obs=("RESID2", "mean"), model=("MODEL_VAR", "mean"),
                                          count=("COUNT_VAR", "mean"), n=("RESID2", "size"))
    return t[t.n >= 5]


def summary(d):
    rows = []
    for f, g in d.groupby("feat", sort=False):
        lo, hi = np.quantile(g.LENGTH, [1 / 3, 2 / 3])
        s, l = g[g.LENGTH <= lo], g[g.LENGTH >= hi]
        r = {"feature": f, "bins": len(g), "+/- (a third)": np.sqrt(2.0 / max(len(g) / 3.0, 1))}
        for c in ("PHI", "TAU2"):
            if c in g.columns and pd.to_numeric(g[c], errors="coerce").notna().any():
                r[c.lower()] = pd.to_numeric(g[c], errors="coerce").median()
        r["obs/model short third"] = s.RESID2.mean() / s.MODEL_VAR.mean()
        r["obs/model long third"] = l.RESID2.mean() / l.MODEL_VAR.mean()
        if g.COUNT_VAR.notna().any():
            r["obs/counting short"] = s.RESID2.mean() / s.COUNT_VAR.mean()
            r["obs/counting long"] = l.RESID2.mean() / l.COUNT_VAR.mean()
        rows.append(r)
    return pd.DataFrame(rows)


def kb_axis(ax, lo, hi):
    """Log x axis in kb with plain tick labels."""
    from matplotlib.ticker import FixedLocator, NullFormatter, FuncFormatter
    ticks = [t for t in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000) if lo <= t <= hi]
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors="#52514e", labelsize=7)
    ax.grid(True, which="major", color="#ecebe8", linewidth=0.6)
    ax.set_axisbelow(True)


def pages(d, title, nbins, pdf, plt):
    feats = list(dict.fromkeys(d.feat))
    ncol = min(4, len(feats))
    nrow = int(np.ceil(len(feats) / ncol))
    rng = np.random.RandomState(0)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.0 * nrow + 0.8), squeeze=False)
    for ax, f in zip(axes.flat, feats):
        g = d[d.feat == f]
        kb = g.LENGTH / 1e3
        show = g.sample(min(len(g), 4000), random_state=rng) if len(g) > 4000 else g
        ax.scatter(show.LENGTH / 1e3, np.maximum(show.RESID2, 1e-12), s=4, color=DOTS, alpha=0.35, linewidths=0, zorder=1,
                   label="squared residual, per bin")
        t = by_length(g, nbins)
        if t["count"].notna().any():
            ax.plot(t.length / 1e3, t["count"], color=COUNT, lw=1.5, ls=(0, (4, 3)), zorder=2, label="counting variance (Poisson / binomial)")
        ax.plot(t.length / 1e3, t.model, color=MODEL, lw=2, marker="o", ms=4, zorder=3, label="model variance")
        ax.plot(t.length / 1e3, t.obs, color=OBS, lw=2, marker="o", ms=4, zorder=4, label="observed variance (mean of residuals)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        top = np.nanmax([t.obs.max(), t.model.max(), t["count"].max() if t["count"].notna().any() else np.nan]) * 30
        low = np.nanmin([t.obs.min(), t.model.min(), t["count"].min() if t["count"].notna().any() else np.nan]) / 30
        ax.set_ylim(low, top)
        ax.set_xlim(kb.min() * 0.9, kb.max() * 1.1)
        kb_axis(ax, kb.min() * 0.9, kb.max() * 1.1)
        ax.set_title(f, fontsize=9, color="#0b0b0b")
        ax.set_xlabel("bin length (kb)", fontsize=8, color="#52514e")
        ax.set_ylabel("variance", fontsize=8, color="#52514e")
        style(ax)
    for ax in list(axes.flat)[len(feats):]:
        ax.axis("off")
    h, l = axes.flat[0].get_legend_handles_labels()
    H = 3.0 * nrow + 0.8
    fig.suptitle(title, fontsize=10, y=1 - 0.12 / H, va="top", color="#0b0b0b")
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1 - 0.38 / H))
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.7 / H))
    pdf.savefig(fig)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    for i, f in enumerate(feats):
        t = by_length(d[d.feat == f], nbins)
        c = cols[i] if i < len(cols) else "#898781"            # more than 8 features: the rest in neutral grey
        ax.plot(t.length / 1e3, t.obs / t.model, color=c, lw=1.5, marker="o", ms=3, label=f if i < len(cols) else None)
    ax.axhline(1.0, color="#52514e", lw=1, ls=(0, (4, 3)))
    ax.set_xscale("log")
    ax.set_yscale("log")
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    lo_r, hi_r = ax.get_ylim()
    ax.yaxis.set_major_locator(FixedLocator([v for v in (0.25, 0.5, 0.67, 0.8, 1, 1.25, 1.5, 2, 4) if lo_r <= v <= hi_r]))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    lens = d.LENGTH / 1e3
    kb_axis(ax, lens.min() * 0.9, lens.max() * 1.1)
    ax.set_xlabel("bin length (kb)", fontsize=8, color="#52514e")
    ax.set_ylabel("observed / model variance", fontsize=8, color="#52514e")
    ax.set_title(f"{title}\nobserved / model variance by bin length (1 = the model fits)", fontsize=9, color="#0b0b0b")
    style(ax)
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("reports", nargs="+", help="hmm_noise.tsv file(s), or run folders (reads calls/hmm_noise.tsv)")
    ap.add_argument("--out", default="hmm_noise.pdf")
    ap.add_argument("--state", default="all", help="all (default), largest, or a state number")
    ap.add_argument("--nbins", type=int, default=12, help="length bins (quantiles) for the averages (default 12)")
    ap.add_argument("--label", action="append", default=[], help="a name per report (default: its folder)")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(a.out) as pdf:
        for i, r in enumerate(a.reports):
            path = os.path.join(r, "calls", "hmm_noise.tsv") if os.path.isdir(r) else r
            if not os.path.isfile(path):
                sys.exit(f"ERROR: no {path}")
            d, note = load(path, a.state)
            parts = os.path.normpath(os.path.abspath(path)).split(os.sep)
            # .../calls/<run>/<phaser>/calls/hmm_noise.tsv -> <run>/<phaser>
            name = a.label[i] if i < len(a.label) else ("/".join(parts[-4:-2]) if parts[-2] == "calls" and len(parts) >= 4 else path)
            model = "counting noise" if "TAU2" in d.columns and pd.to_numeric(d.TAU2, errors="coerce").notna().any() else "legacy (per-state variances)"
            title = f"{name}: {model}, {note}"
            pages(d, title, a.nbins, pdf, plt)
            with pd.option_context("display.width", 200, "display.max_columns", 20):
                print(f"\n== {title}")
                print(summary(d).round(3).to_string(index=False))
    print(f"\n(plots: {a.out})")


if __name__ == "__main__":
    main()
