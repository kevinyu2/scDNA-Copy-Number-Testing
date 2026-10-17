#!/usr/bin/env python3
"""
Draw where bins start and end, one row per bin set, to compare e.g. UGP's and the perfect
phaser's bins by eye. Each bin is a block from START to END (alternating shades, so neighbours
are distinguishable); below, the size of every bin along the chromosome (log scale) for all sets.

Reads any TSV/BED whose first three columns are CHR START END ('#' lines skipped, gzip ok).

Usage: plot_bins.py BINS1 BINS2 [...] --names UGP perfect -o bins_compare.png
"""
import argparse, gzip
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COLORS = [('#1f77b4', '#9ecae1'), ('#d62728', '#fc9272'), ('#2ca02c', '#a1d99b'), ('#9467bd', '#cbc9e2')]


def read_bins(f):
    bins = defaultdict(list)
    with (gzip.open(f, 'rt') if f.endswith('.gz') else open(f)) as i:
        for l in i:
            if l.startswith('#') or not l.strip():
                continue
            p = l.split('\t')
            c = p[0] if p[0].startswith('chr') else 'chr' + p[0]
            bins[c].append((int(p[1]), int(p[2])))
    return {c: sorted(v) for c, v in bins.items()}


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('bins', nargs='+')
ap.add_argument('--names', nargs='+', default=None)
ap.add_argument('-o', '--output', default='bins_compare.png')
a = ap.parse_args()
names = a.names or ['set{}'.format(i + 1) for i in range(len(a.bins))]
assert len(names) == len(a.bins), '--names needs one name per bins file'
sets = [read_bins(f) for f in a.bins]
chrs = sorted(set().union(*sets), key=lambda c: (len(c), c))

fig, axes = plt.subplots(2 * len(chrs), 1, figsize=(16, 4.5 * len(chrs)), squeeze=False,
                         gridspec_kw={'height_ratios': [1, 1.6] * len(chrs)})
for k, c in enumerate(chrs):
    top, low = axes[2 * k, 0], axes[2 * k + 1, 0]
    for r, (name, bins) in enumerate(zip(names, sets)):
        dark, light = COLORS[r % len(COLORS)]
        iv = bins.get(c, [])
        y = len(sets) - 1 - r
        for x, (s, e) in enumerate(iv):
            top.add_patch(Rectangle((s / 1e6, y + 0.1), (e - s) / 1e6, 0.8, color=dark if x % 2 == 0 else light, lw=0))
        for s, e in iv:
            low.hlines((e - s) / 1e3, s / 1e6, e / 1e6, color=dark, lw=2, alpha=0.8)
        low.plot([], [], color=dark, lw=2, label='{} ({} bins)'.format(name, len(iv)))
    end = max(e for b in sets for s, e in b.get(c, [(0, 1)])) / 1e6
    top.set_xlim(0, end)
    top.set_ylim(0, len(sets))
    top.set_yticks([len(sets) - 1 - r + 0.5 for r in range(len(sets))])
    top.set_yticklabels(names)
    top.set_title('{}: bins (alternating shades)'.format(c))
    low.set_xlim(0, end)
    low.set_yscale('log')
    low.set_ylabel('bin size (kb)')
    low.set_xlabel('{} (Mb)'.format(c))
    low.legend(fontsize=8, loc='upper right')
fig.tight_layout()
fig.savefig(a.output, dpi=110)
print('Plot: {}'.format(a.output))
