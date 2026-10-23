#!/usr/bin/env python
"""make_combo.py -- write a small CHISEL combo.tsv (and its truth) from a design you edit below.

    python make_combo.py [OUTDIR]          # writes OUTDIR/combo.tsv and OUTDIR/truth.tsv (default: .)

combo.tsv is what Caller reads (10 columns, no header):
    CHR START END CELL NORM_COUNT COUNT RDR A_COUNT B_COUNT BAF
Counts are drawn the way sequencing makes them: COUNT ~ Poisson(reads per copy x copies), A and B
~ Binomial(SNP reads, haplotype share); RDR = (COUNT / NORM_COUNT) x (normal total / cell total),
as CHISEL's RDREstimator; BAF = B / (A + B), 0.5 without SNP reads. A is the haplotype that is
the minor one in the pseudobulk, as after CHISEL's Combiner (write events as A|B accordingly).
truth.tsv: CHR START END CELL CN_A CN_B.
Works with Python 2.7 and 3 (numpy only). Then:
    chisel_hmmstates OUTDIR/combo.tsv --mink 1 --maxk 6
    python2.7 .../Caller.py OUTDIR/combo.tsv --clustering hmm --hmmmink 1 --hmmmaxk 6 -j 4 --seed 1 > calls.tsv
"""
from __future__ import division, print_function
import os
import sys
import numpy as np

# ------------------------------------------------------------------ the design: edit this part
SEED = 0
CHROMOSOMES = [('chr1', 30), ('chr2', 20)]   # (name, number of bins), in genome order
BIN_SIZE = 1000000
NCELLS = 20
NORMAL_READS = 2000        # reads of the matched normal per bin (diploid)
CELL_READS = 2000          # mean reads per bin of a diploid cell (lower = noisier RDR)
SNP_READS = 40             # mean SNP-covering reads per bin of a diploid cell (lower = noisier BAF)
CELL_DEPTH_SD = 0.2        # spread of sequencing depth between cells (0 = all equal)
# events: (chromosome, first bin, last bin (inclusive), cells, (CN_A, CN_B)); later events win
# cells: a range/list of cell indices, or 'all'
EVENTS = [
    ('chr1', 10, 15, range(0, 10), (1, 2)),    # +1 copy of B in half of the cells
    ('chr1', 22, 23, range(0, 4), (1, 3)),     # small amplification in 4 cells (2 bins)
]
BASE = (1, 1)                                   # copy numbers everywhere else
# --------------------------------------------------------------------------------------------


def main(outdir='.'):
    rng = np.random.RandomState(SEED)
    bins = [(c, i * BIN_SIZE, (i + 1) * BIN_SIZE) for c, n in CHROMOSOMES for i in range(n)]
    index = {(b[0], b[1] // BIN_SIZE): k for k, b in enumerate(bins)}
    cn = np.zeros((len(bins), NCELLS, 2), dtype=int)
    cn[:, :, 0], cn[:, :, 1] = BASE
    for chrom, first, last, cells, state in EVENTS:
        cells = range(NCELLS) if cells == 'all' else list(cells)
        for i in range(first, last + 1):
            for e in cells:
                cn[index[(chrom, i)], e] = state
    tot = cn.sum(axis=2)
    depth = np.exp(rng.normal(0, CELL_DEPTH_SD, NCELLS))            # per-cell depth factor
    norm = rng.poisson(NORMAL_READS, len(bins))
    # reads follow the cell's DNA: a cell with more copies in total spreads its reads over more DNA
    mean_cn = tot.mean(axis=0)
    count = rng.poisson(CELL_READS * depth[None, :] * tot / mean_cn[None, :])
    rdr = (count / np.maximum(norm[:, None], 1)) * (norm.sum() / count.sum(axis=0))[None, :]
    snp = rng.poisson(SNP_READS * depth[None, :] * tot / 2.0)
    share_b = np.where(tot > 0, cn[:, :, 1] / np.maximum(tot, 1), 0.5)
    B = rng.binomial(snp, share_b)
    A = snp - B
    baf = np.where(snp > 0, B / np.maximum(snp, 1), 0.5)

    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    names = ['cell{:02d}'.format(e) for e in range(NCELLS)]
    with open(os.path.join(outdir, 'combo.tsv'), 'w') as f, open(os.path.join(outdir, 'truth.tsv'), 'w') as g:
        g.write('#CHR\tSTART\tEND\tCELL\tCN_A\tCN_B\n')
        for k, (c, s, t) in enumerate(bins):
            for e in range(NCELLS):
                f.write('\t'.join(map(str, [c, s, t, names[e], norm[k], count[k, e], '{:.6g}'.format(rdr[k, e]),
                                            A[k, e], B[k, e], '{:.6g}'.format(baf[k, e])])) + '\n')
                g.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, s, t, names[e], cn[k, e, 0], cn[k, e, 1]))
    print('{} bins x {} cells -> {}/combo.tsv, {}/truth.tsv'.format(len(bins), NCELLS, outdir, outdir))
    print('distinct bin profiles (= the clusters an ideal method finds): {}'.format(
        len(set(tuple(map(tuple, cn[k])) for k in range(len(bins))))))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')