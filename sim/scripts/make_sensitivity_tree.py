#!/usr/bin/env python3
"""Tree for a sensitivity study: how CNA size x number of carrying cells x bin size decide detection.

2000 cells in a fixed hierarchy of clades (every clade is a block of consecutive cells):
    1 x 2000 (root = trunk)  >  2 x 1000  >  4 x 500  >  8 x 200  >  20 x 100
    500 = 200 + 200 + 100,  200 = 100 + 100,  so there are 20 leaves of 100 cells.
Every clade carries events of ONE size S at k locations (k identical events on the same cells, same type and
haplotype). Clades of the same size N are disjoint, so their events never share cells: replicates at one (N, S)
sit on different cells, whereas k > 1 puts identical events on the SAME cells (they can be clustered together,
which is the signal-boost effect). Which clade gets which (S, k) is shuffled, so position in the tree does not
line up with the design. Locations are random on chr1, with a minimum gap between events.

With --batches B (config.sh NUM_BATCHES) the tree is ONE batch: the clade sizes above divided by B (B = 20: 100 cells,
leaves of 5). The pooled data, batch1_cell0 ... batchB_cell<tree cells - 1>, then has the clade sizes above, and the
event positions are the same for any B with the same --seed. first_cell/last_cell in the design table are per batch.

Writes  <out>.tree  (first line = tree, then the CNV table, as scDNA_sim.py reads it; no comments allowed)
        <out>.design.tsv  one row per event: clade, cells (first/last cell index), N, S, k, type, hap, location
        <out>.clades.tsv  one row per clade
"""
import argparse
import random
import sys

K = 1000
M = 1000 * K
# (size, k) of each clade, by number of cells N. The clades of one N are assigned to these in shuffled order.
DESIGN = {
    2000: [(100 * K, 1)],
    1000: [(100 * K, 1), (500 * K, 1)],
    500: [(100 * K, 1), (200 * K, 1), (500 * K, 1), (2 * M, 1)],
    200: [(s, 1) for s in (100 * K, 200 * K, 500 * K, 1 * M, 2 * M, 5 * M)] + [(200 * K, 2), (1 * M, 2)],
    100: [(s, k) for s in (100 * K, 200 * K, 500 * K, 1 * M) for k in (1, 2, 3)]      # 12
         + [(s, 1) for s in (100 * K, 200 * K, 500 * K, 1 * M)]                        # 4 second replicates
         + [(2 * M, 1), (2 * M, 2), (5 * M, 1), (5 * M, 2)],
}
# --profile small: no events on the 1000 / 2000-cell clades, more replicates at 100-500 cells (sizes 100 kb - 2 Mb)
DESIGN_SMALL = {
    2000: [None],
    1000: [None, None],
    500: [(100 * K, 1), (200 * K, 1), (500 * K, 1), (1 * M, 1)],
    200: [(s, 1) for s in (100 * K, 200 * K, 500 * K, 1 * M, 2 * M)] * 1 + [(100 * K, 1), (200 * K, 1), (500 * K, 1)],
    100: [(s, 1) for s in (100 * K, 200 * K, 500 * K, 1 * M)] * 4 + [(2 * M, 1)] * 2 + [(100 * K, 2), (200 * K, 2)],
}
# chr1 regions used, in the simulator's (assembly) coordinates; the centromere / 1q12 region (~110-160 Mb) is left out
USABLE = [(8 * M, 112 * M), (160 * M, 232 * M)]


class Clade:
    def __init__(self, name, cells, children=(), leaf=False):
        self.name, self.cells, self.children, self.leaf = name, cells, list(children), leaf
        self.events = []          # ids
        self.first = None         # first cell index (left-to-right leaf order)

    def string(self):
        ev = "[{}]".format(",".join(map(str, self.events))) if self.events else ""
        if self.leaf:
            return "{}{}".format(self.cells, ev)
        return "(" + ",".join(c.string() for c in self.children) + ")" + ev


def build_clades(B=1):
    """Clade tree for ONE batch: the pooled sizes (2000, 1000, ...) divided by the number of batches B."""
    n = {"L": 0, "K": 0, "J": 0}

    def leaf():
        n["M"] = n.get("M", 0) + 1
        return Clade("M{}".format(n["M"]), 100 // B, leaf=True)

    def l200():
        n["L"] += 1
        return Clade("L{}".format(n["L"]), 200 // B, [leaf(), leaf()])

    def k500():
        n["K"] += 1
        return Clade("K{}".format(n["K"]), 500 // B, [l200(), l200(), leaf()])

    def j1000():
        n["J"] += 1
        return Clade("J{}".format(n["J"]), 1000 // B, [k500(), k500()])

    return Clade("root", 2000 // B, [j1000(), j1000()])


def walk(c):
    yield c
    for ch in c.children:
        yield from walk(ch)


def place(sizes, gap, rng, tries=20000):
    """Random non-overlapping positions (start, end inclusive) for events of these sizes, >= gap apart, inside USABLE."""
    total = sum(e - s for s, e in USABLE)
    need = sum(sizes) + gap * (len(sizes) - 1)
    if need > total:
        raise SystemExit("events need {:.1f} Mb with gaps but only {:.1f} Mb are usable: lower --gap or the design".format(need / M, total / M))
    slack = total - need
    n = len(sizes)
    for _ in range(tries):
        order = list(range(n))
        rng.shuffle(order)
        cuts = sorted(rng.uniform(0, slack) for _ in range(n))
        res, used, ok = {}, 0, True
        for rank, i in enumerate(order):
            x = int(round((cuts[rank] + used) / K)) * K       # position in the concatenated usable space
            used += sizes[i] + gap
            acc = 0
            for a, b in USABLE:
                if x < acc + (b - a):
                    st = a + (x - acc)
                    if st + sizes[i] > b:                      # would straddle the two regions
                        ok = False
                    res[i] = (st, st + sizes[i] - 1)
                    break
                acc += b - a
            else:
                ok = False
            if not ok:
                break
        if ok:
            return [res[i] for i in range(n)]
    raise SystemExit("could not place the events: lower --gap or the design")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", required=True, help="output prefix, e.g. data/trees/sens_del")
    ap.add_argument("--type", default="DEL", choices=["DEL", "AMP"], help="event type (default DEL: one haplotype loses the segment; AMP: one extra copy)")
    ap.add_argument("--hap", default="M", choices=["M", "P"], help="haplotype of every event (default M)")
    ap.add_argument("--gap", type=float, default=2.0, help="minimum gap between events in Mb (default 2)")
    ap.add_argument("--batches", type=int, default=1, help="number of simulation batches (config.sh NUM_BATCHES); the tree is for ONE batch, so its clade sizes are the pooled ones divided by this (default 1). Must divide 100")
    ap.add_argument("--profile", default="spread", choices=["spread", "small"], help="spread: the full design (default); small: no events on the 1000/2000-cell clades, more replicates at 100-500 cells")
    ap.add_argument("--seed", type=int, default=1, help="seed of the clade shuffle and of the positions (default 1)")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    B = a.batches
    if B < 1 or 100 % B:
        sys.exit("--batches must divide 100 (the smallest clade has 100 cells in the pool)")
    root = build_clades(B)
    clades = list(walk(root))
    by_n = {}
    for c in clades:
        by_n.setdefault(c.cells * B, []).append(c)      # keyed by the POOLED number of cells
    design = DESIGN if a.profile == "spread" else DESIGN_SMALL
    for n, cs in by_n.items():
        if len(cs) != len(design[n]):
            sys.exit("design for N={} has {} entries but the tree has {} clades".format(n, len(design[n]), len(cs)))
        d = list(design[n])
        rng.shuffle(d)
        for c, sk in zip(cs, d):
            c.size, c.k = sk if sk else (0, 0)

    # first cell of every clade (leaves in left-to-right order = the simulator's cell numbering)
    nxt = [0]

    def number(c):
        c.first = nxt[0]
        if c.leaf:
            nxt[0] += c.cells
        for ch in c.children:
            number(ch)
    number(root)

    events = []                      # (clade, copy)
    for c in clades:
        for j in range(c.k):
            events.append((c, j + 1))
    locs = place([c.size for c, _ in events], int(a.gap * M), rng)

    rows = []
    for eid, ((c, j), (s, e)) in enumerate(zip(events, locs), start=1):
        c.events.append(eid)
        rows.append((eid, c, j, s, e))

    with open(a.out + ".tree", "w") as f:
        f.write(root.string() + "\n")
        f.write("id\ttype\thap\tvec\tloc\n")
        for eid, c, j, s, e in rows:
            f.write("{}\t{}\t{}\t\tchr1:{}-{}\n".format(eid, a.type, a.hap, s, e))
    with open(a.out + ".design.tsv", "w") as f:
        f.write("event_id\tclade\tclade_cells\ttree_cells\tfirst_cell\tlast_cell\tbatches\tsize_bp\tk_same_cells\tcopy\ttype\thap\tchrom\tstart\tend\n")
        for eid, c, j, s, e in rows:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tchr1\t{}\t{}\n".format(
                eid, c.name, c.cells * B, c.cells, c.first, c.first + c.cells - 1, B, c.size, c.k, j, a.type, a.hap, s, e))
    with open(a.out + ".clades.tsv", "w") as f:
        f.write("clade\tclade_cells\ttree_cells\tfirst_cell\tlast_cell\tbatches\tsize_bp\tk\n")
        for c in clades:
            if c.events:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(c.name, c.cells * B, c.cells, c.first, c.first + c.cells - 1, B, c.size, c.k))
    tot = sum(c.size * c.k for c in clades)
    print("{} events, {:.1f} Mb of CNA, {} cells per batch x {} batches; wrote {}.tree / .design.tsv / .clades.tsv".format(len(rows), tot / M, root.cells, B, a.out))


if __name__ == "__main__":
    main()
