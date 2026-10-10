#!/usr/bin/env python3
"""
Turn bcftools mpileup allele depths at the true het sites into the minimal
VCF HiScanner reads from SCAN2's gatk/hc_raw.mmq60.vcf.gz:
    CHROM POS . REF ALT . . . AD <ref,alt per sample>
REF/ALT are always the true alleles (from perfect_phase.py's sites file), so
every chunk has identical records and `bcftools merge` joins them cleanly.

stdin: bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n'   (one mpileup chunk)
       mpileup's ALT lists the alleles it saw plus <*>; the true ALT's depth is
       taken from its position in that list (0 if it was not seen).
Sites where every sample in the chunk has 0 reads are left out (HiScanner skips
them anyway; after merging they show as missing).
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sites", required=True, help="CHROM POS REF ALT (perfect_phase.py --out-sites, may be .gz)")
    ap.add_argument("--samples", required=True, help="file with one sample name per line, in mpileup order")
    ap.add_argument("--contig-length", type=int, required=True)
    a = ap.parse_args()

    opener = __import__("gzip").open if a.sites.endswith(".gz") else open
    truth = {}
    chrom = None
    with opener(a.sites, "rt") as f:
        for line in f:
            c, pos, ref, alt = line.split()[:4]
            chrom = chrom or c
            truth[(c, int(pos))] = (ref, alt)
    samples = [s.strip() for s in open(a.samples) if s.strip()]
    if not samples:
        sys.exit("ERROR: no samples")

    out = sys.stdout
    out.write("##fileformat=VCFv4.2\n")
    out.write("##source=perfect_ad.py (bcftools mpileup allele depths at the true het SNPs)\n")
    out.write(f"##contig=<ID={chrom},length={a.contig_length}>\n")
    out.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth of REF and ALT">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

    n_in = n_out = n_refmis = 0
    ns = len(samples)
    for line in sys.stdin:
        p = line.rstrip("\n").split("\t")
        n_in += 1
        key = (p[0], int(p[1]))
        t = truth.get(key)
        if t is None:
            continue
        ref, alt = t
        if p[2].upper() != ref:
            n_refmis += 1
            continue
        alts = p[3].split(",")
        ai = alts.index(alt) + 1 if alt in alts else None
        ads = p[4:4 + ns]
        if len(ads) != ns:
            sys.exit(f"ERROR: expected {ns} AD columns at {key}, got {len(ads)}")
        cells, any_reads = [], False
        for x in ads:
            if x == "." or x == "":
                cells.append("0,0")
                continue
            v = x.split(",")
            r = int(v[0]) if v[0] != "." else 0
            k = int(v[ai]) if ai is not None and ai < len(v) and v[ai] != "." else 0
            if r or k:
                any_reads = True
            cells.append(f"{r},{k}")
        if not any_reads:
            continue
        out.write(f"{p[0]}\t{p[1]}\t.\t{ref}\t{alt}\t.\t.\t.\tAD\t" + "\t".join(cells) + "\n")
        n_out += 1
    print(f"[perfect_ad] {n_in:,} mpileup sites in, {n_out:,} written"
          + (f", {n_refmis:,} skipped (REF differs from truth)" if n_refmis else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
