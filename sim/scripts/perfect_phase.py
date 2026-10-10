#!/usr/bin/env python3
"""
"Perfect phasing": the true phased heterozygous SNPs of the simulated genome,
read off the maternal and paternal haplotype assemblies' alignments to hg38.

Input: the minimap2 PAFs (with --cs) of each haplotype against hg38, the same
files eval.sh uses to lift the truth (EVAL_LIFT_DIR).

A position (hg38, 0-based p) becomes a het SNP when
  - each haplotype covers p with exactly one primary alignment (MAPQ >= --min-mapq),
    so duplicated / unaligned sequence is skipped,
  - neither haplotype has an indel within --indel-pad bp (read alignments, and so
    allele counts, are unreliable there),
  - one haplotype carries the hg38 base and the other a different base (A/C/G/T).
    Hom-alt sites and sites where both haplotypes differ from hg38 differently
    (multi-allelic) are skipped.
GT is written maternal|paternal: 1|0 = ALT on the maternal haplotype.

Outputs
  --out-vcf     VCF, one sample (--sample, default phasedgt), FORMAT GT only (uncompressed)
  --out-sites   CHROM <tab> POS <tab> REF <tab> ALT, 1-based (for bcftools -T and perfect_ad.py)
  stats on stderr
"""
import argparse
import re
import sys

import numpy as np

CS_RE = re.compile(r":(\d+)|\*([a-z])([a-z])|\+([a-z]+)|-([a-z]+)|~[a-z]{2}(\d+)[a-z]{2}")


def read_hap(paf, query_chrom, target_chrom, min_mapq, pad):
    """-> (coverage diff array length tlen+1, indel mask diff array, {pos: (ref, alt)}, tlen, stats)"""
    tlen = None
    cov = None
    indel = None
    snps = {}
    n_aln = n_skip = 0
    with open(paf) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or p[5] != target_chrom:
                continue
            if query_chrom and p[0] != query_chrom:
                continue
            tags = {x[:2]: x[5:] for x in p[12:]}
            if tags.get("tp") != "P" or "cs" not in tags or int(p[11]) < min_mapq:
                n_skip += 1
                continue
            if tlen is None:
                tlen = int(p[6])
                cov = np.zeros(tlen + 1, dtype=np.int32)
                indel = np.zeros(tlen + 1, dtype=np.int32)
            n_aln += 1
            # cs is in target (hg38 forward-strand) orientation for both strands
            t = int(p[7])
            bstart = None
            for m in CS_RE.finditer(tags["cs"]):
                match, sref, salt, ins, dele, intron = m.groups()
                if match is not None or sref is not None:
                    if bstart is None:
                        bstart = t
                    if sref is not None:
                        if sref in "acgt" and salt in "acgt" and t not in snps:
                            snps[t] = (sref.upper(), salt.upper())
                        else:                              # N base, or a second alignment here: never call it
                            snps[t] = None
                        t += 1
                    else:
                        t += int(match)
                    continue
                if bstart is not None:                     # close the gap-free block
                    cov[bstart] += 1
                    cov[t] -= 1
                    bstart = None
                if ins is not None:                        # between t-1 and t
                    a, b = max(t - 1 - pad, 0), min(t + pad, tlen)
                elif dele is not None:
                    a, b = max(t - pad, 0), min(t + len(dele) + pad, tlen)
                    t += len(dele)
                else:                                      # intron (not expected for asm alignments)
                    a, b = max(t - pad, 0), min(t + int(intron) + pad, tlen)
                    t += int(intron)
                indel[a] += 1
                indel[b] -= 1
            if bstart is not None:
                cov[bstart] += 1
                cov[t] -= 1
            if t != int(p[8]):
                sys.exit(f"ERROR: cs walk ended at {t}, PAF says target end {p[8]} ({paf}, query {p[0]}:{p[2]})")
    if tlen is None:
        sys.exit(f"ERROR: no primary alignments onto {target_chrom} with cs tags and MAPQ >= {min_mapq} in {paf}")
    return np.cumsum(cov)[:tlen], np.cumsum(indel)[:tlen] > 0, snps, tlen, (n_aln, n_skip)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mat-paf", required=True)
    ap.add_argument("--pat-paf", required=True)
    ap.add_argument("--target-chrom", required=True, help="hg38 chromosome name (e.g. chr1)")
    ap.add_argument("--query-chrom", default="", help="haplotype contig to use (default: any)")
    ap.add_argument("--min-mapq", type=int, default=5)
    ap.add_argument("--indel-pad", type=int, default=10)
    ap.add_argument("--sample", default="phasedgt")
    ap.add_argument("--out-vcf", required=True)
    ap.add_argument("--out-sites", required=True)
    a = ap.parse_args()

    cm, im, sm, tlen_m, st_m = read_hap(a.mat_paf, a.query_chrom, a.target_chrom, a.min_mapq, a.indel_pad)
    cp, ip, sp, tlen_p, st_p = read_hap(a.pat_paf, a.query_chrom, a.target_chrom, a.min_mapq, a.indel_pad)
    if tlen_m != tlen_p:
        sys.exit(f"ERROR: {a.target_chrom} length differs between PAFs ({tlen_m} vs {tlen_p})")
    tlen = tlen_m
    usable = (cm == 1) & (cp == 1) & ~im & ~ip

    stats = dict(candidates=0, het=0, hom_alt=0, multiallelic=0, not_unique=0, near_indel=0, conflict=0)
    rows = []
    for pos in sorted(set(sm) | set(sp)):
        stats["candidates"] += 1
        m, p_ = sm.get(pos, False), sp.get(pos, False)
        if m is None or p_ is None:
            stats["conflict"] += 1
            continue
        if not usable[pos]:
            if cm[pos] != 1 or cp[pos] != 1:
                stats["not_unique"] += 1
            else:
                stats["near_indel"] += 1
            continue
        if m and p_:
            if m[1] == p_[1]:
                stats["hom_alt"] += 1
            else:
                stats["multiallelic"] += 1
            continue
        ref, alt = m if m else p_
        gt = "1|0" if m else "0|1"
        rows.append((pos + 1, ref, alt, gt))
        stats["het"] += 1

    with open(a.out_vcf, "w") as v, open(a.out_sites, "w") as s:
        v.write("##fileformat=VCFv4.2\n")
        v.write("##source=perfect_phase.py (true phase from the HG002 haplotype -> hg38 alignments)\n")
        v.write(f"##contig=<ID={a.target_chrom},length={tlen}>\n")
        v.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="True phased genotype, maternal|paternal">\n')
        v.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{a.sample}\n")
        for pos, ref, alt, gt in rows:
            v.write(f"{a.target_chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
            s.write(f"{a.target_chrom}\t{pos}\t{ref}\t{alt}\n")

    both = int(((cm == 1) & (cp == 1)).sum())
    n_mat = sum(1 for r in rows if r[3] == "1|0")
    print(f"[perfect_phase] alignments used: mat {st_m[0]} (skipped {st_m[1]}), pat {st_p[0]} (skipped {st_p[1]})", file=sys.stderr)
    print(f"[perfect_phase] {a.target_chrom}: {both:,} bp ({both / tlen:.1%}) covered once by both haplotypes", file=sys.stderr)
    print(f"[perfect_phase] SNP positions seen: {stats['candidates']:,}; het: {stats['het']:,} "
          f"(ALT on mat {n_mat:,}, on pat {stats['het'] - n_mat:,}); hom-alt {stats['hom_alt']:,}; "
          f"multi-allelic {stats['multiallelic']:,}; skipped: not covered once by both {stats['not_unique']:,}, "
          f"within {a.indel_pad} bp of an indel {stats['near_indel']:,}, conflicting alignments {stats['conflict']:,}",
          file=sys.stderr)
    if stats["het"] == 0:
        sys.exit("ERROR: no het SNPs found")
    if rows:
        span = (rows[-1][0] - rows[0][0]) / max(len(rows) - 1, 1)
        print(f"[perfect_phase] one het SNP per {span:,.0f} bp on average", file=sys.stderr)


if __name__ == "__main__":
    main()
