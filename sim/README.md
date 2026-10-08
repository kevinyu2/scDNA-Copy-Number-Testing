## Simulated chr1 scDNA benchmark: phasers x CN callers

Simulate ~0.1x single cells + a 30x matched bulk from the HG002 chr1 haplotypes,
phase germline SNPs, call haplotype-specific copy numbers, and score the calls
against the simulator's truth.

### Layout

```
sim/
  config.sh          every setting
  run_sim.sh         stage 0: simulate cells (SLURM array) + bulk
  run_pipeline.sh    phase -> call -> eval for one phaser/caller route
  scripts/           everything the two launchers submit (never called directly)
    common.sh          derived paths, route table, log helpers
    sim_cells.sh, sim_bulk.sh, scDNA_sim.py
    phase_scan2.sh, phase_ugp.sh, snakefile.phasing
    call_hiscanner.sh, chisel_prep.sh, call_chisel.sh
    eval.sh, lift_truth.py, standardize_calls.py, eval_calls.py
```

Always run from `sim/` on the login node; the launchers only submit jobs.
Jobs source `config.sh` when they start, so edits reach queued jobs too.

### Usage

```bash
./run_sim.sh                                          # simulate only
./run_sim.sh -- --phaser ugp --caller chisel          # simulate, then the whole route
./run_pipeline.sh --phaser ugp --caller chisel        # phase -> call -> eval
./run_pipeline.sh --phaser ugp --caller chisel --from call
./run_pipeline.sh --phaser scan2 --caller hiscanner --only eval
./run_pipeline.sh --phaser ugp --to phase             # phasing only, no caller needed
./run_pipeline.sh --list                              # implemented routes
./run_pipeline.sh ... --dry-run                       # print the sbatch commands
```

Implemented routes: `scan2 -> hiscanner`, `ugp -> chisel`. Asking for any other
pair prints why it isn't wired up yet. Starting partway checks that the earlier
stage's output exists (pass `--after JOBID` if it's still being made).

### Outputs (under PROJECT_DIR)

| stage | where |
|---|---|
| simulation | `sim/{bams,bulk,cn_mat,dwgsim}` |
| phase | `scan2_out/` or `ugp/out/phase/phased_het_snps.vcf.gz` |
| call | `<caller>_<binsize bp>/<phaser>/`, e.g. `chisel_5000000/ugp/`, `hiscanner_500000/scan2/` (barcoded BAM shared in `chisel/prep/`) |
| eval | `eval/<phaser>_<caller>_<binsize bp>/`, e.g. `eval/ugp_chisel_5000000/` (lifted truth shared: `eval/truth_hg38.tsv`) |

Bin size comes from `CHISEL_BINSIZE` / `HS_BINSIZE` in `config.sh`, or `--binsize` per run
(`./run_pipeline.sh --phaser ugp --caller chisel --from call --binsize 1Mb`). It is pinned
when you submit, so editing `config.sh` afterwards doesn't affect queued jobs.

Eval outputs: `summary.txt` (accuracy + gain/normal/loss TP/TN/FP/FN), `confusion_segments.pdf`
(+ `section_confusion.tsv`), `per_cell_tracks{,_bin,_ideal}.pdf`, and `calls.tsv` (the
converted calls that were scored).

Every caller's output is converted by `standardize_calls.py` to one table
(`calls.tsv` in the eval folder: cell, chrom, start, end, CN_A, CN_B + extras), so
`eval_calls.py` is the same for every caller. Adding a caller = a `call_<x>.sh`,
a converter in `standardize_calls.py`, and an entry in `ROUTES` in `scripts/common.sh`.

### Setup notes

- Conda envs: `scdna_pipeline` (SCAN2, HiScanner, dwgsim, bcftools, samtools),
  `genotyping-env` (UGP), `chisel` (Python 2.7).
- SCAN2 bug fix: replace `$ENV/lib/scan2/snakefile.phasing` with `scripts/snakefile.phasing`.
- CHISEL needs `Homo_sapiens_assembly38.dict` next to the FASTA (`samtools dict`).
