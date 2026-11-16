# Simulated scDNA benchmark

Simulates low-coverage single cells (about 0.1x) and a 30x bulk from the HG002 chr1 haplotypes along a clone tree with planted copy-number events, then phases SNPs, calls copy numbers and scores the calls against the truth. Every stage is a SLURM job; the launchers only submit.

## Setup

- Conda envs, set in `config.sh`: `scdna_pipeline` (simulation, SCAN2, HiScanner, evaluation), `genotyping-env` (UGP) and the CHISEL envs (Python 2.7; `HMM_CHISEL_ENV` is the [chisel_hmm](https://github.com/kevinyu2/chisel_hmm) fork).
- One config per tree: copy `config.sh`, change `PROJECT_DIR` and `TREE_FILE`, and pass it as `CONFIG_FILE=$PWD/config_<tree>.sh` to every launcher.

## Running

Run from `sim/` on the login node. Each launcher prints its last job ID, for `--after`.

```bash
./run_sim.sh                                                          # simulate
./run_pipeline.sh --phaser ugp --caller chisel --phaserbins 8          # phase -> call -> eval
./run_pipeline.sh --phaser ugp --caller chisel --phaserbins 8 --from call --hmm 24-48 --after <job>
./run_pipeline.sh --phaser ugp --caller chisel --phaserbins 8 --from call --cellhmm 6
EVAL_DETECTION_ONLY=true ./run_pipeline.sh --phaser ugp --caller chisel --phaserbins 8 --only eval
```

| option | what it does |
|---|---|
| `--phaser scan2\|ugp\|perfect` | phasing: SCAN2, UGP, or the true phase |
| `--caller chisel\|hiscanner` | copy-number caller |
| `--from`, `--to`, `--only` | stages: `phase`, `call`, `eval` |
| `--phaserbins N` | the phaser's variable-length bins, N SNP-covering reads per cell per bin |
| `--hmm K` or `--hmm A-B` | re-call with HMM bin clustering, K fixed or chosen in A..B |
| `--cloner` | with `--hmm`: also group cells into clones |
| `--oracle` | re-call with the true bin clusters |
| `--cellhmm K` | re-call with the per-cell HMM |
| `--bincellhmm K` | re-call with RunpengLuo's bin-cell HMM (needs `BINCELLHMM_REPO`) |
| `--after JOBS` | wait for these jobs |
| `--dry-run` | print the jobs without submitting |

The re-call options reuse the `combo.tsv` of the plain run (same command without them), so run that first and pass its job ID with `--after`. Eval settings: `EVAL_DETECTION_ONLY=true` (detection heatmap only, faster), `EVAL_CN_OVERRIDE=raw` (score per-cell calls of a cloned run), `EVAL_SUFFIX=_raw` (eval folder suffix).

## Algorithms

- **Simulation**: every cell's haplotypes are edited along the tree, then sequenced at about 0.1x with allele dropout, which the truth does not include.
- **Phasers**: SCAN2 and UGP phase het SNPs from the bulk (UGP can also make bins); `perfect` uses every true het SNP and its true phase.
- **Variable-length bins**: windows are merged until each bin has N SNP-covering reads per cell on average; smaller N gives shorter, noisier bins.
- **Plain CHISEL**: clusters bins with k-medians across all cells, calls copy numbers per cluster and cell, and corrects them by clone.
- **HMM**: clusters bins with an HMM along the genome whose variances come from read counts; a K range is chosen against a shifted null (BIC with the cloner).
- **Oracle**: the true bin clusters, an upper bound for the clustering step.
- **Cloner**: splits cells top-down on RDR and BAF, and calls every clone from its pooled reads.
- **Per-cell HMMs**: every cell gets its own state in every bin, the cells sharing K (RDR, BAF) states.
- **Scoring**: the truth is lifted to hg38 and onto the bins. Detection is the share of carrier cell-bins whose call shows each planted event (allele-specific), per subclone size x event size.

## Outputs

Under `PROJECT_DIR`: `calls/chisel_pbins<N><suffix>/<phaser>/` (calls) and `eval/<phaser>_chisel_pbins<N><suffix>/` (`summary.txt`, `detection.png`, per-cell plots). Pool detection over trees with `scripts/aggregate_detection.py`; plot a per-cell HMM run with `scripts/plot_cellhmm.py --run <calls folder>`. New sensitivity trees: `scripts/make_sensitivity_tree.py`.
