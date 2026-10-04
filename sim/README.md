## Pipeline for testing simulated chr1 reads for CHISEL, HiScanner, etc

Make one conda env for hiscanner, scan2, dwgsim. CHISEL requires another conda env since the python version is different.s


Small bug fix for SCAN2: replace the existing ```$ENV/lib/scan2/snakefile.phasing``` file with the ```snakefile.phasing``` file in this dir

```config.sh```: all configurations go here

```run_full_sim.sh```: the pipeline from simulation (step 1) to scan2 (step 2) to hiscanner (step3)
Run this with a start and end (```sbatch run_full_sim.sh [START] [END]```, i.e. ```sbatch run_full_sim.sh 1 2```),

```./scripts/ugp_run.sh```: runs Universal Genotyper (replacing step 2)

```./scripts/chisel_prep_run.sh```: preps BAMs for CHISEL