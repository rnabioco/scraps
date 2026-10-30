# SLURM profile

Submits every scraps rule as its own SLURM job through the Snakemake 9 SLURM
executor plugin (`snakemake-executor-plugin-slurm`, installed by
`scraps_conda.yml`). Snakemake itself runs on a login or head node (ideally in
`tmux`/`screen`) and keeps submitting and monitoring jobs.

## Setup

Edit `config.v8+.yaml`:

- `slurm_account`, `slurm_partition`: your allocation and partition.
- `runtime` (minutes): the default and the per-rule overrides under
  `set-resources` (STAR and cutadapt get longer walltimes).
- `resources: total_impact`: each STAR job claims 5, so 50 allows 10 concurrent
  STAR alignments. Lower it if your filesystem struggles.
- `jobs`: maximum number of SLURM jobs queued or running at once.

Each rule's `threads` becomes `--cpus-per-task` and its `mem_mb` becomes
`--mem`, so no per-rule CPU/memory settings are needed here.

STARsolo memory defaults to the genome index size (read from `STAR_INDEX`) plus
8 GB, at least 48 GB: STAR needs the index while mapping, frees it before the
STARsolo counting and BAM sort, and `--limitBAMsortRAM` is set to half the
job's `mem_mb`, the other half kept for STARsolo's per-read data, so the sort
stays inside the allocation. For large libraries
raise it with `STAR_mem_mb` in `config.yaml` (DEFAULTS, a platform, or one
sample) or with `mem_mb` under `set-resources` here. `retries: 1` resubmits a
failed job once, and STAR jobs get twice the memory on the retry. On clusters
that cap memory per core, a large `mem_mb` needs a high-memory partition or
more CPUs (both settable per rule under `set-resources`).

## Run

```bash
conda activate scraps_conda
snakemake -n --profile profiles/slurm --configfile config.yaml   # check first
snakemake --profile profiles/slurm --configfile config.yaml
```

## Alternative: one batch job

Without the profile, the whole pipeline runs locally inside a single
allocation, which needs no plugin:

```bash
#SBATCH --cpus-per-task=24 --mem=96G --time=2-00:00:00
snakemake -j 24 --resources total_impact=5 --keep-going --configfile config.yaml
```

See the plugin documentation for all options:
https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html
