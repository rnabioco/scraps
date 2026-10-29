# AGENTS.md - Developer Guide for scraps

This is the single guide for coding agents (Claude Code, Codex, opencode, ...) and
humans working in the scraps repository: commands, architecture, conventions and
gotchas. Claude Code reads this file directly; there is deliberately no separate
`CLAUDE.md`.

**Project**: scraps - Single Cell RNA PolyA Site Discovery  
**Type**: Snakemake bioinformatics pipeline that extracts mRNA polyadenylation (RT
priming) sites from "TVN"-primed 3' single-cell RNA-seq FASTQs, or from pre-aligned,
barcode-tagged BAMs, at near-nucleotide resolution  
**Primary Languages**: Python 3, R, Snakemake, Shell (zsh)

---

## Quick Start Commands

### Running the Pipeline

```bash
# Dry-run to validate pipeline (recommended before any changes)
snakemake -np --configfile config.yaml

# Run pipeline (--cores/-j is mandatory for a real run; needs the STAR index and
# whitelists from the README Setup section, which are not bundled; each STAR job
# uses total_impact=5, so --resources total_impact=5 runs one STAR job at a time)
snakemake -j 8 --configfile config.yaml --resources total_impact=5 --keep-going

# Generate DAG visualization (the Snakefile's logger output goes to stderr outside
# dry-runs, so stdout is clean DOT)
snakemake --dag --configfile config.yaml | dot -Tpdf > dag.pdf

# Regenerate the README's inst/dag.pdf: rule graph with every arm and discovery on.
# dag_override.yaml contains:
#   GENOME_FASTA: "ref/genome.fa"
#   DISCOVERY: {enabled: true, read: [R1, R2, paired, bam]}
snakemake --rulegraph --configfile config.yaml dag_override.yaml | dot -Tpdf > inst/dag.pdf
```

`--configfile other.yaml` is deep-merged over the hardcoded `config.yaml`
(Snakemake `update_config`), so `SAMPLES` from both files are unioned unless the
override redefines every key. Add `--replace-workflow-config` (Snakemake 9) to use
only the given file instead.

**Cluster (SLURM)**: submit one job per rule with the bundled profile, after setting
the account and partition (see `profiles/slurm/README.md`); `--cluster`/`--drmaa` no
longer exist in Snakemake 8+.

```bash
snakemake --profile profiles/slurm --configfile config.yaml
```

### Testing Changes

```bash
# Always dry-run first to validate Snakemake syntax
snakemake -np --configfile config.yaml

# Test specific rule (target BEFORE --configfile: that flag takes multiple files,
# so a trailing target would be read as another config file)
snakemake -np results/counts/chromiumv2_test_R2_counts.tsv.gz --configfile config.yaml

# List all rules
snakemake --list-rules

# Job reasons are always printed in Snakemake 9 (-r/--reason no longer exist)
```

There is no unit-test suite; CI (`.github/workflows/snakemake_run.yml`) only dry-runs.
Script-level regression checks live in `tests/`:

```bash
tests/validate_synthetic.sh          # fast: 6 cases (R2, R1, paired, bam_same, bam_reverse, bam_paired1); needs bash >= 4, umi_tools CLI, bedtools
tests/validate_kde_peaks.sh          # kde_peaks.py cluster extents: umi_support, unique summits, support floor
Rscript tests/test_geneid.R          # R helpers: delimiter-aware GeneID parsing on the bundled SAFs (needs R + tidyverse)
OLD_REF=<ref> tests/validate_realdata.sh   # slow: full pipeline twice (git worktree of OLD_REF vs working tree); default OLD_REF is the last umi_tools commit
```

Pass criteria the scripts enforce: `count_sites.py` output must be **bit-exact** with
`umi_tools count`; `pileup_sites.py` beds must be **directional** (every legacy position
still present, total signal <= legacy `umi_tools dedup` + `genomecov`), because dedup is
keyed on (cell, UMI, priming position) and ignores the fragmentation end. Keep these
invariants if you touch `count_sites.py`, `pileup_sites.py`, `filter_bam.py` or the
featureCounts flags, and re-run `validate_synthetic.sh`. The per-mode featureCounts and
pileup settings are also hardcoded in that script's `FC_*`/`BED_END`/`GC_FLAG` tables,
so a new mode must be added there too.

### Environment Setup

```bash
# Create conda environment
conda env create -f scraps_conda.yml

# Activate environment
conda activate scraps_conda

# Update environment after changes
conda env update -n scraps_conda -f scraps_conda.yml
```

---

## Project Structure

```
scraps/
├── Snakefile              # Main workflow entry point
├── config.yaml            # Sample and pipeline configuration
├── chemistry.yaml         # Platform-specific chemistry configs
├── scraps_conda.yml       # Conda environment specification
├── rules/                 # Snakemake rule modules
│   ├── cutadapt_star.snake   # Read trimming and alignment
│   ├── count.snake           # Feature counting and quantification
│   ├── qc.snake              # Quality control reports
│   ├── discovery.snake       # De novo RT priming site discovery (optional)
│   └── check_versions.snake  # Dependency version checks
├── inst/scripts/          # Helper scripts called by rules:
│   ├── count_sites.py     #   directional UMI counting per cell x site (replaces umi_tools count)
│   ├── pileup_sites.py    #   strand-aware single-base priming pileup (discovery)
│   ├── filter_bam.py      #   poly(A)/poly(T) soft-clip filter + mate selection
│   ├── cut_paste_fastq.py #   complex-barcode stitching for cutadapt_paired
│   ├── kde_peaks.py, annotate_sites.py   # de novo discovery
│   ├── polyadb/           # polyAdb -> scraps SAF converters (3.2 + v4), standalone
│   ├── R/                 # R analysis functions (downstream, not called by rules)
│   └── bam_check.py, filter_bam_correct*.py, check_longR1/, test/   # legacy / offline
├── tests/                 # validate_synthetic.sh (+ generate_synthetic.py), validate_realdata.sh
├── docs/                  # discovery.md; reports/ (change summaries)
├── ref/                   # polyadb32.{hg38,mm10} + polyadb4.{hg38,mm10}.{main,max} SAFs (converter output); STAR index + whitelists go here
├── sample_data/raw_data/  # Test FASTQs + synthetic_longread.bam (bam-mode fixture)
├── inst/test_output/, inst/dag.pdf, inst/notebooks/, man/figures/   # legacy examples, figures
├── profiles/slurm/        # SLURM executor profile (config.v8+.yaml) + README
├── .github/workflows/snakemake_run.yml   # CI: dry-runs (local + SLURM profile)
└── results/               # Pipeline outputs (generated)
```

---

## Architecture: one DAG arm per alignment mode

Every sample runs one or more of four **alignment modes** (`R1`, `R2`, `paired`, `bam`),
chosen by `alignments` in `config.yaml` and constrained by `valid_modes` in
`chemistry.yaml`. The mode is the `{read}` wildcard. Per-mode rules are suffixed
`_R1`/`_R2`/`_paired`/`_bam` where the body differs (`cutadapt_*` exists for R1/paired,
`starsolo_*` for R1/R2/paired, `filter_*` for R2/bam, `assign_sites_*` for all four);
`count` and `stranded_bed` are single `{read}`-wildcard rules. The strandedness and
priming-end settings must stay consistent across the three rule files, the Snakefile
helpers (`_fc_shift_type`, `_bam_*`) and the tables in `tests/validate_synthetic.sh`:

| mode | trimming | STAR input | featureCounts (`count.snake`) | priming end / strand flip (`discovery.snake`) |
|---|---|---|---|---|
| R2 | none (raw FASTQs) | `starsolo_R2`, then `filter_R2` keeps reads with a soft-clipped A tail (`filter_bam.py`; skipped when `R2_polyA_filter: false`) | `-s 1 --read2pos 3` (+ `--readShiftType downstream --readShiftSize N` when `R2_priming_offset: N`) | 3' end (+ `--shift N`), read strand == transcript strand |
| R1 | `cutadapt_R1` | `starsolo_R1`, barcode mate 1, `STAR_R1` = one `--clip5pNbases` value | `-s 2 --read2pos 5` | 5' end, strand flipped |
| paired | `cutadapt_paired` (plus optional `bc_cut` barcode stitching via `cut_paste_fastq.py`) | `starsolo_paired`, `STAR_paired` = two clip values | `-s 2 -p --read2pos 5` | 5' end, strand flipped |
| bam | none | no STAR: user BAM from the sample's `bam:` key; `filter_bam` drops 0x904 + untagged reads, keeps/unpairs the `bam_mate` mate of proper pairs, optional poly(A)/poly(T) clip filter (`filter_bam.py --end`) | `-L` + `-s 1 --read2pos 3` (`bam_strand: same`) or `-s 2 --read2pos 5` (`reverse`), + shift when `bam_priming_offset: N`; single-threaded | 3' end, identity (`same`) or 5' end, flipped (`reverse`), + `--shift N` |

Per (sample, mode) flow, with `results/` as `RESULTS`:

1. `rules/cutadapt_star.snake`: FASTQs (globbed from `DATA` by `basename*R1*`) ->
   `{sample}/{sample}_{mode}_Aligned.sortedByCoord.out.bam` (STARsolo, CB/UB tags; no cell
   filtering). Not produced for the `bam` mode.
2. `rules/count.snake`: filter to reads with CB and UB (`filter_R2` ->
   `temp/{sample}_R2.bam`; `filter_bam` -> `temp/{sample}_bam.bam`; R1/paired filter inline
   in `assign_sites_*`) -> `featureCounts -R BAM` against `POLYA_SITES` ->
   `temp/{sample}_{mode}_assigned.bam` (XT/XS tags; `temp()`, shared by `count` and
   `stranded_bed`) -> `count_sites.py` (directional UMI collapse per cell x site, parallel
   over contigs) -> `counts/{sample}_{mode}_counts.tsv.gz`. This is the pipeline's primary
   output.
3. `rules/discovery.snake` (only when `DISCOVERY.enabled`): the same assigned BAM ->
   `pileup_sites.py` (strand-aware single-base pileup, kept as
   `discovery/beds/{sample}_{mode}.stranded.bed.gz`) -> `pool_beds` over a **unit** (a set
   of (sample, mode) pairs computed in the Snakefile as `DISC_UNITS`: per sample, per
   `DISCOVERY.groups`, or explicit `DISCOVERY.units`) -> `kde_peaks.py` ->
   `annotate_sites.py` -> `discovery/{unit}_{sites.tsv,sites.bed,denovo.saf}.gz`.
   Details in `docs/discovery.md`.

`rule all` is assembled in the Snakefile from `SAMPLE_OUTS` (+ `DISCOVERY_OUTS`). The
MultiQC report (`rules/qc.snake`) is commented out of `all_outputs`, and `qc_check` expands
over all three FASTQ modes regardless of config, so it only resolves when every sample
ran every mode.

Config lookup for any per-sample item is `_get_config(sample, item)`: sample block ->
`chemistry.yaml[chemistry][platform]` -> `chemistry.yaml[chemistry]` -> `DEFAULTS`, falling
through to `""` (with a one-time `logger.warning` unless the key is in the optional set).
Rules rely on that empty-string fallback (e.g. `cutadapt_*` skips trimming when
`cut_params` is empty), so adding a new chemistry key means checking which rules treat
"missing" as "off". The mode switches `R2_polyA_filter`, `R2_priming_offset` and the
`bam_*` keys go through `_polya_filter` / `_priming_offset` / `_fc_shift_args` /
`_bam_strand` / `_bam_end` / `_bam_mate` (Snakefile; validated at parse time by
`_validate_bam_inputs` / `_validate_offsets`) and `_disc_end` / `_disc_labels`
(`rules/discovery.snake`); see "Configuration Files" below for the key semantics. A
priming offset always moves the counted end *outward* past the read end
(`--readShiftType downstream` for 3'-end modes, `upstream` for 5'-end modes; the same
convention in `pileup_sites.py --shift`), and only `R2_`/`bam_priming_offset` are honoured
(`_OFFSET_MODES`).

---

## Code Style Guidelines

### Python Scripts

**Imports**: Standard library block, blank line, third-party block
```python
import os
import gzip
import argparse
from multiprocessing import Pool

import pysam
```

**Docstrings**: Triple-quoted module docstring as the FIRST statement (as in
`annotate_sites.py`, `kde_peaks.py`, `polyadb/*.py`); several older scripts
(`count_sites.py`, `pileup_sites.py`, `filter_bam.py`, `cut_paste_fastq.py`) put it
after the imports, where it is not a real `__doc__`. Follow the newer form.
```python
"""Call RT priming site summits from a strand-aware priming pileup.
...
"""
```

**Command-line arguments**: Use `argparse` with descriptive help text
```python
parser.add_argument('-i', '--inbam', required=True,
                    help="BAM file to filter")
```

**Naming conventions**:
- Functions: `snake_case` (e.g., `filter_bam_by_A`, `_end_position`, `_pileup_contig`)
- Variables: `snake_case` (e.g., `clip_filter`, `sort_mem`, `fragment_path`)
- Constants: `UPPER_CASE` if truly constant

**File handling**: Use context managers for file operations
```python
with open(file_in) as file, gzip.open(file_out, 'wt') as file2:
    # process files
```

### R Scripts

**Documentation**: Roxygen2-style comments for functions
```r
#' Read scraps output from umi_tools to sparseMatrix
#' 
#' @param file scraps output table
#' @param n_min minimum number of observations
#' @return count matrix
#' @export
```

**Style**: Follow tidyverse conventions
- Use `%>%` pipe operator
- Prefer `dplyr`, `readr`, `stringr`, `tidyr` functions
- Function names: `snake_case`

**Dependencies**: Import packages explicitly
```r
#' @import readr dplyr stringr tidyr
```

### Snakemake Rules

**Shell executable**: Pipeline uses `zsh` (defined in Snakefile line 1)
```python
shell.executable("zsh")
```

**Rule structure**: Include all standard sections
```python
rule rulename:
    input:
        "path/to/input.bam"
    output:
        temp("path/to/output.bam")  # Use temp() for intermediate files
    params:
        job_name = "rulename",
        # Additional parameters
    log:
        "{results}/logs/rulename_{sample}_{read}.txt"
    threads:
        12
    resources:
        mem_mb = 8000
    shell:
        r"""
        exec > {log} 2>&1
        command --arg {input} > {output}
        """
```

**Key conventions**:
- Prefer raw strings `r"""..."""` for shell blocks (count/discovery rules do;
  `cutadapt_star.snake`, `check_versions.snake` and `qc.snake` still use plain `"""`)
- **Log naming**: `{results}/logs/{rule-stem}_{sample}_{mode}.txt`, where
  `{rule-stem}` is the rule name without its `_R1`/`_R2`/`_paired`/`_bam` suffix
  (`filter_R2` -> `filter_{sample}_R2.txt`, `assign_sites_bam` ->
  `assign_sites_{sample}_bam.txt`, `starsolo_paired` -> `starsolo_{sample}_paired.txt`;
  wildcard rules such as `count` and `stranded_bed` use `{read}`), always `.txt`.
  Discovery rules that operate on a unit rather than a single (sample, mode) use
  `{rulename}_{unit}.txt` (the unit id already encodes the mode). Aggregates with
  no per-sample wildcard use a bare descriptive name: `qc_check` -> `multiqc.txt`,
  `check` -> `{data}/logs/check_versions.txt`. Keep this consistent so logs sort by
  rule and are trivially associable with their origin.
- **Capture ALL shell output to the log**: put `exec > {log} 2>&1` as the first
  line of every shell body (preferred over appending `2> {log}` to each command).
  This routes stdout *and* stderr — including bare `echo` status messages like
  cutadapt's `"no trimming"` — into the per-rule log instead of leaking into the
  Snakemake console/master log where they are hard to attribute. Every
  shell-bearing rule, including the discovery rules, does this.
- Mark intermediate files with `temp()`
- Use wildcards in paths: `{sample}`, `{results}`, `{read}`
- Resource specifications: `threads`, `mem_mb`
- Per-(sample, mode) targets for `rule all` are appended inside the
  `for mode in _sample_modes(x)` loop that builds `SAMPLE_OUTS` in the Snakefile;
  `expand()` is used for discovery (`DISCOVERY_OUTS`) and `qc_check`

**Accessing config**: Use helper functions like `_get_config(sample, item)`
```python
def _get_config(sample, item):
    # Hierarchical lookup: sample -> chemistry[platform] -> chemistry -> defaults
```
Diagnostics from the Snakefile (the `_get_config` fallthrough, the parse-time
validators' warnings, the target list) go through **Snakemake's logger**
(`from snakemake.logging import logger`; `logger.warning` / `logger.info`), never
`print()` or `sys.stderr.write`. The logger follows Snakemake's stream choice (stdout in
a dry-run, in order with the job listing; stderr for real runs and for `--dag`, so DOT
output stays clean) and lands in `.snakemake/log/` (verified on Snakemake 9). Two rules:
build the strings with `str.format()` for consistency (Snakemake 7's tokenizer re-spaced
f-string contents; fixed in 8+), and emit each message once via `_diag`
(params lambdas run once per job). `_get_config` stays silent for keys that are
legitimately optional (`_OPTIONAL_KEYS` / `_OPTIONAL_SUFFIXES`: `cutadapt_*`, `bc_cut`,
`bc_length1`, `*_extra_args`, `*_polyA_filter`, `*_priming_offset`, `bam_*`); any other
missing key is a real misconfiguration and is reported once.

---

## Configuration Files

### config.yaml
- `DATA`: Directory containing input FASTQs (and relative `bam:` paths)
- `RESULTS`: Output directory path
- `STAR`: STAR executable (used as `{STAR}` in the starsolo rules)
- `STAR_INDEX`: Path to STAR genome index
- `POLYA_SITES`: PolyA database reference file (SAF format; see SAF spec below; a
  non-`.saf`/`.saf.gz` file is passed to featureCounts as GTF)
- `GENOME_FASTA`: Genome FASTA (with `.fai`); required only when `DISCOVERY.enabled` is true
- `DISCOVERY`: Optional de novo RT priming site discovery block (see `docs/discovery.md`)
- `DEFAULTS`: Lowest-priority values for any per-sample key: `platform`, `chemistry`,
  `alignments`, and the mode switches `R2_polyA_filter`, `R2_priming_offset`,
  `bam_mate`, `bam_strand`, `bam_polyA_filter`, `bam_priming_offset`
- `SAMPLES`: Per-sample configuration (`basename` or `bam`, `platform`, `chemistry`,
  `alignments`, `discovery_read`, the `bam_*`/`R2_*` switches, and any
  chemistry/platform key such as `STAR_R2` or `STAR_R2_extra_args` as an override)
- `report_section_order`: MultiQC layout, used only by the disabled `qc_check` rule

Config errors (unknown chemistry/platform: `_validate_config`; unsupported
alignment mode: `_validate_modes`; bam-mode keys: `_validate_bam_inputs`;
`DISCOVERY.enabled` without `GENOME_FASTA`) are raised at Snakefile parse time, so
`snakemake -n` catches them.

### chemistry.yaml
Platform-specific configurations organized hierarchically:
```yaml
chemistry_name:
  bc_whitelist: path/to/whitelist
  platform_name:
    valid_modes: [R1, R2, paired]       # alignment modes this platform supports
    cutadapt_R1: "trimming parameters"  # omit a cutadapt_* key to skip trimming
    cutadapt_paired: "trimming parameters"
    STAR_R1: "alignment parameters"     # single-mate R1 rule: --clip5pNbases = 1 value
    STAR_paired: "alignment parameters" # two-mate paired rule: --clip5pNbases = 2 values (<R1clip> 0)
    STAR_R2: "alignment parameters"
    R2_polyA_filter: true               # optional; false skips filter_bam.py in filter_R2
    R2_priming_offset: 0                # optional; bp shift of the R2 3' end toward transcript 3'
  generic:                              # pre-aligned BAM input, no STAR/cutadapt keys
    valid_modes: [bam]
```

`bc_whitelist`, `bc_cut` and `bc_length1` are chemistry-level scalars and must stay in
the `meta_keys` set of `_validate_config` so they are not mistaken for platforms.
`bc_cut` is honoured by `cutadapt_paired` only. STAR keys reach STAR as
`--soloCBwhitelist <bc_whitelist> <STAR_{mode}> <STAR_{mode}_extra_args>`
(`_get_chem_version` in `rules/cutadapt_star.snake`).

**R2-arm switches** (`R2_polyA_filter`, `R2_priming_offset`): resolved per sample by
`_polya_filter(sample, mode, default)` / `_priming_offset(sample, mode)` /
`_fc_shift_args(sample, mode)` in the `Snakefile`, with defaults in `config.yaml`
`DEFAULTS`. An empty `_get_config` fallthrough (configs that predate the keys) means
"filter on, offset 0". A non-zero offset renders
`--readShiftType downstream --readShiftSize N` in `assign_sites_R2` and `--shift N` in
discovery's `stranded_bed` (`pileup_sites.py`), so the count table and the discovery bed
agree on the site coordinate. The Ultima platforms use them: `ultima_cram` (CRAM-derived
single long read, R1 mode, no cutadapt) vs `ultima_fastq` (delivered pseudo-paired FASTQs,
R2 mode, `R2_polyA_filter: false`, `R2_priming_offset: 5`; reads must NOT be required
to carry a poly(A) anywhere: Ultima trims the homopolymer, `filter_R2` then keeps every
CB/UB-tagged read, and `starsolo_R2`'s `--clipAdapterType CellRanger4` was checked on
simulated Ultima reads and leaves their 3' ends exact, including reads ending just before
or a few bases into a genomic A-run; it only clips reads ending in a ~30-A tract, which
Ultima trimming does not produce). The offset key is per mode
(`{read}_priming_offset`), but only R2 and bam honour it (`_OFFSET_MODES`): the R1/paired
count rules apply no shift, so `_priming_offset` ignores (with a warning) any
`R1_`/`paired_priming_offset`, keeping discovery beds and count tables in sync.

**`bam` mode / `generic` platform**: pre-aligned, CB/UB-tagged BAM input with no STAR
step (`rules/count.snake` `filter_bam` -> `assign_sites_bam` -> `count`). The platform
block is just `valid_modes: [bam]`; the sample supplies `bam:` (path) and the switches
`bam_strand` (`same` -> `-s 1 --read2pos 3`, 3' end, like R2; `reverse` -> `-s 2
--read2pos 5`, 5' end, like R1/paired), `bam_polyA_filter` (poly(A)/poly(T) soft clip at
the counted end via `filter_bam.py --end`), `bam_mate` (0 single-end; 1/2 = the mate of a
paired BAM carrying the priming site, which is kept from proper pairs and unpaired so
positions are per fragment) and `bam_priming_offset`. Helpers: `_bam_strand`, `_bam_end`,
`_bam_mate`, `_fc_shift_type` (Snakefile), `_disc_end` / `_disc_labels`
(`rules/discovery.snake`). featureCounts runs with `-L` (long-read mode, single thread).
Assumptions (documented, not enforced): full-length molecules, CB/UB tags, strand
accurate for the declared `bam_strand`.

**clip5pNbases mate-count rule**: STAR requires one `--clip5pNbases` value per
mate. `STAR_R1` feeds `starsolo_R1` (one `--readFilesIn` mate) so it must carry
exactly one value; `STAR_paired` feeds `starsolo_paired` (two mates) so it
carries two (`<R1clip> 0`, R2/cDNA unclipped). These are distinct keys — do NOT
reuse a two-value `STAR_R1` for both, or R1-only alignment fails with
"--clip5pNbases has to contain 1 values to match the number of mates".
`STAR_paired` is fully self-contained (clip/solo args plus the paired-only
`--alignEndsProtrude 58 ConcordantPair`); it is no longer split across a
separate DEFAULTS key.

### SAF format (POLYA_SITES)

Tab-separated, header `GeneID  Chr  Start  End  Strand`. The **GeneID has 7
fields**:

```
gene_symbol{D}refseq_gene_id{D}ensembl_id{D}chrom{D}pos{D}strand{D}pas_type
```

- Delimiter `{D}`: `;` human, `_` mouse (matches `ref/polyadb32.{hg38,mm10}.saf.gz`;
  the de novo SAF written by `annotate_sites.py` is always `;`-joined, also for mouse;
  region SAFs built by `inst/scripts/R/expand_region_saf.R` reuse the input's delimiter
  with `pos`=`NA` and are not usable as a discovery `--polya` reference).
- `pos` = field 5 (1-based cleavage position); `pas_type` = field 7 (class).
  Missing values are `NA`. `pos` may occasionally appear in float/scientific
  notation (e.g. `7.7e+07` from a float round-trip, as in one legacy row of the
  bundled `ref/polyadb32.mm10.saf.gz`, since repaired). Parsers must coerce it
  with `int(round(float(pos)))` rather than `int(pos)` so such rows are recovered
  instead of silently skipped — both `annotate_sites.py` (`load_polya`) and
  `inst/scripts/polyadb/validate_saf.py` do this. The converters in
  `inst/scripts/polyadb/` always emit plain integers.
- `Start`/`End` = strand-specific 16-base window (`End - Start = 15`; transcript
  -10/+5): `+` → `[pos-10, pos+5]`, `-` → `[pos-5, pos+10]`. This is what the
  converters (`saf_common.site_window`) emit and what the bundled
  `ref/polyadb32.{hg38,mm10}.saf.gz` references use; `validate_saf.py` enforces it.
- The `_` (mouse) delimiter appears inside RefSeq IDs and scaffold chrom names;
  parse by anchoring on the SAF `Chr`/`Strand` columns and reading trailing
  fields from the right (see `annotate_sites.py` `_split_geneid`).
- PAS type vocab is release-specific and passed through verbatim (3.2 vs v4
  differ, e.g. `intergenic` vs `Intergenic`).

---

## Common Development Tasks

### Adding a New Rule

1. Create rule in appropriate file under `rules/`
2. Follow naming convention: `{step}_{mode}` for per-mode variants (`filter_R2`,
   `assign_sites_bam`), or one `{read}`-wildcard rule when the body is mode-agnostic
   (`count`, `stranded_bed`)
3. Add to workflow: per-(sample, mode) targets go inside the `for mode in
   _sample_modes(x)` loop that builds `SAMPLE_OUTS` (Snakefile; note the `bam` mode
   yields no `*_Aligned.sortedByCoord.out.bam`); discovery targets go in `DISCOVERY_OUTS`
4. Test with dry-run: `snakemake -np`

### Modifying Chemistry Configuration

1. Edit `chemistry.yaml`
2. Ensure the `STAR_*` field for every mode in `valid_modes` is present (`cutadapt_*` is
   optional: omitting it skips trimming). The `bam` mode needs neither: its platform block
   is just `valid_modes: [bam]`
3. Optional fields: `bc_whitelist`, `bc_cut`, `bc_length1` (chemistry level),
   `STAR_{R1,R2,paired}_extra_args`, `R2_polyA_filter`, `R2_priming_offset`, and the
   `bam_*` switches (platform or sample level)
4. Test with dry-run to validate YAML syntax

### Adding Python Helper Script

1. Place in `inst/scripts/`
2. Use argparse for CLI interface
3. Include a module docstring explaining purpose (first statement)
4. Call from Snakemake rule with `python3 inst/scripts/script.py` (no shebang or
   executable bit needed)

### De novo RT Priming Site Discovery

Optional module (off by default); see `docs/discovery.md` for full details.

- Rules: `rules/discovery.snake` (`stranded_bed`, `pool_beds`, `kde_peaks`,
  `annotate_sites`). Wired into `all_outputs` in `Snakefile` only when
  `DISCOVERY.enabled` is true; requires `GENOME_FASTA`.
- Scripts: `inst/scripts/pileup_sites.py` (UMI-deduplicated strand-aware
  single-base pileup, rule `stranded_bed`), `inst/scripts/kde_peaks.py` (KDE peak
  calling), `inst/scripts/annotate_sites.py` (polyAdb match + sequence-context
  classification + de novo SAF).
- Intermediates: `discovery/beds/{sample}_{mode}.stranded.bed.gz` is kept (5 columns:
  chrom, start, end, UMI count, transcript strand); `discovery/beds/{unit}.pooled.bed.gz`
  and `discovery/{unit}_peaks.tsv.gz` are `temp()`.
- Annotation classes/SAF `class` values: `known_pas`,
  `suspected_internal_priming` (downstream A-rich: internal priming or a genuine
  site in an A-tract, not resolvable from sequence; not to be dropped by
  default), `potential_novel_pas`, `potential_novel_no_pas` (no A-rich
  signature, no PAS; lower confidence). `inst/scripts/R/discovery.R` reads the
  legacy `likely_internal_priming` label as `suspected_internal_priming`. The de novo SAF reuses the
  standard `gene;genbank;id;chrom;pos;strand;class` GeneID encoding so it can be
  set as `POLYA_SITES` to re-quantify (filter artifacts first).
- Outputs under `{results}/discovery/`: `<unit>_sites.tsv.gz`,
  `<unit>_sites.bed.gz`, `<unit>_denovo.saf.gz` (unit = group name or sample).
- Validate with `snakemake -np` after setting `DISCOVERY.enabled: true`.

### Generating/converting polyAdb references

Standalone converters in `inst/scripts/polyadb/` (stdlib only + UCSC liftOver,
auto-downloaded). Not part of the Snakemake DAG.

- `pas32_to_saf.py`: polyAdb 3.2 `*.PAS.txt` (both tables need liftOver: human
  hg19→hg38, mouse mm9→mm10 despite the "mm10" label in the file). Resolves
  columns by header name.
- `pasv4_to_saf.py`: polyAdb 4 `*.PAS.{main,max}.tsv` (human/mouse; main/max
  auto-detected; optional liftOver for future builds).
- `saf_common.py`: shared 7-field GeneID encoding, window, dedup/sort, liftOver.
- `validate_saf.py`: format checks + optional content comparison to a reference.
- `--species` sets the delimiter (`;` human / `_` mouse; `--delim` overrides).
  `--source-build`/`--target-build` are required on both converters (liftOver is
  skipped when equal). Output `.gz` writes gzipped, droppable straight into `ref/`.
- The auto-downloaded liftOver binary is the macOS arm64 build; on Linux pass
  `--liftover /path/to/liftOver` (e.g. from conda `ucsc-liftover`).
- `validate_saf.py` prints PASS/FAIL but always exits 0.
- PAS type passed through verbatim (3.2 vs v4 class vocab differs).
- Verify: run on a head-slice, then `validate_saf.py generated.saf ref/...saf.gz`
  (human 3.2 should match the bundled reference 100% at shared positions).

---

## Error Handling and Debugging

**Log files**: All rules write logs to `{results}/logs/` (the version-check rule to
`{data}/logs/`), named `{rule-stem}_{sample}_{mode}.txt` (or `{rulename}_{unit}.txt`
for discovery, `multiqc.txt` / `check_versions.txt` for aggregates; see Log naming
above).
- Check logs for detailed error messages
- Logs capture BOTH stdout and stderr (`exec > {log} 2>&1`), so status echoes
  (e.g. cutadapt's `"no trimming"`) and tool output land in the per-rule log,
  not the Snakemake console/master log.

**Common issues**:
- Missing conda dependencies → check `scraps_conda.yml`
- YAML syntax errors → validate with `snakemake -np`
- Missing input files → check `DATA` path in config.yaml
- Resource exhaustion → adjust `mem_mb` or `threads` in rules

**Debugging Snakemake**:
```bash
# Show detailed execution plan
snakemake -np --verbose

# Print shell commands without execution
snakemake -np --printshellcmds

# Force re-run specific rule
snakemake -j 8 --forcerun rulename   # preview first with -n
```

---

## Dependencies and Tools

**Core requirements** (installed via conda):
- Python >= 3.11 (required by Snakemake 9)
- Snakemake >= 9, < 10 (`min_version("9.0")` in the Snakefile; Snakemake 8 removed
  `-r/--reason`, `--cluster` and `--drmaa`)
- snakemake-executor-plugin-slurm (only for `profiles/slurm/`)
- STAR >= 2.7.9a (RNA-seq aligner)
- cutadapt >= 3.4 (adapter trimming; zsh `=(...)` process substitution, so zsh is required)
- samtools >= 1.15 (BAM manipulation; >= 1.16 for `[CB]!="-"` to drop untagged reads)
- subread >= 2.0.1 (featureCounts; `--readShiftType` and `-L` are used)
- pysam >= 0.16.0 (Python BAM interface) and UMI-tools >= 1.1.2 (used as a Python
  library by `count_sites.py`; the CLI is only used by `tests/validate_synthetic.sh`)
- numpy (`kde_peaks.py`) and pandas (imported at parse time by `rules/cutadapt_star.snake`);
  scipy and scikit-learn are in `scraps_conda.yml` but currently unused
- MultiQC >= 1.6 (only for the disabled `qc_check` rule)
- bedtools >= 2.30.0 (no longer called by any rule; `tests/validate_synthetic.sh` only)
- zsh

**Version checking**: `snakemake -j 1 results/check_versions.txt --configfile config.yaml`
(rule `check`, not part of `rule all`; checks python, snakemake, pysam/umi_tools import,
cutadapt, samtools, featureCounts, STAR, multiqc)

---

## Notes for coding agents

- **Always dry-run first**: Use `snakemake -np` before any pipeline changes
- **Respect shell choice**: Pipeline explicitly uses `zsh`, not bash (the cutadapt rules
  use zsh `=(...)` process substitution)
- **Preserve temp files**: Snakemake manages cleanup via `temp()` directive
- **Follow hierarchical config**: sample → chemistry.platform → chemistry → DEFAULTS → `""`
- **Log everything**: capture stdout AND stderr with `exec > {log} 2>&1` as the first
  shell line
- **Resource awareness**: Bioinformatics tools are memory/CPU intensive
- **Tests**: see "Testing Changes" above; re-run `tests/validate_synthetic.sh` after
  touching `count_sites.py`, `pileup_sites.py`, `filter_bam.py` or the featureCounts flags.
- **Legacy / offline code, not used by any rule**: `inst/scripts/*_old.py`,
  `bam_check.py`, `filter_bam_correct*.py`, `inst/scripts/check_longR1/`,
  `inst/scripts/test/`, `inst/test_output/` (outputs of the removed bed arm) and
  `inst/dag.pdf` (predates the current DAG). The R helpers in `inst/scripts/R/` are
  downstream analysis code, also not called by rules.
- **Untracked local artifacts**: the working tree may carry large run artifacts
  (`AEG1_*` multi-GB BAM/bed, `MousePas/`, `results/`, `old_results/`, `new_results/`,
  `ref/genome.fa*`, `*.out` job logs, `config41.yaml`, `chemistry41.yaml`, `.snakemake/`,
  `__pycache__/`, `.validate_realdata/` from the real-data test). They are not covered by
  `.gitignore`; never `git add -A`, stage explicit paths.
- **Fixture to keep tracked**: `sample_data/raw_data/synthetic_longread.bam` (+ `.bai`) is
  referenced by the bam-mode samples in `config.yaml`; without it the dry-run fails
- **Git**: default branch for PRs is `master`. Background for the pysam/discovery work:
  `docs/reports/pipeline_changes_since_94e9782.md` (tracked) summarizes every pipeline
  change since the discovery module landed; `.opencode/plans/*.md` (local, untracked) hold
  the design notes behind the bed-generation, shared-dedup-BAM and confident-peaks changes.
