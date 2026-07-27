# AGENTS.md - Developer Guide for scraps

This guide provides conventions and commands for AI coding agents working in the scraps repository.

**Project**: scraps - Single Cell RNA PolyA Site Discovery  
**Type**: Snakemake bioinformatics pipeline for analyzing mRNA polyadenylation sites from single-cell RNA-seq data  
**Primary Languages**: Python 3, R, Snakemake, Shell (zsh)

---

## Quick Start Commands

### Running the Pipeline

```bash
# Dry-run to validate pipeline (recommended before any changes)
snakemake -npr --configfile config.yaml

# Run pipeline with test data
snakemake --snakefile Snakefile \
  --configfile config.yaml \
  --resources total_impact=5 \
  --keep-going

# Run with specific number of cores
snakemake -j 8 --configfile config.yaml

# Generate DAG visualization
snakemake --dag | dot -Tpdf > dag.pdf
```

### Testing Changes

```bash
# Always dry-run first to validate Snakemake syntax
snakemake -npr --configfile config.yaml

# Test specific rule
snakemake -npr --configfile config.yaml results/counts/chromiumv2_test_R2_counts.tsv.gz

# List all rules
snakemake --list

# Show reason for rule execution
snakemake -npr --reason --configfile config.yaml
```

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
├── inst/scripts/          # Helper scripts
│   ├── *.py               # Python utilities (BAM filtering, KDE, annotation)
│   ├── polyadb/           # polyAdb -> scraps SAF converters (3.2 + v4)
│   └── R/                 # R analysis functions
├── docs/                  # Long-form docs (discovery.md, etc.)
├── ref/                   # Reference files (polyA_DB, etc.)
├── sample_data/           # Test data location
└── results/               # Pipeline outputs (generated)
```

---

## Code Style Guidelines

### Python Scripts

**Imports**: Standard library → Third party → Local, grouped and sorted
```python
import os
import re
import argparse

import pysam
import pandas as pd
import numpy as np
```

**Docstrings**: Triple-quoted strings describing script/function purpose
```python
""" Filter BAM files to only reads with soft-clipped A tail,
suitable for cellranger and starsolo output
"""
```

**Command-line arguments**: Use `argparse` with descriptive help text
```python
parser.add_argument('-i', '--inbam',
                    help="Bam file to correct",
                    required=True)
```

**Naming conventions**:
- Functions: `snake_case` (e.g., `filter_bam_by_A`, `correct_bam_read1`)
- Variables: `snake_case` (e.g., `target_len`, `filter_cut`, `single_end`)
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
- Use raw strings `r"""..."""` for shell blocks
- **Log naming**: `{results}/logs/{rulename}_{sample}_{read}.txt`, i.e.
  `{rulename}` first, then sample, then the alignment mode (`R1`/`R2`/`paired`),
  and always the `.txt` extension. Discovery rules that operate on a unit rather
  than a single (sample, mode) use `{rulename}_{unit}.txt` (the unit id already
  encodes the mode). Aggregate rules with no per-sample wildcard use a bare
  descriptive name (`multiqc.txt`, `check_versions.txt`). Keep this consistent so
  logs sort by rule and are trivially associable with their origin.
- **Capture ALL shell output to the log**: put `exec > {log} 2>&1` as the first
  line of every shell body (preferred over appending `2> {log}` to each command).
  This routes stdout *and* stderr — including bare `echo` status messages like
  cutadapt's `"no trimming"` — into the per-rule log instead of leaking into the
  Snakemake console/master log where they are hard to attribute. All shell-bearing
  rules do this; discovery rules use the inline `> {log} 2>&1` / `2> {log}`
  equivalents on their single pipeline.
- Mark intermediate files with `temp()`
- Use wildcards in paths: `{sample}`, `{results}`, `{read}`
- Resource specifications: `threads`, `mem_mb`
- Use `expand()` for generating multiple outputs

**Accessing config**: Use helper functions like `_get_config(sample, item)`
```python
def _get_config(sample, item):
    # Hierarchical lookup: sample -> chemistry[platform] -> chemistry -> defaults
```
Diagnostic messages from `_get_config` (missing chemistry/platform, or an item
falling through to the empty-string default) are written to **`sys.stderr`**, not
`print()`/stdout. Snakemake buffers stdout independently of its own stderr log
stream, so `print()` here lands detached (often bunched at the end of the run);
`sys.stderr.write(...)` keeps the messages in order at their point of origin.

---

## Configuration Files

### config.yaml
- `DATA`: Directory containing input FASTQs
- `RESULTS`: Output directory path
- `STAR_INDEX`: Path to STAR genome index
- `POLYA_SITES`: PolyA database reference file (SAF format; see SAF spec below)
- `GENOME_FASTA`: Genome FASTA (with `.fai`); required only when `DISCOVERY.enabled` is true
- `DISCOVERY`: Optional de novo RT priming site discovery block (see `docs/discovery.md`)
- `DEFAULTS`: Default chemistry and platform settings
- `SAMPLES`: Per-sample configuration (basename, chemistry, alignments)

### chemistry.yaml
Platform-specific configurations organized hierarchically:
```yaml
chemistry_name:
  bc_whitelist: path/to/whitelist
  platform_name:
    cutadapt_R1: "trimming parameters"
    STAR_R1: "alignment parameters"     # single-mate R1 rule: --clip5pNbases = 1 value
    STAR_paired: "alignment parameters" # two-mate paired rule: --clip5pNbases = 2 values (<R1clip> 0)
    STAR_R2: "alignment parameters"
```

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

- Delimiter `{D}`: `;` human, `_` mouse (matches `ref/polyadb32.{hg38,mm10}.saf.gz`).
- `pos` = field 5 (1-based cleavage position); `pas_type` = field 7 (class).
  Missing values are `NA`. `pos` may occasionally appear in float/scientific
  notation (e.g. `7.7e+07` from a float round-trip, as in one legacy row of the
  bundled `ref/polyadb32.mm10.saf.gz`, since repaired). Parsers must coerce it
  with `int(round(float(pos)))` rather than `int(pos)` so such rows are recovered
  instead of silently skipped — both `annotate_sites.py` (`load_polya`) and
  `inst/scripts/polyadb/validate_saf.py` do this. The converters in
  `inst/scripts/polyadb/` always emit plain integers.
- `Start`/`End` = strand-specific 15 bp window (transcript -10/+5):
  `+` → `[pos-10, pos+5]`, `-` → `[pos-5, pos+10]`.
- The `_` (mouse) delimiter appears inside RefSeq IDs and scaffold chrom names;
  parse by anchoring on the SAF `Chr`/`Strand` columns and reading trailing
  fields from the right (see `annotate_sites.py` `_split_geneid`).
- PAS type vocab is release-specific and passed through verbatim (3.2 vs v4
  differ, e.g. `intergenic` vs `Intergenic`).

---

## Common Development Tasks

### Adding a New Rule

1. Create rule in appropriate file under `rules/`
2. Follow naming convention: `verb_target` (e.g., `assign_sites_R1`)
3. Add to workflow by including outputs in `SAMPLE_OUTS` (Snakefile)
4. Test with dry-run: `snakemake -npr`

### Modifying Chemistry Configuration

1. Edit `chemistry.yaml`
2. Ensure all required fields present: `cutadapt_*`, `STAR_*`
3. Optional fields: `bc_whitelist`, `bc_cut`, `bc_length1`
4. Test with dry-run to validate YAML syntax

### Adding Python Helper Script

1. Place in `inst/scripts/`
2. Use argparse for CLI interface
3. Include docstring explaining purpose
4. Make executable: `chmod +x script.py`
5. Call from Snakemake rule with `python3 inst/scripts/script.py`

### De novo RT Priming Site Discovery

Optional module (off by default); see `docs/discovery.md` for full details.

- Rules: `rules/discovery.snake` (`stranded_bed`, `pool_beds`, `kde_peaks`,
  `annotate_sites`). Wired into `all_outputs` in `Snakefile` only when
  `DISCOVERY.enabled` is true; requires `GENOME_FASTA`.
- Scripts: `inst/scripts/kde_peaks.py` (KDE peak calling),
  `inst/scripts/annotate_sites.py` (polyAdb match + sequence-context
  classification + de novo SAF).
- Annotation classes/SAF `class` values: `known_pas`,
  `likely_internal_priming`, `potential_novel_pas`. The de novo SAF reuses the
  standard `gene;genbank;id;chrom;pos;strand;class` GeneID encoding so it can be
  set as `POLYA_SITES` to re-quantify (filter artifacts first).
- Outputs under `{results}/discovery/`: `<unit>_sites.tsv.gz`,
  `<unit>_sites.bed.gz`, `<unit>_denovo.saf.gz` (unit = group name or sample).
- Validate with `snakemake -npr` after setting `DISCOVERY.enabled: true`.

### Generating/converting polyAdb references

Standalone converters in `inst/scripts/polyadb/` (stdlib only + UCSC liftOver,
auto-downloaded). Not part of the Snakemake DAG.

- `pas32_to_saf.py`: polyAdb 3.2 `*.PAS.txt` (human hg19→hg38 via liftOver;
  mouse mm10 no liftOver). Resolves columns by header name.
- `pasv4_to_saf.py`: polyAdb 4 `*.PAS.{main,max}.tsv` (human/mouse; main/max
  auto-detected; optional liftOver for future builds).
- `saf_common.py`: shared 7-field GeneID encoding, window, dedup/sort, liftOver.
- `validate_saf.py`: format checks + optional content comparison to a reference.
- `--species` sets the delimiter (`;` human / `_` mouse). Output `.gz` writes
  gzipped, droppable straight into `ref/`.
- PAS type passed through verbatim (3.2 vs v4 class vocab differs).
- Verify: run on a head-slice, then `validate_saf.py generated.saf ref/...saf.gz`
  (human 3.2 should match the bundled reference 100% at shared positions).

---

## Error Handling and Debugging

**Log files**: All rules write logs to `{results}/logs/`, named
`{rulename}_{sample}_{read}.txt` (or `{rulename}_{unit}.txt` for discovery,
`multiqc.txt` / `check_versions.txt` for aggregates).
- Check logs for detailed error messages
- Logs capture BOTH stdout and stderr (`exec > {log} 2>&1`), so status echoes
  (e.g. cutadapt's `"no trimming"`) and tool output land in the per-rule log,
  not the Snakemake console/master log.

**Common issues**:
- Missing conda dependencies → check `scraps_conda.yml`
- YAML syntax errors → validate with `snakemake -npr`
- Missing input files → check `DATA` path in config.yaml
- Resource exhaustion → adjust `mem_mb` or `threads` in rules

**Debugging Snakemake**:
```bash
# Show detailed execution plan
snakemake -npr --verbose

# Print shell commands without execution
snakemake -np --printshellcmds

# Force re-run specific rule
snakemake --forcerun rulename
```

---

## Dependencies and Tools

**Core requirements** (installed via conda):
- Python >= 3.7
- Snakemake >= 5.3.0, < 8
- STAR >= 2.7.9a (RNA-seq aligner)
- UMI-tools >= 1.1.2 (UMI handling)
- cutadapt >= 3.4 (adapter trimming)
- samtools >= 1.15 (BAM manipulation)
- bedtools >= 2.30.0 (genomic intervals)
- subread >= 2.0.1 (featureCounts)
- MultiQC >= 1.6 (report generation)
- pysam >= 0.16.0 (Python BAM interface)
- numpy, pandas, scipy, scikit-learn (de novo discovery: KDE + annotation)

**Version checking**: Run `snakemake --configfile config.yaml` to trigger version checks

---

## Notes for AI Agents

- **Always dry-run first**: Use `snakemake -npr` before any pipeline changes
- **Respect shell choice**: Pipeline explicitly uses `zsh`, not bash
- **Preserve temp files**: Snakemake manages cleanup via `temp()` directive
- **Follow hierarchical config**: Sample → Chemistry/Platform → Defaults
- **Log everything**: Redirect stderr to log files for debugging
- **Resource awareness**: Bioinformatics tools are memory/CPU intensive
- **No traditional tests**: Validation is via successful Snakemake dry-run
