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
│   └── check_versions.snake  # Dependency version checks
├── inst/scripts/          # Helper scripts
│   ├── *.py               # Python utilities (BAM filtering, etc.)
│   └── R/                 # R analysis functions
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
        "{results}/logs/{sample}_rulename.txt"
    threads:
        12
    resources:
        mem_mb = 8000
    shell:
        r"""
        command --arg {input} > {output} 2> {log}
        """
```

**Key conventions**:
- Use raw strings `r"""..."""` for shell blocks
- Redirect stderr to log files: `2> {log}`
- Mark intermediate files with `temp()`
- Use wildcards in paths: `{sample}`, `{results}`, `{read}`
- Resource specifications: `threads`, `mem_mb`
- Use `expand()` for generating multiple outputs

**Accessing config**: Use helper functions like `_get_config(sample, item)`
```python
def _get_config(sample, item):
    # Hierarchical lookup: sample -> chemistry[platform] -> chemistry -> defaults
```

---

## Configuration Files

### config.yaml
- `DATA`: Directory containing input FASTQs
- `RESULTS`: Output directory path
- `STAR_INDEX`: Path to STAR genome index
- `POLYA_SITES`: PolyA database reference file (SAF format)
- `DEFAULTS`: Default chemistry and platform settings
- `SAMPLES`: Per-sample configuration (basename, chemistry, alignments)

### chemistry.yaml
Platform-specific configurations organized hierarchically:
```yaml
chemistry_name:
  bc_whitelist: path/to/whitelist
  platform_name:
    cutadapt_R1: "trimming parameters"
    STAR_R1: "alignment parameters"
    STAR_R2: "alignment parameters"
```

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

---

## Error Handling and Debugging

**Log files**: All rules write logs to `{results}/logs/`
- Check logs for detailed error messages
- Logs include stderr from all commands

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
