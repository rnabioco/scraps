# scraps <img src="man/figures/logo.png" align="right">

<!-- badges: start -->
![Snakemake Check](https://github.com/rnabioco/scraps/actions/workflows/snakemake_run.yml/badge.svg)
![Github Clones](https://img.shields.io/badge/dynamic/json.svg?label=Downloads&url=https://raw.githubusercontent.com/raysinensis/clone_counts_public/main/scraps.json&query=downloads&colorB=brightgreen)
<!-- badges: end -->

### scraps extracts mRNA polyadenylation sites from "TVN"-primed single-cell RNA-seq libraries at near-nucleotide resolution.

scraps (**S**ingle **C**ell **R**N**A** **P**olyA **S**ite Discovery) is currently implemented as a <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> pipeline for 
10X Genomics 3' end v2/3 libraries (and other platforms with similar library structure, including Drop-seq, 
Microwell-seq, and BD Rhapsody). If long Read1 is available (estimated ~6% of SRA-deposited data, or now planning new experiments), positional information will be calculated from paired realignment; otherwise, the less optimal anchored Read2 approach is used. scraps will eventually be expanded for analyzing a range of RNA processing 
changes in single-cell RNA-seq data.

For additional discussions and usage cases, please see [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2022.08.22.504859v1).

<img src="inst/flow.png" width="400" align="right">

---
-   [Example usage](#example-usage)
-   [Configuration](#configuration)
-   [Supported scRNA-seq platforms](#supported-scrna-seq-platforms)
-   [Output](#output)
-   [Setup](#setup)
-   [Dependencies](#dependencies)
-   [For Developers](#for-developers)
-   [Extended function](#bonus-function)
    
## Example usage

scraps requires the following as input (defined in config.yaml):

  - 10X Genomics 3' v2/3 single-cell FASTQs or other platforms (with names "_R1.fastq.gz"" and "_R2.fastq.gz"")
  - A STAR genome index (must be generated with STAR 2.7.4a and above)
  - Whitelist for cell barcodes (optional but recommended to speed up run time)
  - A featureCounts reference (SAF-formatted <a href="http://exon.umdnj.edu/polya_db/">polya_db</a>, hg38 and mm10 files are included in [ref](https://github.com/rnabioco/scraps/tree/master/ref) subdirectory)

### Quick Start

1. Set up conda environment:
```bash
conda env create -f scraps_conda.yml
conda activate scraps_conda
```

2. Configure your samples in `config.yaml` under the `SAMPLES` section

3. Run the pipeline:
```bash
snakemake --configfile config.yaml --resources total_impact=5 --keep-going
```

### Detailed Usage

To run test data, simply execute:
```
snakemake --snakefile Snakefile \
  --configfile config.yaml \
  --resources total_impact=5 \
  --keep-going
```
[DAG steps illustration](https://raw.githack.com/rnabioco/scraps/master/inst/dag.pdf)

[submit jobs in cluster mode](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)

Notes: `total_impact` is set to 5 for each sample, change this to control how many samples are processed in parallel

---

## Configuration

scraps uses two main configuration files for flexible pipeline setup:

### config.yaml

Main pipeline configuration file containing:

- **DATA**: Directory containing input FASTQ files
- **RESULTS**: Output directory for pipeline results
- **STAR_INDEX**: Path to STAR genome index directory
- **POLYA_SITES**: PolyA database reference file (SAF format, provided in `ref/`)
- **DEFAULTS**: Default chemistry and platform settings
  - `platform`: Default sequencing platform (e.g., illumina, element, ultima)
  - `chemistry`: Default chemistry type (e.g., chromiumV3, chromiumV2, dropseq)
  - `alignments`: Which alignment modes to run ([R1, R2, paired])
- **SAMPLES**: Per-sample configuration

#### Sample Configuration Example:
```yaml
SAMPLES:
  sample_name:
    basename: sample-           # FASTQ file prefix
    platform: illumina          # Sequencing platform
    chemistry: chromiumV3       # Platform chemistry
    alignments:                 # Optional: override default alignments
      - R2
      - paired
```

### chemistry.yaml

Platform and chemistry-specific parameters organized hierarchically. Each chemistry type (chromiumV3, chromiumV2, dropseq, microwellseq, bd, indrop) contains:

- **bc_whitelist**: Path to barcode whitelist file (optional)
- **bc_cut**: Adapter sequences for complex barcode extraction (optional)
- **Platform-specific settings**: Nested configurations for different sequencing platforms
  - `cutadapt_R1` / `cutadapt_paired`: Adapter trimming parameters
  - `STAR_R1` / `STAR_R2`: STAR alignment parameters (UMI/barcode positions)

### Configuration Hierarchy

The pipeline uses hierarchical configuration lookup to determine parameters for each sample:

```
┌─────────────────────────────────────────────────────────┐
│  1. Sample-specific settings (config.yaml SAMPLES)     │
│     Highest priority - overrides everything             │
└────────────────────┬────────────────────────────────────┘
                     │ If not found ↓
┌─────────────────────────────────────────────────────────┐
│  2. Chemistry + Platform (chemistry.yaml)               │
│     e.g., chromiumV3 → illumina → STAR_R1               │
└────────────────────┬────────────────────────────────────┘
                     │ If not found ↓
┌─────────────────────────────────────────────────────────┐
│  3. Chemistry defaults (chemistry.yaml)                 │
│     e.g., chromiumV3 → bc_whitelist                     │
└────────────────────┬────────────────────────────────────┘
                     │ If not found ↓
┌─────────────────────────────────────────────────────────┐
│  4. Global defaults (config.yaml DEFAULTS)              │
│     Lowest priority - fallback values                   │
└─────────────────────────────────────────────────────────┘
```

This allows platform-specific customization (e.g., Illumina vs Ultima Genomics) while maintaining chemistry-specific defaults.

---

## Supported scRNA-seq platforms
| Platform | Library (BC+UMI+A) | Setting | Test data |
| :--------|:------------| :------------| :---------|
| 10x Chromium V3 | [16 + 12 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | chromiumV3 | ✓ |
| 10x V3 - Ultima Genomics | [adapter + 16 + 9 + 3 ignored + 8](https://www.nature.com/articles/s41587-022-01452-6/figures/1) | chromiumV3UG | |
| 10x Chromium V2 | [16 + 10 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | chromiumV2 | ✓ |
| 10x Chromium Visium | [16 + 10 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | visium | |
| Drop-seq | [12 + 8 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/Drop-seq.html) | dropseq | ✓ |
| Microwell-seq | [6x3 + 6 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html) | microwellseq | ✓ |
| BD Rhapsody | [9x3 + 8 + 18](https://teichlab.github.io/scg_lib_structs/methods_html/BD_Rhapsody.html) | bd | |
| inDrop | [8 + 6 + 18](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html) | indrop | |

Custom chemistry supported, by editing [chemistry.yaml](https://github.com/rnabioco/scraps/tree/master/chemistry.yaml). Also see synthetic FASTQ [tool](https://github.com/raysinensis/bc_umi_gen).

---

## Output
1. bedgraph : TVN-priming site pileup
```
chr11   215106  215107  1
chr11   689216  689217  1
chr11   812862  812863  1
chr11   812870  812871  2
chr11   812871  812872  2
```
2. count table : +-10 around PolyA_DB sites, by cell barcode
```
gene    cell    count
AC135178.2_NA_ENSG00000263809_chr17_8377523_-_Intron,RPL26_6154_ENSG00000161970_chr17_8377523_-_3'UTR(M)        AACTCCCGTTCCTCCA        1
AC135178.2_NA_ENSG00000263809_chr17_8377523_-_Intron,RPL26_6154_ENSG00000161970_chr17_8377523_-_3'UTR(M)        CCCATACGTTAAAGAC        1
AC135178.2_NA_ENSG00000263809_chr17_8377523_-_Intron,RPL26_6154_ENSG00000161970_chr17_8377523_-_3'UTR(M)        CGTCCATTCGACAGCC        1
ACTG1_71_ENSG00000184009_chr17_81509999_-_3'UTR(M)      ACATCAGGTGATGTCT        1
ADRM1_11047_ENSG00000130706_chr20_62308862_+_3'UTR(M)   CAGCGACTCTGCCCTA        1
```
3. [html report](https://raw.githack.com/rnabioco/scraps/master/inst/test_output/report/multiqc_report.html) : various metrics from steps in the pipeline

[R functions](https://github.com/rnabioco/scraps/tree/master/inst/scripts/R) available for importing results into Seurat object, and finding differential PA site usage. Alternatively, a package of the same functions can be installed with [`remotes::install_github("rnabioco/scrapR")`](https://github.com/rnabioco/scrapR)

<p float="left">
  <img src="man/figures/example.png" height="300"/>
  <img src="man/figures/pa_setx.png" height="300"/>
</p>

___

## Setup

### 1. Clone repository
```bash
git clone https://github.com/rnabioco/scraps
cd scraps
```

### 2. Set up conda environment (recommended)
```bash
conda env create -f scraps_conda.yml
conda activate scraps_conda
```

Alternatively, ensure all [dependencies](#dependencies) are installed and available in your PATH.

### 3. Prepare reference files

#### STAR genome index
Place STAR index in the `ref/` directory or specify custom path in config.yaml (`STAR_INDEX`)

Download link (extract after download):
- [GRCh38 index](https://scrapsaccessory.s3.us-west-2.amazonaws.com/GRCh38_cr2020A_star.tar.gz)

#### Barcode whitelists (optional but recommended)
Whitelist paths are configured per chemistry in `chemistry.yaml`. Place downloaded whitelists in the `ref/` directory.

Download links (extract after download):
- [10x V2 barcodes](https://scrapsaccessory.s3.us-west-2.amazonaws.com/737K-august-2016.txt.gz) → `ref/737K-august-2016.txt`
- [10x V3 barcodes](https://scrapsaccessory.s3.us-west-2.amazonaws.com/3M-february-2018.txt.gz) → `ref/3M-february-2018.txt`

Update `chemistry.yaml` with the correct paths:
```yaml
chromiumV3:
  bc_whitelist: ref/3M-february-2018.txt
chromiumV2:
  bc_whitelist: ref/737K-august-2016.txt
```

### 4. Configure your samples

Edit `config.yaml` to specify:
- **DATA**: Path to directory containing FASTQ files (with naming pattern: `*_R1.fastq.gz`, `*_R2.fastq.gz`)
- **RESULTS**: Output directory path
- **STAR_INDEX**: Path to STAR genome index
- **POLYA_SITES**: PolyA database reference (provided: `ref/polyadb32.hg38.saf.gz` or `ref/polyadb32.mm10.saf.gz`)
- **SAMPLES**: Define each sample with:
  - `basename`: FASTQ filename prefix
  - `chemistry`: Platform chemistry type (chromiumV2, chromiumV3, dropseq, microwellseq, bd, indrop)
  - `platform`: Sequencing platform (illumina, element, ultima)
  - `alignments`: Optional list of alignment modes to run

Example:
```yaml
SAMPLES:
  my_sample:
    basename: SRR9887775_        # Matches SRR9887775_R1.fastq.gz, SRR9887775_R2.fastq.gz
    chemistry: chromiumV3
    platform: illumina
```

**Note**: SRA accessions (e.g., `SRR9887775`) can be used directly as basenames for automatic download.

### 5. Run the pipeline

```bash
# Dry-run to check configuration
snakemake -npr --configfile config.yaml

# Run pipeline
snakemake --configfile config.yaml --resources total_impact=5 --keep-going

# Or with specific core count
snakemake -j 8 --configfile config.yaml
```

Sample test results can be found at [inst/test_output/](inst/test_output/)

___

## Dependencies

scraps requires the following executables in your PATH:

  - <a href="https://www.python.org">Python 3</a> (>= 3.7)
  - <a href="https://snakemake.readthedocs.io">Snakemake</a> (>= 5.3.0, < 8.0)
  - <a href="https://github.com/CGATOxford/UMI-tools">UMI-tools</a> (>= 1.1.2)
  - <a href="https://cutadapt.readthedocs.io">cutadapt</a> (>= 3.4)
  - <a href="https://github.com/alexdobin/STAR">STAR</a> (>= 2.7.9a)
  - <a href="https://www.htslib.org">Samtools</a> (>= 1.15)
  - <a href="https://bedtools.readthedocs.io/en/latest">Bedtools</a> (>= 2.30.0)
  - <a href="http://subread.sourceforge.net">Subread</a> (>= 2.0.1)
  - <a href="https://multiqc.info">MultiQC</a> (>= 1.6)
  - <a href="https://pysam.readthedocs.io">pysam</a> (>= 0.16.0)
  - <a href="https://github.com/drmaa-python/drmaa-python">drmaa</a> (>= 0.7.9, for cluster execution)
  - zsh (shell used for pipeline execution)
  
**Recommended**: Use [Conda](https://docs.conda.io/en/latest/) to manage these dependencies:
```bash
conda env create -f scraps_conda.yml
conda activate scraps_conda
```
All required dependencies (including zsh) will be installed automatically.

Docker image for automated deployment can also be found at https://hub.docker.com/r/rnabioco/scraps.

Please also see the <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> documentation
for general information on executing and manipulating snakemake pipelines.

---

## For Developers

For detailed development guidelines including code style conventions, testing procedures, and instructions for adding new rules or chemistry configurations, see [AGENTS.md](AGENTS.md).

Key resources:
- **Code style**: Python, R, and Snakemake conventions
- **Testing**: Using Snakemake dry-run for validation
- **Adding rules**: How to extend the pipeline
- **Chemistry config**: Adding new platform support
- **Debugging**: Common issues and solutions

---

## Extended function

<img src="inst/apop.png" height="250" align="left">

**1) Measuring internal priming as indicator of apoptotic cytoplasmic poly(A) RNA decay**

(Based on widespread RNA decay during apoptosis: [Liu and Fu et al.](https://www.sciencedirect.com/science/article/pii/S0092867418305105))
Use SAF (hg38 version provided in [ref](https://github.com/rnabioco/scraps/tree/master/ref) subdirectory) file marking all gene regions (5'UTR, intron, CDS, 3'UTR), and helper [R functions](https://github.com/rnabioco/scraps/tree/master/inst/scripts/R/scraps_priming_region.R) to process output. Please see [Rmarkdown notebook](https://github.com/rnabioco/scraps/blob/master/inst/notebooks/scraps_internalpriming_5fu_2023.Rmd) for more.

---

**2) Accurate intron/exon quantification for RNA velocity**

(See discussions on quantification approaches and pitfalls: [Soneson et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008585))

| Consideration | scraps |
|:--------|:------------|
| Avoid feature double-counting | ✓ |
| Take strandedness into account | ✓ |
| Avoid count substraction | ✓ |
| Resolve spliced vs unspliced target | ✓ |
| Speed | ✓ |
