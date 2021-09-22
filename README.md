# scraps <img src="man/figures/logo.png" align="right">

<!-- badges: start -->

![](https://img.shields.io/badge/dynamic/json.svg?label=Downloads&url=https://raw.githubusercontent.com/raysinensis/clone_counts_public/main/scraps.json&query=downloads&colorB=brightgreen)
<!-- badges: end -->

### scraps extracts mRNA polyadenylation sites from "TVN"-primed single-cell RNA-seq libraries.

scraps (**S**ingle **C**ell **R**N**A** **P**rocessing **S**oftware) is currently implemented as a <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> pipeline for 
10X Genomics 3' end v2/3 libraries (and other platforms with similar library structure, including Drop-seq, 
Microwell-seq, and BD Rhapsody). It will eventually be expanded for analyzing a range of RNA processing 
changes in single-cell RNA-seq data.

<img src="inst/flow.png" width="230" align="right">

---

## Example usage

scraps requires the following as input (defined in config.yaml and sample_fastqs.tsv):

  - 10X Genomics 3' v2/3 single-cell FASTQs or other platforms (with names "_R1.fastq.gz"" and "_R2fastq.gz"")
  - A STAR genome index (must be generated with STAR 2.7.4a and above)
  - Whitelist for cell barcodes (optional but recommended to speed up run time)
  - A featureCounts reference (SAF-formatted <a href="http://exon.umdnj.edu/polya_db/">polya_db</a>, hg38 and mm10 files are included in ref subdirectory)

To run test data, simply execute:
```
Snakemake --snakefile Snakefile \
  --configfile config.yaml \
  --resources total_impact=5 \
  --keep-going
```

---

## Supported scRNA-seq platforms
| Platform | Lib diagram | Setting name | Test data |
| :--------|:------------| :------------| :---------|
| 10x Chromium V3 | [16 + 12 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | chromiumV3 | ✓ |
| 10x Chromium V2 | [16 + 10 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | chromiumV2 | ✓ |
| 10x Chromium Visium | [16 + 10 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) | visium | |
| Drop-seq | [12 + 8 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/Drop-seq.html) | dropseq | ✓ |
| Microwell-seq | [6x3 + 6 + 30](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html) | microwellseq | ✓ |
| BD Rhapsody | [9x3 + 8 + 18](https://teichlab.github.io/scg_lib_structs/methods_html/BD_Rhapsody.html) | bd | |
| inDrop | [8 + 6 + 18](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html) | indrop | |

---

## Setup
1. Clone repository:
`
git clone https://github.com/rnabioco/scraps
`
2. Place appropriate STAR index in `index/` folder, and barcode whitelists in `ref/` <br>Download links: [GRCh38 index](https://scrapsaccessory.s3.us-west-2.amazonaws.com/GRCh38_cr2020A_star.tar.gz); [10x V2 barcodes](https://scrapsaccessory.s3.us-west-2.amazonaws.com/737K-august-2016.txt.gz); [10x V3 barcodes](https://scrapsaccessory.s3.us-west-2.amazonaws.com/3M-february-2018.txt.gz)
3. Edit settings in `config.yaml`
4. List files in `sample_fastqs.tsv`
5. Run!

___

## Output
1. bed : TVN-priming site pileup
2. count table : +-10 around PolyA_DB sites, by cell barcode

[R functions](https://github.com/rnabioco/scraps/tree/master/inst/scripts/R) available for importing results into Seurat object, and finding differential PA site usage

___

## Dependencies

scraps requires the following executables in your PATH:

  - <a href="https://www.python.org">Python 3</a> (developed with version 3.8.5)
  - <a href="https://bitbucket.org/snakemake/snakemake/src/master/">Snakemake</a> (developed with version 3.11.2)
  - <a href="https://github.com/CGATOxford/UMI-tools">UMI-tools</a> (developed with version 1.1.1)
  - <a href="https://cutadapt.readthedocs.io">cutadapt</a> (developed with version 3.4)
  - <a href="https://github.com/alexdobin/STAR">STAR</a> (developed with version 2.7.9a)
  - <a href="http://subread.sourceforge.net">Subread</a> (developed with version 1.6.2)
  - <a href="https://multiqc.info">MultiQC</a> (developed with version 1.9)

See the <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> documentation
for general information on executing and manipulating snakemake pipelines.
