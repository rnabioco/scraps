# scraps <img src="man/figures/logo.png" align="right">

scraps (**S**ingle **C**ell **R**NA **P**rocessing **S**oftware) extracts
mRNA polyadenylation sites from "TVN"-primed single-cell RNA-seq
libraries.

Currently implemented as a snakemake pipeline for 10X Genomics
3' end v2/3 libraries, scraps will eventually be expanded into 
a standalone package for analyzing a range of RNA processing
changes in single-cell RNA-seq data.

For more information, see our manuscript [IN PREPARATION].

## Dependencies

scraps requires the following executables in your PATH:

  - <a href="https://www.python.org/">Python 3</a> (developed with version 3.6.3)
  - <a href="https://bitbucket.org/snakemake/snakemake/src/master/">Snakemake</a> (developed with version 3.11.2)
  - <a href="https://github.com/CGATOxford/UMI-tools">UMI-tools</a> (developed with version 1.0.0)
  - <a href="http://samtools.sourceforge.net/">Samtools</a> (developed with version 1.5)
  - <a href="http://subread.sourceforge.net/">Subread</a> (developed with version 1.6.1)

## Example usage

scraps requires the following as input (defined in config.yaml):

  - 10X Genomics 3' v2/3 single-cell FASTQs 
  - Optional cell barcode whitelists for each sample, in FASTQ folder  (can be auto-generated in UMI-tools format, single column lists also compatible)
  - A STAR genome index (compatible with cellranger indices <a href="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz">GRCh38</a> and <a href="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz">mm10</a>)
  - A featureCounts reference (SAF-formatted <a href="http://exon.umdnj.edu/polya_db/">polya_db</a>, hg38 and mm10 files are included in ref subdirectory)

See the <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> documentation
for general information on executing and manipulating snakemake pipelines.

<https://github.com/rnabioco/scraps>
