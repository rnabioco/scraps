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
  - <a href="https://github.com/CGATOxford/UMI-tools">UMI-tools<a> (developed with version 1.0.0)
  - <a href="http://samtools.sourceforge.net/">Samtools</a> (developed with version 1.5)
  - <a href="http://subread.sourceforge.net/">Subread<a> (developed with version 1.6.1)

## Example usage

scraps requires the following as input (defined in config.yaml):

  - 10X Genomics 3' v2/3 single-cell FASTQs 
  - Cell barcode whitelists for each sample (UMI-tools format)
  - A STAR genome index
  - A featureCounts reference (an SAF-formatted <a href="http://exon.umdnj.edu/polya_db/">pola_db</a> example is included)

See the <a href="https://snakemake.readthedocs.io/en/stable/">Snakemake</a> documentation
for general information on executing and manipulating snakemake pipelines.
