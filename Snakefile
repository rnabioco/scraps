shell.executable("/usr/bin/env bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq 3' end counting """

configfile: "config.yaml"
  
DATA = config["DATA"]
SAMPLE = config["SAMPLES"]
GENOME_DIR = config["GENOME_DIR"]
POLYA_SITES = config["POLYA_SITES"]
BC_PATTERN = config["BC_PATTERN"]
READS = config["READS"]

rule all:
  input:
    # Generates both read 1 (positional) and read 2 (polyA trimmed) counts by default;
    expand("{data}/counts/{sample}_{read}_counts.tsv.gz", data = DATA, sample = SAMPLE, read = READS),
    "{data}/multiqc_report.html"

include: "rules/check_versions.snake"
include: "rules/preprocess.snake"
include: "rules/star_two_pass.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
