shell.executable("/usr/bin/env bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq 3' end counting """

configfile: "config.yaml"
  
DATA = config["DATA"]
SAMPLE = config["SAMPLES"]
GENOME_DIR = config["GENOME_DIR"]
POLYA_SITES = config["POLYA_SITES"]
#POLYA_FORMAT = config["POLYA_FORMAT"]

rule all:
  input:
    # Generates both read 1 (positional) and read 2 (polyA trimmed) counts by default;
    # Remove "R1" or "R2" to use only one method.
    expand("{data}/counts/{sample}_{read}_counts.tsv.gz", data = DATA, sample = SAMPLE, read = ["R1","R2"])

include: "rules/preprocess.snake"
include: "rules/star_two_pass.snake"
include: "rules/count.snake"
