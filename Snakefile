shell.executable("/usr/bin/env bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq 3' end counting """

configfile: "config.yaml"
  
DATA = config["DATA"]
RESULTS = config["RESULTS"]
R1_SAMPLES = config["R1_SAMPLES"]
R2_SAMPLES = config["R2_SAMPLES"]
STAR_INDEX = config["STAR_INDEX"]
POLYA_SITES = config["POLYA_SITES"]
FASTQS = config["FASTQS"]
CAPTURES = config["CAPTURES"]
WHITELIST_V2 = config["WHITELIST_V2"]
WHITELIST_V3 = config["WHITELIST_V3"]
STAR = config["STAR"]


rule all:
  input:
    # Generates read 1 (positional) counts;
    expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R1_SAMPLES, read = "R1"),
    # Generates read 2 (trimmed0 counts;
    expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R2_SAMPLES, read = "R2"),
    expand("{data}/multiqc_report.html", data = DATA)

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
