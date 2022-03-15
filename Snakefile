shell.executable("/bin/bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq 3' end counting """

configfile: "config.yaml"

DATA = config["DATA"]
RESULTS = config["RESULTS"]
R1_SAMPLES = config["R1_SAMPLES"]
R2_SAMPLES = config["R2_SAMPLES"]
STAR_INDEX = config["STAR_INDEX"]
POLYA_SITES = config["POLYA_SITES"]
FASTQS = config["FASTQS"]
WHITELIST_V2 = config["WHITELIST_V2"]
WHITELIST_V3 = config["WHITELIST_V3"]
STAR = config["STAR"]
READS = ["R1", "R2"]

import json
with open('chemistry.json') as fp:
   chemistry = json.load(fp)

rule all:
  input:
    # Generate read 1 & 2 BAMs;
    expand("{results}/{sample}/{sample}_{read}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = R1_SAMPLES, read = "R1"),
    # Generate read 2 BAMS;
    expand("{results}/{sample}/{sample}_{read}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = R2_SAMPLES, read = "R2"),
    # Generates read 1 (positional) counts;
    expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R1_SAMPLES, read = "R1"),
    # Generates read 2 (trimmed) counts;
    expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R2_SAMPLES, read = "R2"),
    expand("{results}/report/multiqc_report.html", results = RESULTS),
    # Generates bed files;
    expand("{results}/bed/{sample}_{read}.bed.gz", results = RESULTS, sample = R1_SAMPLES, read = "R1"),
    expand("{results}/bed/{sample}_{read}.bed.gz", results = RESULTS, sample = R2_SAMPLES, read = "R2")

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
