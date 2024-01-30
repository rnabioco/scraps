shell.executable("/bin/bash")

""" Snakemake pipeline for 10X Genomics 3' single-cell RNA-seq 3' end counting """

configfile: "config.yaml"

DATA = config["DATA"]
RESULTS = config["RESULTS"]
PAIRED_SAMPLES = config["PAIRED_SAMPLES"]
R1_SAMPLES = config["R1_SAMPLES"]
R2_SAMPLES = config["R2_SAMPLES"]
STAR_INDEX = config["STAR_INDEX"]
POLYA_SITES = config["POLYA_SITES"]
FASTQS = config["FASTQS"]
WHITELIST_V2 = config["WHITELIST_V2"]
WHITELIST_V3 = config["WHITELIST_V3"]
STAR = config["STAR"]
READS = ["paired", "R1", "R2"]

import json
with open('chemistry.json') as fp:
   chemistry = json.load(fp)

# paired positional approach
PAIRs = []
if PAIRED_SAMPLES:
  PAIRs.extend(expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = PAIRED_SAMPLES, read = "paired"))
  PAIRs.extend(expand("{results}/{sample}/{sample}_{read}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = PAIRED_SAMPLES, read = "paired"))
  PAIRs.extend(expand("{results}/bed/{sample}_{read}.bed.gz", results = RESULTS, sample = PAIRED_SAMPLES, read = "paired"))  
# read 1 positional approach
R1s = []
if R1_SAMPLES:
  R1s.extend(expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R1_SAMPLES, read = "R1"))
  R1s.extend(expand("{results}/{sample}/{sample}_{read}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = R1_SAMPLES, read = "R1"))
  R1s.extend(expand("{results}/bed/{sample}_{read}.bed.gz", results = RESULTS, sample = R1_SAMPLES, read = "R1"))  
# read 2 trimmed approach
R2s = []
if R2_SAMPLES:
  R2s.extend(expand("{results}/counts/{sample}_{read}_counts.tsv.gz", results = RESULTS, sample = R2_SAMPLES, read = "R2"))
  R2s.extend(expand("{results}/{sample}/{sample}_{read}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = R2_SAMPLES, read = "R2"))
  R2s.extend(expand("{results}/bed/{sample}_{read}.bed.gz", results = RESULTS, sample = R2_SAMPLES, read = "R2"))
# combine
all_outputs = PAIRs + R1s + R2s #+ expand("{results}/report/multiqc_report.html", results = RESULTS)
print(all_outputs)

rule all:
  input:
    all_outputs = all_outputs

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
