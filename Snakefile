shell.executable("/bin/bash")

""" Snakemake pipeline for single-cell RNA-seq 3' end counting """

configfile: "config.yaml"

DATA = config["DATA"]
RESULTS = config["RESULTS"]
STAR_INDEX = config["STAR_INDEX"]
POLYA_SITES = config["POLYA_SITES"]
WHITELIST_V2 = config["WHITELIST_V2"]
WHITELIST_V3 = config["WHITELIST_V3"]
STAR = config["STAR"]
DEFAULTS = config["DEFAULTS"]
SAMPLES = config["SAMPLES"]
READS = ["R1", "R2", "paired"]

import json
with open('chemistry.json') as fp:
   chemistry = json.load(fp)
   
# with open('platform.json') as fp:
#    platform = json.load(fp)

def _get_config(sample, item):
  try:
    return SAMPLES[sample][item]
  except KeyError:
    return DEFAULTS[item]
  
# assemble outputs for rule all
SAMPLE_OUTS = []
for x in SAMPLES:
  SAMPLE_OUTS.extend(expand("{results}/counts/{sample}_{alignments}_counts.tsv.gz", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))
  SAMPLE_OUTS.extend(expand("{results}/{sample}/{sample}_{alignments}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))
  SAMPLE_OUTS.extend(expand("{results}/bed/{sample}_{alignments}.bed.gz", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))  

# optionally add multiqc
all_outputs = SAMPLE_OUTS #+ expand("{results}/report/multiqc_report.html", results = RESULTS)
print(all_outputs)

rule all:
  input:
    all_outputs = all_outputs

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
