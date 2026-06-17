shell.executable("zsh")

""" Snakemake pipeline for single-cell RNA-seq 3' end counting """

configfile: "config.yaml"

DATA = config["DATA"]
RESULTS = config["RESULTS"]
STAR_INDEX = config["STAR_INDEX"]
POLYA_SITES = config["POLYA_SITES"]
STAR = config["STAR"]
DEFAULTS = config["DEFAULTS"]
SAMPLES = config["SAMPLES"]
READS = ["R1", "R2", "paired"]
GENOME_FASTA = config.get("GENOME_FASTA")
DISCOVERY = config.get("DISCOVERY", {}) or {}

import yaml
with open('chemistry.yaml') as fp:
   CHEMISTRY = yaml.safe_load(fp)

def _get_config(sample, item):
  try:
    return SAMPLES[sample][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[SAMPLES[sample]["chemistry"]][SAMPLES[sample]["platform"]][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[SAMPLES[sample]["chemistry"]][DEFAULTS["platform"]][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[SAMPLES[sample]["chemistry"]][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[DEFAULTS["chemistry"]][SAMPLES[sample]["platform"]][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[DEFAULTS["chemistry"]][DEFAULTS["platform"]][item]
  except KeyError:
    pass
  try:
    return CHEMISTRY[DEFAULTS["chemistry"]][item]
  except KeyError:
    return DEFAULTS[item]

# assemble outputs for rule all
SAMPLE_OUTS = []
for x in SAMPLES:
  SAMPLE_OUTS.extend(expand("{results}/counts/{sample}_{alignments}_counts.tsv.gz", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))
  SAMPLE_OUTS.extend(expand("{results}/{sample}/{sample}_{alignments}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))
  SAMPLE_OUTS.extend(expand("{results}/bed/{sample}_{alignments}.bed.gz", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))  

# optionally add de novo RT priming site discovery outputs
DISCOVERY_OUTS = []
if DISCOVERY.get("enabled", False):
  if GENOME_FASTA is None:
    raise ValueError("DISCOVERY.enabled is true but GENOME_FASTA is not set in config.yaml")
  disc_read = DISCOVERY.get("read", "R2")
  # samples that actually run the chosen discovery read
  disc_samples = [x for x in SAMPLES if disc_read in _get_config(x, "alignments")]
  groups = DISCOVERY.get("groups", {}) or {}
  grouped_samples = set()
  for members in groups.values():
    grouped_samples.update(members)
  # one discovery track per named group + per ungrouped sample
  disc_units = list(groups.keys()) + [x for x in disc_samples if x not in grouped_samples]
  for unit in disc_units:
    DISCOVERY_OUTS.extend(expand("{results}/discovery/{unit}_sites.tsv.gz", results = RESULTS, unit = unit))
    DISCOVERY_OUTS.extend(expand("{results}/discovery/{unit}_sites.bed.gz", results = RESULTS, unit = unit))
    DISCOVERY_OUTS.extend(expand("{results}/discovery/{unit}_denovo.saf.gz", results = RESULTS, unit = unit))

# optionally add multiqc
all_outputs = SAMPLE_OUTS + DISCOVERY_OUTS #+ expand("{results}/report/multiqc_report.html", results = RESULTS)
print(all_outputs)

rule all:
  input:
    all_outputs = all_outputs

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
include: "rules/discovery.snake"
