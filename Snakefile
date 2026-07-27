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

import sys
import yaml
with open('chemistry.yaml') as fp:
   CHEMISTRY = yaml.safe_load(fp)

def _validate_config():
  """Fail fast on invalid chemistry/platform before _get_config's empty-string
  fallback silently propagates empty args into cutadapt/STAR rules."""
  errors = []
  meta_keys = {"bc_whitelist", "bc_cut", "bc_length1"}  # scalar chemistry-level keys, not platforms
  for sample, cfg in SAMPLES.items():
    cfg = cfg or {}
    chemistry = cfg.get("chemistry", DEFAULTS.get("chemistry"))
    platform  = cfg.get("platform",  DEFAULTS.get("platform"))

    if chemistry not in CHEMISTRY:
      errors.append(
        f"sample '{sample}': chemistry '{chemistry}' not found in chemistry.yaml. "
        f"Valid options: {sorted(CHEMISTRY)}"
      )
      continue  # platform check meaningless without valid chemistry

    platforms = [k for k, v in CHEMISTRY[chemistry].items()
                 if isinstance(v, dict) and k not in meta_keys]
    if platform not in platforms:
      errors.append(
        f"sample '{sample}': platform '{platform}' not found under chemistry "
        f"'{chemistry}'. Valid options: {sorted(platforms)}"
      )

  if errors:
    raise ValueError(
      "Invalid chemistry/platform configuration:\n  - " + "\n  - ".join(errors)
    )

_validate_config()

def _validate_modes():
  """Fail fast when a sample requests alignment modes its chemistry/platform
  does not support (valid_modes in chemistry.yaml). Runs after
  _validate_config() so chemistry/platform are known-valid."""
  errors = []
  for sample, cfg in SAMPLES.items():
    cfg = cfg or {}
    chemistry = cfg.get("chemistry", DEFAULTS.get("chemistry"))
    platform  = cfg.get("platform",  DEFAULTS.get("platform"))
    plat_cfg = CHEMISTRY.get(chemistry, {}).get(platform, {})
    valid = plat_cfg.get("valid_modes")
    if not valid:
      errors.append(
        f"sample '{sample}': chemistry '{chemistry}' platform '{platform}' "
        f"defines no valid_modes in chemistry.yaml"
      )
      continue
    modes = cfg.get("alignments", DEFAULTS.get("alignments", []))
    if isinstance(modes, str):
      modes = [modes]
    bad = [m for m in modes if m not in valid]
    if bad:
      errors.append(
        f"sample '{sample}': alignment mode(s) {bad} not supported by "
        f"chemistry '{chemistry}' platform '{platform}'. Valid: {sorted(valid)}"
      )
  if errors:
    raise ValueError(
      "Invalid alignment mode configuration:\n  - " + "\n  - ".join(errors)
    )

_validate_modes()

def _get_config(sample, item):
  sample_cfg = SAMPLES[sample]

  try:
    chemistry = sample_cfg.get("chemistry", DEFAULTS["chemistry"])
  except KeyError:
    sys.stderr.write("Error: Chemistry must be defined per sample and/or in Defaults.\n")
    raise

  try:
    platform = sample_cfg.get("platform", DEFAULTS["platform"])
  except KeyError:
    sys.stderr.write("Error: Platform must be defined per sample and/or in Defaults.\n")
    raise

  if item in sample_cfg:
    return sample_cfg[item]

  if item in CHEMISTRY.get(chemistry, {}).get(platform, {}):
    return CHEMISTRY[chemistry][platform][item]

  if item in CHEMISTRY.get(chemistry, {}):
    return CHEMISTRY[chemistry][item]

  if item in DEFAULTS:
    return DEFAULTS[item]

  sys.stderr.write(
    f"Message: {item} not found in config or chemistry "
    f"for sample {sample}. Returning an empty string.\n"
  )
  return ""

# assemble outputs for rule all
SAMPLE_OUTS = []
for x in SAMPLES:
  SAMPLE_OUTS.extend(expand("{results}/counts/{sample}_{alignments}_counts.tsv.gz", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))
  SAMPLE_OUTS.extend(expand("{results}/{sample}/{sample}_{alignments}_Aligned.sortedByCoord.out.bam", results = RESULTS, sample = x, alignments = _get_config(x, "alignments")))

# --- discovery unit resolution (shared by output assembly and rules) ---------
# A discovery "unit" is an aggregation of one or more (sample, mode) stranded
# priming beds. DISC_UNITS maps unit_id -> [(sample, mode), ...]. Units come from
# three sources, all mode-suffixed and mutually unambiguous:
#   - per (sample, mode)         : id = "{sample}_{mode}"   (default)
#   - per (group, mode)          : id = "{group}_{mode}"    (DISCOVERY.groups)
#   - explicit DISCOVERY.units   : id = user name, members = [[sample, mode],...]
# Modes are resolved per sample (per-sample override -> global list -> the
# sample's own `alignments`) and always intersected with what the sample was
# actually aligned in, so only runnable (sample, mode) beds are requested.
def _as_mode_list(val):
  """Coerce a scalar-or-list mode spec into a list."""
  if val is None:
    return None
  return [val] if isinstance(val, str) else list(val)

def _sample_disc_modes(sample):
  """Discovery modes for one sample: per-sample override -> global read ->
  the sample's alignments; then intersect with the sample's alignments."""
  aln = _as_mode_list(_get_config(sample, "alignments")) or []
  override = _as_mode_list(SAMPLES[sample].get("discovery_read")
                           if isinstance(SAMPLES[sample], dict) else None)
  glob = _as_mode_list(DISCOVERY.get("read"))
  want = override if override is not None else (glob if glob is not None else aln)
  # preserve requested order, keep only modes the sample was aligned in
  return [m for m in want if m in aln]

def _discovery_units():
  """Return {unit_id: [(sample, mode), ...]} for all discovery tracks."""
  units = {}
  groups = DISCOVERY.get("groups", {}) or {}
  explicit = DISCOVERY.get("units", {}) or {}
  grouped_samples = set()
  for members in groups.values():
    grouped_samples.update(members)

  # explicit units: arbitrary (sample, mode) aggregation (cross-sample/mode)
  for name, members in explicit.items():
    pairs = []
    for m in members:
      s, mode = (m[0], m[1]) if not isinstance(m, str) else (m, None)
      modes = [mode] if mode else _sample_disc_modes(s)
      for md in modes:
        if md in _sample_disc_modes(s) or mode:
          pairs.append((s, md))
    if pairs:
      units[name] = pairs

  # groups: pool member samples within each shared mode -> "{group}_{mode}"
  for gname, members in groups.items():
    # modes shared by all members (intersection, order from first member)
    per = [ _sample_disc_modes(s) for s in members ]
    if not per:
      continue
    shared = [m for m in per[0] if all(m in p for p in per[1:])]
    for md in shared:
      units["{}_{}".format(gname, md)] = [(s, md) for s in members]

  # ungrouped samples: one unit per (sample, mode) -> "{sample}_{mode}"
  for s in SAMPLES:
    if s in grouped_samples:
      continue
    for md in _sample_disc_modes(s):
      units["{}_{}".format(s, md)] = [(s, md)]

  return units

DISC_UNITS = {}
DISCOVERY_OUTS = []
if DISCOVERY.get("enabled", False):
  if GENOME_FASTA is None:
    raise ValueError("DISCOVERY.enabled is true but GENOME_FASTA is not set in config.yaml")
  DISC_UNITS = _discovery_units()
  for unit in DISC_UNITS:
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
