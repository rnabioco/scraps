shell.executable("zsh")

""" Snakemake pipeline for single-cell RNA-seq 3' end counting """

from snakemake.utils import min_version
min_version("9.0")

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

import os
import sys
import yaml
from snakemake.logging import logger
with open('chemistry.yaml') as fp:
   CHEMISTRY = yaml.safe_load(fp)

# Diagnostics go through Snakemake's logger (verified on Snakemake 9): in a
# dry-run that is stdout, in order with the job listing; otherwise stderr plus
# .snakemake/log/, so `--dag`/`--rulegraph` stdout stays clean DOT. Strings are
# built with str.format() for consistency (Snakemake 7's tokenizer re-spaced
# f-string contents; fixed in 8+).
_diag_seen = set()

def _diag(key, msg):
  """Emit `msg` once per `key` (params lambdas are evaluated once per job)."""
  if key in _diag_seen:
    return
  _diag_seen.add(key)
  logger.warning(msg)

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
        "sample '{}': chemistry '{}' not found in chemistry.yaml. "
        "Valid options: {}".format(sample, chemistry, sorted(CHEMISTRY))
      )
      continue  # platform check meaningless without valid chemistry

    platforms = [k for k, v in CHEMISTRY[chemistry].items()
                 if isinstance(v, dict) and k not in meta_keys]
    if platform not in platforms:
      errors.append(
        "sample '{}': platform '{}' not found under chemistry '{}'. "
        "Valid options: {}".format(sample, platform, chemistry, sorted(platforms))
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
        "sample '{}': chemistry '{}' platform '{}' defines no valid_modes in "
        "chemistry.yaml".format(sample, chemistry, platform)
      )
      continue
    modes = cfg.get("alignments", DEFAULTS.get("alignments", []))
    if isinstance(modes, str):
      modes = [modes]
    bad = [m for m in modes if m not in valid]
    if bad:
      errors.append(
        "sample '{}': alignment mode(s) {} not supported by chemistry '{}' "
        "platform '{}'. Valid: {}".format(sample, bad, chemistry, platform,
                                          sorted(valid))
      )
  if errors:
    raise ValueError(
      "Invalid alignment mode configuration:\n  - " + "\n  - ".join(errors)
    )

_validate_modes()

# keys that are legitimately absent (rules treat "" as "off"); no diagnostic
_OPTIONAL_KEYS = {"cutadapt_R1", "cutadapt_paired", "bc_cut", "bc_length1",
                  "bam_strand", "bam_mate", "STAR_mem_mb"}
_OPTIONAL_SUFFIXES = ("_extra_args", "_polyA_filter", "_priming_offset")

def _get_config(sample, item):
  """Hierarchical lookup: sample -> chemistry.platform -> chemistry -> DEFAULTS
  -> "" (with a one-time warning unless the key is optional)."""
  sample_cfg = SAMPLES[sample] or {}

  # per-sample value wins; DEFAULTS is consulted only when the sample has none
  chemistry = sample_cfg.get("chemistry") or DEFAULTS.get("chemistry")
  if chemistry is None:
    raise ValueError(
      "sample '{}': chemistry must be defined per sample or in DEFAULTS".format(sample))

  platform = sample_cfg.get("platform") or DEFAULTS.get("platform")
  if platform is None:
    raise ValueError(
      "sample '{}': platform must be defined per sample or in DEFAULTS".format(sample))

  if item in sample_cfg:
    return sample_cfg[item]

  if item in CHEMISTRY.get(chemistry, {}).get(platform, {}):
    return CHEMISTRY[chemistry][platform][item]

  if item in CHEMISTRY.get(chemistry, {}):
    return CHEMISTRY[chemistry][item]

  if item in DEFAULTS:
    return DEFAULTS[item]

  if item not in _OPTIONAL_KEYS and not item.endswith(_OPTIONAL_SUFFIXES):
    _diag(("missing", sample, item),
          "Warning: '{}' not found in config.yaml or chemistry.yaml for sample "
          "'{}' (chemistry '{}', platform '{}'); using an empty "
          "string.".format(item, sample, chemistry, platform))
  return ""

# --- STARsolo memory ----------------------------------------------------------
# STAR's peak RAM is the larger of two phases (STAR source: STAR.cpp,
# bamSortByCoordinate.cpp; --genomeLoad NoSharedMemory, STAR's default):
#   1. mapping: genome index + per-thread I/O and BAM-bin buffers (a few GB);
#   2. after mapping the index is freed (genomeMain.freeMemory()), then STARsolo
#      builds its per-read structures (scale with read count; kept through the
#      sort for the CB/UB tags) and the BAM sort loads bins up to
#      --limitBAMsortRAM (a cap on in-flight bins, not an allocation).
# STAR_INDEX_MB is the index size on disk (32000, a human index, when STAR_INDEX
# is absent, e.g. a CI dry-run). STAR_mem_mb (sample/platform/chemistry/
# DEFAULTS) sets the job memory; unset, it covers phase 1 with 8 GB of buffers,
# floored at 48 GB. The starsolo rules multiply it by the attempt number
# (profile `retries`) and derive --limitBAMsortRAM from the job's final mem_mb
# (_star_sort_ram), reserving a quarter (at least 8 GB) for STARsolo and process
# overhead, so a set-resources override keeps the two in sync.
def _dir_size_mb(path):
  try:
    files = [os.path.join(path, f) for f in os.listdir(path)]
  except OSError:
    return None
  return sum(os.path.getsize(f) for f in files if os.path.isfile(f)) // 1000000

STAR_INDEX_MB = _dir_size_mb(STAR_INDEX) or 32000

def _star_mem_mb(wildcards, attempt):
  mem = _get_config(wildcards.sample, "STAR_mem_mb")
  if mem in ("", None):
    mem = max(STAR_INDEX_MB + 8000, 48000)
  return int(mem) * attempt

def _star_sort_ram(wildcards, resources):
  """--limitBAMsortRAM in bytes: job memory minus the STARsolo/overhead reserve."""
  if resources.mem_mb < STAR_INDEX_MB + 4000:
    _diag(("star_mem", wildcards.sample),
          "Warning: sample '{}': STARsolo mem_mb {} is under the {} MB genome "
          "index + 4 GB needed for mapping; raise STAR_mem_mb.".format(
            wildcards.sample, resources.mem_mb, STAR_INDEX_MB))
  sort_mb = resources.mem_mb - max(8000, resources.mem_mb // 4)
  return max(sort_mb, 2000) * 1000000

# --- per-mode behaviour switches ----------------------------------------------
# Per-mode keys ({mode} = R2 or bam), settable per sample / platform / chemistry
# with defaults in config.yaml DEFAULTS:
#   {mode}_polyA_filter   (bool) - require a poly(A)/poly(T) soft clip at the
#                                  counted read end (filter_bam.py). Default
#                                  true for R2, false for bam.
#   {mode}_priming_offset (int)  - bp the priming site lies BEYOND the counted
#                                  read end (transcript-3' direction), applied in
#                                  featureCounts (--readShiftType) and the
#                                  discovery pileup (pileup_sites.py --shift).
#                                  Ultima FASTQs: 5. Default 0.
# bam-mode only (pre-aligned BAM input, see rules/count.snake filter_bam):
#   bam        (path) - input BAM (absolute, or relative to DATA). Required.
#   bam_mate   (0|1|2, default 0) - 0: single-end BAM; 1/2: paired-end BAM, the
#                                   mate carrying the priming site (kept from
#                                   proper pairs, unpaired; the other mate dropped)
#   bam_strand (same|reverse, default same) - read (or selected mate) strand vs
#                                   transcript: same -> count the read 3' end
#                                   (featureCounts -s 1 --read2pos 3, like R2);
#                                   reverse -> count the read 5' end (-s 2
#                                   --read2pos 5, like R1/paired)
# The _get_config fallthrough ("" for configs predating these keys) must resolve
# to the defaults, and YAML booleans must not reach the shell as "False".
def _as_bool(val, default):
  """YAML/CLI-robust truthiness: bool, '' (-> default), true/false, yes/no,
  on/off, 1/0. Anything else is a config error."""
  if isinstance(val, bool):
    return val
  if val is None:  # YAML null (key present, no value)
    return default
  s = str(val).strip().lower()
  if s == "":
    return default
  if s in ("true", "yes", "on", "1"):
    return True
  if s in ("false", "no", "off", "0"):
    return False
  raise ValueError("expected a boolean, got " + repr(val))

def _polya_filter(sample, mode, default):
  """'true' / 'false' for the shell test in filter_{mode}: whether a
  poly(A)/poly(T) soft clip is required at the counted read end."""
  val = _get_config(sample, mode + "_polyA_filter")
  return "true" if _as_bool(val, default) else "false"

def _bam_strand(sample):
  """'same' | 'reverse': strand of the counted read (or selected mate) relative
  to the transcript, bam mode only. '' / None -> 'same'."""
  val = _get_config(sample, "bam_strand")
  s = "" if val is None else str(val).strip().lower()
  if s == "":
    return "same"
  if s in ("same", "reverse"):
    return s
  raise ValueError(
    "sample '{}': bam_strand must be 'same' or 'reverse', got {}".format(sample, val)
  )

def _bam_end(sample):
  """Counted read end for bam mode: '3' when the read is on the transcript
  strand (same), '5' when it is antisense (reverse)."""
  return "3" if _bam_strand(sample) == "same" else "5"

def _bam_mate(sample):
  """0 (single-end BAM) or the mate (1|2) carrying the priming site in a
  paired-end BAM. '' / None -> 0."""
  val = _get_config(sample, "bam_mate")
  if val is None or str(val).strip() == "":
    return 0
  try:
    m = int(val)
  except (TypeError, ValueError):
    m = -1
  if m not in (0, 1, 2):
    raise ValueError(
      "sample '{}': bam_mate must be 0, 1 or 2, got {}".format(sample, val)
    )
  return m

# modes whose site-assignment rule applies the priming offset in featureCounts.
# Other modes (R1, paired) have no shift wired into assign_sites_*, so their
# offset is ignored everywhere (incl. discovery) to keep counts and beds in sync.
_OFFSET_MODES = ("R2", "bam")

def _priming_offset(sample, mode):
  """Non-negative int bp shift of the {mode} priming end toward transcript 3'.
  Looks up '{mode}_priming_offset'; '' (unset) -> 0. Only R2/bam are honoured."""
  val = _get_config(sample, mode + "_priming_offset")
  if val is None or str(val).strip() == "":
    return 0
  if mode not in _OFFSET_MODES:
    return 0  # warned once at parse time by _validate_offsets()
  n = int(val)
  if n < 0:
    raise ValueError(
      "sample '{}': {}_priming_offset must be >= 0, got {}".format(sample, mode, n)
    )
  return n

# featureCounts --readShiftType is read-strand-aware. A priming offset always
# moves the counted end OUTWARD (past the read end): 'downstream' (toward the
# read 3') for modes that count the 3' end (R2, bam/same) and 'upstream' for
# modes that count the 5' end (R1, paired, bam/reverse). pileup_sites.py --shift
# uses the same outward convention. Only R2 and bam honour an offset (see
# _OFFSET_MODES); the discovery stranded_bed rule uses the same helper.
def _fc_shift_type(sample, mode):
  if mode == "bam":
    return "downstream" if _bam_strand(sample) == "same" else "upstream"
  return {"R2": "downstream", "R1": "upstream", "paired": "upstream"}[mode]

def _fc_shift_args(sample, mode):
  """featureCounts shift flags for {mode}, or '' when the offset is 0."""
  n = _priming_offset(sample, mode)
  if n == 0:
    return ""
  return "--readShiftType {} --readShiftSize {}".format(_fc_shift_type(sample, mode), n)

def _sample_modes(sample):
  modes = _get_config(sample, "alignments") or []
  return [modes] if isinstance(modes, str) else list(modes)

def _validate_bam_inputs():
  """Fail fast for bam-mode samples: the input `bam` key must be present and the
  bam_strand / bam_mate switches must be valid (the helpers raise)."""
  errors = []
  for sample in SAMPLES:
    if "bam" not in _sample_modes(sample):
      continue
    cfg = SAMPLES[sample] or {}
    if not cfg.get("bam"):
      errors.append(
        "sample '{}': alignment mode 'bam' requires a per-sample 'bam' key "
        "(input BAM path, absolute or relative to DATA)".format(sample)
      )
    for fn in (_bam_strand, _bam_mate):
      try:
        fn(sample)
      except ValueError as e:
        errors.append(str(e))
  if errors:
    raise ValueError(
      "Invalid bam-mode configuration:\n  - " + "\n  - ".join(errors)
    )

_validate_bam_inputs()

def _validate_offsets():
  """Warn once, at parse time, when a sample sets a priming offset for a mode
  whose assign_sites rule applies none (lazy params are only evaluated for
  scheduled jobs, so this is the one reliable place to say so)."""
  for sample in SAMPLES:
    for mode in ("R1", "paired"):
      val = _get_config(sample, mode + "_priming_offset")
      if val is not None and str(val).strip() not in ("", "0"):
        logger.warning(
          "Warning: {}_priming_offset is set for sample '{}' but only {} "
          "support a priming offset; it is ignored.".format(
            mode, sample, "/".join(_OFFSET_MODES)))

_validate_offsets()

# assemble outputs for rule all. The bam mode consumes a user-supplied BAM, so it
# yields no *_Aligned.sortedByCoord.out.bam of its own.
SAMPLE_OUTS = []
for x in SAMPLES:
  for mode in _sample_modes(x):
    SAMPLE_OUTS.append("{}/counts/{}_{}_counts.tsv.gz".format(RESULTS, x, mode))
    if mode != "bam":
      SAMPLE_OUTS.append("{}/{}/{}_{}_Aligned.sortedByCoord.out.bam".format(RESULTS, x, x, mode))

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
  """Return {unit_id: [(sample, mode), ...]} for all discovery tracks.

  Unit ids must be unique across the three sources (explicit DISCOVERY.units,
  "{group}_{mode}", "{sample}_{mode}"); a collision is a config error rather
  than a silent overwrite. Explicit (sample, mode) members must name a mode the
  sample is actually aligned in (its `alignments`), otherwise discovery would
  schedule an alignment the user never asked for.
  """
  units = {}
  origin = {}
  errors = []
  groups = DISCOVERY.get("groups", {}) or {}
  explicit = DISCOVERY.get("units", {}) or {}
  grouped_samples = set()
  for members in groups.values():
    grouped_samples.update(members)

  def _add(name, pairs, source):
    if name in units:
      errors.append(
        "discovery unit id '{}' from {} collides with the {} unit of the same "
        "name; rename the explicit unit or group".format(name, source, origin[name]))
      return
    units[name] = pairs
    origin[name] = source

  # explicit units: arbitrary (sample, mode) aggregation (cross-sample/mode)
  for name, members in explicit.items():
    pairs = []
    for m in members:
      s, mode = (m[0], m[1]) if not isinstance(m, str) else (m, None)
      if s not in SAMPLES:
        errors.append("DISCOVERY.units '{}': unknown sample '{}'".format(name, s))
        continue
      aligned = _as_mode_list(_get_config(s, "alignments")) or []
      if mode:
        if mode not in aligned:
          errors.append(
            "DISCOVERY.units '{}': sample '{}' is not aligned in mode '{}' "
            "(alignments: {})".format(name, s, mode, aligned))
          continue
        pairs.append((s, mode))
      else:
        pairs.extend((s, md) for md in _sample_disc_modes(s))
    if pairs:
      _add(name, pairs, "DISCOVERY.units")

  # groups: pool member samples within each shared mode -> "{group}_{mode}"
  for gname, members in groups.items():
    unknown = [s for s in members if s not in SAMPLES]
    if unknown:
      errors.append("DISCOVERY.groups '{}': unknown sample(s) {}".format(gname, unknown))
      continue
    # modes shared by all members (intersection, order from first member)
    per = [ _sample_disc_modes(s) for s in members ]
    if not per:
      continue
    shared = [m for m in per[0] if all(m in p for p in per[1:])]
    for md in shared:
      _add("{}_{}".format(gname, md), [(s, md) for s in members], "DISCOVERY.groups")

  # ungrouped samples: one unit per (sample, mode) -> "{sample}_{mode}"
  for s in SAMPLES:
    if s in grouped_samples:
      continue
    for md in _sample_disc_modes(s):
      _add("{}_{}".format(s, md), [(s, md)], "per-sample")

  if errors:
    raise ValueError(
      "Invalid DISCOVERY configuration:\n  - " + "\n  - ".join(errors))
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
# via the logger (never print()): stdout must stay clean for `snakemake --dag | dot`
logger.info("targets: " + " ".join(all_outputs))

rule all:
  input:
    all_outputs = all_outputs

include: "rules/check_versions.snake"
include: "rules/cutadapt_star.snake"
include: "rules/count.snake"
include: "rules/qc.snake"
include: "rules/discovery.snake"
