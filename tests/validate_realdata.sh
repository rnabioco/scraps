#!/usr/bin/env bash
# End-to-end validation on REAL data. Runs the pipeline twice from the same
# FASTQs -- once with the original umi_tools code path, once with the new pysam
# scripts -- into separate results dirs, then compares count tables (bit-exact)
# and discovery stranded beds (directional: every legacy position still present,
# total signal <= old). Count-arm bed tracks no longer exist and are not compared.
#
# The original (umi_tools) path is obtained from git: the old ref is checked out
# into a scratch git worktree (the current working tree is untouched) and run
# there; the new code runs in the current tree. Both runs use the snakemake on
# PATH (Snakemake 9; the legacy commit d335350 dry-runs fine under it). The
# umi_tools CLI must be installed for the legacy side. Assumptions:
#   - a configured config.yaml with a valid STAR_INDEX and DATA, and `snakemake`
#     + tools on PATH;
#   - DISCOVERY.enabled: true in that config, otherwise no stranded beds exist
#     and the bed comparison silently matches nothing;
#   - no bam-mode / Ultima-FASTQ / data-less fixture samples in SAMPLES (the old
#     ref cannot run them); the default config.yaml contains such fixtures.
#
# Usage:
#   OLD_REF=<git-ref-before-change> tests/validate_realdata.sh
# where OLD_REF is a commit/branch/tag containing the umi_tools rules. Default:
# d335350, the last commit before the pysam replacement (05cbc3f).
# RESULTS_old / RESULTS_new are created under a scratch dir.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"

OLD_REF="${OLD_REF:-d335350}"
CONFIG="${CONFIG:-config.yaml}"
CORES="${CORES:-8}"
# fresh result directories per invocation: the count/pileup helper scripts are
# invoked from shell commands, not declared as Snakemake inputs, so reusing a
# results dir would let an edited script "pass" on stale outputs
RUN_ID="$(date +%Y%m%d_%H%M%S)"
SCRATCH="${SCRATCH:-$REPO/.validate_realdata}/$RUN_ID"
OLD_DIR="$SCRATCH/results_old"
NEW_DIR="$SCRATCH/results_new"
OLD_TREE="$SCRATCH/worktree_old"
RESOLVED_CONFIG="$SCRATCH/config.resolved.yaml"

mkdir -p "$SCRATCH"
echo ">>> scratch: $SCRATCH"

# Resolve the config's resource paths (DATA, STAR_INDEX, POLYA_SITES,
# GENOME_FASTA) to absolute paths so both runs see the same files even though
# the legacy run executes inside the old worktree. chemistry.yaml whitelist
# paths are relative (ref/...), so missing ref/ entries are linked into the
# worktree from this repository as well.
python3 - "$REPO" "$CONFIG" "$RESOLVED_CONFIG" <<'EOF'
import os, sys, yaml
repo, cfg_in, cfg_out = sys.argv[1:4]
cfg = yaml.safe_load(open(os.path.join(repo, cfg_in) if not os.path.isabs(cfg_in) else cfg_in))
for key in ("DATA", "STAR_INDEX", "POLYA_SITES", "GENOME_FASTA"):
    val = cfg.get(key)
    if isinstance(val, str) and val and not os.path.isabs(val):
        cfg[key] = os.path.abspath(os.path.join(repo, val))
yaml.safe_dump(cfg, open(cfg_out, "w"), sort_keys=False)
EOF

echo ">>> OLD path from git ref: $OLD_REF"
git worktree remove --force "$OLD_TREE" 2>/dev/null || true
git worktree add --force "$OLD_TREE" "$OLD_REF"
mkdir -p "$OLD_TREE/ref"
for f in "$REPO"/ref/*; do
  [ -e "$OLD_TREE/ref/$(basename "$f")" ] || ln -s "$f" "$OLD_TREE/ref/$(basename "$f")"
done
( cd "$OLD_TREE" && snakemake -j "$CORES" --configfile "$RESOLVED_CONFIG" \
    --config RESULTS="$OLD_DIR" --resources total_impact=5 --keep-going )

echo ">>> NEW path (current working tree)"
snakemake -j "$CORES" --configfile "$RESOLVED_CONFIG" \
    --config RESULTS="$NEW_DIR" --resources total_impact=5 --keep-going

echo ">>> DIFF"
fail=0
# count tables: must be BIT-EXACT vs umi_tools count (directional reproduced).
compare_exact() {  # $1 = relative glob under results dir
  for oldf in "$OLD_DIR"/$1; do
    [ -e "$oldf" ] || continue
    rel="${oldf#$OLD_DIR/}"
    newf="$NEW_DIR/$rel"
    if [ ! -e "$newf" ]; then echo "  MISSING new: $rel"; fail=1; continue; fi
    if diff <(zcat < "$oldf" | sort) <(zcat < "$newf" | sort) >/dev/null; then
      echo "  OK(exact)   $rel"
    else
      echo "  DIFF        $rel"; fail=1
    fi
  done
}

# beds: DIRECTIONAL vs legacy umi_tools dedup + genomecov. The new pileup keys
# dedup on (cell, UMI, priming position); umi_tools keyed on the mapping
# coordinate and kept ONE read per (cell, UMI, start), whose priming end may
# differ from its duplicates'. So every legacy priming position must still be
# present in the new output (only-old == 0), while the new output may contain
# positions the legacy one lacked and per-position counts are not comparable;
# total signal must not grow (fragmentation-split duplicates collapse).
# Keys are (chrom, pos, strand) for stranded (5-col) beds.
compare_directional() {  # $1 = relative glob
  for oldf in "$OLD_DIR"/$1; do
    [ -e "$oldf" ] || continue
    rel="${oldf#$OLD_DIR/}"
    newf="$NEW_DIR/$rel"
    if [ ! -e "$newf" ]; then echo "  MISSING new: $rel"; fail=1; continue; fi
    only_old=$(comm -23 \
      <(zcat < "$oldf" | awk -v OFS='\t' '{print $1,$2,(NF>=5?$5:"")}' | sort -u) \
      <(zcat < "$newf" | awk -v OFS='\t' '{print $1,$2,(NF>=5?$5:"")}' | sort -u) | wc -l | tr -d ' ')
    old_sum=$(zcat < "$oldf" | awk '{s+=$4}END{print s+0}')
    new_sum=$(zcat < "$newf" | awk '{s+=$4}END{print s+0}')
    if [ "$only_old" -eq 0 ] && [ "$new_sum" -le "$old_sum" ]; then
      echo "  OK(dir)     $rel  (only_old=0, new=$new_sum <= old=$old_sum)"
    else
      echo "  VIOLATION   $rel  (only_old=$only_old, new=$new_sum, old=$old_sum)"; fail=1
    fi
  done
}

compare_exact "counts/*.tsv.gz"
# count-arm bed/*.bed.gz outputs were removed (unstranded pileups are less
# informative for polyA data than the strand-aware discovery beds).
compare_directional "discovery/beds/*.stranded.bed.gz"

git worktree remove --force "$OLD_TREE" || true

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "REAL-DATA VALIDATION PASSED (counts exact; beds directional)"
else
  echo "REAL-DATA VALIDATION FAILED"; exit 1
fi
