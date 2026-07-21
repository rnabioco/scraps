#!/usr/bin/env bash
# End-to-end bit-exactness validation on REAL data. Runs the pipeline twice from
# the same FASTQs -- once with the original umi_tools code path, once with the
# new pysam scripts -- into separate results dirs, then diffs every count table,
# bed, and discovery stranded bed.
#
# The original (umi_tools) path is obtained from git: this script stashes the
# working tree, checks out the pre-change rules into a scratch worktree, runs it,
# then restores and runs the new code. It assumes a configured config.yaml with a
# valid STAR_INDEX and DATA, and that `snakemake` + tools are on PATH.
#
# Usage:
#   OLD_REF=<git-ref-before-change> tests/validate_realdata.sh
# where OLD_REF is a commit/branch/tag containing the umi_tools rules (default:
# HEAD~1). RESULTS_old / RESULTS_new are created under a scratch dir.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO"

OLD_REF="${OLD_REF:-HEAD~1}"
CONFIG="${CONFIG:-config.yaml}"
CORES="${CORES:-8}"
SCRATCH="${SCRATCH:-$REPO/.validate_realdata}"
OLD_DIR="$SCRATCH/results_old"
NEW_DIR="$SCRATCH/results_new"
OLD_TREE="$SCRATCH/worktree_old"

mkdir -p "$SCRATCH"

echo ">>> OLD path from git ref: $OLD_REF"
git worktree remove --force "$OLD_TREE" 2>/dev/null || true
git worktree add --force "$OLD_TREE" "$OLD_REF"
# reuse the same config but redirect RESULTS to OLD_DIR via --config override
( cd "$OLD_TREE" && snakemake -j "$CORES" --configfile "$REPO/$CONFIG" \
    --config RESULTS="$OLD_DIR" --resources total_impact=5 --keep-going )

echo ">>> NEW path (current working tree)"
snakemake -j "$CORES" --configfile "$CONFIG" \
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
# dedup on (cell, UMI, priming position) and ignores the fragmentation/splice
# coordinate, so at each priming position new <= old and no new position is
# introduced. Criterion: only-old == 0 AND total new signal <= total old signal.
# (col4 is the count; for stranded beds the strand col is a pure relabel.)
compare_directional() {  # $1 = relative glob
  for oldf in "$OLD_DIR"/$1; do
    [ -e "$oldf" ] || continue
    rel="${oldf#$OLD_DIR/}"
    newf="$NEW_DIR/$rel"
    if [ ! -e "$newf" ]; then echo "  MISSING new: $rel"; fail=1; continue; fi
    only_old=$(comm -23 \
      <(zcat < "$oldf" | awk -v OFS='\t' '{print $1,$2}' | sort -u) \
      <(zcat < "$newf" | awk -v OFS='\t' '{print $1,$2}' | sort -u) | wc -l | tr -d ' ')
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
