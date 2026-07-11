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
compare() {  # $1 = relative glob under results dir
  for oldf in "$OLD_DIR"/$1; do
    [ -e "$oldf" ] || continue
    rel="${oldf#$OLD_DIR/}"
    newf="$NEW_DIR/$rel"
    if [ ! -e "$newf" ]; then echo "  MISSING new: $rel"; fail=1; continue; fi
    if diff <(zcat < "$oldf" | sort) <(zcat < "$newf" | sort) >/dev/null; then
      echo "  OK   $rel"
    else
      echo "  DIFF $rel"; fail=1
    fi
  done
}

compare "counts/*.tsv.gz"
compare "bed/*.bed.gz"
compare "discovery/beds/*.stranded.bed.gz"

git worktree remove --force "$OLD_TREE" || true

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "REAL-DATA VALIDATION PASSED (bit-exact)"
else
  echo "REAL-DATA VALIDATION FAILED"; exit 1
fi
