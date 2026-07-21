#!/usr/bin/env bash
# Validation for the pysam replacements of umi_tools, on synthetic aligned BAMs
# (all read modes). For each mode:
#   1. run REAL featureCounts (as assign_sites_* does) -> assigned BAM
#   2. count path:  umi_tools count            vs  count_sites.py   [bit-exact]
#   3. bed path:    umi_tools dedup + genomecov vs  pileup_sites.py
#   4. stranded bed (discovery)                vs  pileup_sites.py --strand-split
#
# Correctness criteria differ by path, matching the agreed semantics:
#   - count (all modes): BIT-EXACT vs umi_tools (directional collapse reproduced).
#   - bed + sbed (ALL modes R1/R2/paired): DIRECTIONAL. The new pileup keys dedup
#     on (cell, UMI, priming position) and ignores the opposite-end / fragment /
#     splice coordinate (correct for post-amplification fragmentation, consistent
#     with the count table), whereas umi_tools keys dedup on the read mapping
#     coordinate. So at each priming position the new output must (a) introduce
#     NO position absent from old (only-old == 0) and (b) have total signal <=
#     old total signal (fragmentation-split PCR duplicates collapse to one).
# Exits non-zero on any violation. Requires the scraps_conda tools on PATH.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SAF="${SAF:-$REPO/ref/polyadb32.hg38.saf.gz}"
WORK="${WORK:-$(mktemp -d)}"
echo "workdir: $WORK"
echo "saf:     $SAF"

python3 "$REPO/tests/generate_synthetic.py" --saf "$SAF" --outdir "$WORK"

fail=0

# Directional check for a bed: every OLD position must appear in NEW
# (only-old == 0) and total NEW signal must not exceed total OLD signal (new
# collapses fragmentation-split PCR duplicates). Accepts an optional label.
# Args: old.bed.gz new.bed.gz LABEL
check_directional_bed() {
  local old="$1" new="$2" label="$3"
  local only_old old_sum new_sum
  only_old=$(comm -23 \
    <(zcat < "$old" | awk -v OFS='\t' '{print $1,$2}' | sort -u) \
    <(zcat < "$new" | awk -v OFS='\t' '{print $1,$2}' | sort -u) | wc -l | tr -d ' ')
  old_sum=$(zcat < "$old" | awk '{s+=$4}END{print s+0}')
  new_sum=$(zcat < "$new" | awk '{s+=$4}END{print s+0}')
  if [ "$only_old" -eq 0 ] && [ "$new_sum" -le "$old_sum" ]; then
    echo "  $label:  OK (directional: only_old=0, new_sum=$new_sum <= old_sum=$old_sum)"
    return 0
  fi
  echo "  $label:  VIOLATION (only_old=$only_old, new_sum=$new_sum, old_sum=$old_sum)"
  return 1
}

# featureCounts strandedness / end per mode (matches count.snake assign_sites_*)
declare -A FC_S=( [R1]=2 [R2]=1 [paired]=2 )
declare -A FC_POS=( [R1]=5 [R2]=3 [paired]=5 )
declare -A FC_EXTRA=( [R1]="" [R2]="" [paired]="-p" )
# bed end / mode
declare -A BED_END=( [R1]=5 [R2]=3 [paired]=5 )
# genomecov flag + samtools prefilter for the OLD bed path
declare -A GC_FLAG=( [R1]="-5" [R2]="-3" [paired]="-5" )

for mode in R2 R1 paired; do
  echo "=================== mode: $mode ==================="
  aln="$WORK/synth_${mode}_aligned.bam"

  # 1. real featureCounts -> assigned BAM
  featureCounts -s "${FC_S[$mode]}" -Q 10 -O ${FC_EXTRA[$mode]} \
    --read2pos "${FC_POS[$mode]}" \
    -o "$WORK/${mode}_fc" -F SAF -a "$SAF" -R BAM -T 4 "$aln" \
    2> "$WORK/${mode}_fc.log"
  samtools sort "$aln.featureCounts.bam" -o "$WORK/${mode}_assigned.bam"
  samtools index "$WORK/${mode}_assigned.bam"
  assigned="$WORK/${mode}_assigned.bam"

  # 2. COUNT path
  umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell \
    --extract-umi-method=tag --umi-tag=UB --cell-tag=CB \
    -I "$assigned" -S "$WORK/${mode}_count_old.tsv.gz" 2> "$WORK/${mode}_count_old.log"
  python3 "$REPO/inst/scripts/count_sites.py" \
    -i "$assigned" -o "$WORK/${mode}_count_new.tsv.gz" -t 4

  if diff <(zcat < "$WORK/${mode}_count_old.tsv.gz" | sort) \
          <(zcat < "$WORK/${mode}_count_new.tsv.gz" | sort) > "$WORK/${mode}_count.diff"; then
    echo "  COUNT: IDENTICAL"
  else
    echo "  COUNT: DIFFERENT -> $WORK/${mode}_count.diff"; fail=1
  fi

  # 3. BED path
  paired_flag=""; [ "$mode" = "R1" ] && paired_flag="--paired"
  [ "$mode" = "paired" ] && paired_flag="--paired"
  umi_tools dedup --extract-umi-method=tag --umi-tag=UB --per-cell --cell-tag=CB \
    $paired_flag --method=unique \
    -I "$assigned" -S "$WORK/${mode}_dedup.bam" 2> "$WORK/${mode}_dedup.log"
  samtools index "$WORK/${mode}_dedup.bam"

  # 4-COL UNIT CHECK: the count arm no longer emits per-sample bed tracks, but
  # pileup_sites.py retains its non-strand-split (4-col) code path. Exercise it
  # here against the legacy genomecov 4-col output so the branch stays covered.
  # This is a unit check of the script, not a pipeline-output check.
  case "$mode" in
    R1)     samtools view -F 4 -b "$WORK/${mode}_dedup.bam" > "$WORK/${mode}_pre.bam" ;;
    paired) samtools view -f 0x42 -b "$WORK/${mode}_dedup.bam" > "$WORK/${mode}_pre.bam" ;;
    R2)     cp "$WORK/${mode}_dedup.bam" "$WORK/${mode}_pre.bam" ;;
  esac
  bedtools genomecov -ibam "$WORK/${mode}_pre.bam" ${GC_FLAG[$mode]} -dz \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' | gzip -c > "$WORK/${mode}_bed_old.bed.gz"

  python3 "$REPO/inst/scripts/pileup_sites.py" \
    -i "$assigned" -o "$WORK/${mode}_bed_new.bed.gz" \
    --end "${BED_END[$mode]}" --mode "$mode"

  check_directional_bed "$WORK/${mode}_bed_old.bed.gz" "$WORK/${mode}_bed_new.bed.gz" "BED(4col unit)" || fail=1

  # 4. STRANDED bed (discovery stranded_bed path). Uses labels +/- here; the
  # real rule may relabel to transcript strand, but that is a pure per-strand
  # rename that does not affect the dedup/pileup equivalence being tested.
  ( bedtools genomecov -ibam "$WORK/${mode}_pre.bam" ${GC_FLAG[$mode]} -dz -strand + \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3,"+"}'
    bedtools genomecov -ibam "$WORK/${mode}_pre.bam" ${GC_FLAG[$mode]} -dz -strand - \
      | awk 'BEGIN{OFS="\t"}{print $1,$2,$2+1,$3,"-"}' ) \
    | gzip -c > "$WORK/${mode}_sbed_old.bed.gz"
  python3 "$REPO/inst/scripts/pileup_sites.py" \
    -i "$assigned" -o "$WORK/${mode}_sbed_new.bed.gz" \
    --end "${BED_END[$mode]}" --mode "$mode" \
    --strand-split --label-plus + --label-minus -

  # compare per (chrom,pos,count) — strand column (5) is a pure relabel.
  # write to real temp files (the check reads its inputs twice; process
  # substitution FIFOs cannot be re-read).
  zcat < "$WORK/${mode}_sbed_old.bed.gz" | cut -f1-4 | gzip -c > "$WORK/${mode}_sbed_old4.bed.gz"
  zcat < "$WORK/${mode}_sbed_new.bed.gz" | cut -f1-4 | gzip -c > "$WORK/${mode}_sbed_new4.bed.gz"
  check_directional_bed "$WORK/${mode}_sbed_old4.bed.gz" "$WORK/${mode}_sbed_new4.bed.gz" "SBED" || fail=1
done

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "ALL SYNTHETIC CHECKS PASSED (count bit-exact; bed/sbed directional)"
else
  echo "SOME CHECKS FAILED"; exit 1
fi
