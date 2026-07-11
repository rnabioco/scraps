#!/usr/bin/env bash
# Bit-exactness validation for the pysam replacements of umi_tools, on synthetic
# aligned BAMs (all read modes). For each mode:
#   1. run REAL featureCounts (as assign_sites_* does) -> assigned BAM
#   2. count path:  umi_tools count            vs  count_sites.py
#   3. bed path:    umi_tools dedup + genomecov vs  pileup_sites.py
# Exits non-zero on any diff. Requires the scraps_conda tools on PATH.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SAF="${SAF:-$REPO/ref/polyadb32.hg38.saf.gz}"
WORK="${WORK:-$(mktemp -d)}"
echo "workdir: $WORK"
echo "saf:     $SAF"

python3 "$REPO/tests/generate_synthetic.py" --saf "$SAF" --outdir "$WORK"

fail=0

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

  # old bed path: replicate the exact samtools prefilter + genomecov + awk
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

  if diff <(zcat < "$WORK/${mode}_bed_old.bed.gz" | sort) \
          <(zcat < "$WORK/${mode}_bed_new.bed.gz" | sort) > "$WORK/${mode}_bed.diff"; then
    echo "  BED:   IDENTICAL"
  else
    echo "  BED:   DIFFERENT -> $WORK/${mode}_bed.diff"; fail=1
  fi

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

  if diff <(zcat < "$WORK/${mode}_sbed_old.bed.gz" | sort) \
          <(zcat < "$WORK/${mode}_sbed_new.bed.gz" | sort) > "$WORK/${mode}_sbed.diff"; then
    echo "  SBED:  IDENTICAL"
  else
    echo "  SBED:  DIFFERENT -> $WORK/${mode}_sbed.diff"; fail=1
  fi
done

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "ALL SYNTHETIC CHECKS PASSED (bit-exact)"
else
  echo "SOME CHECKS FAILED"; exit 1
fi
