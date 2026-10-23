#!/usr/bin/env bash
# Validation for the pysam replacements of umi_tools (and the bam-mode / priming-
# offset / soft-clip machinery), on synthetic aligned BAMs from
# tests/generate_synthetic.py. For each of the 6 cases
# (R2, R1, paired, bam_same, bam_reverse, bam_paired1):
#   0. bam_* only: reproduce rule filter_bam (samtools -F 0x904 [+ mate flags] +
#      CB/UB, filter_bam.py --mate for paired input) and assert one unpaired
#      record per fragment, no secondary/supplementary survivors
#   1. run REAL featureCounts (as assign_sites_* does; -L for bam_*) -> assigned BAM
#   2. count path:  umi_tools count            vs  count_sites.py   [bit-exact]
#   3. 4-col bed:   umi_tools dedup + genomecov vs  pileup_sites.py (unit check of
#      the non-strand-split branch; the pipeline no longer emits 4-col beds)
#   4. stranded bed (discovery)                vs  pileup_sites.py --strand-split
# Then:
#   5. priming offset on R2: featureCounts --readShiftType downstream vs a
#      fixture with known expected assignments (+ --readExtension3 INFO), and
#      pileup_sites.py --shift as a per-strand bijection (5-col; 4-col non-empty)
#   6. the same under featureCounts -L for bam_same / bam_reverse, incl. a
#      body-overlap case proving --read2pos is honoured under -L
#   7. pileup_sites.py --end 5 --shift (outward = read 5' direction)
#   8. filter_bam.py --end 3 / --end 5 soft-clip filter: kept set must be exact
#
# Correctness criteria differ by path, matching the agreed semantics:
#   - count (all cases): BIT-EXACT vs umi_tools (directional collapse reproduced).
#   - bed + sbed (all cases): DIRECTIONAL. The new pileup keys dedup on
#     (cell, UMI, priming position) and ignores the opposite-end / fragment /
#     splice coordinate (correct for post-amplification fragmentation, keyed like
#     the count table), whereas umi_tools keys dedup on the read mapping
#     coordinate. So every position present in the OLD output must still be
#     present in the new one (only-old == 0) and the new total signal must be <=
#     the old total signal (fragmentation-split PCR duplicates collapse to one).
# Exits non-zero on any violation. Requires bash >= 4 (associative arrays;
# macOS /bin/bash 3.2 will not do) and the scraps_conda tools on PATH, including
# the umi_tools CLI and bedtools (test-only dependencies).
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

# featureCounts strandedness / end per mode (matches count.snake assign_sites_*).
# bam_* are the pre-aligned `bam` mode: bam_same = class 1 (3' end), bam_reverse
# = class 3 (5' end), bam_paired1 = class 3 on mate 1 of a paired BAM. They use
# featureCounts -L (long-read mode, single thread) exactly like assign_sites_bam.
declare -A FC_S=( [R1]=2 [R2]=1 [paired]=2 [bam_same]=1 [bam_reverse]=2 [bam_paired1]=2 )
declare -A FC_POS=( [R1]=5 [R2]=3 [paired]=5 [bam_same]=3 [bam_reverse]=5 [bam_paired1]=5 )
declare -A FC_EXTRA=( [R1]="" [R2]="" [paired]="-p" [bam_same]="-L" [bam_reverse]="-L" [bam_paired1]="-L" )
declare -A FC_T=( [R1]=4 [R2]=4 [paired]=4 [bam_same]=1 [bam_reverse]=1 [bam_paired1]=1 )
# bed end / mode, and the --mode passed to pileup_sites.py
declare -A BED_END=( [R1]=5 [R2]=3 [paired]=5 [bam_same]=3 [bam_reverse]=5 [bam_paired1]=5 )
declare -A PILEUP_MODE=( [R1]=R1 [R2]=R2 [paired]=paired [bam_same]=bam [bam_reverse]=bam [bam_paired1]=bam )
# filter_bam prefilter for bam_* modes: mate selection (0 = single-end)
declare -A BAM_MATE=( [bam_same]=0 [bam_reverse]=0 [bam_paired1]=1 )
# genomecov flag + samtools prefilter for the OLD bed path
declare -A GC_FLAG=( [R1]="-5" [R2]="-3" [paired]="-5" [bam_same]="-3" [bam_reverse]="-5" [bam_paired1]="-5" )

for mode in R2 R1 paired bam_same bam_reverse bam_paired1; do
  echo "=================== mode: $mode ==================="
  aln="$WORK/synth_${mode}_aligned.bam"

  # 0. bam modes: reproduce rule filter_bam (samtools -F 0x904 [+ mate flags] +
  #    CB/UB, then filter_bam.py --mate M --no-clip-filter for paired input) and
  #    feed the result to BOTH paths. Asserts the per-fragment reduction.
  case "$mode" in
    bam_*)
      mate="${BAM_MATE[$mode]}"
      mflags=""; [ "$mate" = 1 ] && mflags="-f 0x42"; [ "$mate" = 2 ] && mflags="-f 0x82"
      samtools view -b -F 0x904 $mflags -e '[CB]!="-" && [UB]!="-"' "$aln" > "$WORK/${mode}_pre.bam"
      if [ "$mate" != 0 ]; then
        python3 "$REPO/inst/scripts/filter_bam.py" --mate "$mate" --end "${BED_END[$mode]}" \
          --no-clip-filter -i "$WORK/${mode}_pre.bam" -o "$WORK/${mode}_filtered.bam" > "$WORK/${mode}_filter.log"
        # one UNPAIRED record per proper-pair fragment of the selected mate
        n_paired=$(samtools view -c -f 0x1 "$WORK/${mode}_filtered.bam")
        if [ "$n_paired" -ne 0 ]; then echo "  FILTER(mate $mate): FAIL ($n_paired records still flagged paired)"; fail=1; fi
        if diff <(samtools view "$WORK/${mode}_filtered.bam" | cut -f1 | sort) \
                <(samtools view -F 0x904 $mflags "$aln" | cut -f1 | sort) > /dev/null; then
          echo "  FILTER(mate $mate): OK (one single-end record per proper-pair fragment)"
        else
          echo "  FILTER(mate $mate): FAIL (qname set differs from proper-pair mate-$mate set)"; fail=1
        fi
      else
        mv "$WORK/${mode}_pre.bam" "$WORK/${mode}_filtered.bam"
      fi
      n_bad=$(samtools view -c -f 0x900 "$WORK/${mode}_filtered.bam")
      if [ "$n_bad" -ne 0 ]; then echo "  FILTER: FAIL ($n_bad secondary/supplementary records survived)"; fail=1; fi
      samtools sort -o "$WORK/${mode}_filtered.sorted.bam" "$WORK/${mode}_filtered.bam"
      samtools index "$WORK/${mode}_filtered.sorted.bam"
      aln="$WORK/${mode}_filtered.sorted.bam"
      ;;
  esac

  # 1. real featureCounts -> assigned BAM
  featureCounts -s "${FC_S[$mode]}" -Q 10 -O ${FC_EXTRA[$mode]} \
    --read2pos "${FC_POS[$mode]}" \
    -o "$WORK/${mode}_fc" -F SAF -a "$SAF" -R BAM -T "${FC_T[$mode]}" "$aln" \
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
  [ "$mode" = "paired" ] && paired_flag="--paired"   # bam_*: single-end after filter_bam
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
    *)      cp "$WORK/${mode}_dedup.bam" "$WORK/${mode}_pre.bam" ;;   # R2, bam_* (already single-end)
  esac
  bedtools genomecov -ibam "$WORK/${mode}_pre.bam" ${GC_FLAG[$mode]} -dz \
    | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' | gzip -c > "$WORK/${mode}_bed_old.bed.gz"

  python3 "$REPO/inst/scripts/pileup_sites.py" \
    -i "$assigned" -o "$WORK/${mode}_bed_new.bed.gz" \
    --end "${BED_END[$mode]}" --mode "${PILEUP_MODE[$mode]}"

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
    --end "${BED_END[$mode]}" --mode "${PILEUP_MODE[$mode]}" \
    --strand-split --label-plus + --label-minus -

  # compare per (chrom,pos,count) — strand column (5) is a pure relabel.
  # write to real temp files (the check reads its inputs twice; process
  # substitution FIFOs cannot be re-read).
  zcat < "$WORK/${mode}_sbed_old.bed.gz" | cut -f1-4 | gzip -c > "$WORK/${mode}_sbed_old4.bed.gz"
  zcat < "$WORK/${mode}_sbed_new.bed.gz" | cut -f1-4 | gzip -c > "$WORK/${mode}_sbed_new4.bed.gz"
  check_directional_bed "$WORK/${mode}_sbed_old4.bed.gz" "$WORK/${mode}_sbed_new4.bed.gz" "SBED" || fail=1
done

echo "=================== priming offset (Ultima FASTQ: R2 3' end + 5) ==================="
# The ultima_fastq platform sets R2_priming_offset=5: assign_sites_R2 adds
# `--readShiftType downstream --readShiftSize 5` to featureCounts and stranded_bed
# passes `--shift 5` to pileup_sites.py. Both are checked here.
#   1. featureCounts: synthetic reads at known offsets from a SAF-subset of
#      windows must flip assignment exactly as the qname prefix predicts (see
#      tests/generate_synthetic.py build_shift). This also pins down the shift
#      vs --read2pos order of operations empirically.
#   2. pileup_sites.py --shift 5 must equal the unshifted stranded bed with every
#      position moved +5 (+ strand) / -5 (- strand): a per-strand bijection, so
#      counts must match exactly.
SHIFT=5
sbam="$WORK/synth_shift_R2_aligned.bam"
ssaf="$WORK/shift_windows.saf"

# run featureCounts (R2-arm flags + extras) and emit "<prefix>\t<XS>" per read.
# $sbam, FC_BASE (strand/end flags) and FC_THREADS are set by the caller.
FC_BASE="-s 1 --read2pos 3"; FC_THREADS=4
fc_assign() {  # $1 = label, remaining args = extra featureCounts flags
  local label="$1"; shift
  featureCounts $FC_BASE -Q 10 -O "$@" \
    -o "$WORK/shift_${label}_fc" -F SAF -a "$ssaf" -R BAM -T "$FC_THREADS" "$sbam" \
    2> "$WORK/shift_${label}_fc.log" || return 1
  samtools view "$sbam.featureCounts.bam" \
    | awk -v OFS='\t' '{xs="NA"; for(i=12;i<=NF;i++) if($i ~ /^XS:Z:/) xs=substr($i,6);
                        split($1,a,"_"); print a[1], xs}' \
    | sort > "$WORK/shift_${label}.tsv"
}

# assert every read with a given qname prefix has the expected status
# Args: tsv label prefix=Assigned|Unassigned ...
check_shift_expect() {
  local tsv="$1" label="$2"; shift 2
  local ok=1 spec prefix want n bad
  for spec in "$@"; do
    prefix="${spec%%=*}"; want="${spec#*=}"
    n=$(awk -v p="$prefix" '$1==p' "$tsv" | wc -l | tr -d ' ')
    if [ "$want" = "Assigned" ]; then
      bad=$(awk -v p="$prefix" '$1==p && $2!="Assigned"' "$tsv" | wc -l | tr -d ' ')
    else
      bad=$(awk -v p="$prefix" '$1==p && $2=="Assigned"' "$tsv" | wc -l | tr -d ' ')
    fi
    if [ "$n" -eq 0 ] || [ "$bad" -ne 0 ]; then
      echo "  $label $prefix: FAIL (n=$n, violating=$bad, expected $want)"; ok=0
    else
      echo "  $label $prefix: OK (n=$n, all $want)"
    fi
  done
  [ "$ok" -eq 1 ]
}

# 1a. control: no shift
fc_assign noshift
check_shift_expect "$WORK/shift_noshift.tsv" "FC(no shift)" \
  shiftin=Unassigned inside=Assigned edge=Assigned far=Unassigned || fail=1

# 1b. --readShiftType downstream (what assign_sites_R2 emits)
if fc_assign readshift --readShiftType downstream --readShiftSize "$SHIFT"; then
  check_shift_expect "$WORK/shift_readshift.tsv" "FC(readShift $SHIFT)" \
    shiftin=Assigned inside=Assigned edge=Unassigned far=Unassigned || fail=1
else
  echo "  FC(readShift): FAIL - featureCounts rejected --readShiftType/--readShiftSize" \
       "(subread too old?). Fallback: --readExtension3 $SHIFT (only valid if extension" \
       "precedes --read2pos reduction; see INFO line below)."; fail=1
fi

# 1c. INFO: does the documented --readExtension3 behave identically? (not a failure)
if fc_assign readext --readExtension3 "$SHIFT"; then
  if diff -q "$WORK/shift_readshift.tsv" "$WORK/shift_readext.tsv" > /dev/null 2>&1; then
    echo "  INFO: --readExtension3 $SHIFT gives the SAME assignments as --readShiftType downstream"
  else
    echo "  INFO: --readExtension3 $SHIFT DIFFERS from --readShiftType downstream (not a usable fallback)"
  fi
fi

# 2. pileup_sites.py --shift: per-strand bijection of the unshifted stranded bed
python3 "$REPO/inst/scripts/pileup_sites.py" \
  -i "$WORK/R2_assigned.bam" -o "$WORK/R2_sbed_shift.bed.gz" \
  --end 3 --mode R2 --strand-split --label-plus + --label-minus - --shift "$SHIFT"
zcat < "$WORK/R2_sbed_new.bed.gz" \
  | awk -v s="$SHIFT" 'BEGIN{OFS="\t"}{p=($5=="+")?$2+s:$2-s; print $1,p,p+1,$4,$5}' \
  | sort -k1,1 -k2,2n -k5,5 > "$WORK/R2_sbed_shift_expected.bed"
zcat < "$WORK/R2_sbed_shift.bed.gz" | sort -k1,1 -k2,2n -k5,5 > "$WORK/R2_sbed_shift_got.bed"
if diff "$WORK/R2_sbed_shift_expected.bed" "$WORK/R2_sbed_shift_got.bed" > "$WORK/R2_sbed_shift.diff"; then
  echo "  PILEUP(--shift $SHIFT): IDENTICAL to unshifted bed moved +/-$SHIFT by strand"
else
  echo "  PILEUP(--shift $SHIFT): DIFFERENT -> $WORK/R2_sbed_shift.diff"; fail=1
fi
# 4-col path with --shift (branch coverage only)
python3 "$REPO/inst/scripts/pileup_sites.py" \
  -i "$WORK/R2_assigned.bam" -o "$WORK/R2_bed_shift.bed.gz" --end 3 --mode R2 --shift "$SHIFT"
if [ "$(zcat < "$WORK/R2_bed_shift.bed.gz" | wc -l | tr -d ' ')" -gt 0 ]; then
  echo "  PILEUP(4col --shift $SHIFT): OK (non-empty)"
else
  echo "  PILEUP(4col --shift $SHIFT): FAIL (empty output)"; fail=1
fi

# 3. bam mode under featureCounts -L: long multi-block reads, both orientations.
#    Proves -L + --read2pos + --readShiftType compose (lrbody_ reads overlap the
#    window with their body but their counted end is 200 bp past it).
echo "--- bam mode (-L): same orientation (3' end, shift downstream)"
sbam="$WORK/synth_shift_bam_aligned.bam"; FC_BASE="-L -s 1 --read2pos 3"; FC_THREADS=1
fc_assign lr_noshift
check_shift_expect "$WORK/shift_lr_noshift.tsv" "FC(-L same, no shift)" \
  lrshiftin=Unassigned lrinside=Assigned lredge=Assigned lrfar=Unassigned lrbody=Unassigned || fail=1
if fc_assign lr_readshift --readShiftType downstream --readShiftSize "$SHIFT"; then
  check_shift_expect "$WORK/shift_lr_readshift.tsv" "FC(-L same, readShift $SHIFT)" \
    lrshiftin=Assigned lrinside=Assigned lredge=Unassigned lrfar=Unassigned lrbody=Unassigned || fail=1
else
  echo "  FC(-L same, readShift): FAIL - featureCounts rejected the flags under -L"; fail=1
fi
echo "--- bam mode (-L): reverse orientation (5' end, shift upstream)"
sbam="$WORK/synth_shift_bamrev_aligned.bam"; FC_BASE="-L -s 2 --read2pos 5"; FC_THREADS=1
fc_assign lrrev_noshift
check_shift_expect "$WORK/shift_lrrev_noshift.tsv" "FC(-L reverse, no shift)" \
  lrshiftin=Unassigned lrinside=Assigned lredge=Assigned lrfar=Unassigned lrbody=Unassigned || fail=1
if fc_assign lrrev_readshift --readShiftType upstream --readShiftSize "$SHIFT"; then
  check_shift_expect "$WORK/shift_lrrev_readshift.tsv" "FC(-L reverse, readShift $SHIFT)" \
    lrshiftin=Assigned lrinside=Assigned lredge=Unassigned lrfar=Unassigned lrbody=Unassigned || fail=1
else
  echo "  FC(-L reverse, readShift): FAIL - featureCounts rejected the flags under -L"; fail=1
fi

# 4. pileup --end 5 --shift (outward = read 5' direction): bam_reverse assigned BAM
python3 "$REPO/inst/scripts/pileup_sites.py" \
  -i "$WORK/bam_reverse_assigned.bam" -o "$WORK/bamrev_sbed_shift.bed.gz" \
  --end 5 --mode bam --strand-split --label-plus + --label-minus - --shift "$SHIFT"
zcat < "$WORK/bam_reverse_sbed_new.bed.gz" \
  | awk -v s="$SHIFT" 'BEGIN{OFS="\t"}{p=($5=="+")?$2-s:$2+s; print $1,p,p+1,$4,$5}' \
  | sort -k1,1 -k2,2n -k5,5 > "$WORK/bamrev_sbed_shift_expected.bed"
zcat < "$WORK/bamrev_sbed_shift.bed.gz" | sort -k1,1 -k2,2n -k5,5 > "$WORK/bamrev_sbed_shift_got.bed"
if diff "$WORK/bamrev_sbed_shift_expected.bed" "$WORK/bamrev_sbed_shift_got.bed" > "$WORK/bamrev_sbed_shift.diff"; then
  echo "  PILEUP(--end 5 --shift $SHIFT): IDENTICAL to unshifted bed moved -/+$SHIFT by strand (outward)"
else
  echo "  PILEUP(--end 5 --shift $SHIFT): DIFFERENT -> $WORK/bamrev_sbed_shift.diff"; fail=1
fi

echo "=================== soft-clip filter (filter_bam.py --end 3 / --end 5) ==================="
# qname = "<keep3|drop3>-<keep5|drop5>_<case>_<rep>"; the kept set must match exactly
cbam="$WORK/synth_clip_aligned.bam"
for e in 3 5; do
  python3 "$REPO/inst/scripts/filter_bam.py" --end "$e" -i "$cbam" -o "$WORK/clip_end${e}.bam" > "$WORK/clip_end${e}.log"
  samtools view "$WORK/clip_end${e}.bam" | cut -f1 | sort > "$WORK/clip_end${e}_kept.txt"
  if [ "$e" = 3 ]; then
    samtools view "$cbam" | cut -f1 | grep '^keep3-' | sort > "$WORK/clip_end${e}_expected.txt"
  else
    samtools view "$cbam" | cut -f1 | grep -- '-keep5_' | sort > "$WORK/clip_end${e}_expected.txt"
  fi
  if diff "$WORK/clip_end${e}_expected.txt" "$WORK/clip_end${e}_kept.txt" > "$WORK/clip_end${e}.diff"; then
    echo "  FILTER(--end $e): OK ($(wc -l < "$WORK/clip_end${e}_kept.txt" | tr -d ' ') kept, exactly the expected set)"
  else
    echo "  FILTER(--end $e): FAIL -> $WORK/clip_end${e}.diff"; fail=1
  fi
done

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "ALL SYNTHETIC CHECKS PASSED (count bit-exact; bed/sbed directional; priming offset; bam mode; clip filter)"
else
  echo "SOME CHECKS FAILED"; exit 1
fi
