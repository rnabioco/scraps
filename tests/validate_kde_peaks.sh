#!/usr/bin/env bash
# Regression checks for inst/scripts/kde_peaks.py cluster semantics:
#   A. umi_support is the UMI sum over the whole cluster extent (not only the
#      span between density maxima): a single-maximum 88-UMI pileup reports 88,
#      and an 83-UMI pileup passes the default --min-umi 10 instead of being
#      dropped.
#   B. two clusters joined by sparse single-UMI bases get DISTINCT summits, each
#      inside its own cluster, under both --summit-method values (no duplicate
#      (chrom, summit, strand) rows); their supports add up to the track total.
#   C. an isolated 5-UMI bump is dropped by the support floor and reported in
#      the stderr summary.
# Exits non-zero on any violation. Requires python3 + numpy (scraps_conda).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${WORK:-$(mktemp -d)}"
mkdir -p "$WORK"
KDE="$REPO/inst/scripts/kde_peaks.py"
echo "workdir: $WORK"
fail=0

python3 - "$WORK" <<'EOF'
import gzip, sys
work = sys.argv[1]
def bed(path, rows):
    with gzip.open(path, "wt") as fh:
        for chrom, pos, count, strand in rows:
            fh.write(f"{chrom}\t{pos}\t{pos+1}\t{count}\t{strand}\n")
a1 = [1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10, 12, 9, 7, 5, 4, 3, 1, 1, 1]        # 88
a3 = [1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 9, 8, 9, 7, 5, 4, 3, 1, 1, 1]          # 83
bridge = list(zip([998, 999, 1000, 1001, 1020, 1040, 1058, 1059, 1060, 1061],
                  [5, 10, 20, 5, 1, 1, 10, 30, 60, 10]))                  # 152
isolated = [(4999, 1), (5000, 3), (5001, 1)]                               # 5
for strand in "+-":
    bed(f"{work}/a1_{strand}.bed.gz", [("chr1", 1000 + i, c, strand) for i, c in enumerate(a1)])
    bed(f"{work}/a3_{strand}.bed.gz", [("chr1", 1000 + i, c, strand) for i, c in enumerate(a3)])
    bed(f"{work}/b_{strand}.bed.gz", [("chr1", p, c, strand) for p, c in bridge])
bed(f"{work}/c.bed.gz", [("chr1", 1000 + i, c, "+") for i, c in enumerate(a1)]
                       + [("chr1", p, c, "+") for p, c in isolated])
EOF

run() {  # $1 = bed, $2 = out tsv, rest = extra args; stderr -> $2.log
  python3 "$KDE" -i "$1" -o "$2" "${@:3}" 2> "$2.log"
}
rows() { zcat < "$1" | tail -n +2; }

echo "--- A. umi_support covers the whole cluster"
for strand in + -; do
  for method in downstream kde; do
    run "$WORK/a1_$strand.bed.gz" "$WORK/a1_${strand}_$method.tsv.gz" --summit-method "$method"
    n=$(rows "$WORK/a1_${strand}_$method.tsv.gz" | wc -l | tr -d ' ')
    sup=$(rows "$WORK/a1_${strand}_$method.tsv.gz" | cut -f5)
    summit=$(rows "$WORK/a1_${strand}_$method.tsv.gz" | cut -f2)
    if [ "$n" = 1 ] && [ "$sup" = 88 ]; then
      echo "  A1 $strand $method: OK (1 peak, umi_support=88, summit=$summit)"
    else
      echo "  A1 $strand $method: FAIL (n=$n, umi_support=$sup)"; fail=1
    fi
    # downstream '+': most-downstream base with count >= 0.75*12 = 9 -> 1012; '-': lowest -> 1010
    if [ "$method" = downstream ]; then
      want=1012; [ "$strand" = "-" ] && want=1010
      [ "$summit" = "$want" ] && echo "  A1 $strand summit: OK ($summit)" || { echo "  A1 $strand summit: FAIL ($summit != $want)"; fail=1; }
    fi
    run "$WORK/a3_$strand.bed.gz" "$WORK/a3_${strand}_$method.tsv.gz" --summit-method "$method"
    n=$(rows "$WORK/a3_${strand}_$method.tsv.gz" | wc -l | tr -d ' ')
    sup=$(rows "$WORK/a3_${strand}_$method.tsv.gz" | cut -f5)
    if [ "$n" = 1 ] && [ "$sup" = 83 ]; then
      echo "  A3 $strand $method: OK (83-UMI peak retained)"
    else
      echo "  A3 $strand $method: FAIL (n=$n, umi_support=$sup)"; fail=1
    fi
  done
done

echo "--- B. bridged clusters keep distinct summits and disjoint support"
for strand in + -; do
  for method in downstream kde; do
    run "$WORK/b_$strand.bed.gz" "$WORK/b_${strand}_$method.tsv.gz" --summit-method "$method"
    n=$(rows "$WORK/b_${strand}_$method.tsv.gz" | wc -l | tr -d ' ')
    nuniq=$(rows "$WORK/b_${strand}_$method.tsv.gz" | cut -f1-3 | sort -u | wc -l | tr -d ' ')
    s1=$(rows "$WORK/b_${strand}_$method.tsv.gz" | sed -n 1p | cut -f2)
    s2=$(rows "$WORK/b_${strand}_$method.tsv.gz" | sed -n 2p | cut -f2)
    total=$(rows "$WORK/b_${strand}_$method.tsv.gz" | awk '{s+=$5} END{print s+0}')
    ok=1
    [ "$n" = 2 ] && [ "$nuniq" = 2 ] || ok=0
    [ "$s1" -ge 998 ] && [ "$s1" -le 1001 ] || ok=0     # first cluster's summit inside 998-1001
    [ "$s2" -ge 1058 ] && [ "$s2" -le 1061 ] || ok=0    # second cluster's summit inside 1058-1061
    [ "$total" = 152 ] || ok=0                           # extents partition the island
    if [ "$method" = downstream ]; then
      [ "$s1" = 1000 ] && [ "$s2" = 1060 ] || ok=0       # thr 0.75*20 -> 1000; 0.75*60 -> 1060
    fi
    if [ "$ok" = 1 ]; then
      echo "  B $strand $method: OK (summits $s1,$s2; supports sum to $total)"
    else
      echo "  B $strand $method: FAIL (n=$n uniq=$nuniq summits=$s1,$s2 total=$total)"; fail=1
      rows "$WORK/b_${strand}_$method.tsv.gz" | sed 's/^/    /'
    fi
  done
done

echo "--- C. support floor drops an isolated 5-UMI bump and reports it"
run "$WORK/c.bed.gz" "$WORK/c.tsv.gz"
n=$(rows "$WORK/c.tsv.gz" | wc -l | tr -d ' ')
if [ "$n" = 1 ] && grep -q '1 clusters dropped by the support floor' "$WORK/c.tsv.gz.log"; then
  echo "  C: OK (1 peak kept, 1 dropped; summary: $(cat "$WORK/c.tsv.gz.log" | cut -c1-80)...)"
else
  echo "  C: FAIL (n=$n; log: $(cat "$WORK/c.tsv.gz.log"))"; fail=1
fi

echo "================================================="
if [ "$fail" -eq 0 ]; then
  echo "ALL KDE_PEAKS CHECKS PASSED"
else
  echo "SOME KDE_PEAKS CHECKS FAILED"; exit 1
fi
