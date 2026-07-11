import gzip
import argparse
from collections import defaultdict

import pysam

""" UMI-deduplicated single-base priming pileup from a featureCounts assigned
BAM, replacing `umi_tools dedup` + `bedtools genomecov`.

Counts distinct (cell-barcode, UMI) pairs at each single-base read end position.
This reproduces `umi_tools dedup --method=unique` followed by
`bedtools genomecov -5/-3 -dz` bit-for-bit: umi_tools dedup groups reads by
alignment position and, with method=unique, keeps one read per distinct UMI at a
position; genomecov then piles those up per base. Counting unique (CB, UB) per
end position yields the same per-base counts in a single pass, with no
intermediate deduplicated BAM.

End selection mirrors the original rules:
  R2     : read 3' end  (bed_R2:     genomecov -3)
  R1     : read 5' end, mapped reads only        (bed_R1:     -F 4, genomecov -5)
  paired : read 5' end, first-in-pair mate only  (bed_paired: -f 0x42, -5)

Strand handling:
  --strand-split emits a 5-column stranded bed (chrom start end count strand)
  used by discovery. The printed strand is the TRANSCRIPT strand, supplied via
  --label-plus/--label-minus for reads aligned to the + / - strand respectively
  (genomecov -strand selects by aligned strand; the label is the corresponding
  transcript strand). Without --strand-split a 4-column bed
  (chrom start end count) is emitted, matching the bed_* rules.
"""


def _end_position(read, end):
    """Single-base end coordinate (0-based) matching bedtools genomecov -5/-3.

    genomecov reports the 5' or 3' most base of the alignment in transcript-read
    orientation: for a forward read the 5' end is reference_start and the 3' end
    is reference_end-1; for a reverse read these swap.
    """
    if end == "5":
        return read.reference_end - 1 if read.is_reverse else read.reference_start
    else:  # "3"
        return read.reference_start if read.is_reverse else read.reference_end - 1


def _passes_mode(read, mode):
    """Read-mode mate/flag filter matching the samtools prefilters in the rules."""
    if mode == "R1":
        # bed_R1: samtools view -F 4 (mapped only)
        return not read.is_unmapped
    if mode == "paired":
        # bed_paired: samtools view -f 0x42 (read paired + first in pair)
        return read.is_paired and read.is_read1
    # R2: single-end, no extra filter (unmapped dropped below anyway)
    return not read.is_unmapped


def pileup(bam_path, out_path, end, mode, strand_split, label_plus, label_minus):
    # key -> set of (cell, umi). key includes strand only when splitting.
    counts = defaultdict(set)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            if not _passes_mode(read, mode):
                continue
            try:
                cb = read.get_tag("CB")
                ub = read.get_tag("UB")
            except KeyError:
                continue
            if cb == "-" or ub == "-":
                continue
            pos = _end_position(read, end)
            if strand_split:
                strand = label_minus if read.is_reverse else label_plus
                key = (read.reference_name, pos, strand)
            else:
                key = (read.reference_name, pos)
            counts[key].add((cb, ub))

    with gzip.open(out_path, "wt") as out:
        if strand_split:
            # discovery stranded_bed awk: print $1, $2, $2 + 1, $3, strand on a
            # 0-based genomecov -dz position -> interval [pos, pos+1)
            for (chrom, pos, strand) in sorted(counts):
                n = len(counts[(chrom, pos, strand)])
                out.write(f"{chrom}\t{pos}\t{pos + 1}\t{n}\t{strand}\n")
        else:
            # bed_* awk: print $1, $2 - 1, $2, $3 on a 0-based genomecov -dz
            # position -> interval [pos-1, pos)
            for (chrom, pos) in sorted(counts):
                n = len(counts[(chrom, pos)])
                out.write(f"{chrom}\t{pos - 1}\t{pos}\t{n}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Unique-UMI single-base priming pileup from an assigned BAM"
    )
    parser.add_argument("-i", "--inbam", required=True,
                        help="featureCounts assigned BAM (CB/UB tags)")
    parser.add_argument("-o", "--out", required=True,
                        help="output gzipped bed")
    parser.add_argument("--end", choices=["5", "3"], required=True,
                        help="read end to pile up (5' or 3')")
    parser.add_argument("--mode", choices=["R1", "R2", "paired"], required=True,
                        help="alignment mode (sets mate/flag filter)")
    parser.add_argument("--strand-split", action="store_true",
                        help="emit 5-col stranded bed (discovery)")
    parser.add_argument("--label-plus", default="+",
                        help="transcript strand label for +-aligned reads")
    parser.add_argument("--label-minus", default="-",
                        help="transcript strand label for --aligned reads")
    args = parser.parse_args()

    pileup(args.inbam, args.out, args.end, args.mode,
           args.strand_split, args.label_plus, args.label_minus)


if __name__ == "__main__":
    main()
