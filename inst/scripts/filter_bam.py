import re
import argparse

import pysam

""" Filter BAM files to reads carrying a poly(A)/poly(T) soft clip at the
counted (priming) end, suitable for cellranger / starsolo output and for
pre-aligned BAMs (scraps `bam` mode). Optionally selects one mate of a
paired-end BAM and rewrites it as a single-end record.

Which soft clip is inspected depends on the counted end and the read's aligned
strand (the priming site is the transcript 3' end; the clip is the untemplated
poly(A) tail or, for antisense reads, its poly(T) complement):

    --end 3 (read is on the transcript strand; poly(A) at the read 3' end)
        forward read : trailing S, count A      reverse read : leading S, count T
    --end 5 (read is antisense; poly(T) at the read 5' end)
        forward read : leading S, count T       reverse read : trailing S, count A

--end 3 with --mate 0 is the historical behaviour (filter_R2). Only a soft clip
directly adjacent to an M block is recognised: CIGARs using =/X operations,
hard clips (pbmm2, minimap2 --eqx), or an I/D/N operation between the clip and
the first/last M never match and such reads are dropped.
"""

# (end, is_reverse) -> (clip regex, slice of the read sequence, base to count)
_CLIP_RULES = {
    ("3", False): (re.compile(r"M([0-9]+)S$"), "tail", "A"),
    ("3", True):  (re.compile(r"^([0-9]+)S[0-9]+M"), "head", "T"),
    ("5", False): (re.compile(r"^([0-9]+)S[0-9]+M"), "head", "T"),
    ("5", True):  (re.compile(r"M([0-9]+)S$"), "tail", "A"),
}

# flag bits cleared when a selected mate is rewritten as a single-end record:
# paired, proper pair, mate unmapped, mate reverse, read1, read2
_PAIR_BITS = 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80


def _clip_passes(read, end, nmin, fraq_min):
    """True if the read carries a soft clip > nmin bases at the counted end
    with a poly(A)/poly(T) base fraction >= fraq_min."""
    pattern, where, base = _CLIP_RULES[(end, read.is_reverse)]
    cigstring = read.cigarstring
    if cigstring is None:
        return False
    m = pattern.search(cigstring)
    if m is None:
        return False  # no soft clip at the counted end
    clipn = int(m.group(1))
    if clipn <= nmin:
        return False  # too short
    readseq = read.query_sequence
    clipseq = readseq[-clipn:] if where == "tail" else readseq[:clipn]
    return clipseq.count(base) / clipn >= fraq_min


def _unpair(read):
    """Rewrite a selected mate as a single-end record in place."""
    read.flag &= ~_PAIR_BITS
    read.next_reference_id = -1
    read.next_reference_start = -1
    read.template_length = 0
    return read


def filter_bam_by_A(bam, outbam, nmin, fraq_min, end="3", mate=0,
                    clip_filter=True):
    """Stream `bam` -> `outbam`, keeping mapped reads that (a) are the selected
    mate of a proper pair when mate != 0 (rewritten single-end) and (b) pass the
    soft-clip test at `end` when clip_filter is True. Returns (kept, total)."""
    samfile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template=samfile)
    i = 0
    k = 0
    for read in samfile.fetch(until_eof=True):
        i += 1
        if read.is_unmapped:
            continue
        if mate:
            if not (read.is_paired and read.is_proper_pair):
                continue
            if mate == 1 and not read.is_read1:
                continue
            if mate == 2 and not read.is_read2:
                continue
        if clip_filter and not _clip_passes(read, end, nmin, fraq_min):
            continue
        if mate:
            _unpair(read)
        k += 1
        outfile.write(read)
    outfile.close()
    samfile.close()
    return (k, i)


def main():
    parser = argparse.ArgumentParser(description="""
        Filter a BAM to reads with a poly(A)/poly(T) soft clip at the counted
        end; optionally select and unpair one mate of a paired-end BAM.
        """)
    parser.add_argument('-i', '--inbam', required=True,
                        help="BAM file to filter")
    parser.add_argument('-o', '--outbam', required=True,
                        help="output BAM file after filtering")
    parser.add_argument('-n', '--ntail', default=5, type=int, required=False,
                        help="minimum soft-clip length (bases); clips <= this "
                             "are rejected (default 5)")
    parser.add_argument('-f', '--fractionmin', default=0.8, type=float,
                        required=False,
                        help="minimum fraction of A (or T) in the soft clip "
                             "(default 0.8)")
    parser.add_argument('-e', '--end', choices=["3", "5"], default="3",
                        help="counted read end at which the clip is checked: "
                             "3 = poly(A) at the read 3' end (read on the "
                             "transcript strand; default), 5 = poly(T) at the "
                             "read 5' end (antisense read)")
    parser.add_argument('--mate', type=int, choices=[0, 1, 2], default=0,
                        help="0 = single-end BAM (default); 1 or 2 = paired-end "
                             "BAM: keep only that mate of proper pairs and "
                             "rewrite it as a single-end record")
    parser.add_argument('--no-clip-filter', action="store_true",
                        help="skip the soft-clip test (e.g. only select/unpair "
                             "a mate)")
    args = parser.parse_args()

    print("settings- end:", args.end, "; mate:", args.mate,
          "; clip_filter:", not args.no_clip_filter,
          "; min_soft_clip:", args.ntail, "; min_fraction_A/T:", args.fractionmin)
    k, i = filter_bam_by_A(args.inbam, args.outbam, nmin=args.ntail,
                           fraq_min=args.fractionmin, end=args.end,
                           mate=args.mate, clip_filter=not args.no_clip_filter)
    print("number of reads kept: ", k)
    print("fraction of reads kept: ", k / i if i else 0.0)


if __name__ == '__main__':
    main()
