import pysam
import re
import argparse

""" Filter BAM files to only reads with soft-clipped A tail,
suitable for cellranger and starsolo output
"""

def filter_bam_by_A(bam, outbam, nmin, fraq_min):
    samfile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template = samfile)
    #samfile.references
    i = 0
    k = 0
    for read in samfile.fetch(until_eof = True):
        i += 1
        if read.is_unmapped:
            # not mapped, toss
            continue

        cigstring = read.cigarstring
        readseq = read.seq

        if not read.is_reverse:
            pattern = "M([0-9]{1,})S$"
            try:
                clipn = int(re.search(pattern, cigstring).group(1))
                if clipn <= nmin:
                    # too short, toss
                    clipn = 0
                    continue
            except AttributeError:
                # no soft clipping at the 3'end, toss
                clipn = 0
                continue
            clipseq = readseq[-clipn: ]
            afrac = clipseq.count("A") / clipn

        else:
            pattern = "^([0-9]{1,})S[0-9]+M"
            try:
                clipn = int(re.search(pattern, cigstring).group(1))
                if clipn <= nmin:
                    # too short, toss
                    clipn = 0
                    continue
            except AttributeError:
                # no soft clipping at the 3'end, toss
                clipn = 0
                continue
            clipseq = readseq[ :clipn]
            afrac = clipseq.count("T") / clipn

        if afrac >= fraq_min:
            k += 1
            outfile.write(read)
    return(k, i)

def main():
    parser = argparse.ArgumentParser(description = """
        Utility to filter BAM files according to soft-clipped A tail length
        """)

    parser.add_argument('-i',
                        '--inbam',
                        help = """
                        Bam file to filter
                        """,
                        required = True)

    parser.add_argument('-o',
                        '--outbam',
                        help ="""
                        output Bam file after filtering
                        """,
                        default = ".",
                        required = False)

    parser.add_argument('-n',
                        '--ntail',
                        help ="""number of reads to stream from each file for calculations""",
                        default = 5,
                        required = False)

    parser.add_argument('-f',
                        '--fractionmin',
                        help ="""number of reads to stream from each file for calculations""",
                        default = 0.8,
                        required = False)

    args=parser.parse_args()

    in_bam = args.inbam
    out_bam = args.outbam
    n_tail = int(args.ntail)
    fraction_min = float(args.fractionmin)

    print("settings- min_soft_clip: ", n_tail, "; min_fraction_A: ", fraction_min)
    k, i = filter_bam_by_A(in_bam, out_bam, nmin = n_tail, fraq_min = fraction_min)
    print("number of reads kept: ", k)
    print("fraction of reads kept: ", k/i)



if __name__ == '__main__': main()
