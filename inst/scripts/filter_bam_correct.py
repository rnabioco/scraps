import pysam
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from functools import reduce

""" Filter BAM files to only reads with soft-clipped A tail,
suitable for cellranger and starsolo output
"""

def correct_bam_read1(bam, outbam, target_len, filter_cut, single_end):

    samfile = pysam.AlignmentFile(bam, "rb")
    outfile = pysam.AlignmentFile(outbam, "wb", template = samfile)

    # total numbers
    n = reduce(lambda x, y: x + y, [ int(chrom.split("\t")[2]) for chrom in pysam.idxstats(bam).split("\n")[:-1] ])

    #samfile.references
    i = 0
    j = 0
    k = 0
    for read in samfile.fetch(until_eof = True):

        if read.is_unmapped:
            # not mapped, toss
            continue
        
        if not single_end:
            if not read.is_proper_pair:
            # not properly paired, toss
                continue

        if not single_end:
            if not read.is_read1:
            # not read1, toss
                continue

        j += 1

        cigstring = read.cigarstring

        if not read.is_reverse:
            pattern = "^([0-9]{1,})S"
            try:
                clipn = int(re.search(pattern, cigstring).group(1))
                diffn = clipn - target_len
                if abs(diffn) > filter_cut:
                    #outlier, toss
                    continue
            except AttributeError:
                # no soft clipping at the 5'end, toss
                continue

            read.reference_start = read.reference_start - diffn
            if read.reference_start > 0:
                outfile.write(read)
                i += 1
                if not diffn == 0:
                    k += 1

        else:
            pattern = "M([0-9]{1,})S$"
            try:
                clipn = int(re.search(pattern, cigstring).group(1))
                diffn = clipn - target_len
                if abs(diffn) > filter_cut:
                    #outlier, toss
                    continue
            except AttributeError:
                # no soft clipping at the 5'end, toss
                continue
            read.reference_start = read.reference_start + diffn
            if read.reference_start > 0:
                outfile.write(read)
                i += 1
                if not diffn == 0:
                    k += 1

        if int(n / i * 100) == 0:
            print(i)
            #break

    return(k, i, j)

def main():
    parser = argparse.ArgumentParser(description = """
        Utility to correct Read1 position according to soft-clipping
        """)

    parser.add_argument('-i',
                        '--inbam',
                        help = """
                        Bam file to correct
                        """,
                        required = True)

    parser.add_argument('-o',
                        '--outbam',
                        help ="""
                        output Bam file after correction
                        """,
                        default = ".",
                        required = False)

    parser.add_argument('-l',
                        '--targetlen',
                        help ="""number of nt softclipping expected""",
                        default = 58,
                        required = False)

    parser.add_argument('-c',
                        '--filtercut',
                        help ="""maximum nts off of target_len to keep""",
                        default = 10,
                        required = False)
    parser.add_argument('-s',
                        '--single',
                        action="store_true",
                        required = False)
    args=parser.parse_args()

    in_bam = args.inbam
    out_bam = args.outbam
    target_len = int(args.targetlen)
    filter_cut = float(args.filtercut)
    single_end = args.single
    print(single_end)
    print("settings- target_len: ", target_len, "; filter_cut: ", filter_cut)
    k,i,j = correct_bam_read1(in_bam, out_bam, target_len = target_len, filter_cut = filter_cut, single_end = single_end)
    print("fraction of reads kept: ", i/j)
    print("fraction of reads corrected: ", k/i)


if __name__ == '__main__': main()
