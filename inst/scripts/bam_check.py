import pysam
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from functools import reduce
from  collections import defaultdict

""" Filter BAM files to only reads with soft-clipped A tail,
suitable for cellranger and starsolo output
"""

def qc_bam_read1(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    # total numbers
    n = reduce(lambda x, y: x + y, [ int(chrom.split("\t")[2]) for chrom in pysam.idxstats(bam).split("\n")[:-1] ])
    i = 0
    j = 0
    k = 0

    for read in samfile.fetch(until_eof = True):
        if read.is_unmapped:
            # not mapped, toss
            continue
        if not read.is_proper_pair:
            # not properly paired, toss
            continue
        if not read.is_read1:
            # not read1, toss
            continue

        j += 1
        if j == 1:
            readlen = len(read.query_qualities)
            d = defaultdict(lambda: [0]*readlen)
            d2 = defaultdict(lambda: [0]*readlen)
            phred = []

        if not read.is_reverse:
            i += 1
            for x, char in enumerate(read.seq):
                d[char][x] += 1
            if i <= 1000000:
                phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
                phredfill2 = np.concatenate([read.query_qualities,np.array(phredfill)])
                phred.append(phredfill2)
            # if i == 1000000:
            #     print("finished " + "1000000" + " reads for phred score")
            if i == 10000000:
                # print("finished " + "10000000" + " reads for nucleotide composition")
                break

        else:
            i += 1
            for x, char in enumerate(read.get_forward_sequence()):
                d[char][x] += 1
            if i <= 1000000:
                phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
                phredfill2 = np.concatenate([np.flip(read.query_qualities),np.array(phredfill)])
                phred.append(phredfill2)
            # if i == 1000000:
            #     print("finished " + "1000000" + " reads for phred score")
            if i == 10000000:
                # print("finished " + "10000000" + " reads for nucleotide composition")
                break

    return(d, phred)

def qc_bam_unmapped_read1(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    # total numbers
    n = reduce(lambda x, y: x + y, [ int(chrom.split("\t")[2]) for chrom in pysam.idxstats(bam).split("\n")[:-1] ])
    i = 0
    j = 0
    k = 0

    for read in samfile.fetch(until_eof = True):
        if not read.is_unmapped:
            # mapped, toss
            continue
        if not read.is_read1:
            # not read1, toss
            continue

        j += 1
        if j == 1:
            readlen = len(read.query_qualities)
            d = defaultdict(lambda: [0]*readlen)
            d2 = defaultdict(lambda: [0]*readlen)
            phred = []

        if not read.is_reverse:
            i += 1
            for x, char in enumerate(read.seq):
                d[char][x] += 1
            if i <= 1000000:
               phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
               phredfill2 = np.concatenate([read.query_qualities,np.array(phredfill)])
               phred.append(phredfill2)
            # if i % 1000000 == 0:
            #     print("finished " + "1000000" + " reads for phred score")
            if i % 10000000 == 0:
                # print("finished " + "10000000" + " reads for nucleotide composition")
                break

        else:
            i += 1
            for x, char in enumerate(read.get_forward_sequence()):
                d[char][x] += 1
            if i <= 1000000:
               phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
               phredfill2 = np.concatenate([np.flip(read.query_qualities),np.array(phredfill)])
               phred.append(phredfill2)
            # if i % 1000000 == 0:
            #    print("finished " + "1000000" + " reads for phred score")
            if i % 10000000 == 0:
            #     print("finished " + "10000000" + " reads for nucleotide composition")
                 break
                
    return(d, phred)

def plot_line(df, plotname):
    plot = df.plot.line()
    fig = plot.get_figure()
    fig.set_size_inches(0.25 * len(df), 10, forward=True)
    fig.set_dpi(300)
    fig.savefig(plotname)

def plot_box(df, plotname):
    fig, ax = plt.subplots()
    fig.set_size_inches(0.5 * len(df.columns), 10, forward=True)
    fig.set_dpi(300)
    df.boxplot(ax = ax, showfliers = False)
    ax.set_ylim(bottom = 0)
    fig.savefig(plotname)

def do_plot_nt(d, out_name):
    df_atcg = pd.DataFrame.from_dict(d)
    try:
        df_atcg = df_atcg[["A", "T", "C", "G", "N"]]
    except KeyError:
        df_atcg = df_atcg[["A", "T", "C", "G"]]

    plot_line(df_atcg/df_atcg.iloc[1].sum(), out_name + "_ntcomp.pdf")

def do_plot_box(phred, out_name):
    df = pd.DataFrame(phred)
    plot_box(df, out_name + "_phred.pdf")

def main():
    parser = argparse.ArgumentParser(description = """
        Utility to render reports on Read1 quality
        """)

    parser.add_argument('-i',
                        '--inbam',
                        help = """
                        Bam file to generate quality control plots
                        """,
                        required = True)

    parser.add_argument('-o',
                        '--outname',
                        help ="""
                        output prefix for pdf files
                        """,
                        default = "qc",
                        required = False)

    parser.add_argument('-m',
                        '--qcmode',
                        help ="""mode, default to check mapped read1, can be set to unmapped""",
                        default = "mapped",
                        required = False)

    args=parser.parse_args()

    in_bam = args.inbam
    out_name = args.outname
    qc_mode = args.qcmode

    if qc_mode == "mapped":
        print("checking paired and mapped read1")
        d,phred = qc_bam_read1(in_bam)
    elif qc_mode == "unmapped":
        print("checking unmapped read1")
        d,phred = qc_bam_unmapped_read1(in_bam)

    do_plot_nt(d, out_name)
    do_plot_box(phred, out_name)

if __name__ == '__main__': main()
#
# bam = "/Users/rf/scraps072321/qc/PT46L_miseq_R1_Aligned.sortedByCoord.out.bam"
# outbam = "/Users/rf/scraps072321/DX12_r1_corrected2.bam"
# samfile = pysam.AlignmentFile(bam, "rb")
# outfile = pysam.AlignmentFile(outbam, "w", template = samfile)
# target_len = 58
# i = 0
# j = 0
# k = 0
# filter_cut = 5
# d = defaultdict(lambda: [0]*151)
# d2 = defaultdict(lambda: [0]*151)
# phred = []
# np.empty((0, 151), int)
# for read in samfile.fetch(until_eof = True):
#     if read.is_unmapped:
#         # not mapped, toss
#         continue
#
#     if not read.is_proper_pair:
#         # not properly paired, toss
#         continue
#
#     if not read.is_read1:
#         # not read1, toss
#         continue
#
#     i += 1
#     if i == 1:
#         readlen = len(read.query_qualities)
#
#     cigstring = read.cigarstring
#
#     if not read.is_reverse:
#         # pattern = "^([0-9]{1,})S"
#         # try:
#         #     clipn = int(re.search(pattern, cigstring).group(1))
#         #     diffn = clipn - target_len
#         #     if abs(diffn) > filter_cut:
#         #         #outlier, toss
#         #         continue
#         # except AttributeError:
#         #     # no soft clipping at the 5'end, toss
#         #     continue
#         # read.reference_start = read.reference_start - diffn
#         #outfile.write(read)
#         k += 1
#         for x, char in enumerate(read.seq):
#             d[char][x] += 1
#         if k <= 1000000:
#            phredfill = [np.nan for _ in range(151 - len(read.query_qualities))]
#            phredfill2 = np.concatenate([read.query_qualities,np.array(phredfill)])
#            phred.append(phredfill2)
#
#     else:
#         pattern = "M([0-9]{1,})S$"
#         try:
#             clipn = int(re.search(pattern, cigstring).group(1))
#             diffn = clipn - target_len
#             if abs(diffn) > filter_cut:
#                 #outlier, toss
#                 continue
#         except AttributeError:
#             # no soft clipping at the 5'end, toss
#             continue
#         # read.reference_start = read.reference_start + diffn
#         #outfile.write(read)
#         k += 1
#         for x, char in enumerate(read.get_forward_sequence()):
#             d[char][x] += 1
#         if k <= 1000000:
#            phredfill = [np.nan for _ in range(151 - len(read.query_qualities))]
#            phredfill2 = np.concatenate([np.flip(read.query_qualities),np.array(phredfill)])
#            phred.append(phredfill2)
#
#     if k % 1000000 == 0:
#         print(k)
#     if k % 10000000 == 0:
#         print(k)
#         break

# df_atcg = pd.DataFrame.from_dict(d)
# try:
#     df_atcg = df_atcg[["A", "T", "C", "G", "N"]]
# except KeyError:
#     df_atcg = df_atcg[["A", "T", "C", "G"]]
#     # df_atcg = df_atcg.reindex(sorted(df_atcg.columns), axis=1)
# df_atcg/df_atcg.iloc[1].sum()
#
# plot_line(df_atcg/df_atcg.iloc[1].sum(), "PT46L_miseq_r1_ntcomp.pdf")

# for read in samfile.fetch(until_eof = True):
#     if not read.is_read1:
#         # not read1, toss
#         continue
#     print(read.seq)
#     print(len(read.query_qualities))
#     break
#
#
# df = pd.DataFrame(phred)
# plot_box(df, "PT46L_miseq_r1_phred.pdf")
#
#
# pysam.idxstats(bam)
# int(n/n * 100)
