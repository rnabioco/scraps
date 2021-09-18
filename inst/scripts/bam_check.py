import pysam
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from functools import reduce
from  collections import defaultdict
import statistics

""" QC for various metrics
"""

def qc_bam_read1(bam, nreads, readlen_default = None):
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

        # parse read name
        id,run,cell,lane,tile,x,y =read.query_name.split(":")
        tile = lane + ":" + tile
        try:
            CB = read.get_tag('CR')
        except KeyError:
            continue

        j += 1

        if j == 1:
            if readlen_default is None:
                readlen = len(read.query_qualities)
            else:
                readlen = readlen_default
            d = defaultdict(lambda: [0]*readlen)
            d2 = defaultdict(lambda: [0]*readlen)
            phred = []
            cs = defaultdict(int)
            cp = defaultdict(int)
            ts = defaultdict(int)
            tp = defaultdict(int)

        if not read.is_reverse:
            i += 1
            for x, char in enumerate(read.seq):
                d[char][x] += 1
            if i <= 1000000:
                phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
                phredfill2 = np.concatenate([read.query_qualities,np.array(phredfill)])
                phred.append(phredfill2)

        else:
            i += 1
            for x, char in enumerate(read.get_forward_sequence()):
                d[char][x] += 1
            if i <= 1000000:
                phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
                phredfill2 = np.concatenate([np.flip(read.query_qualities),np.array(phredfill)])
                phred.append(phredfill2)

        phred_part1 = statistics.mean(phredfill2[58:(58 + 40)])
        phred_part2 = statistics.mean(phredfill2[(58-20):58])
        if not np.isnan(phred_part1):
            cs[CB] += 1
            cp[CB] = cp[CB] + phred_part1
            ts[tile] += 1
            tp[tile] = tp[tile] + phred_part1

        if i > nreads:
            break

    return(d, phred, cs, cp, ts, tp)

def qc_bam_unmapped_read1(bam, nreads):
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

        # parse read name
        id,run,cell,lane,tile,x,y =read.query_name.split(":")
        tile = lane + ":" + tile
        try:
            CB = read.get_tag('CR')
        except KeyError:
            continue

        j += 1

        if j == 1:
            if readlen_default is None:
                readlen = len(read.query_qualities)
            else:
                readlen = readlen_default
            d = defaultdict(lambda: [0]*readlen)
            d2 = defaultdict(lambda: [0]*readlen)
            phred = []
            cs = defaultdict(int)
            cp = defaultdict(int)
            ts = defaultdict(int)
            tp = defaultdict(int)

        if not read.is_reverse:
            i += 1
            for x, char in enumerate(read.seq):
                d[char][x] += 1
            if i <= 1000000:
               phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
               phredfill2 = np.concatenate([read.query_qualities,np.array(phredfill)])
               phred.append(phredfill2)

        else:
            i += 1
            for x, char in enumerate(read.get_forward_sequence()):
                d[char][x] += 1
            if i <= 1000000:
               phredfill = [np.nan for _ in range(readlen - len(read.query_qualities))]
               phredfill2 = np.concatenate([np.flip(read.query_qualities),np.array(phredfill)])
               phred.append(phredfill2)

        phred_part1 = statistics.mean(phredfill2[58:(58 + 40)])
        phred_part2 = statistics.mean(phredfill2[(58-20):58])
        if not np.isnan(phred_part1):
            cs[CB] += 1
            cp[CB] = cp[CB] + phred_part1
            ts[tile] += 1
            tp[tile] = tp[tile] + phred_part1

        if i > nreads:
            break

    return(d, phred, cs, cp, ts, tp)

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

def dicts_to_csv(s, p, out_name, name2):
    df_count = pd.DataFrame.from_dict([s])
    df_phred = pd.DataFrame.from_dict([p])
    df2 = df_count.append(df_phred, ignore_index=True)
    df2.iloc[1] = df2.iloc[1] / df2.iloc[0]
    df2.T.to_csv(out_name + name2 + ".csv")

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
    parser.add_argument('-w',
                        '--writecsv',
                        help ="""whether to output csv qc metrics""",
                        action="store_true")
    parser.add_argument('-n',
                        '--numberreads',
                        help ="""upper limit for reads to check""",
                        default = 10000000,
                        required = False)
    parser.add_argument('-l',
                        '--readlength',
                        help ="""default read length, otherwise inferred""",
                        default = None,
                        required = False)
    args=parser.parse_args()

    in_bam = args.inbam
    out_name = args.outname
    qc_mode = args.qcmode
    wr = args.writecsv
    nreads = int(args.numberreads)
    if args.readlength is None:
        lread = None
    else:
        lread = int(args.readlength)

    if qc_mode == "mapped":
        print("checking paired and mapped read1")
        d,phred,cs,cp,ts,tp = qc_bam_read1(in_bam, nreads, lread)
    elif qc_mode == "unmapped":
        print("checking unmapped read1")
        d,phred,cs,cp,ts,tp = qc_bam_unmapped_read1(in_bam, nreads, lread)

    do_plot_nt(d, out_name)
    do_plot_box(phred, out_name)

    if wr:
        dicts_to_csv(cs, cp, out_name, "_barcode")
        dicts_to_csv(ts, tp, out_name, "_tile")

if __name__ == '__main__': main()
