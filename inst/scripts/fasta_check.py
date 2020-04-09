#! /usr/bin/env python3

import sys
import errno
import os
import re
import gzip
import threading
import requests
from  collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import urllib.request
import urllib.error

""" For given project identifier (i.e SRX.... or PRJN....),
stream data from ENA database to check for compatibility
with scraps workflow, ie read length, qc score, base comp
"""

def get_study_metadata(study_id, logfile):
    baseurl = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
    fields_url = "&result=read_run&fields=run_accession,sample_alias,fastq_ftp"
    full_url = baseurl + study_id + fields_url

  # get study info
    try:
        fp = urllib.request.urlopen(full_url)
    except urllib.error.HTTPError as e:
        sys.exit("unable to retrieve study from ENA: error {}".format(e.code))
    except urllib.error.URLError as e:
        sys.exit("unable to retrieve study from ENA: error {}".format(e.reason))

    mdata = []
    for line in fp:
        line = line.decode("utf8")
        logfile.write(line)
        mdata.append(line)

    fp.close()

    if len(mdata) < 2:
      # ENA will return only headers with some input strings
        sys.exit("unable to retrieve study from ENA")

    return mdata
# m = get_study_metadata("PRJNA428979", open("testlog.txt", 'w'))
# m2 = get_fastq_links(m)
# m2[5]


def get_fastq_links(metadata, fq_ids):
    dl_cmds = []
    header = metadata.pop(0)

    if fq_ids:
        filter_fqs = True
    else:
        filter_fqs = False

    for line in metadata:
        accession, title, ftp_url = line.rstrip().split("\t")

        if filter_fqs:
            if accession not in fq_ids:
                continue

        ftp_urls = ftp_url.split(";")
        for i,furl in enumerate(ftp_urls):
            if furl.startswith('ftp://'):
                pass
            else:
                ftp_urls[i] =  "ftp://" + furl

        if len(ftp_urls) == 2:
            libtype = "paired_end"
        elif len(ftp_urls) == 1:
            libtype = "single_end"
        else:
            sys.exit("unknown ftp urls in metadata")
        dl_cmds.append(ftp_urls)

    return dl_cmds


# if not sys.stdin.isatty():
#     input_stream = sys.stdin

#filepath1 = "/Users/rf/Downloads/tinygex_1_S1_L001_R1_001.fastq.gz"
#filepath2 = "/Users/rf/Downloads/tinygex_1_S1_L001_R2_001.fastq.gz"

# from 10x
# dictionary of instrument id regex: [platform(s)]
InstrumentIDs = {"HWI-M[0-9]{4}$" : ["MiSeq"],
        "HWUSI" : ["Genome Analyzer IIx"],
        "M[0-9]{5}$" : ["MiSeq"],
        "HWI-C[0-9]{5}$" : ["HiSeq 1500"],
        "C[0-9]{5}$" : ["HiSeq 1500"],
        "HWI-D[0-9]{5}$" : ["HiSeq 2500"],
        "D[0-9]{5}$" : ["HiSeq 2500"],
        "J[0-9]{5}$" : ["HiSeq 3000"],
        "K[0-9]{5}$" : ["HiSeq 3000","HiSeq 4000"],
        "E[0-9]{5}$" : ["HiSeq X"],
        "NB[0-9]{6}$": ["NextSeq"],
        "NS[0-9]{6}$" : ["NextSeq"],
        "MN[0-9]{5}$" : ["MiniSeq"]}

# dictionary of flow cell id regex: ([platform(s)], flow cell version and yeild)
FCIDs = {"C[A-Z,0-9]{4}ANXX$" : (["HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], "High Output (8-lane) v4 flow cell"),
         "C[A-Z,0-9]{4}ACXX$" : (["HiSeq 1000", "HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], "High Output (8-lane) v3 flow cell"),
         "H[A-Z,0-9]{4}ADXX$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}BCXX$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
         "H[A-Z,0-9]{4}BCXY$" : (["HiSeq 1500", "HiSeq 2500"], "Rapid Run (2-lane) v2 flow cell"),
         "H[A-Z,0-9]{4}BBXX$" : (["HiSeq 4000"], "(8-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}BBXY$" : (["HiSeq 4000"], "(8-lane) v1 flow cell"),
         "H[A-Z,0-9]{4}CCXX$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}CCXY$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}ALXX$" : (["HiSeq X"], "(8-lane) flow cell"),
         "H[A-Z,0-9]{4}BGXX$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}BGXY$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}BGX2$" : (["NextSeq"], "High output flow cell"),
         "H[A-Z,0-9]{4}AFXX$" : (["NextSeq"], "Mid output flow cell"),
         "A[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq flow cell"),
         "B[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq flow cell"),
         "D[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq nano flow cell"),
         "G[A-Z,0-9]{4}$" : (["MiSeq"], "MiSeq micro flow cell"),
         "H[A-Z,0-9]{4}DMXX$" : (["NovaSeq"], "S2 flow cell")}


# check lengths
def get_lengths(filepath, nentries = 100, verbose = False):
    vec = []
    nlines = nentries * 4

    with gzip.open(filepath, 'r') as f:
        i = 0
        for line in f:
            i += 1
            if i == 1:
                header = line.decode('utf8').rstrip('\n').split(":")
                if len(header) < 3:
                    print("strange header : " + "/".join(header))
                    iid = "unknown"
                    fcid = "unknown"
                else:
                    iid = header[0][1:]
                    if " " in iid:
                        iid = iid.split(" ")[1]
                    fcid = header[2]
            if i%4 == 2:
                line_de = line.decode('utf8').rstrip('\n') #decode('utf8').
                vec.append(len(line_de))
            if i >= nlines:
                break
    f.close()
    # print(iid)
    # print(1)
    # print(fcid)
    sequencer = infer_sequencer(iid, fcid)
    return vec, sequencer
# get_lengths(filepath1)

# infer sequencer from ids from single fastq
def infer_sequencer(iid, fcid):
    seq_by_iid = []
    for key in InstrumentIDs:
        if re.search(key,iid):
            seq_by_iid += InstrumentIDs[key]

    seq_by_fcid = []
    for key in FCIDs:
        if re.search(key,fcid):
            seq_by_fcid += FCIDs[key][0]

    if seq_by_iid and seq_by_fcid:
        sequencer = list(set(seq_by_iid) & set(seq_by_fcid))
        if sequencer:
            return " or ".join(sequencer)
        else:
            return str(seq_by_iid) + " or " + str(seq_by_fcid)
    elif seq_by_iid or seq_by_fcid:
        sequencer = seq_by_iid + seq_by_fcid
        return " or ".join(sequencer)
    else:
        return "unknown"
# infer_sequencer("A00228", "HFWFVDMXX")
# infer_sequencer("unknown", "unknown")

# def nog_get_lengths(filepath, nentries = 10, verbose = False):
#     vec = []
#     nlines = nentries * 4
#     with open(filepath, 'r') as f:
#         i = 0
#         while i <= 10:
#             i += 1
#             line_de = f.readline().rstrip('\n')
#             if i%4 == 2:
#                 print(line_de)
#                 vec.append(len(line_de))
#     f.close()
#     return vec

def check_lengths(lengths_vec, verbose = False):
    lengths_set = list(set(lengths_vec))
    if len(lengths_set) > 1:
        if verbose == True:
            print("inconsistent lengths from file")
        return -1
    else:
        if verbose == True:
            print("uniform length of " + str(lengths_set[0]))
        return lengths_set[0]

def get_lines(filepath, length, nentries = 100000, verbose = False):
    # set up counting
    vec = []
    max_len = length
    d = defaultdict(lambda: [0]*max_len)

    # file
    nlines = nentries * 4
    with gzip.open(filepath, 'r') as f:
        i = 0
        for line in f:
            i += 1
            if i%4 == 0:
                line_de = line.decode('utf8').rstrip('\n')
                if verbose == True:
                    print(line_de)
                vec.append([ord(c) - 33 for c in line_de])
            if i%4 == 2:
                line_de = line.decode('utf8').rstrip('\n')
                for j, char in enumerate(line_de):
                        d[char][j] += 1
            if i >= nlines:
                break
    f.close()
    df_qc = pd.DataFrame(vec)
    df_qc.columns = list(range(1, len(df_qc.columns) + 1))
    df_atcg = pd.DataFrame.from_dict(d)
    df_atcg = df_atcg.reindex(sorted(df_atcg.columns), axis=1)
    return df_qc, df_atcg/df_atcg.iloc[1].sum()

# def get_qcs(filepath, nentries = 100000, verbose = False):
#     vec = []
#     nlines = nentries * 4
#     with gzip.open(filepath, 'r') as f:
#         i = 0
#         for line in f:
#             i += 1
#             if i%4 == 0:
#                 line_de = line.decode('utf8').rstrip('\n')
#                 if verbose == True:
#                     print(line_de)
#                 vec.append([ord(c) - 33 for c in line_de])
#             if i >= nlines:
#                 break
#     f.close()
#     df = pd.DataFrame(vec)
#     df.columns = list(range(1, len(df.columns) + 1))
#     return df
# res = get_qcs(filepath2)
# check_lengths(get_lengths(filepath2))
# sys.getsizeof(get_qcs(filepath2))/1048576

# def ATCG_comp(filepath, length, nentries = 100000, verbose = False):
#     max_len = length
#     d = defaultdict(lambda: [0]*max_len)  # d[char] = [pos0, pos12, ...]
#     nlines = nentries * 4
#     with gzip.open(filepath, 'r') as f:
#         i = 0
#         for line in f:
#             i += 1
#             if i%4 == 2:
#                 line_de = line.decode('utf8').rstrip('\n')
#                 #print(line_de)
#                 for j, char in enumerate(line_de):
#                         d[char][j] += 1
#             if verbose == True:
#                 if i%40 == 0:
#                     print(i/4)
#             if i >= nlines:
#                 break
#     f.close()
#     df = pd.DataFrame.from_dict(d)
#     df = df.reindex(sorted(df.columns), axis=1)
#     return (df/df.iloc[1].sum())
#res2 = ATCG_comp(filepath1, verbose = True)

def plot_box(df, plotname):
    fig, ax = plt.subplots()
    fig.set_size_inches(0.5 * len(df.columns), 10, forward=True)
    fig.set_dpi(300)
    df.boxplot(ax = ax, showfliers = False)
    ax.set_ylim(bottom = 0)
    fig.savefig(plotname)
# plot_box(res, "read2.pdf")

def plot_line(df, plotname):
    plot = df.plot.line()
    fig = plot.get_figure()
    fig.set_size_inches(0.25 * len(df), 10, forward=True)
    fig.set_dpi(300)
    fig.savefig(plotname)
# plot_line(ATCG_comp(filepath2), "ntcomp.pdf")

def stream_file(url, path):
    os.mkfifo(path)
    def _target():
        try:
            t.result = urllib.request.urlretrieve(url, path)
        except Exception as exc:
            t.failure = 1
    t = threading.Thread(target=_target, daemon=True)
    t.start()
        # except IOError as e:
        #     if e.errno == errno.EPIPE:
        #         ignore = True
    # return path
# stream_file("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/005/SRR7617315/SRR7617315_2.fastq.gz", "pipe")
# get_lengths("pipe", verbose = True)
# res2 = ATCG_comp("pipe", 150, nentries = 100000, verbose = False)
# plot_line(res2, "SRR7617315_2_comp.pdf")

def wrap_stream_analysis(link, output_dir, readn = 100000, cutoff = 28):
    if os.path.exists("pipe"):
        os.remove("pipe")
    print("for " + link + " : " + "check read lengths")
    stream_file(link, "pipe")
    res_l, sequencer = get_lengths("pipe")
    os.remove("pipe")
    seq_length = check_lengths(res_l, verbose = True)
    print("sequenced with " + sequencer)
    if seq_length == -1:
        print("inconsistent read length, aborted")
        return()
    elif seq_length <= cutoff:
        print("read length too short, aborted")
        return()

    print("for " + link + " : " + "streaming " + str(readn) + " reads")
    stream_file(link, "pipe")
    res_qc, res_atcg = get_lines("pipe", seq_length, nentries = readn)
    os.remove("pipe")

    print("for " + link + " : " + "plotting")
    plot_box(res_qc, os.path.join(output_dir, link.split("/")[-1].split(".")[0] + "_qc.pdf"))
    plot_line(res_atcg, os.path.join(output_dir, link.split("/")[-1].split(".")[0] + "_base.pdf"))
# wrap_stream_analysis("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/005/SRR7617315/SRR7617315_2.fastq.gz")

def wrap_links_stream(links, output_dir, nreads = 100000, cutoff = 28):
    for linkset in links:
        if len(linkset) == 1:
            print("skipping single end sample")
            continue
        elif len(linkset) == 2:
            for link in linkset:
                wrap_stream_analysis(link, output_dir, nreads, cutoff)

def main():
    parser = argparse.ArgumentParser(description = """
        Utility to stream and qc fastqs from European Nucleotide Archive
        """)

    parser.add_argument('-s',
                        '--study',
                        help = """
                        Study accessions numbers(ERP, SRP, DRP, PRJ prefixes)
                        e.g. SRX017289 or PRJNA342320
                        """,
                        required = True)

    parser.add_argument('-i',
                        '--ids',
                        help = """
                        Specific fastqs to qc, defaults to all fastqs
                        in study
                        """,
                        required = False, nargs = "+")

    parser.add_argument('-o',
                        '--outputdir',
                        help ="""output directory to place fastqs, defaults
                        to '.' """,
                        default = ".",
                        required = False)

    parser.add_argument('-n',
                        '--nreads',
                        help ="""number of reads to stream from each file for calculations""",
                        default = 100000,
                        required = False)

    parser.add_argument('-c',
                        '--cutoff',
                        help ="""lower limit of read length to keep processing""",
                        default = 28,
                        required = False)

    args=parser.parse_args()

    study_id = args.study
    fq_ids = args.ids
    output_dir = args.outputdir
    n_reads = int(args.nreads)
    cut_off = int(args.cutoff)
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = "."
    logfile = os.path.join(output_dir, study_id + "_meta_log.txt")
    log_fp = open(logfile, 'w')

    mdata = get_study_metadata(study_id, log_fp)
    links = get_fastq_links(mdata, fq_ids)
    wrap_links_stream(links, output_dir, n_reads, cut_off)

    log_fp.close()
    print("all finished")

if __name__ == '__main__': main()

# wrap_links_stream(m2, "test")
