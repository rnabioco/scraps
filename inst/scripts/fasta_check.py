#! /usr/bin/env python3

import sys
import errno
import os
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

# check lengths
def get_lengths(filepath, nentries = 10, verbose = False):
    vec = []
    nlines = nentries * 4

    with gzip.open(filepath, 'r') as f:
        i = 0
        for line in f:
            i += 1
            if i%4 == 2:
                line_de = line.decode('utf8').rstrip('\n') #decode('utf8').
                if verbose == True:
                    print(line_de)
                vec.append(len(line_de))
            if i >= nlines:
                break
    f.close()
    return vec

def nog_get_lengths(filepath, nentries = 10, verbose = False):
    vec = []
    nlines = nentries * 4
    with open(filepath, 'r') as f:
        i = 0
        while i <= 10:
            i += 1
            line_de = f.readline().rstrip('\n')
            if i%4 == 2:
                print(line_de)
                vec.append(len(line_de))
    f.close()
    return vec

def check_lengths(lengths_vec, verbose = False):
    lengths_set = list(set(lengths_vec))
    if len(lengths_set) > 1:
        if verbose == True:
            print("inconsistent lengths from file")
        return(-1)
    else:
        if verbose == True:
            print("uniform length of " + str(lengths_set[0]))
        return(lengths_set[0])

def get_qcs(filepath, nentries = 100000, verbose = False):
    vec = []
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
            if i >= nlines:
                break
    f.close()
    df = pd.DataFrame(vec)
    df.columns = list(range(1, len(df.columns) + 1))
    return df
# res = get_qcs(filepath2)
# check_lengths(get_lengths(filepath2))
# sys.getsizeof(get_qcs(filepath2))/1048576

def ATCG_comp(filepath, length, nentries = 100000, verbose = False):
    max_len = length
    d = defaultdict(lambda: [0]*max_len)  # d[char] = [pos0, pos12, ...]
    nlines = nentries * 4
    with gzip.open(filepath, 'r') as f:
        i = 0
        for line in f:
            i += 1
            if i%4 == 2:
                line_de = line.decode('utf8').rstrip('\n')
                #print(line_de)
                for j, char in enumerate(line_de):
                        d[char][j] += 1
            if verbose == True:
                if i%40 == 0:
                    print(i/4)
            if i >= nlines:
                break
    f.close()
    df = pd.DataFrame.from_dict(d)
    df = df.reindex(sorted(df.columns), axis=1)
    return (df/df.iloc[1].sum())
#res2 = ATCG_comp(filepath1, verbose = True)
#res2
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
    # ax.set_ylim(bottom = 0)
    fig.savefig(plotname)
# plot_line(ATCG_comp(filepath2), "ntcomp.pdf")

# def setup():
#     filepipe = 'pipe'
#     os.mkfifo(filepipe)
#
# def cleanup():
#     filepipe = 'pipe'
#     os.remove(filepipe)

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

def wrap_stream_analysis(link, output_dir):
    if os.path.exists("pipe"):
        os.remove("pipe")
    print("for " + link + " : " + "check read lengths")
    stream_file(link, "pipe")
    res_l = get_lengths("pipe")
    os.remove("pipe")
    seq_length = check_lengths(res_l, verbose = True)
    if res_l == -1:
        return()

    print("for " + link + " : " + "check qc scores")
    stream_file(link, "pipe")
    res_qc = get_qcs("pipe")
    os.remove("pipe")
    plot_box(res_qc, os.path.join(output_dir, link.split("/")[-1].split(".")[0] + "_qc.pdf"))

    print("for " + link + " : " + "check ATCG content")
    stream_file(link, "pipe")
    res_atcg = ATCG_comp("pipe", seq_length)
    os.remove("pipe")
    plot_line(res_atcg, os.path.join(output_dir, link.split("/")[-1].split(".")[0] + "_base.pdf"))
# wrap_stream_analysis("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/005/SRR7617315/SRR7617315_2.fastq.gz")

def wrap_links_stream(links, output_dir):
    for linkset in links:
        if len(linkset) == 1:
            continue
        elif len(linkset) == 2:
            for link in linkset:
                wrap_stream_analysis(link, output_dir)

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

    args=parser.parse_args()

    study_id = args.study
    fq_ids = args.ids
    output_dir = args.outputdir
    if output_dir:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    else:
        output_dir = "."
    logfile = os.path.join(output_dir, study_id + "_download_log.txt")
    log_fp = open(logfile, 'w')

    mdata = get_study_metadata(study_id, log_fp)
    links = get_fastq_links(mdata, fq_ids)
    wrap_links_stream(links, output_dir)

    log_fp.close()

if __name__ == '__main__': main()

# wrap_links_stream(m2, "test")
