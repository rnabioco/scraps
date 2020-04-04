import sys
import atexit
import os
import subprocess
import gzip
import threading
import requests
from  collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import urllib.request
if not sys.stdin.isatty():
    input_stream = sys.stdin

#filepath1 = "/Users/rf/Downloads/tinygex_1_S1_L001_R1_001.fastq.gz"
#filepath2 = "/Users/rf/Downloads/tinygex_1_S1_L001_R2_001.fastq.gz"

# check lengths
def get_lengths(filepath, nentries = 10, verbose = False):
    vec = []
    nlines = nentries * 4

    # if not sys.stdin.isatty():
    #     f = input_stream
    #     i = 0
    #     for line in f:
    #         i += 1
    #         if i%4 == 2:
    #             try:
    #                 line_de = line.decode('utf8').rstrip('\n')
    #             except:
    #                 line_de = line.rstrip('\n')
    #             if verbose == True:
    #                 print(line_de)
    #             vec.append(len(line_de))
    #         if i >= nlines:
    #             break
    #     return vec

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
    fig.set_size_inches(0.5 * len(res.columns), 10, forward=True)
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

def setup():
    filepipe = 'pipe'
    os.mkfifo(filepipe)

def cleanup():
    filepipe = 'pipe'
    os.remove(filepipe)

# atexit.register(cleanup)
# setup()
# cleanup()
# filepipe='pipe'
# cmd = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/004/SRR7617314/SRR7617314_1.fastq.gz"
# p = subprocess.Popen(["curl", "-Ls", cmd ], stdout = subprocess.PIPE)
# p = subprocess.Popen(["gunzip", "-c", "/Users/rf/Downloads/tinygex_1_S1_L001_R2_001.fastq.gz"], stdout = filepipe)
# p = subprocess.call(["less", "/Users/rf/newselect201905/sq_gen/sq_gen/IBAdown.txt"])
#
# stdout = p.communicate()
# stdout
#
# res2 = ATCG_comp(filepipe, verbose = True)
# nog_get_lengths(stdout, verbose = True)
#
# with open(filepipe, 'r') as read_fifo:
#     read_fifo.readline()
# filepipe
#
# nog_get_lengths("/Users/rf/newselect201905/sq_gen/sq_gen/IBAdown.txt", verbose = True)
#
# response = requests.get("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/004/SRR7617314/SRR7617314_1.fastq.gz", stream=True)
def stream_file(url, path):
    os.remove(path)
    os.mkfifo(path)
    t = threading.Thread(target=urllib.request.urlretrieve, args=(url, path), daemon=True)
    t.start()
    return path
stream_file("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR761/005/SRR7617315/SRR7617315_2.fastq.gz", "pipe")
get_lengths("pipe", verbose = True)
res2 = ATCG_comp("pipe", 150, nentries = 100000, verbose = False)

plot_line(res2, "SRR7617315_2_comp.pdf")
