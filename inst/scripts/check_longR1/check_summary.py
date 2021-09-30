import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

""" Summarize and highlight fastq_check results
"""

# file1 = "/Users/rf/git_libs/scraps/inst/scripts/panglaodb_human_summary.csv"
# df = pd.read_csv(file1)
# df['judgment'].value_counts().plot(kind='barh')
#
# fig, ax = plt.subplots()
# fig.set_size_inches(8, 8, forward=True)
# fig.set_dpi(300)
# df['judgment'].value_counts().plot(kind='barh')
# fig.savefig("test.pdf")

# samples = df[1].tolist()
# txt = os.path.join(samples[1], samples[1] + "_meta_log.txt")
# csv = os.path.join(samples[1], samples[1] + ".csv")
# os.path.exists(csv)
# df_csv = pd.read_csv("/Users/rf/git_libs/scraps/inst/scripts/SRS2493625/SRS2493625.csv")
# df_csv['read'] = df_csv['link'].str.slice(-10,-9)
# res = df_csv.groupby('read')['seq_length'].max().reset_index()
# res["seq_length"].tolist() > 0
# all(x > 1 for x in res["seq_length"].tolist())
# ["longR1"] + res["seq_length"].tolist()
#


def check_files(sample):
    txt = os.path.join(sample, sample + "_meta_log.txt")
    csv = os.path.join(sample, sample + ".csv")
    if not os.path.exists(txt):
        return([sample, "fail", 0, 0])
    elif not os.path.exists(csv):
        return([sample, "no_fastq", 0, 0])
    else:
        df_csv = pd.read_csv(csv)
        if df_csv.shape[0] == 0:
            return([sample, "not_paired", 0, 0])
        else:
            df_csv['read'] = df_csv['link'].str.slice(-10,-9)
            df_group = df_csv.groupby('read')['seq_length'].max().reset_index()
            if all(x > 40 for x in df_group["seq_length"].tolist()):
                return([sample, "longR1"] + df_group["seq_length"].tolist())
            else:
                return([sample, "shortR1"] + df_group["seq_length"].tolist())



def main():
    parser = argparse.ArgumentParser(description = """
        Utility to stream and qc fastqs from European Nucleotide Archive
        """)

    parser.add_argument('-f',
                        '--file',
                        help = """
                        Tsv file with sample info, used previously in scrape_all.sh
                        """,
                        required = True)

    args=parser.parse_args()

    file1 = args.file
    df = pd.read_csv(file1, header = None)
    samples = df[0].tolist()
    df_final = pd.DataFrame(columns = ["sample", "judgment", "R1", "R2"])

    for sample in samples:
        check = check_files(sample)
        df_final = df_final.append(pd.Series(check, index = df_final.columns), ignore_index=True)

    df_final.to_csv(file1.replace(".txt", "_summary.csv"))

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6, forward=True)
    fig.set_dpi(300)
    df_final['judgment'].value_counts().plot(kind='bar', ax = ax)
    fig.savefig(file1.replace(".txt", ".pdf"))

if __name__ == '__main__': main()
