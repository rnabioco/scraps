import gzip
import argparse

""" Paste Cutadapt internal sequence removal back into FASTQ form
"""

def paste_fastq(file_in, file_out, length1):
    print(length1)
    with open(file_in) as file, gzip.open(file_out, 'wt', compresslevel = 1) as file2:
        i = 0
        for line in file:
            i += 1
            elements = line.split("\t")
            if len(elements) < 9:
                #print("wrong length detected at line", i)
                #print(line)
                continue
            if ";1" in elements[7]:
                if length1:
                    seq1 = elements[4][0:length1]
                    qual1 = elements[8][0:length1]
                else:
                    seq1 = elements[4]
                    qual1 = elements[8]
            if ";2" in elements[7]:
                seq = "@" + elements[0] + "\n" + seq1 + elements[4] + elements[6] + "\n+\n" + qual1 + elements[8] + elements[10] + "\n"
                file2.write(seq)
def main():
    parser = argparse.ArgumentParser(description = """
        Utility to recreate FASTQ from Cutadapt info file, needed for removal of internal adapters/sequences
        """)
    parser.add_argument('-l',
                        '--length1',
                        help = """
                        trim first fragment to length
                        """,
                        required = False)
    parser.add_argument('-i',
                        '--file_in',
                        help = """
                        Cutadapt info file
                        """,
                        required = True)

    parser.add_argument('-o',
                        '--file_out',
                        help ="""
                        output fastq
                        """,
                        required = True)

    args=parser.parse_args()
    length1 = args.length1
    if length1:
        length1 = int(length1)
    file_in = args.file_in
    file_out = args.file_out
    print(length1)
    paste_fastq(file_in, file_out, length1)

if __name__ == '__main__': main()
