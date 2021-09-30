awk -F "\"*,\"*" '{system("python3 fastq_check.py -s " $1 " -o " $1)}' 10x_studies.txt
python3 check_summary.py -f 10x_studies.txt