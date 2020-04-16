awk -F "\"*,\"*" '{system("python3 fastq_check.py -s " $2 " -o " $2)}' panglaodb_human.csv
