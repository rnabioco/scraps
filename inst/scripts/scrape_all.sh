awk -F "\"*,\"*" '{system("python3 fastq_check.py -s " $2 " -o " $2)}' panglaodb_mouse.csv
python3 check_summary.py -f panglaodb_mouse.csv