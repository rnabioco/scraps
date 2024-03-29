from collections import defaultdict
import json
import pandas
chromiumV3 = {"length" : "58", "bc_cut" : "", "R1" : '[WHITELIST_V3,"--soloUMIlen 12 --clip5pNbases 58 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17"]', "R2" : '[WHITELIST_V3,"--soloUMIlen 12"]'}

chromiumV3UG = {"length": "58", "bc_cut": "", "R1": '[WHITELIST_V3,"--soloUMIlen 9 --clip5pNbases 58 --soloCBstart 23 --soloCBlen 16 --soloUMIstart 39"]', "R2": '[WHITELIST_V3,"--soloUMIlen 9"]'}

chromiumV2 = {"length" : "56", "bc_cut" : "", "R1" : '[WHITELIST_V2,"--soloUMIlen 10 --clip5pNbases 56 0 --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17"]', "R2" : '[WHITELIST_V2,"--soloUMIlen 10"]'}

dropseq = {"length" : "50", "bc_cut" : "", "R1" : '["None --soloUMIlen 8 --clip5pNbases 50 0 --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13"]', "R2" : '["None --soloUMIlen 8 --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13"]'}

microwellseq = {"length" : "54", "bc_cut" : "CGACTCACTACAGGG...TCGGTGACACGATCG", "R1" : '["None --soloUMIlen 6 --clip5pNbases 54 0 --soloCBstart 1 --soloCBlen 18 --soloUMIstart 19"]', "R2" : '["None --soloUMIlen 6 --soloCBstart 1 --soloCBlen 18 --soloUMIstart 19"]'}

bd = {"length" : "53", "bc_cut" : "ACTGGCCTGCGA...GGTAGCGGTGACA", "R1" : '["None --soloUMIlen 8 --clip5pNbases 53 0 --soloCBstart 1 --soloCBlen 27 --soloUMIstart 28"]', "R2" : '["None --soloUMIlen 8 --soloCBstart 1 --soloCBlen 27 --soloUMIstart 28"]'}

indrop = {"length" : "32", "bc_cut" : "", "R1" : '["None --soloUMIlen 6 --clip5pNbases 32 0 --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9"]', "R2" : '["None --soloUMIlen 6 --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9"]'}

chemistry = {"chromiumV3" : chromiumV3, "chromiumV2" : chromiumV2, "dropseq" : dropseq, "microwellseq" :  microwellseq, "bd" : bd, "indrop" : indrop}
chemjson = json.dumps(chemistry, indent = 4)

with open('chemistry.json', 'w') as fp:
    fp.write(chemjson)

with open('chemistry.json') as fp:
   chemistry = json.load(fp)

chem_version = "chromiumV3"
eval(chemistry[chem_version]["R1"])
