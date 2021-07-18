# -*- coding: utf-8 -*-
Description = '''
This scripts can calculate the duplicated sequence caused by Tn5 insertion
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-05-07.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import Levenshtein

def fuzz_align(query_seq,match_seq,mismatch):
    dist = Levenshtein.distance(query_seq, match_seq)
    if dist <= mismatch:  # find first then break
        return True

def run(options):
    
    for n in range(options.TEST):
        
        print("-->>>Starting assumption that Tn5 introduces {}bp duplicates and allow {} mismatch".format(n,options.MISMATCH))
        bed_pos_dict = {}
        bed_neg_dict = {}
        
        pos = 0
        neg = 0
        
        with open(options.BED) as bed:
            for line in bed:
                line = line.strip().split("\t")
                strand = line[5]
                if strand == "+":
                    pos += 1
                    key = line[0] + "\t" + str(int(line[1])+ n)
                    if key not in bed_pos_dict.keys():
                        bed_pos_dict[key] = line[4] + "\t" + str(1)
                    else:
                        value = bed_pos_dict[key].split("\t")
                        num = int(value[1]) + 1
                        bed_pos_dict[key] = value[0] + "\t" + str(num)
                        
                if strand == "-":
                    neg += 1
                    key = line[0] + "\t" + line[2]
                    if key not in bed_neg_dict.keys():
                        bed_neg_dict[key] = line[4] + "\t" + str(1)
                    else:
                        value = bed_neg_dict[key].split("\t")
                        num = int(value[1]) + 1
                        bed_neg_dict[key] = value[0] + "\t" + str(num)  
                        
        print("Totally {} sequences in forward (+) strand".format(pos))
        print("Totally {} sequences in reverse (-) strand".format(neg))
        
        O = open(str(n) + "bp_duplicates.txt","w+")
        pair = 0
        match = 0
        for key in bed_pos_dict.keys():
            if key in bed_neg_dict.keys(): 
                pos_v = bed_pos_dict[key].split("\t")[0]
                pos_n = bed_pos_dict[key].split("\t")[1]
                
                neg_v = bed_neg_dict[key].split("\t")[0]
                neg_n = bed_neg_dict[key].split("\t")[1]
                times = max(int(pos_n),int(neg_n))
                pair += times
                
                dup = neg_v[-n:]
                if pos_v.startswith(dup) or fuzz_align(pos_v[n],dup,options.MISMATCH):
                    O.write(dup + "\n")
                    match += times
        O.close()
        print("Totally {} reads were paired around Tn5 cutting sites".format(pair))
        print("Totally {}({}%) paired reads were found to be {}bp duplicated\n".format(match,match*100/pos,n))

def main():   
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--BED',action = 'store', type=str, dest = 'BED',required=True,
                        help = 'BED file, better deduplicated!')
    parser.add_argument('-m', '--MISMATCH',action = 'store', type=int, dest = 'MISMATCH',required=True,
                        help = 'Allow mismatches when calculate overlaped reads. [1 is for recommendation]')
    parser.add_argument('-n', '--TEST',action = 'store', type=int, dest = 'TEST',required=True,
                        help = 'How many duplicates you think it will be?')
    
    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    run(options)

if __name__ == '__main__':
    main()



