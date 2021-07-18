# -*- coding: utf-8 -*-
Description = '''
This script can calculate the duplicated sequence and location caused by Tn5 insertion given bed file.
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-05-27.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import argparse

def run(options):

    bed_pos_dict = {}
    bed_neg_dict = {}
    
    n = 9 
    with open(options.BED) as bed:
        for line in bed:
            line = line.strip().split("\t")
            strand = line[5]
            if strand == "+":
                key = line[0] + "\t" + str(int(line[1])+ n)
                if key not in bed_pos_dict.keys():
                    bed_pos_dict[key] = line[4] + "\t" + str(1)
                else:
                    value = bed_pos_dict[key].split("\t")
                    num = int(value[1]) + 1
                    bed_pos_dict[key] = value[0] + "\t" + str(num)
                    
            if strand == "-":
                key = line[0] + "\t" + line[2]
                if key not in bed_neg_dict.keys():
                    bed_neg_dict[key] = line[4] + "\t" + str(1)
                else:
                    value = bed_neg_dict[key].split("\t")
                    num = int(value[1]) + 1
                    bed_neg_dict[key] = value[0] + "\t" + str(num)  
    
    bed = open(options.BED + "_9bp_duplicates.bed","w+")
    pair = 0
    
    for key in bed_pos_dict.keys():
        if key in bed_neg_dict.keys():
            pos_v = bed_pos_dict[key].split("\t")[0]
            pos_n = bed_pos_dict[key].split("\t")[1]
                
            neg_v = bed_neg_dict[key].split("\t")[0]
            neg_n = bed_neg_dict[key].split("\t")[1]
            times = max(int(pos_n),int(neg_n))
        
            dup = neg_v[-9:]
            if pos_v.startswith(dup):
                pair += times
                #Write bed files
                st = int(key.split("\t")[1]) - 9
                ed = int(key.split("\t")[1]) 
                bed.write(key.split("\t")[0] + "\t" + str(st) + "\t" + str(ed) + "\t" + str(times) + "\n")
    bed.close()
    
    print("Totally {} reads were paired around Tn5 cutting sites".format(pair))  

def main():   
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--BED',action = 'store', type=str, dest = 'BED',required=True,
                        help = 'BED file, better deduplicated!')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    run(options)

if __name__ == '__main__':
    main()

