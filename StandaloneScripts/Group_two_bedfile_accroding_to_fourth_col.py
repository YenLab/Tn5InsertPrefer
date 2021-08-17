# -*- coding: utf-8 -*-

Description = '''
Given two bedgraph files, this script can group them to different groups according to score range
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2020-06-24.
Issue report on Hughiez047@gmail.com
Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
'''.lstrip()

import sys, argparse

def group_bg(options):
    
    out_prefix = options.bed1.strip(".bed") + "_" + options.bed2.strip(".bed") + "_"
    f0 = open(out_prefix + "Allzeros.bed","w+")
    f1 = open(out_prefix + options.bed1 + "_bigger.bed","w+")
    f2 = open(out_prefix + options.bed2 + "_bigger.bed","w+")
    fe = open(out_prefix + "Similar.bed","w+")
    
    with open(options.bed1) as textfile1, open(options.bed2) as textfile2: 
        for x, y in zip(textfile1, textfile2):
            feature1 = x.strip().split("\t")
            feature2 = y.strip().split("\t")
            v1 = float(feature1[3])
            v2 = float(feature2[3])
            range_value = v1*options.Similarity
            out = feature1[0] + "\t" + str(feature1[1]) + "\t" + str(feature1[2]) + "\t" + str(feature1[3])  + "\t" + str(feature2[3] + "\n")  
            
            if (v1 == 0 and v2 == 0):
                f0.write(out)
            if (v1 - range_value < v2 < v1 + range_value):
                fe.write(out)
            if (v1 + range_value < v2):
                f2.write(out)
            if (v1 - range_value > v2):
                f1.write(out)            
        
    print("Finished") 
        
def main():   
    parser = argparse.ArgumentParser(prog = sys.argv[0], usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-B', '--bed1',action = 'store', type=str, 
                        help = 'A bedgraph file as reference file')
    parser.add_argument('-b', '--bed2',action = 'store', type=str, 
                        help = 'A bedgraph file for sort')
    parser.add_argument('-s', '--Similarity',action = 'store', type=float, 
                        help = 'Define a percentage for regard feature in two bed file as a group')

    options = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    group_bg(options)

if __name__ == '__main__':
    main()