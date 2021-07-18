# -*- coding: utf-8 -*-
Description = '''
This scripts can calculate the GC frequency of input fasta file.
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-05-21.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import argparse
import numpy as np
from Bio import SeqIO  

def run(options):
    
    if options.SAMPLE:
        os.system("reformat.sh in="+ options.FASTA + " out=" + options.FASTA + "_sample.fa fastawrap=2000 ow=t samplerate=" + options.SAMPLE)
        FASTA = options.FASTA + "_sample.fa"
    else:
        FASTA = options.FASTA

    file_line = 0
    length = []
    with open(FASTA) as f:
        for line in f:
            length.append(len(line.strip()))
            file_line +=1
    
    seq_len = max(length)
    del length[:]
    
    print("--->>> Start parsing {} Fasta sequence...".format(round(file_line/2)))
    fasta_list = np.chararray((round(file_line/2), seq_len),unicode=True)
    index = 0
    for record in SeqIO.parse(FASTA, "fasta"):
        seq = list(record.seq)
        #print(seq)
        if len(seq) == seq_len:
                fasta_list[index,] = seq
                index += 1
        
    #print(index)
    print("--->>> Start calculating GC content...")
    res = np.zeros((4,seq_len))
    for i in range(seq_len):
        seq = fasta_list[:,i]
        for index,char in enumerate("ACTG"):
            res[index,i] = sum(seq.count(char))/(len(seq) - sum(seq.count("N")))
    
    np.savetxt(options.OUT + "_nuc.txt", res, fmt='%.4f', delimiter="\t")
    
    if options.SAMPLE:
        os.system("rm -f " + FASTA)
        
def main():   
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-f', '--FASTA',action = 'store', type=str, dest = 'FASTA',required=True,
                        help = 'Whole genome fasta file')
    parser.add_argument('-s', '--SAMPLE',action = 'store', type=str, dest = 'SAMPLE',
                        help = 'Sample fraction from fasta file')
    parser.add_argument('-o', '--OUT',action = 'store', type=str, dest = 'OUT',required=True,
                        help = 'Destination of result file prefix')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    run(options)
    print("Successful!!!")
    
if __name__ == '__main__':
    main()