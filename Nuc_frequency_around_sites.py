# -*- coding: utf-8 -*-
Description = '''
This scripts can calculate the single-base sites enrichment in given genomic features
'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-04-30.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import argparse
import pybedtools 
import numpy as np
#import matplotlib.pyplot as plt
'''
def check_bed_boundary(options,bedfile):
    chrom_dict ={}
    with open(options.GENOME) as cs:
        for line in cs:
            line = line.strip().split("\t")
            chrom_dict[line[0]] = line[1]

    s = ""
    with open(bedfile.fn) as bf:
        for line in bf:
            line = line.strip().split()
            if int(line[1]) <= int(chrom_dict[line[0]]):
                s = s + "\t".join(line) + "\n"
    
    a = pybedtools.BedTool(s, from_string=True)
    return a
'''

def run(options):
    '''
    This is main function for run 
    '''
    #Read sites bed file
    sites = pybedtools.BedTool(options.BED_SITES)
    print("Totally {} sites for processing...".format(sites.count()))
    
    #Read genome fasta file for getfasta
    fasta = pybedtools.BedTool(options.FASTA)
    
    #shift file by given length
    subset = sites.filter(lambda b: b.chrom != "MtDNA")
    shifted_sites = subset.slop(l=options.EXT_LEFT, r=options.EXT_RIGHT, genome=options.GENOME)
    #a = len(shifted_sites)
    
    #check bad lines
    print("Removing bad lines and line of chrMT in extended bed file...")
    shifted_sites.remove_invalid().truncate_to_chrom(options.GENOME)
    
    '''
    shifted_sites = check_bed_boundary(shifted_sites)
    b = len(shifted_sites)
    print("{} sites removed, and {} sites for processing...".format(b-a,b))
    '''
    #Get fasta file for extended sites region
    print("Calculating nucleotide frequency...")
    shifted_sites_fasta = shifted_sites.getfasta(fi=fasta,fo=options.OUT + ".fa")
    
    #Convert fasta file into a big numpy chararray  
    total_base = options.EXT_LEFT + options.EXT_RIGHT + 1
    fasta_list = np.chararray((shifted_sites_fasta.count(), total_base),unicode=True)
    
    index = 0
    for seq in open(shifted_sites_fasta.seqfn):
        if not seq.startswith(">"):
            seq = list(seq.strip())
            if len(seq) == total_base:
                fasta_list[index,] = seq
                index += 1
    
    fasta_list = fasta_list[0:index,]
    #Calculate ATCG frequency
    res = np.zeros((4,total_base))
    for i in range(total_base):
        seq = fasta_list[:,i]
        for index,char in enumerate("ACTG"):
            res[index,i] = sum(seq.count(char))/len(seq)
    
    np.savetxt(options.OUT + "_nuc.txt", res, fmt='%.4f', delimiter="\t")
    
    #Plot frequency
'''
    plt.plot(, res[1,:], label='A')
    plt.plot(x, res[2,:], label='C')
    plt.plot(x, res[3,:], label='T')
    plt.plot(x, res[4,:], label='G')
    plt.xlabel('Position releative to Sites')
    plt.ylabel('Nucleotide frequency')

    plt.legend()

    plt.show()
'''
    
def main():   
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-b', '--BED_SITES',action = 'store', type=str, dest = 'BED_SITES',required=True,
                        help = 'Sites for research in bed format')
    parser.add_argument('-g', '--GENOME',action = 'store', type=str, dest = 'GENOME',required=True,
                        help = 'Genome for this data [mm10,hg19,ce11...]')
    parser.add_argument('-f', '--FASTA',action = 'store', type=str, dest = 'FASTA',required=True,
                        help = 'Whole genome fasta file')
    parser.add_argument('-L', '--EXT_LEFT',action = 'store', type=int, dest = 'EXT_LEFT',required=True,
                        help = 'Extend bps to the left of sites')
    parser.add_argument('-R', '--EXT_RIGHT',action = 'store', type=int, dest = 'EXT_RIGHT',required=True,
                        help = 'Extend bps to the right of sites')
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
