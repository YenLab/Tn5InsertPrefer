# -*- coding: utf-8 -*-
Description = '''
This scripts can check whether your bed file is valid (start > 0; start < end).

I strongly recommend you extend bed files using bedtools slop, 
however, if inevasible unvalid bed files occur, use this will help you!

'''.lstrip()

Copyright = '''
***
Created by Houyu Zhang on 2019-06-03.
Issue report on Hughiez047@gmail.com
Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import argparse

def run(options):
    
    print("processing {} ...".format(options.BED))
    
    total = 0
    start_st0 = 0
    start_lt_end = 0
    
    OUT = open(options.BED + ".alt","w+")
    with open (options.BED) as bed:
        for line in bed:
            total += 1
            lines = line.strip().split("\t")
            if int(lines[1]) < 0:
                start_st0 += 1
                continue
            if int(lines[1]) >= int(lines[2]):
                start_lt_end += 1
                continue
            OUT.write(line)
    removed = total - start_st0 - start_lt_end
    print("Totally {} sequences as input \n {} were removed because of start < 0\n \
    {} were removed because of start >= end\nTotally {} sequences as output".format(total,start_st0,start_lt_end,removed))
    
    OUT.close()
    
def main():   
    parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]', 
                                     description = Description, epilog = Copyright,
                                     formatter_class = argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--BED',action = 'store', type=str, dest = 'BED',required=True,
                        help = 'BED file')

    options = parser.parse_args()
    
    if not options:
        parser.print_help()
        sys.exit(3)

    run(options)

if __name__ == '__main__':
    main()
    
