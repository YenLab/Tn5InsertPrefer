# -*- coding: utf-8 -*-
import getopt,sys

opts, args = getopt.getopt(sys.argv[1:], 'hi:o:l:r:t:T:')
for op, value in opts:
    if op == '-i':
        input_file = value
    if op == '-o':
        output_file = value
    if op == '-l':
        ext_left = int(value)
    if op == '-r':
        ext_right = int(value)
    if op == '-t':
        thread = int(value)
    if op == '-T':
        THREAD = int(value)        

sites_dict = {}
with open(input_file) as f:
    for line in f:
        line = line.strip().split()
        if "\t".join(line[0:2]) + "\t" + line[5] in sites_dict.keys():
            sites_dict["\t".join(line[0:2]) + "\t" + line[5]] += 1
        else:
            sites_dict["\t".join(line[0:2]) + "\t" + line[5]] = 1

out = open(output_file,"w+")

for k in sorted(sites_dict.keys()):
    #print(k)
    if(THREAD > sites_dict[k] > thread):
        new = k.strip().split("\t")
        st = str(int(new[1]) - ext_left)
        ed = str(int(new[1]) + ext_right)
        out.write(new[0] + "\t"+ st + "\t" + ed + "\t" + str(sites_dict[k]) +"\t" + new[2] + "\n")
out.close()
