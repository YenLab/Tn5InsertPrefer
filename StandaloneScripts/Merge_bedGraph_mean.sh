#!/bin/bash
#This script if for merging two bedGraph files by averaging the 4th value
#Usage bash this.sh file1 file2 finalfile

file1=$1
file2=$2
finalfile=$3

cat ${file1} ${file2} | awk '!a[$1$2$3]++' | sort -k1,1 -k2,2n | cut -f 1-3 > ${file1}_tmp0

bedtools intersect -a ${file1}_tmp0 -b ${file1} -wao | cut -f 1-3,7 | sed 's/\.$/0/g' > ${file1}_tmp1
bedtools intersect -a ${file1}_tmp0 -b ${file2} -wao | cut -f 1-3,7 | sed 's/\.$/0/g' > ${file1}_tmp2
paste tmp1 <(cut -f 4 tmp2) | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,($4+$5)/2}' > ${finalfile}

rm -f ${file1}_tmp0 ${file1}_tmp1 ${file1}_tmp2
echo "-->> Done for" $1 $2
