#!/bin/bash

#This function recieve bw/wig/bedgraph files as input, and liftover them to correct bedgraph/bw files
#Usage: convert2corbw originalFile prefix hg38

source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh
originalFile=$1
finalFile=$2
Configuration_info $3
remove_pattern="KI|GL|^GL|^JH|^KI|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg"

#passed test
if [[ ${originalFile} == *bw || ${originalFile} == *bigWig ]];then
	echo "--->>> bigwig file ${originalFile} was detected! bigWigToBedGraph..."
	bigWigToBedGraph ${originalFile} ${originalFile}_tmp1

#passed test
elif [[ ${originalFile} == *bg || ${originalFile} == *bedGraph ]];then
	echo "--->>> bedGraph file ${originalFile} was detected! Remove header..."
	grep -v "track" ${originalFile} > ${originalFile}_tmp1

elif [[ ${originalFile} == *wig ]];then
	echo "--->>> Wig file ${originalFile} was detected! Convert2bed..."
	convert2bed --input=wig --output=bed < ${originalFile} | cut -f 1-3,5 > ${originalFile}_tmp1

#elif [[ ${originalFile} == *bed$ ]];then
#	echo "--->>> bed file ${originalFile} was detected!"
#	echo "--->>> liftover using ${Liftover_file}..."
#	liftOver ${originalFile} ${Liftover_file} ${originalFile}_tmp2 unlifted.bed
#	grep -vE ${remove_pattern} ${originalFile}_tmp2 | sort -k1,1 -k2,2n > ${finalFile}.bed
fi

echo "--->>> liftover using ${Liftover_file}..."
liftOver ${originalFile}_tmp1 ${Liftover_file} ${originalFile}_tmp2 unlifted.bed
echo "--->>> Remove patches, overlapping regions and sort to a bedGraph file..."
grep -vE ${remove_pattern} ${originalFile}_tmp2 \
| sort -k1,1 -k2,2n | bedmap --echo --count - | grep "|1$" | sed "s/|1//" > ${finalFile}.bedGraph
echo "--->>> bedGraphToBigWig..."
bedGraphToBigWig ${finalFile}.bedGraph ${GENOME_SIZE_chrM} ${finalFile}.bw

#Clear 
rm -f ${originalFile}_tmp1 ${originalFile}_tmp2 unlifted.bed

