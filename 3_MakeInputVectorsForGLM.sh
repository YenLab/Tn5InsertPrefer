#!/bin/bash
##########################################################################
#This script is for preparing input data for GLM 
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh
dir=rawdataML2
[ ! -d ${dir} ] && mkdir ${dir}
cd ${dir}

#ln -sf /public/home/zhy/Tn5_bias/pre-processing/mapping/all_sites/*t0Tinf_l0r1.bed .

for i in *t0Tinf_l0r1.bed
do
	base=`basename ${i} _t0Tinf_l0r1.bed`
	
	#Get species information for current sample
	id=`echo ${base} | sed "s/_.*//"`
	Configuration_info ${id} 
	echo "--->>> Start Process ${i} for ${GENOME_NAME} ..."

	sample_size=10000
	dl_sample_size=`expr ${sample_size} \* 2`
	tr_sample_size=`expr ${sample_size} \* 3`

	Prefix=${base}_shuf${sample_size}
	dl_Prefix=${base}_shuf${dl_sample_size}
	tr_Prefix=${base}_shuf${tr_sample_size}

	#Skip current sample if has been processed
	if [[ -f ${ft_Prefix}_ext25bp.bed ]]; then
		continue
	fi

	bedtools random -l 51 -n 10000000 -g ${GENOME_SIZE_chrM} > ${base}_random.bed

	if [[ ! $base == *input* ]]; then
		echo "Your input file is for Chromatin, will only use Tn5 cut sites in open peaks!!!"
		
		PeakFile=/public/home/zhy/Tn5_bias/pre-processing/PeakCalling/${base}/${base}_peaks.broadPeak 
		bedtools intersect -nonamecheck -a ${i} -b ${PeakFile} \
		| awk -v FS="\t" -v OFS="\t" '{if(($3-$2)==1) {print $0}}' > ${i}_Peak && mv ${i}_Peak ${i}
		bedtools intersect -nonamecheck -a ${base}_9bp_duplicates.bed -b ${PeakFile} \
		| awk -v FS="\t" -v OFS="\t" '{if(($3-$2)==9) {print $0}}' > ${base}_dup_Peak && mv ${base}_dup_Peak ${base}_9bp_duplicates.bed
		
		#grep -vE ${CONSENSUS_REMOVE} ${PeakFile} | bedtools slop -b 1000 -i - -g ${GENOME_SIZE_chrM} \
		#| bedtools intersect -nonamecheck -a <(bedtools random -l 51 -n 10000000 -g ${GENOME_SIZE_chrM}) -b - > ${base}_random.bed
	fi

	#Handle on negative sites file
	bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_random.bed \
	| awk -v RS=">" '!/N/{printf $0RT}' | grep ">" | sed 's/[:-]/\t/g' | sed 's/>//' \
	| awk -v FS="\t" -v OFS="\t" '{if(($3-$2) == 51) {print $0}}' \
	| awk '!a[$1$2$3]++' | head -${tr_sample_size} | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,"0","."}' > ${tr_Prefix}_uncut_ext25bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${tr_Prefix}_uncut_ext25bp.bed -fo ${tr_Prefix}_uncut_ext25bp.fa

	rm -f ${base}_random.bed
	#Remove blacklist because I saw some sites are extremely Tn5 preferred, which fall within blacklist, may bias our classical analysis
	#There are some patches and scaffolds need to be remove avoid misleading

	if [[ ${blacklist} ]]; then
		echo "Blacklist is found for ${sp}, remove them for better explanation..."
		bedtools intersect -nonamecheck -a ${i} -b ${blacklist} -v > tmp && mv tmp ${i}
		bedtools intersect -nonamecheck -a ${base}_9bp_duplicates.bed -b ${blacklist} -v > tmp && mv tmp ${base}_9bp_duplicates.bed
	else
		echo "No blacklist is available for ${sp}, skip filtering!"
	fi

	#Handle on cut sites file
	#remove chrM, patches and scaffolds
	grep -vE ${CONSENSUS_REMOVE} ${i} | grep -E "^chr" | grep -v "chrM" \
	| bedtools slop -b 25 -i - -g ${GENOME_SIZE_chrM} | bedtools getfasta -fi ${FASTA_FILE} -bed - \
	| awk -v RS=">" '!/N/{printf $0RT}' | grep '>' | sed 's/[:-]/\t/g' | sed 's/>//' \
	| awk -v FS="\t" -v OFS="\t" '{if(($3-$2) == 51) {print $0}}' | awk '!a[$1$2$3]++' \
	| shuf | head -${dl_sample_size} | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,"1","."}' > ${dl_Prefix}_all_ext25bp.bed
	bedtools getfasta -fi ${FASTA_FILE} -bed ${dl_Prefix}_all_ext25bp.bed -fo ${dl_Prefix}_all_ext25bp.fa

	#Merge cut and uncut files
	cat ${dl_Prefix}_all_ext25bp.bed ${tr_Prefix}_uncut_ext25bp.bed > ${base}_trainning_ext25bp.bed
	cat ${dl_Prefix}_all_ext25bp.fa ${tr_Prefix}_uncut_ext25bp.fa > ${base}_trainning_ext25bp.fa
	nohup python ~/Packages/BiasAway/BiasAway.py m -f ${base}_trainning_ext25bp.fa | sed "s/ >/>/" > ${base}_trainning_ext25bp_shuffled.fa &
done

ln -sf ~/Tn5_bias/scripts/generate_shapes.R .

for i in *_trainning_ext25bp.fa
do
	nohup Rscript generate_shapes.R -f ${i} &
done