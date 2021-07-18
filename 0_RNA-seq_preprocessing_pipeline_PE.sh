#!/bin/bash
##########################################################################
#This script takes a directory of Pair-end fastq RNA-seq reads, 
#Creat by Houyu Zhang 
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
#USAGE EXAMPLE (Note using absolute path!!!)
#bash this.sh /public/home/apoptosis/ATAC_input/ce/rawdata/test_pipeline/rawdata/ /public/home/apoptosis/ATAC_input/ce/rawdata/test_pipeline/mapping/ mm10

if [ -n "$1" ]; then
    echo "--->>> You provided the arguments: $@, Let's do it..."
else
    echo "USEAGE: ATAC-seq read trimming, mapping, filtering, adjusting, trimming to insert sites:" 1>&2
    echo "This program accepts the original PE fastq files and takes care of everything else up until Tn5 cut sites files generated." 1>&2
    echo "$(basename $0) [FASTQDir] [BAMDir] [Genome]"  1>&2
    echo "[FASTQDir]: absolute directory path containing FASTQ files (with format: *_R1.fastq and *_R2.fastq)" 1>&2
    echo "[BAMDir]: absolute directory path to put result files" 1>&2
    echo "[Genome]: genome name [ce11, mm10, hg38, dm6, tair10, danRer11, pfa2, zm3] or alias" 1>&2
    exit 1
fi

# =======================================================================================
# 0. Set up user local enviroment
# =======================================================================================
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info ${3}

RUN_THREASHOLD=30
trimmomatic_path="/public/home/zhy/Packages/trimmomatic/bin"

FASTQ_DIR=`echo $1 | sed 's:/$::g'`
BAM_DIR=`echo $2 | sed 's:/$::g'`

# make directories if not exist
[ ! -d ${FASTQ_DIR}/Raw ] && mkdir ${FASTQ_DIR}/Raw
[ ! -d ${FASTQ_DIR}/fastqc ] && mkdir ${FASTQ_DIR}/fastqc
[ ! -d ${BAM_DIR}/ ] && mkdir ${BAM_DIR}/
[ ! -d ${BAM_DIR}/metrics ] && mkdir ${BAM_DIR}/metrics

echo "--->>> Local environment set done..."
# =======================================================================================
# 1. Trim adapter sequences, mapping and post-mapping process
# =======================================================================================
cd ${FASTQ_DIR}

for FIRST_FQ in *R1.fastq.gz
do
	# Get file names
	FILE_PREFIX=`basename ${FIRST_FQ} _R1.fastq.gz`
	echo "--->>> Start processing ${FILE_PREFIX}..."
	
	# =============
	#quality check before trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz
	
	#Trimming adapters
	echo "Trimming adapters by Trimmomatic..."
	java -jar ${trimmomatic_path}/trimmomatic.jar PE -threads ${RUN_THREASHOLD} -phred33 ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz \
	${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_trimmed.fastq ${FILE_PREFIX}_R2_unpaired.fastq \
	ILLUMINACLIP:${trimmomatic_path}/adapters/TruSeq3-PE-2.fa:2:30:10:1:true MINLEN:36  
	
	rm -f ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_unpaired.fastq
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R2_trimmed.fastq
	
	#quality check after trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz

	# Map with STAR
	echo "Aligning by STAR..."
	STAR --runThreadN ${RUN_THREASHOLD} --genomeDir ${STAR_INDEXES} --sjdbGTFfile ${GTF_FILE} --limitBAMsortRAM 300000000 \
	--readFilesCommand zcat --readFilesIn ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz \
	--outFileNamePrefix ${FILE_PREFIX}. --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20
	
	echo "Removing un-paired and poor quality reads by samtools..."  
    samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b ${FILE_PREFIX}.Aligned.sortedByCoord.out.bam \
    | samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam
    samtools index ${FILE_PREFIX}.filtered.bam
    
    tail -n +5 ${FILE_PREFIX}.ReadsPerGene.out.tab | cut -f 1,2 | sort -k1,1 > ${FILE_PREFIX}_Rawcount.txt

	mv *zip *html ${FASTQ_DIR}/fastqc
	mv ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz ${FASTQ_DIR}/Raw
    mv ${FILE_PREFIX}.Log.progress.out ${FILE_PREFIX}.Log.out ${BAM_DIR}/metrics 
done
