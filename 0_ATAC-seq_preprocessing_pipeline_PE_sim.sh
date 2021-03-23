#!/bin/bash
##########################################################################
#This script takes a directory of Pair-end fastq ATAC-seq reads, 
#1. Fastq quality control using [fastqc]
#2. Trims adapters using [Trmmomatic]
#3. Mapping using [Bowtie2]
#4. Mark PCR duplicates and generate fragment size using [Picard]
#5. Filter bam reads using [bamtools]
#6. Adjusts for the Tn5 binding sites, and creates bed and bam files of just the insert site (i.e. the 5' end of each read)
#7. Format converting for downstream analysis [bedtools, deeptools]

#Created by Houyu Zhang
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2019 __YenLab@SCUT__. All rights reserved.
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
#get genome information of species
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
[ ! -d ${BAM_DIR}/InsertSizeMetrics ] && mkdir ${BAM_DIR}/InsertSizeMetrics

echo "--->>> Local environment set done..."

# =======================================================================================
# 1. Trim adapter sequences, mapping and post-mapping process
# =======================================================================================
cd ${FASTQ_DIR}

for FIRST_FQ in *R1.fastq.gz
do
	# Get file prefixes
	FILE_PREFIX=`basename ${FIRST_FQ} _R1.fastq.gz`
	echo "--->>> Start processing ${FILE_PREFIX}..."

	# =============
	#echo "quality check before trimming adapters..."
	#fastqc -t 2 -q ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz
	
	echo "Trimming adapters by Trimmomatic..."
	java -jar ${trimmomatic_path}/trimmomatic.jar PE -threads ${RUN_THREASHOLD} -phred33 ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz \
	${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_trimmed.fastq ${FILE_PREFIX}_R2_unpaired.fastq \
	ILLUMINACLIP:${trimmomatic_path}/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:20  
	
	rm -f ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_unpaired.fastq
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R2_trimmed.fastq
	
	#echo "quality check after trimming adapters..."
	#fastqc -t 2 -q ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz
	#mv *zip *html ${FASTQ_DIR}/fastqc

	# =============
	echo "Aligning by bowtie2..."
	bowtie2 --threads ${RUN_THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${FILE_PREFIX}.alignMetrics.txt -q \
	-x ${BOWTIE2_INDEXES} -1 ${FILE_PREFIX}_R1_trimmed.fastq.gz -2 ${FILE_PREFIX}_R2_trimmed.fastq.gz \
	| samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b - \
	| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam

    echo "Marking PCR duplicates by picard..." 
    picard MarkDuplicates INPUT=${FILE_PREFIX}.filtered.bam OUTPUT=${FILE_PREFIX}.filtered.dupmark.bam METRICS_FILE=${FILE_PREFIX}.filtered.deduplicate.qc \
    VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
    samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dupmark.bam
    
    echo "Removing duplicates, Index final position sorted BAM file..."
    # remove read is PCR or optical duplicate (1024)
    samtools view -@ ${RUN_THREASHOLD} -F 1804 -b ${FILE_PREFIX}.filtered.dupmark.bam > ${FILE_PREFIX}.filtered.dedup.bam
    samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam
    samtools flagstat -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam > ${FILE_PREFIX}.filtered.dedup.flagstat.qc
    
    echo "Cleaning up FASTQ dir..."
    rm -f ${FILE_PREFIX}.filtered.bam 
	mv ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz ${FASTQ_DIR}/Raw
	mv ${FILE_PREFIX}.filtered.dedup.insertSizes.txt ${FILE_PREFIX}.filtered.dedup.insertSizes.pdf ${BAM_DIR}/InsertSizeMetrics
    mv ${FILE_PREFIX}.alignMetrics.txt ${FILE_PREFIX}.filtered.dedup.flagstat.qc ${FILE_PREFIX}.filtered.deduplicate.qc ${BAM_DIR}/metrics
done
