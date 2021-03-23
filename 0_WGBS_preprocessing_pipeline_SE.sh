#!/bin/bash
##########################################################################
#This script takes a directory of Single-end fastq WGBS reads, 
#1. Fastq quality control using [fastqc]
#2. Trims adapters using [Trmmomatic]
#3. Mapping and methylation extraction using [Bismark2]

#Creat by Houyu Zhang
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
#USAGE EXAMPLE (Note using absolute path!!!)
#bash this.sh /public/home/apoptosis/ATAC_input/ce/rawdata/test_pipeline/rawdata/ /public/home/apoptosis/ATAC_input/ce/rawdata/test_pipeline/mapping/ mm10

if [ -n "$1" ]; then
    echo "--->>> You provided the arguments: $@, Let's do it..."
else
    echo "USEAGE: WGBS read trimming, mapping, filtering, methylation extraction:" 1>&2
    echo "$(basename $0) [FASTQDir] [BAMDir] [Genome]"  1>&2
    echo "[FASTQDir]: absolute directory path containing FASTQ files (with format: *_R1.fastq and *_R2.fastq)" 1>&2
    echo "[BAMDir]: absolute directory path to put result files" 1>&2
    echo "[Genome]: genome name [ce11, mm10, hg38, dm6, tair10, danRer11, pfa2, zm3] or alias" 1>&2
    exit 1
fi

# =======================================================================================
# 0. These need to be customized for each user's environment
# =======================================================================================
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh
Configuration_info ${3}

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

for FIRST_FQ in *.fastq.gz
do
	# Get file names
	FILE_PREFIX=`basename ${FIRST_FQ} .fastq.gz`
	echo "--->>> Start processing ${FILE_PREFIX}..."

	# =============
	#quality check before trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}.fastq.gz
	
	#Trimming adapters
	echo "Trimming adapters by Trimmomatic..."
	java -jar ${trimmomatic_path}/bin/trimmomatic.jar SE -threads ${RUN_THREASHOLD} -phred33 ${FILE_PREFIX}.fastq.gz \
	${FILE_PREFIX}_trimmed.fastq ILLUMINACLIP:${trimmomatic_path}/adapters/TruSeq3-SE.fa:2:30:10:1 MINLEN:36 
	rm -f ${FILE_PREFIX}_unpaired.fastq
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_trimmed.fastq
	
	#quality check after trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}_trimmed.fastq.gz
	mv *zip *html ${FASTQ_DIR}/fastqc

	# =============
	# Map with BisMark2
	bismark -p 8 --multicore 1 --single_end ${FILE_PREFIX}_trimmed.fastq.gz --non_directional --nucleotide_coverage --output_dir ./ \
	--genome ${sp_dir}

	deduplicate_bismark --bam --single  --output_dir ./ ${FILE_PREFIX}_trimmed_bismark_bt2.bam

	# --parallel 10 is likely to use around 30 cores of system resources.
	bismark_methylation_extractor --single-end --no_overlap --gzip --multicore 10 --buffer_size 30G --ignore 10 \
	--cytosine_report --output ./ --genome_folder ${sp_dir} ${FILE_PREFIX}_trimmed_bismark_bt2.deduplicated.bam

	bismark2report --dir ./
	bismark2summary ./

    echo "Cleaning up FASTQ dir..."
	mv ${FILE_PREFIX}.fastq.gz ${FASTQ_DIR}/Raw
	mv *html *txt ${BAM_DIR}/metrics
	mv *bismark_bt2* ${BAM_DIR}
done
