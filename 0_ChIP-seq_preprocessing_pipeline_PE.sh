#!/bin/bash
##########################################################################
#This script takes a directory of Pair-end fastq ChIP-seq reads, 
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
	ILLUMINACLIP:${trimmomatic_path}/adapters/TruSeq3-PE-2.fa:2:30:10:1:true MINLEN:20 
	rm -f ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_unpaired.fastq
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R2_trimmed.fastq
	
	#quality check after trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz

	bowtie2 --threads ${RUN_THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${FILE_PREFIX}.alignMetrics.txt -q \
	-x ${BOWTIE2_INDEXES} -1 ${FILE_PREFIX}_R1_trimmed.fastq.gz -2 ${FILE_PREFIX}_R2_trimmed.fastq.gz \
	| samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b - \
	| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam

    # =============
    echo "Marking PCR duplicates by picard..." 
    picard MarkDuplicates INPUT=${FILE_PREFIX}.filtered.bam OUTPUT=${FILE_PREFIX}.filtered.dupmark.bam METRICS_FILE=${FILE_PREFIX}.filtered.deduplicate.qc \
    VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
    samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dupmark.bam
    
    echo "Removing duplicates, Index final position sorted BAM file..."
    # remove read is PCR or optical duplicate (1024)
    samtools view -@ ${RUN_THREASHOLD} -F 1804 -b ${FILE_PREFIX}.filtered.dupmark.bam > ${FILE_PREFIX}.filtered.dedup.bam
    samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam
    samtools flagstat -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam > ${FILE_PREFIX}.filtered.dedup.flagstat.qc

    #Convert to bigwig format
	bamCoverage --numberOfProcessors ${RUN_THREASHOLD} --bam ${FILE_PREFIX}.filtered.dedup.bam --outFileFormat bigwig \
	--outFileName ${FILE_PREFIX}.filtered.dedup.bw --binSize 1 \
	--normalizeUsing CPM --effectiveGenomeSize ${GENOME_SIZE_NUM} --blackListFileName ${blacklist} 

    mv ${FILE_PREFIX}_R1.fastq.gz ${FILE_PREFIX}_R2.fastq.gz ${FASTQ_DIR}/Raw
    mv *zip *html ${FASTQ_DIR}/fastqc
    rm -f ${FILE_PREFIX}.filtered.bam ${FILE_PREFIX}.filtered.dupmark.bam* 
    mv ${FILE_PREFIX}.alignMetrics.txt ${FILE_PREFIX}.filtered.deduplicate.qc ${FILE_PREFIX}.filtered.dedup.flagstat.qc ${BAM_DIR}/metrics
done
