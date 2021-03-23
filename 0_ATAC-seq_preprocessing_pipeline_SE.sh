#!/bin/bash
##########################################################################
#This script takes a directory of SE fastq ATAC-seq reads, 
#1. Fastq quality control using [fastqc]
#2. Trims adapters using [Trmmomatic]
#3. Mapping using [Bowtie2]
#4. Mark PCR duplicates and generate fragment size using [Picard]
#5. Filter bam reads using [bamtools]
#6. Adjusts for the Tn5 binding footprint, and creates bed and bam files of just the insert site (i.e. the 5' end of each read)
#7. Format converting for downstream analysis [bedtools, deeptools]

#Created by Houyu Zhang on 2020-09-25.
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
    echo "[FASTQDir]: absolute directory path containing FASTQ files (with format: *.fastq.gz )" 1>&2
    echo "[BAMDir]: absolute directory path to put result files" 1>&2
    echo "[Genome]: genome name [ce11, mm10, hg38, dm6, tair10, danRer11, pfa2, zm3] or alias" 1>&2
    exit 1
fi

# =======================================================================================
# 0. Set up user local enviroment
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
[ ! -d ${BAM_DIR}/InsertSizeMetrics ] && mkdir ${BAM_DIR}/InsertSizeMetrics

echo "--->>> Local environment set done..."

#function HandleBed {
#	grep -vE <(echo ${CONSENSUS_REMOVE}) $1 | sort -k1,1 -k2,2n 
#}

# =======================================================================================
# 1. Trim adapter sequences, mapping and post-mapping process
# =======================================================================================
cd ${FASTQ_DIR}

for FIRST_FQ in *fastq.gz
do
	# Get file names
	FILE_PREFIX=`basename ${FIRST_FQ} .fastq.gz`
	echo "--->>> Start processing ${FILE_PREFIX}..."

	# =============
	#quality check before trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}.fastq.gz 
	
	#Trimming adapters
	echo "Trimming adapters by Trimmomatic..."
	java -jar ${trimmomatic_path}/trimmomatic.jar SE -threads ${RUN_THREASHOLD} -phred33 ${FILE_PREFIX}.fastq.gz \
	${FILE_PREFIX}_trimmed.fastq ILLUMINACLIP:${trimmomatic_path}/adapters/NexteraPE-PE.fa:2:30:10:1:true MINLEN:30  
	
	rm -f ${FILE_PREFIX}_unpaired.fastq 
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_trimmed.fastq 
	
	#quality check after trimming adapters
	fastqc -t 2 -q ${FILE_PREFIX}_trimmed.fastq.gz
	mv *zip *html ${FASTQ_DIR}/fastqc

	# =============
	# Map with bowite2
	echo "Aligning by bowtie2..."
	bowtie2 --threads ${RUN_THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${FILE_PREFIX}.alignMetrics.txt -q \
	-x ${BOWTIE2_INDEXES} -U ${FILE_PREFIX}_trimmed.fastq.gz > ${FILE_PREFIX}.sam

	# =============
	echo "Removing read pairs with bad CIGAR strings and sort by position..." 
	awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && $6!="*" { cigar=$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) print $1 ; }' ${FILE_PREFIX}.sam \
	| sort | uniq > badCigars.temp

	if [[ $(cat badCigars.temp | wc -l) -gt 0 ]]
	then
		cat ${FILE_PREFIX}.sam | grep -v -F -f badCigars.temp | samtools view -@ ${RUN_THREASHOLD} -Sb - | samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.bam
	else
		samtools view -@ ${RUN_THREASHOLD} -Sb ${FILE_PREFIX}.sam | samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.bam
	fi

	rm -f badCigars.temp ${FILE_PREFIX}.sam
	samtools flagstat -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.bam > ${FILE_PREFIX}.flagstat.qc

	# =============
	echo "Removing un-paired and poor quality reads by samtools..."  
	# Remove read unmapped (4), mate unmapped (8), not primary alignment (256), read fails platform/vendor quality checks (512)
    # Remove low MAPQ reads (q < 30) and obtain name sorted BAM file

    samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b ${FILE_PREFIX}.bam | samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam
    rm -f ${FILE_PREFIX}.bam

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
    
    echo "Getting Insert Size metrics..."
    picard CollectInsertSizeMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS INPUT=${FILE_PREFIX}.filtered.dedup.bam \
    OUTPUT=${FILE_PREFIX}.filtered.dedup.insertSizes.txt HISTOGRAM_FILE=${FILE_PREFIX}.filtered.dedup.insertSizes.pdf 
    
    # =============
    echo "Getting insert sites by shifting reads start sites..."
    bedtools bamtobed -i ${FILE_PREFIX}.filtered.dedup.bam > temp.bed

    #This is for find paired cut sites
  	samtools view -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam | cut -f 1,10 > insert_cutsites.temp
  	paste temp.bed insert_cutsites.temp | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$8,$6}' \
  	| grep -vE ${CONSENSUS_REMOVE} | sort --parallel ${RUN_THREASHOLD} -k1,1 -k2,2n > ${FILE_PREFIX}_forfindcutsites.bed

	samtools view -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam | cut -f 1,9 > insert.temp
	paste temp.bed insert.temp | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$8,$6}' > ${FILE_PREFIX}_tmp.bed
	#Adjust 4 or 5bp for the Tn5 binding site, and finally trim to the 5prime end of the read.
	awk -v FS="\t" -v OFS="\t" '{if ($6 == "+") { $2 = $2 + 4 } else if ($6 == "-") {$3 = $3 - 5} print $0}' ${FILE_PREFIX}_tmp.bed \
	| awk -v OFS="\t" '{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}' \
	| grep -vE ${CONSENSUS_REMOVE} | sort -k1,1 -k2,2n > ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bed	

	rm -f temp.bed insert.temp insert_cutsites.temp ${FILE_PREFIX}_tmp.bed
	
	# =============
	echo "Converting insertion sites into bam/bigwig format..."
	bedtools bedtobam -i ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bed -g ${GENOME_SIZE_chrM} > ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam
	samtools sort -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam > ${FILE_PREFIX}_tmp \
	&& mv ${FILE_PREFIX}_tmp ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam
    samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam

   	if [[ -n ${blacklist} ]]
   	then
		bamCoverage --numberOfProcessors ${RUN_THREASHOLD} --bam ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam --outFileFormat bigwig \
		--outFileName ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bw --binSize 1 --normalizeUsing CPM \
		--effectiveGenomeSize ${GENOME_SIZE_NUM} --blackListFileName ${blacklist}    
	else
		bamCoverage --numberOfProcessors ${RUN_THREASHOLD} --bam ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bam --outFileFormat bigwig \
		--outFileName ${FILE_PREFIX}.filtered.dedup.shifted.insertSites.bw --binSize 1 --normalizeUsing CPM \
		--effectiveGenomeSize ${GENOME_SIZE_NUM} 
	fi  

	# ============= 
	#This part is for generating Tn5 cut sites with duplicates
	dup_cutsites="false"
	if [[ ${dup_cutsites} == "True" ]]
	then
		bamToBed -i ${FILE_PREFIX}.filtered.dupmark.bam > temp.bed
		samtools view -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dupmark.bam | cut -f 1,9 > insert.temp
		paste temp.bed insert.temp | awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3,$4,$8,$6}' > ${FILE_PREFIX}_tmp.bed
		
		awk -v FS="\t" -v OFS="\t" '{if ($6 == "+") { $2 = $2 + 4 } else if ($6 == "-") {$3 = $3 - 5} print $0}' ${FILE_PREFIX}_tmp.bed \
		| awk -v OFS="\t" '{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}' | \
		grep -vE ${CONSENSUS_REMOVE} | sort -k1,1 -k2,2n > ${FILE_PREFIX}.filtered.dupmark.shifted.insertSites.bed

		bedtools bedtobam -i ${FILE_PREFIX}.filtered.dupmark.shifted.insertSites.bed -g ${GENOME_SIZE_chrM} > ${FILE_PREFIX}.filtered.dupmark.shifted.insertSites.bam
		samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dupmark.shifted.insertSites.bam

		rm -f temp.bed insert.temp ${FILE_PREFIX}_tmp.bed
	fi

    echo "Cleaning up FASTQ dir..."
	mv ${FILE_PREFIX}.fastq.gz ${FASTQ_DIR}/Raw
	mv ${FILE_PREFIX}.filtered.dedup.insertSizes.txt ${FILE_PREFIX}.filtered.dedup.insertSizes.pdf  ${BAM_DIR}/InsertSizeMetrics
	mv *bam *bam.bai *forfindcutsites.bed *shifted.insertSites* ${BAM_DIR}
    mv ${FILE_PREFIX}.alignMetrics.txt ${FILE_PREFIX}.flagstat.qc ${FILE_PREFIX}.filtered.dedup.flagstat.qc ${FILE_PREFIX}.filtered.deduplicate.qc ${BAM_DIR}/metrics
done
