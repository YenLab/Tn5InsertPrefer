#!/bin/bash
##########################################################################
#This pipeline is for ATAC-seq analysis, featured on Tn5 bias-correction.
#All dependency can be installed through conda using BiasFreeATAC.yaml 

#Created by Houyu Zhang
#Version 1.0
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2021 __YenLab@SCUT__. All rights reserved.
##########################################################################
set -ue
set -o pipefail

export LC_ALL=C

# help message
usage() {
	echo -e "Usage: $0 [options]"
	echo -e ""
	echo -e "-r <string>		If start from raw fastq file, please provide Read1 (_R1.fastq.gz) with this parameter."
	echo -e "-b <string>		If start from mapped bam file, please provide with this parameter and you can miss the -r/-m/-i parameter."
	echo -e "-m <string>		Indicate which mapper you wish to use (Bowtie2/BWA)"	
	echo -e "-i <string>		Location of mapping index."
	echo -e "-g <string>		Genome size for each chromosome."
	echo -e "-f <string>		Genome reference fasta file."
	echo -e "-t <string>		tallymer mappability file [optional]."
	echo -e "-l <string>		Blacklist regions [optional]."
	echo -e "-p <string>		Thresholds for parallelly run this pipeline."
	echo -e "-o <string>		Working directory, all output files will be generated here."
	exit 1
}

Tallymer=""
FqFile=""
BamFile=""

while getopts ":r:b:m:i:g:f:t:l:p:o:" op; do
	case $op in
		r) FqFile=${OPTARG} ;;
		b) BamFile=${OPTARG} ;;
		m) MAPPER=${OPTARG} ;;
		i) MAPPING_INDEXES=${OPTARG} ;;
		g) GENOMESIZE=${OPTARG} ;;
		f) FASTA_FILE=${OPTARG} ;;
		t) Tallymer=${OPTARG} ;;
		l) blacklist=${OPTARG} ;;
		p) RUN_THREASHOLD=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# make directories if not exist
[ ! -d ${WDIR} ] && mkdir ${WDIR}
cd ${WDIR}

#=========================================================================================
# Pipeline main body
#=========================================================================================
if [[ -n ${FqFile} ]]; then
	# Get file prefixes
	FILE_PATH=`dirname ${FqFile}`
	FILE_PREFIX=`basename ${FqFile} _R1.fastq.gz`
	echo "You provied fastq file (${FqFile}), we will work from reads quality check ......"
	echo "------------------------------->>> Start processing ${FILE_PREFIX} <<<-------------------------------"
	
	echo "Step1. Checking quality of raw data ......"
	if [[ ! -f ${FILE_PATH}/${FILE_PREFIX}_R2.fastq.gz ]];then
		echo "Sorry, BiasFreeATAC can't find you fastq Read2 file (${FILE_PATH}/${FILE_PREFIX}_R2.fastq.gz), please check!!!"
		exit -1
	fi
	fastqc --threads ${RUN_THREASHOLD} --quiet --outdir ./ ${FILE_PATH}/${FILE_PREFIX}_R1.fastq.gz ${FILE_PATH}/${FILE_PREFIX}_R2.fastq.gz 
	
	echo "Step2. Trimming Tn5 adapters ......"
	
#Output the Tn5 adapter sequence
cat > NexteraPE-PE.fa <<__EOF__
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
__EOF__

	trimmomatic PE -threads ${RUN_THREASHOLD} -phred33 \
	${FILE_PATH}/${FILE_PREFIX}_R1.fastq.gz ${FILE_PATH}/${FILE_PREFIX}_R2.fastq.gz \
	${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R1_unpaired.fastq \
	${FILE_PREFIX}_R2_trimmed.fastq ${FILE_PREFIX}_R2_unpaired.fastq \
	ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:1:true MINLEN:20 

	rm -f ${FILE_PREFIX}_R1_unpaired.fastq ${FILE_PREFIX}_R2_unpaired.fastq
	pigz -p ${RUN_THREASHOLD} ${FILE_PREFIX}_R1_trimmed.fastq ${FILE_PREFIX}_R2_trimmed.fastq
	
	echo "Step3. Checking quality of trimmed data ......"
	fastqc --threads ${RUN_THREASHOLD} --quiet --outdir ./ ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz
	
	echo "Step4. Aligning trimmed reads to reference genome and do initial filtering ......"
	# Remove read unmapped (4), mate unmapped (8), not primary alignment (256), read fails platform/vendor quality checks (512)
	# Remove low MAPQ reads (q < 30) and obtain sorted BAM file
	if [[ ${MAPPER} == "BWA" ]];then
		bwa mem ${MAPPING_INDEXES} -t ${RUN_THREASHOLD} \
		<(zcat ${FILE_PREFIX}_R1_trimmed.fastq.gz) <(zcat ${FILE_PREFIX}_R2_trimmed.fastq.gz) \
		| samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b - \
		| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam
	
	elif [[ ${MAPPER} == "Bowtie2" ]];then
		bowtie2 --threads ${RUN_THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${FILE_PREFIX}.alignMetrics.txt -q \
		-x ${MAPPING_INDEXES} -1 ${FILE_PREFIX}_R1_trimmed.fastq.gz -2 ${FILE_PREFIX}_R2_trimmed.fastq.gz \
		| samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b - \
		| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam
	fi
	
	samtools flagstat -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.bam > ${FILE_PREFIX}.filtered.flagstat.qc
	
	echo "Step5. Getting Insert Size distribution ......"
	picard CollectInsertSizeMetrics --METRIC_ACCUMULATION_LEVEL ALL_READS --INPUT ${FILE_PREFIX}.filtered.bam \
	--OUTPUT ${FILE_PREFIX}.filtered.insertSizes.txt --Histogram_FILE ${FILE_PREFIX}.filtered.insertSizes.pdf
	
	echo "Step6. Marking PCR duplicates ......" 
	picard MarkDuplicates -INPUT ${FILE_PREFIX}.filtered.bam --OUTPUT ${FILE_PREFIX}.filtered.dupmark.bam \
	--METRICS_FILE ${FILE_PREFIX}.filtered.deduplicate.qc --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false
	samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dupmark.bam
	
	echo "Step7. Removing duplicates, Index final position sorted BAM file ......"
	# remove read is PCR or optical duplicate (1024)
	samtools view -@ ${RUN_THREASHOLD} -F 1024 -b ${FILE_PREFIX}.filtered.dupmark.bam > ${FILE_PREFIX}.filtered.dedup.bam
	samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam
	samtools flagstat -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.filtered.dedup.bam > ${FILE_PREFIX}.filtered.dedup.flagstat.qc
	
	rm -f ${FILE_PREFIX}.filtered.bam* ${FILE_PREFIX}.filtered.dupmark.bam*
	
	echo "Step8. Shift the reads in bam file to get Tn5 insertion sites..."
	#Adjust 4 or 5bp for the Tn5 binding site and get cut sites bam file
	bedtools bamtobed -i ${FILE_PREFIX}.filtered.dedup.bam \
	| awk -v FS="\t" -v OFS="\t" '{ if ($6 == "+") {print $1,$2+4,$2+4+1,$4,$5,$6} else if ($6 == "-") {print $1,$3-5-1,$3-5,$4,$5,$6} }' \
	| bedtools bedtobam -i - -g ${GENOMESIZE} \
	| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.InsertSites.bam
	
	samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.InsertSites.bam

	[ ! -d ${WDIR}/metrics ] && mkdir -p ${WDIR}/metrics
	mv *zip *html *insertSizes.txt *insertSizes.pdf *alignMetrics.txt *.flagstat.qc *deduplicate.qc ./metrics
fi

if [[ -n ${BamFile} ]]; then
	BAM=${BamFile}
	echo "You provied bam file (${BamFile}), we will direct work on bias correction ......"
	FILE_PREFIX=`basename ${BamFile} .bam`
	echo "------------------------------->>> Start processing ${FILE_PREFIX} <<<-------------------------------"
else
	BAM=${FILE_PREFIX}.filtered.dedup.bam
fi

echo "Step9. Correcting Tn5 bias using nucleotide dependency information ......"
#Since the Tn5 can recognize a total 19bp region, we use --kmer-size=19,
#The cut sites is 5bp offset the 19bp regions, so we use --plus-offset=5 --minus-offset=5
#The mask scheme is slightly adapted from (Martins et al., 2017).
#After correction, we shifted the read position to count Tn5 insertion center following (Buenrostro et al., 2013)
#The read length was used to calculate genome-wide mappability, 36 generally perform well (Derrien et al., 2012).
ReadLen=36
#Returned uncorrected ATAC-seq signals
seqOutBias ${FASTA_FILE} ${BAM} --skip-bed --read-size=${ReadLen} --kmer-size=19 --plus-offset=5 --minus-offset=5 \
--strand-specific --custom-shift=4,5 --no-scale --bw=${FILE_PREFIX}_uncorrected.bigWig 

bigWigToBedGraph ${FILE_PREFIX}_uncorrected.bigWig /dev/stdout \
| grep "^chr" | sort --parallel=${RUN_THREASHOLD} -k1,1 -k2,2n > ${FILE_PREFIX}_uncorrected.bedGraph

#Returned corrected ATAC-seq signals
kmer_mask=XNXXXCXXNNXNNNXXNNX
seqOutBias ${FASTA_FILE} ${BAM} --skip-bed --read-size=${ReadLen} --kmer-size=19 --plus-offset=5 --minus-offset=5 \
--strand-specific --custom-shift=4,5 --kmer-mask ${kmer_mask} --bw=${FILE_PREFIX}_corrected.bigWig 

bigWigToBedGraph ${FILE_PREFIX}_corrected.bigWig /dev/stdout \
| grep "^chr" | sort --parallel=${RUN_THREASHOLD} -k1,1 -k2,2n > ${FILE_PREFIX}_corrected.bedGraph

echo "Step10. Peakcalling using uncorrected and corrected signals ......"
GENOME_SIZE_NUM=`awk '{SUM+=$2}END{print SUM}' ${GENOMESIZE}`

macs2 callpeak --broad --treatment ${FILE_PREFIX}_uncorrected.bedGraph --format BED --gsize ${GENOME_SIZE_NUM} --qvalue 0.01 --broad-cutoff 0.01 \
--outdir ${FILE_PREFIX}_uncorrected_peaks --name ${FILE_PREFIX}_uncorrected --bdg --nomodel --max-gap 100 --shift -100 --extsize 200 

macs2 callpeak --broad --treatment ${FILE_PREFIX}_corrected.bedGraph --format BED --gsize ${GENOME_SIZE_NUM} --qvalue 0.01 --broad-cutoff 0.01 \
--outdir ${FILE_PREFIX}_corrected_peaks --name ${FILE_PREFIX}_corrected --bdg --nomodel --max-gap 100 --shift -100 --extsize 200 

uncorrectedBroadPeak=${FILE_PREFIX}_uncorrected_peaks/${FILE_PREFIX}_uncorrected_peaks.broadPeak
correctedBroadPeak=${FILE_PREFIX}_corrected_peaks/${FILE_PREFIX}_corrected_peaks.broadPeak

bedtools intersect -a ${uncorrectedBroadPeak} -b ${correctedBroadPeak} -wa -v > ${FILE_PREFIX}_uncorrected_specific_peaks.bed
bedtools intersect -a ${uncorrectedBroadPeak} -b ${correctedBroadPeak} -wa    > ${FILE_PREFIX}_uncorrected_shared_peaks.bed
bedtools intersect -b ${uncorrectedBroadPeak} -a ${correctedBroadPeak} -v -wa > ${FILE_PREFIX}_corrected_specific_peaks.bed
bedtools intersect -b ${uncorrectedBroadPeak} -a ${correctedBroadPeak} -wa    > ${FILE_PREFIX}_corrected_share_peaks.bed

echo "------------------------------->>> Done processing ${FILE_PREFIX} <<<-------------------------------"