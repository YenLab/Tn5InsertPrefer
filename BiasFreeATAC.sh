#!/bin/bash
##########################################################################
#This pipeline is for ATAC-seq analysis, featured on Tn5 bias-correction.
#All dependency can be installed through conda using BiasFreeATAC.yaml 

#Created by Houyu Zhang
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
	echo -e "-r <string>		ATAC-seq Read1 fastq file (_R1.fastq.gz)"
	echo -e "-i <string>		Location of Bowtie2 mapping index"
	echo -e "-g <string>		Genome size for each chromosome"
	echo -e "-f <string>		Genome reference fasta file"
	echo -e "-l <string>		Blacklist regions (provide if avaliable)"
	echo -e "-p <string>		Thresholds for parallelly run this pipeline"
	echo -e "-o <string>		Working directory, all output files will be generated here"
	exit 1
}

while getopts ":r:i:g:f:l:p:o:" op; do
	case $op in
		r) FqFile=${OPTARG} ;;
		i) BOWTIE2_INDEXES=${OPTARG} ;;
		g) GENOMESIZE=${OPTARG} ;;
		f) FASTA_FILE=${OPTARG} ;;
		l) blacklist=${OPTARG} ;;
		p) RUN_THREASHOLD=${OPTARG} ;;
		o) WDIR=${OPTARG} ;;
		*) usage ;;
	esac
done
shift $((OPTIND-1))

# check parameters
if [[ -z $FqFile ]] || [[ -z $BOWTIE2_INDEXES ]] || [[ -z $GENOMESIZE ]] || [[ -z $FASTA_FILE ]] || [[ -z $RUN_THREASHOLD ]]|| [[ -z $WDIR ]]; then
	echo -e "Input parameters missing, please specify each clearly..."
	usage
	exit -1
fi

# make directories if not exist
[ ! -d ${WDIR} ] && mkdir ${WDIR}
[ ! -d ${WDIR}/metrics ] && mkdir -p ${WDIR}/metrics

cd ${WDIR}
# Get file prefixes
FILE_PATH=`dirname ${FqFile}`
FILE_PREFIX=`basename ${FqFile} _R1.fastq.gz`

echo "------------------>>> Start processing ${FILE_PREFIX} <<<------------------"

echo "Step1. Checking quality of raw data ......"
fastqc -t 2 -q ${FILE_PATH}/${FILE_PREFIX}_R1.fastq.gz ${FILE_PATH}/${FILE_PREFIX}_R2.fastq.gz

echo "Step2. Trimming Tn5 adapters ......"

cat >> NexteraPE-PE.fa <<__EOF__
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
fastqc -t 2 -q ${FILE_PREFIX}_R1_trimmed.fastq.gz ${FILE_PREFIX}_R2_trimmed.fastq.gz

echo "Step4. Aligning trimmed reads to reference genome and do initial filtering ......"
# Remove read unmapped (4), mate unmapped (8), not primary alignment (256), read fails platform/vendor quality checks (512)
# Remove low MAPQ reads (q < 30) and obtain name sorted BAM file
bowtie2 --threads ${RUN_THREASHOLD} --end-to-end --no-mixed --maxins 2000 --met 1 --met-file ${FILE_PREFIX}.alignMetrics.txt -q \
-x ${BOWTIE2_INDEXES} -1 ${FILE_PREFIX}_R1_trimmed.fastq.gz -2 ${FILE_PREFIX}_R2_trimmed.fastq.gz \
| samtools view -@ ${RUN_THREASHOLD} -F 780 -q 30 -b - \
| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.filtered.bam

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
#Adjust 4 or 5bp for the Tn5 binding site, and finally trim to the 5prime end of the read.
bedtools bamtobed -i ${FILE_PREFIX}.filtered.dedup.bam \
| awk -v FS="\t" -v OFS="\t" '{ if ($6 == "+") {print $1,$2+4,$2+4+1,$4,$5,$6} else if ($6 == "-") {print $1,$3-5-1,$3-5,$4,$5,$6} }' \
| bedtools bedtobam -i - -g ${GENOMESIZE} \
| samtools sort -@ ${RUN_THREASHOLD} - > ${FILE_PREFIX}.InsertSites.bam

samtools index -@ ${RUN_THREASHOLD} ${FILE_PREFIX}.InsertSites.bam

echo "Step9. Correcting Tn5 bias (DNA motif and DNA shape) ......"
#This mask schme is from (Martins et al., NAR, 2017)
plus_mask=NXNXXXCXXNNXNNNXXNN
minus_mask=NNXXNNNXNNXXCXXXNXN

#Split bam file and do bias correction on each strand
samtools view -@ ${RUN_THREASHOLD} -bh -F 20 ${FILE_PREFIX}.filtered.dedup.bam > ${FILE_PREFIX}_plus.bam
samtools view -@ ${RUN_THREASHOLD} -bh -f 0x10 ${FILE_PREFIX}.filtered.dedup.bam > ${FILE_PREFIX}_minus.bam

ReadLen=`samtools view ${FILE_PREFIX}_plus.bam | cut -f 6 | sort --parallel=${RUN_THREASHOLD} \
| uniq -c | sort -k1nr | head -1 | cut -d ' ' -f4 | sed "s/M//"`

#Returned uncorrected ATAC-seq signals
seqOutBias ${FASTA_FILE} ${FILE_PREFIX}_minus.bam --no-scale --bw=${FILE_PREFIX}_uncorrected_minus.bigWig --shift-counts --read-size=${ReadLen}
seqOutBias ${FASTA_FILE} ${FILE_PREFIX}_plus.bam --no-scale --bw=${FILE_PREFIX}_uncorrected_plus.bigWig --shift-counts --read-size=${ReadLen}
bigWigMerge ${FILE_PREFIX}_uncorrected_minus.bigWig ${FILE_PREFIX}_uncorrected_plus.bigWig /dev/stdout \
sort --parallel=${RUN_THREASHOLD} -k1,1 -k2,2n > ${FILE_PREFIX}_uncorrected.bedGraph
bedGraphToBigWig ${FILE_PREFIX}_uncorrected.bedGraph ${FILE_PREFIX}_uncorrected.bigWig

rm -f ${FILE_PREFIX}_uncorrected_minus.bigWig ${FILE_PREFIX}_uncorrected_plus.bigWig 

#Returned corrected ATAC-seq signals
seqOutBias ${FASTA_FILE} ${FILE_PREFIX}_plus.bam --kmer-mask ${plus_mask} --bw=${FILE_PREFIX}_corrected_plus.bigWig --shift-counts --read-size=${ReadLen}
seqOutBias ${FASTA_FILE} ${FILE_PREFIX}_minus.bam --kmer-mask ${minus_mask} --bw=${FILE_PREFIX}_corrected_minus.bigWig --shift-counts --read-size=${ReadLen}
bigWigMerge ${FILE_PREFIX}_corrected_plus.bigWig ${FILE_PREFIX}_corrected_minus.bigWig /dev/stdout \
sort --parallel=${RUN_THREASHOLD} -k1,1 -k2,2n > ${FILE_PREFIX}_corrected.bedGraph
bedGraphToBigWig ${FILE_PREFIX}_corrected.bedGraph ${FILE_PREFIX}_corrected.bigWig

rm -f ${FILE_PREFIX}_corrected_minus.bigWig ${FILE_PREFIX}_corrected_plus.bigWig

echo "Step10. Peakcalling using uncorrected and corrected signals ......"
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

echo "Step11. Arranging file locations Tn5 bias ......"
mv *zip *html *insertSizes.txt *insertSizes.pdf *alignMetrics.txt *.flagstat.qc ${WDIR}/metrics

echo "------------------>>> Done processing ${FILE_PREFIX} <<<------------------"
