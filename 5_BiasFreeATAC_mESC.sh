#!/bin/bash

#=========================================================================================
# Correct genome-wide Tn5 bias using seqOutBias
#=========================================================================================
#Since the Tn5 can recognize a total 19bp region, we use --kmer-size=19,
#The cut sites is 5bp offset the 19bp regions, so we use --plus-offset=5 --minus-offset=5
#The mask scheme is slightly adapted from (Martins et al., 2017).
#After correction, we shifted the read position to count Tn5 insertion center following (Buenrostro et al., 2013)
#The read length was used to calculate genome-wide mappability, 36 generally perform well (Derrien et al., 2012).
plus_mask=NXNXXXCXXNNXNNNXXN
minus_mask=NXXNNNXNNXXCXXXNXN
name=Mouse_ESC_rep1_ATACseq_2017Gary.filtered.dedup

source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10

samtools view -bh -F 20 ${name}.bam > ${name}_plus.bam
samtools view -bh -f 0x10 ${name}.bam > ${name}_minus.bam

seqOutBias ${FASTA_FILE} ${name}_plus.bam --kmer-mask ${plus_mask} --bw=${name}_plus_${plus_mask}-mer.bigWig --shift-counts --read-size=51
seqOutBias ${FASTA_FILE} ${name}_minus.bam --kmer-mask ${minus_mask} --bw=${name}_minus_${minus_mask}-mer.bigWig --shift-counts --read-size=51

seqOutBias ${FASTA_FILE} ${name}_minus.bam --no-scale --bw=${name}_no_scale_minus.bigWig --shift-counts --read-size=51
seqOutBias ${FASTA_FILE} ${name}_plus.bam --no-scale --bw=${name}_no_scale_plus.bigWig --shift-counts --read-size=51

bigWigMerge ${name}_plus_${plus_mask}-mer.bigWig ${name}_minus_${minus_mask}-mer.bigWig ${name}_${plus_mask}_${minus_mask}_merged.bedGraph
bigWigMerge ${name}_no_scale_plus.bigWig ${name}_no_scale_minus.bigWig ${name}_no_scale_merged.bedGraph

sort -k1,1 -k2,2n ${name}_${plus_mask}_${minus_mask}_merged.bedGraph | grep -vE "GL|JH" > ${name}_${plus_mask}_${minus_mask}_merged.sorted.bedGraph
sort -k1,1 -k2,2n ${name}_no_scale_merged.bedGraph | grep -vE "GL|JH" > ${name}_no_scale_merged.sorted.bedGraph

bedGraphToBigWig ${name}_${plus_mask}_${minus_mask}_merged.sorted.bedGraph ${GENOME_SIZE_chrM} ${name}_${plus_mask}_${minus_mask}_merged.bigWig
bedGraphToBigWig ${name}_no_scale_merged.sorted.bedGraph ${GENOME_SIZE_chrM} ${name}_no_scale_merged.bigWig

#=========================================================================================
# Peakcalling using MACS2
#=========================================================================================
for i in *.bedGraph
do
	base=`basename ${i} .bedGraph` 
	source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10
	cu=0.01
	
	#standard NFR calling + merge regions using max-gap 100 to joint linker regions
	macs2 callpeak --broad --treatment ${i} --format BED --gsize ${GENOME_SIZE_NUM} --qvalue ${cu} --broad-cutoff ${cu}  \
	--outdir MACS_${base}_${cu} --name MACS_${base}_${cu} --bdg --nomodel --max-gap 100 --shift -100 --extsize 200
done

bedtools intersect -a MACS_ESC_uncorrected_peaks.broadPeak -b MACS_ESC_corrected_peaks.broadPeak -wa -v > SOB_uncorrected_specific.bed
bedtools intersect -a MACS_ESC_uncorrected_peaks.broadPeak -b MACS_ESC_corrected_peaks.broadPeak -wa  > SOB_uncorrected_shared.bed
bedtools intersect -b MACS_ESC_uncorrected_peaks.broadPeak -a MACS_ESC_corrected_peaks.broadPeak -v -wa > SOB_corrected_specific.bed
bedtools intersect -b MACS_ESC_uncorrected_peaks.broadPeak -a MACS_ESC_corrected_peaks.broadPeak -wa > SOB_corrected_share.bed

#=========================================================================================
# Transcription factors enrichment analysis
#=========================================================================================
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10
bash Negetive_sequence_matched_length_GC.sh SOB_corrected_share.bed ${GENOME_SIZE_chrM} ${FASTA_FILE} 100 
bash Negetive_sequence_matched_length_GC.sh SOB_corrected_specific.bed ${GENOME_SIZE_chrM} ${FASTA_FILE} 100 
bash Negetive_sequence_matched_length_GC.sh SOB_uncorrected_shared.bed ${GENOME_SIZE_chrM} ${FASTA_FILE} 100 
bash Negetive_sequence_matched_length_GC.sh SOB_uncorrected_specific.bed ${GENOME_SIZE_chrM} ${FASTA_FILE} 100 

for i in *bed;do
	bgzip -@ 50 ${i}
done

for i in *.bed_matched.bed
do
	base=`basename ${i} .bed_matched.bed`
	sort -k1,1 -k2,2n ${i} | bgzip -c > ${base}_matched.bed.gz
done

#Build giggle reference on all TFs, filtered using the RIKEN TF list 
#(http://genome.gsc.riken.jp/TFdb/htdocs/tf_list.html)
tail -n +2 ../Riken_filteredTFs.tsv | cut -f 6 | while read ipline
do 
	cp *${ipline}* ../TFs_all
done
giggle index -i "TFs_all/*.gz" -o TFs_all_index -f -s

#Do giggle searching
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10
for i in *bed.gz
do
	base=`basename ${i} .bed.gz`
	giggleindex=/public/home/zhy/Genomes_info/Mus_musculus/Cistrome/ESC_factors/TFs_all_index/
	zcat ${i} | shuf | head -4000 | sort -k1,1 -k2,2n | bgzip -c > ${base}_equal.bed.gz

	nohup giggle search -i ${giggleindex} -q ${base}_equal.bed.gz -s -g ${GENOME_SIZE_NUM} > ${base}_equal_giggle.tsv  &
	nohup giggle search -i ${giggleindex} -q ${i} -s -g ${GENOME_SIZE_NUM} > ${base}_TFs_giggle.tsv  &

	giggleindex=/public/home/zhy/Genomes_info/Mus_musculus/Cistrome/ESC_factors/TFs_unique_index/
	nohup giggle search -i ${giggleindex} -q ${base}_equal.bed.gz -s -g ${GENOME_SIZE_NUM} > ${base}_equal_unique_giggle.tsv  &
	nohup giggle search -i ${giggleindex} -q ${i} -s -g ${GENOME_SIZE_NUM} > ${base}_TFs_unique_giggle.tsv  &
done

#=========================================================================================
# Pre-creat files for seqOutBias
#=========================================================================================
for i in mm10 hg38 dm6 danRer11 tair10 ce11 
do
	#mkdir ${i} && 
	cd ${i}
	source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info ${i}

	#1. Maked gtTxt 
	seqOutBias tallymer ${FASTA_FILE} 36

	faspre=`basename ${FASTA_FILE} .fa`
	#2. Get TBL file
	BAM=../Mouse_ESC_ATACseq.filtered.dedup.bam
	nohup seqOutBias ${FASTA_FILE} ${BAM} --skip-bed --read-size=36 --strand-specific --custom-shift=4,5 \
	--kmer-mask XNXXXCXXNNXNNNXXNNX --bw=${FILE_PREFIX}_corrected.bigWig --tallymer=${faspre}.tal_36.gtTxt.gz &
	cd ../
done

./BiasFreeATAC \
-r ./Mouse_ESC_ATACseq_R1.fastq.gz \
-t ~/Tn5_bias/pre-processing/PeakCalling/ESC/testpipeline/Tallymer \
-m Bowtie2 \
-i ~/Genomes_info/Mus_musculus/bowtie2_index/mm10_index \
-g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-l ~/Genomes_info/Mus_musculus/mm10-blacklist.v2.bed \
-p 30 \
-o ./test5

nohup ./BiasFreeATAC \
-r ~/Tn5_bias/pre-processing/PeakCalling/ESC/testpipeline/Mouse_ZHBTC4_rep1_ATACseq_2017King_R1.fastq.gz \
-t ~/Tn5_bias/pre-processing/PeakCalling/ESC/testpipeline/Tallymer \
-m Bowtie2 \
-i ~/Genomes_info/Mus_musculus/bowtie2_index/mm10_index \
-g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-l ~/Genomes_info/Mus_musculus/mm10-blacklist.v2.bed \
-p 30 \
-o ./ESC &> ESC_BiasFreeATAC.log &

nohup ./BiasFreeATAC -b ESC/Mouse_ZHBTC4_rep1_ATACseq_2017King.filtered.dedup.bam \
-t ~/Tn5_bias/pre-processing/PeakCalling/ESC/testpipeline/Tallymer -m Bowtie2 \
-i ~/Genomes_info/Mus_musculus/bowtie2_index/mm10_index -g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-l ~/Genomes_info/Mus_musculus/mm10-blacklist.v2.bed -p 30 -o ./ESC &>> ESC_BiasFreeATAC.log &