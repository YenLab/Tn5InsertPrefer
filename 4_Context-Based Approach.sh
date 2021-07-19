#!/bin/bash
##########################################################################
#This script is used for dissecting DNA methylation effect on
#Tn5 insertion. Since the DNA sequence (motif and shape) affect Tn5 insertion
#a lot, to explore the true effect of DNA methylation, 
#we followed the kmer-based approach (Lazarovici et al., 2013). Because the Tn5 
#core motif is 9bp, 9mers-based approach should largely reflect the Tn5 preference.

#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################

# =======================================================================================
# Make whole matrix of 9mer-based information
# =======================================================================================
#Make tiling windows along the Mouse genome (This take many resources!)
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10
bedtools makewindows -g ${GENOME_SIZE_chrM} -w 9 -s 1 | bedtools getfasta -fi ${FASTA_FILE} -bed - \
| awk -v RS=">" '!/N/{printf $0RT}' | paste - - | sed 's/[:-]/\t/g' | sed 's/>//' > Mouse.9mers.bedfa

target_file=Mouse.9mers_collection.bedfa
cat <(echo -e "#Chr\tStart\tEnd\tSeq") Mouse.9mers.bedfa > ${target_file}

#Get DNA methylation values in each 9mer 
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10
for i in Mouse_E14_rep2_5mC_WGBS_2019Schule.bedGraph Mouse_Germ_5mC_BSseq_2014Pastor.bedGraph
do
	base=`basename ${i} .bedGraph`
	bedtools intersect -a Mouse.9mers.bedfa -b ${i} -wao -sorted -g ${GENOME_SIZE_chrM} \
	| cut -f 1-4,8 | sed 's/\.$/0/g' | bedtools groupby -i - -g 1,2,3 -c 5 -o sum -full | cut -f 1-4,6 > Mouse.9mers_${base}.bedfa
	paste ${target_file} <(cat <(echo ${base}) <(cut -f 5 Mouse.9mers_${base}.bedfa)) > ${target_file}_tmp && mv ${target_file}_tmp ${target_file}
done

#Get Tn5 cut values in each 9mer
for i in *.filtered.dedup.shifted.insertSites.bed
do
	base=`basename ${i} .filtered.dedup.shifted.insertSites.bed`
	bedtools intersect -a Mouse.9mers.bedfa -b ${i} -c -sorted -g ${GENOME_SIZE_chrM} > Mouse.9mers_${base}.bedfa
	paste ${target_file} <(cat <(echo ${base}) <(cut -f 5 Mouse.9mers_${base}.bedfa)) > ${target_file}_tmp && mv ${target_file}_tmp ${target_file}
done

# =======================================================================================
# Context-dependent approach on Mouse ESC and Germ cell
# =======================================================================================
# naked DNA
cut -f 4,5,6,13,14 Mouse.9mers_collection.bedfa | grep -v "Seq" | awk -v FS="\t" -v OFS="\t" '{
	if($2 != 0 && $3 == 0)  {
		print $0 >> "E14_higher.bedfa"
	} else if($2 == 0 && $3 != 0) {
		print $0 >> "Germ_higher.bedfa"
	} else if($2 == 0 && $3 == 0) {
		print $0 >> "All_zeros.bedfa"
	} else if($2 > 0 && $3 > 0) {
		print $0 >> "All_high.bedfa"
	}
}'
nohup sort -k1 E14_higher.bedfa  | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > E14_higher_unique.bedfa &
nohup sort -k1 Germ_higher.bedfa | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > Germ_higher_unique.bedfa &
nohup sort -k1 All_zeros.bedfa   | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > All_zeros_unique.bedfa &
nohup sort -k1 All_high.bedfa    | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > All_high_unique.bedfa &

# Chromatin
cut -f 4,5,6,12,15 Mouse.9mers_collection.bedfa | grep -v "Seq" | awk -v FS="\t" -v OFS="\t" '{
	if($2 != 0 && $3 == 0)  {
		print $0 >> "E14_higher_peaks.bedfa"
	} else if($2 == 0 && $3 != 0) {
		print $0 >> "Germ_higher_peaks.bedfa"
	} else if($2 == 0 && $3 == 0) {
		print $0 >> "All_zeros_peaks.bedfa"
	} else if($2 > 0 && $3 > 0) {
		print $0 >> "All_high_peaks.bedfa"
	}
}'
nohup sort -k1 E14_higher_peaks.bedfa  | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > E14_higher_peaks_unique.bedfa &
nohup sort -k1 Germ_higher_peaks.bedfa | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > Germ_higher_peaks_unique.bedfa &
nohup sort -k1 All_zeros_peaks.bedfa   | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > All_zeros_peaks_unique.bedfa &
nohup sort -k1 All_high_peaks.bedfa    | bedtools groupby -i - -g 1 -c 2,3,4,5 -o mean > All_high_peaks_unique.bedfa &

# =======================================================================================
# Deciles-based analysis 
# =======================================================================================
res=refined10_E14
cut -f 4,5,14 Mouse.9mers_collection.bedfa | awk -v FS="\t" -v OFS="\t" '{
	if($2 == 0)  {
	print $1,"0",$3  } else if($2 > 0 && $2 <= 10)   {
	print $1,"1",$3  } else if($2 > 10 && $2 <= 20)  {
	print $1,"2",$3  } else if($2 > 20 && $2 <= 30)  {
	print $1,"3",$3  } else if($2 > 30 && $2 <= 40)  {
	print $1,"4",$3  } else if($2 > 40 && $2 <= 50)  {
	print $1,"5",$3  } else if($2 > 50 && $2 <= 60)  {
	print $1,"6",$3  } else if($2 > 60 && $2 <= 70)  {
	print $1,"7",$3  } else if($2 > 70 && $2 <= 80)  {
	print $1,"8",$3  } else if($2 > 80 && $2 <= 90)  {
	print $1,"9",$3  } else if($2 > 90 && $2 <= 100) {
	print $1,"10",$3 } else if($2 > 100) {
	print $1,"11",$3 } 
}' | grep -v "Seq" > ${res}.kmc
nohup sort --parallel 30 -k1,1 -k2,2n ${res}.kmc | bedtools groupby -i - -g 1,2 -c 3 -o mean > ${res}_collapsed.txt &

res=refined10_E14chromatin
cut -f 4,5,10 Mouse.9mers_collection.bedfa | awk -v FS="\t" -v OFS="\t" '{
	if($2 == 0)  {
	print $1,"0",$3  } else if($2 > 0 && $2 <= 10)   {
	print $1,"1",$3  } else if($2 > 10 && $2 <= 20)  {
	print $1,"2",$3  } else if($2 > 20 && $2 <= 30)  {
	print $1,"3",$3  } else if($2 > 30 && $2 <= 40)  {
	print $1,"4",$3  } else if($2 > 40 && $2 <= 50)  {
	print $1,"5",$3  } else if($2 > 50 && $2 <= 60)  {
	print $1,"6",$3  } else if($2 > 60 && $2 <= 70)  {
	print $1,"7",$3  } else if($2 > 70 && $2 <= 80)  {
	print $1,"8",$3  } else if($2 > 80 && $2 <= 90)  {
	print $1,"9",$3  } else if($2 > 90 && $2 <= 100) {
	print $1,"10",$3 } else if($2 > 100) {
	print $1,"11",$3 } 
}' | grep -v "Seq" > ${res}.kmc
nohup sort --parallel 30 -k1,1 -k2,2n ${res}.kmc | bedtools groupby -i - -g 1,2 -c 3 -o mean > ${res}_collapsed.txt &

#=========================================================================================
# 9mers-based Approach on open chromatin regions (ESC Naked DNA)
#=========================================================================================
peakfile="Mouse_ESC_rep1_ATACseq_2017Gary_peaks.broadPeak"
#MNase="Mouse_ESC_MNaseseq_2017Dieuleveult.filtered.dedup.bedGraph"
bonefile2="Mouse.9mers_collection_peaks_Naked_DNAme.bedfa"
bonefile2_kmc="Mouse.9mers_collection_peaks_Naked_DNAme.kmc"
bonefile2_collapsed="Mouse.9mers_collection_peaks_Naked_DNAme_collapsed.txt"

bedtk flt ${peakfile} <(cut -f 1-5,14 ../Mouse.9mers_collection.bedfa) | cut -f 4,5,6 > ${bonefile2}
awk -v FS="\t" -v OFS="\t" '{
	if($2 == 0)  {
	print $1,"0",$3  } else if($2 > 0 && $2 <= 10)   {
	print $1,"1",$3  } else if($2 > 10 && $2 <= 20)  {
	print $1,"2",$3  } else if($2 > 20 && $2 <= 30)  {
	print $1,"3",$3  } else if($2 > 30 && $2 <= 40)  {
	print $1,"4",$3  } else if($2 > 40 && $2 <= 50)  {
	print $1,"5",$3  } else if($2 > 50 && $2 <= 60)  {
	print $1,"6",$3  } else if($2 > 60 && $2 <= 70)  {
	print $1,"7",$3  } else if($2 > 70 && $2 <= 80)  {
	print $1,"8",$3  } else if($2 > 80 && $2 <= 90)  {
	print $1,"9",$3  } else if($2 > 90 && $2 <= 100) {
	print $1,"10",$3 } else if($2 > 100) {
	print $1,"11",$3 } 
}' ${bonefile2} > ${bonefile2_kmc}

nohup sort --parallel 30 -k1,1 -k2,2n ${bonefile2_kmc} | bedtools groupby -i - -g 1,2 -c 3 -o mean > ${bonefile2_collapsed} &

peakfile="Mouse_ESC_rep1_ATACseq_2017Gary_peaks.broadPeak"
MNase="Mouse_ESC_MNaseseq_2017Dieuleveult.filtered.dedup.bedGraph"
bonefile1="Mouse.9mers_collection_peaks_chromatin_DNAme_MNase.bedfa"
bonefile1_kmc="Mouse.9mers_collection_peaks_chromatin_DNAme_MNase.kmc"
bonefile1_collapsed="Mouse.9mers_collection_peaks_chromatin_DNAme_MNase_collapsed.txt"
bonefile1_deMNase_collapsed="Mouse.9mers_collection_peaks_chromatin_DNAme_deMNase_collapsed.txt"

#15,16,17
bedtk flt ${peakfile} <(cut -f 1-5,15 ../Mouse.9mers_collection.bedfa) | \
bedtools intersect -a - -b ${MNase} -wao | cut -f 1-6,10 | sed 's/\.$/0/g' | \
bedtools groupby -i - -g 1,2,3 -c 7 -o sum -full | cut -f 4,5,6,8 > ${bonefile1}

#Decile the DNA methylation level
awk -v FS="\t" -v OFS="\t" '{
	if($2 == 0)  {
	print $1,"0",$3,$4  } else if($2 > 0 && $2 <= 10)   {
	print $1,"1",$3,$4  } else if($2 > 10 && $2 <= 20)  {
	print $1,"2",$3,$4  } else if($2 > 20 && $2 <= 30)  {
	print $1,"3",$3,$4  } else if($2 > 30 && $2 <= 40)  {
	print $1,"4",$3,$4  } else if($2 > 40 && $2 <= 50)  {
	print $1,"5",$3,$4  } else if($2 > 50 && $2 <= 60)  {
	print $1,"6",$3,$4  } else if($2 > 60 && $2 <= 70)  {
	print $1,"7",$3,$4  } else if($2 > 70 && $2 <= 80)  {
	print $1,"8",$3,$4  } else if($2 > 80 && $2 <= 90)  {
	print $1,"9",$3,$4  } else if($2 > 90 && $2 <= 100) {
	print $1,"10",$3,$4 } else if($2 > 100) {
	print $1,"11",$3,$4 } 
}' ${bonefile1} > ${bonefile1_kmc}

#Average all 9mers
nohup sort --parallel 30 -k1,1 -k2,2n ${bonefile1_kmc} | bedtools groupby -i - -g 1,2 -c 3,4 -o mean > ${bonefile1_collapsed} &
