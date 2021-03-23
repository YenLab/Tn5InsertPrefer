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
# 9mers-based Approach on Mouse naked genome
# =======================================================================================
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh && Configuration_info mm10

#Make tiling windows along the Human genome (This take many resources!)
bedtools makewindows -g ${GENOME_SIZE_chrM} -w 9 -s 1 | bedtools getfasta -fi ${FASTA_FILE} -bed - \
| awk -v RS=">" '!/N/{printf $0RT}' | paste - - | sed 's/[:-]/\t/g' | sed 's/>//' > Mouse.9mers.bedfa

backbone_file=Mouse.9mers
target_file=${backbone_file}_collection.bedfa
cat <(echo -e "#Chr\tStart\tEnd\tSeq") ${backbone_file}.bedfa > ${target_file}

#Get DNA me values in each 9mer 
for i in Mouse_E14_rep1_5mC_WGBS_2019Schule.bedGraph Mouse_Germ_5mC_BSseq_2014Pastor.bedGraph
do
	base=`basename ${i} .bedGraph`
	bedtools intersect -a ${backbone_file}.bedfa -b ${i} -wao -sorted -g ${GENOME_SIZE_chrM} \
	| cut -f 1-4,8 | sed 's/\.$/0/g' | bedtools groupby -i - -g 1,2,3 -c 5 -o sum -full | cut -f 1-4,6 > ${backbone_file}_${base}.bedfa
	paste ${target_file} <(cat <(echo ${base}) <(cut -f 5 ${backbone_file}_${base}.bedfa)) > ${target_file}_tmp && mv ${target_file}_tmp ${target_file}
done

#Get Tn5 cut values in each 9mer
for i in *.filtered.dedup.shifted.insertSites.bed
do
	base=`basename ${i} .filtered.dedup.shifted.insertSites.bed`
	bedtools intersect -a ${backbone_file}.bedfa -b ${i} -c -sorted -g ${GENOME_SIZE_chrM} > ${backbone_file}_${base}.bedfa
	paste ${target_file} <(cat <(echo ${base}) <(cut -f 5 ${backbone_file}_${base}.bedfa)) > ${target_file}_tmp && mv ${target_file}_tmp ${target_file}
done

#Chr    Start   End     Seq     Mouse_E14_rep1_5mC_WGBS_2019Schule
#Mouse_Germ_5mC_BSseq_2014Pastor
#
#Mouse_ESC_input_rep1_ATACseq_2017Gary   
#Mouse_ESC_input_rep2_ATACseq_2017Gary  
#Mouse_ESC_rep1_ATACseq_2017Gary 
#Mouse_ESC_rep2_ATACseq_2017Gary 
#Mouse_Germ_ATACseq_2014Pastor   
#Mouse_Germ_input_ATACseq_2014Pastor

# 9mers-based Approach for dissecting DNA methylation effect on Germ and ESC 
file=Mouse.9mers_collection_DNAme.bedfa
awk -v FS="\t" -v OFS="\t" '{
	if($5 != 0 && $6 == 0)  {
		print $0 >> "E14_higher.bedfa"
	} else if($5 == 0 && $6 != 0) {
		print $0 >> "Germ_higher.bedfa"
	} else if($5 == 0 && $6 == 0) {
		print $0 >> "All_zeros.bedfa"
	} else if($5 > 0 && $6 > 0) {
		print $0 >> "All_high.bedfa"
	}
}' ${file}

nohup sort -k4 E14_higher.bedfa  | bedtools groupby -i - -g 4 -c 5,6,7,8 -o mean > E14_higher_unique.bedfa &
nohup sort -k4 Germ_higher.bedfa | bedtools groupby -i - -g 4 -c 5,6,7,8 -o mean > Germ_higher_unique.bedfa &
nohup sort -k4 All_zeros.bedfa   | bedtools groupby -i - -g 4 -c 5,6,7,8 -o mean > All_zeros_unique.bedfa &
nohup sort -k4 All_high.bedfa    | bedtools groupby -i - -g 4 -c 5,6,7,8 -o mean > All_high_unique.bedfa &

# 9mers-based Approach for dissecting quantitative DNA methylation effect on ESC
file=Mouse.9mers_collection_DNAme.bedfa
res=refined10_ESC
cut -f 4,5,7 ${file} | awk -v FS="\t" -v OFS="\t" '{
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

#all=`awk '{s+=$3; ss+=$3^2} END{print m=s/NR, sqrt(ss/NR-m^2)}' ${res}.kmc`
#awk -v FS="\t" -v OFS="\t" '{print $1,$2,($3-0.27155)/0.59517}' ${res}.kmc > ${res}_zscaled.kmc
nohup sort --parallel 30 -k1,1 -k2,2n ${res}.kmc | bedtools groupby -i - -g 1,2 -c 3 -o mean > ${res}_collapsed.txt &

#=========================================================================================
# 9mers-based Approach on ESC chromatin 
#=========================================================================================
bedtools intersect -a Mouse_ESC_rep1_ATACseq_2017Gary_peaks.broadPeak -b Mouse_Germ_ATACseq_2014Pastor_peaks.broadPeak -nonamecheck \
grep -vE "GL|JH" > Mouse_ESC_Germ_peaks.broadPeak
bedtk flt Mouse_ESC_Germ_peaks.broadPeak <(cut -f 1-6,11-13 ../Mouse.9mers_collection.bedfa) > Mouse.9mers_collection_peaks_DNAme.bedfa

backbone_file=Mouse.9mers_collection_peaks_DNAme_peaks.bedfa
res=refined10_ESC_peaks

cut -f 4,5,7 ${backbone_file} | awk -v FS="\t" -v OFS="\t" '{
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

all=`awk '{s+=$3; ss+=$3^2} END{print m=s/NR, sqrt(ss/NR-m^2)}' ${res}.kmc`
awk -v FS="\t" -v OFS="\t" '{print $1,$2,($3-0.27155)/0.59517}' ${res}.kmc > ${res}_zscaled.kmc
nohup sort --parallel 30 -k1,1 -k2,2n ${res}.kmc | bedtools groupby -i - -g 1,2 -c 3 -o mean > ${res}_collapsed.txt &

###210303
bedtk flt Mouse_ESC_rep1_ATACseq_2017Gary_peaks.broadPeak <(cut -f 1-5,10 ../Mouse.9mers_collection.bedfa) > Mouse.9mers_collection_DNAme_peaks.bedfa
nohup bedtools intersect -a Mouse.9mers_collection_DNAme_peaks.bedfa -b Mouse_ESC_MNaseseq_2017Dieuleveult.filtered.dedup.bedGraph -wao \
| cut -f 1-6,10 | sed 's/\.$/0/g' | bedtools groupby -i - -g 1,2,3 -c 7 -o sum -full | cut -f 4,5,6,8 > Mouse.9mers_collection_peaks_DNAme_MNase_peaks.bedfa &

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
}' Mouse.9mers_collection_peaks_DNAme_MNase_peaks.bedfa > Mouse.9mers_collection_peaks_DNAme_MNase_peaks.kmc
nohup sort --parallel 30 -k1,1 -k2,2n Mouse.9mers_collection_peaks_DNAme_MNase_peaks.kmc \
| bedtools groupby -i - -g 1,2 -c 3,4 -o mean > Mouse.9mers_collection_peaks_DNAme_MNase_peaks_collapsed.txt &

nohup awk -v FS="\t" -v OFS="\t" '{if($4 == 0) {print $0}}' Mouse.9mers_collection_peaks_DNAme_MNase_peaks.kmc \
| sort --parallel 30 -k1,1 -k2,2n | bedtools groupby -i - -g 1,2 -c 3,4 -o mean > Mouse.9mers_collection_peaks_DNAme_deMNase_peaks_collapsed.txt &
