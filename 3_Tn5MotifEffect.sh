#!/bin/bash
#This script is used for searching potential Tn5 motifs on all species

#=========================================================================================
# Search genome-wide potential motifs using FIMO, with MEME returned PWM
#=========================================================================================
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh
echo -ne "sample\tTn5_cut_sites\tmotif_explained\ttotal_motifs\tmotif_used\n" >> motif_usage_0.05.txt
echo -ne "sample\tTn5_cut_sites\tmotif_explained\ttotal_motifs\tmotif_used\n" >> motif_usage_0.01.txt
echo -ne "sample\tTn5_cut_sites\tmotif_explained\ttotal_motifs\tmotif_used\n" >> motif_usage_0.005.txt
echo -ne "sample\tTn5_cut_sites\tmotif_explained\ttotal_motifs\tmotif_used\n" >> motif_usage_0.001.txt
echo -ne "sample\tTn5_cut_sites\tmotif_explained\ttotal_motifs\tmotif_used\n" >> motif_usage_0.0001.txt

for motif in *Tn5_motif_full.txt
do
	base_sp=`basename ${motif} _Tn5_motif_full.txt`
  id=`echo ${motif} | sed "s/_.*//"`
  echo "Processing ${id} ..."
  Configuration_info ${id}

  [ ! -d ${sp} ] && mkdir ${sp}
  mv ${motif} ${sp} && cd ${sp}
  ln -sf /public/home/zhy/Tn5_bias/${sp}/Cutsites/arch/all_sites/*t0Tinf_l0r1.bed .

  for thread in 0.05 0.01 0.005 0.001 0.0001
  do
    fimo --verbosity 4 --max-stored-scores 50000000 --thresh ${thread} --oc ${base_sp}_${thread} ${motif} ${FASTA_FILE}

    awk '!a[$2$3$4]++' ${base_sp}_${thread}/fimo.tsv | awk -v FS="\t" -v OFS="\t" '{print $3,$4-1,$5,"Potential_site_"NR-1,$7,".",$8,$9,$10}' | \
    sed '1d' | head -n -4 > ${base_sp}_potential_motifs_${thread}.bed

    tmp=`cat ${base_sp}_potential_motifs_${thread}.bed | wc -l`
    
    #Overlap FIMO motifs and cut sites (Motif usage analysis)
    for bed in $(find -L . -xtype l)
    do
      base=`basename ${bed} .bed`
      Tn5_cut_sites=`cat ${bed} | wc -l`
      bedtools intersect -nonamecheck -a ${base_sp}_potential_motifs_${thread}.bed -b ${bed} -wo > ${base_sp}_potential_motifs_${thread}_${base}
      motif_explained=`awk -v FS="\t" -v OFS="\t" '{print $10,$11,$12}' ${base_sp}_potential_motifs_${thread}_${base} | wc -l`
      motif_used=`awk -v FS="\t" -v OFS="\t" '{print $4}' ${base_sp}_potential_motifs_${thread}_${base} | sort | uniq | wc -l`
  
      echo -ne ${base}"\t"${Tn5_cut_sites}"\t"${motif_explained}"\t"${tmp}"\t"${motif_used}"\n" >> ../motif_usage_${thread}.txt
    done

  done

  cd ../
done

#=========================================================================================
# Dissect FIMO results, to four groups + 3 negative controls
#=========================================================================================
#1. motif been cut 
awk -v FS="\t" -v OFS="\t" '{print $1,$2,$3}' mm_potential_motifs_0.01_mm_input_rep1_t0Tinf_l0r1 | sort | uniq > motif_cut.bed
#2. motif been uncut 
bedtools intersect -a mm_potential_motifs_0.01.bed -b motif_cut.bed -v > motif_uncut.bed

#All Tn5 cut sites
#mm_input_rep1_t0Tinf_l0r1.bed 

#3. fall within motifs
awk -v FS="\t" -v OFS="\t" '{print $10,$11,$12}' mm_potential_motifs_0.01_mm_input_rep1_t0Tinf_l0r1 | sort | uniq > cut_motif.bed 
#4. fall within unmotifs
bedtools intersect -a mm_input_rep1_t0Tinf_l0r1.bed -b cut_motif.bed -v > cut_unmotif.bed

#run_motifs_select.sh
NUM=100000
for i in motif_cut.bed motif_uncut.bed
do
  base=`basename ${i} .bed`
  grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" ${i} | shuf | head -${NUM} \
  | bedtools slop -b 91 -i - -g ${GENOME_SIZE} > ${base}_+-100bp_shuf${NUM}.bed
    bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_+-100bp_shuf${NUM}.bed -fo ${base}_+-100bp_shuf${NUM}.fa
    python ~/Packages/BiasAway/BiasAway.py m -f ${base}_+-100bp_shuf${NUM}.fa | sed "s/ >/>/" > ${base}_+-100bp_shuf${NUM}_shuffled.fa
done

for i in cut_motif.bed cut_unmotif.bed
do
  base=`basename ${i} .bed`
  grep -vE "GL|JH|KI|GL|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg" ${i} | shuf | head -${NUM} \
  | bedtools slop -b 100 -i - -g ${GENOME_SIZE} > ${base}_+-100bp_shuf${NUM}.bed
    bedtools getfasta -fi ${FASTA_FILE} -bed ${base}_+-100bp_shuf${NUM}.bed -fo ${base}_+-100bp_shuf${NUM}.fa
    python ~/Packages/BiasAway/BiasAway.py m -f ${base}_+-100bp_shuf${NUM}.fa | sed "s/ >/>/" > ${base}_+-100bp_shuf${NUM}_shuffled.fa
done

bedtools random -l 201 -n ${NUM} -g ${GENOME_SIZE} > mm10.random_201bp_shuf${NUM}.bed
bedtools getfasta -fi ${FASTA_FILE} -bed mm10.random_201bp_shuf${NUM}.bed -fo mm10.random_201bp_shuf${NUM}.fa
python ~/Packages/BiasAway/BiasAway.py m -f mm10.random_201bp_shuf${NUM}.fa | sed "s/ >/>/" > mm10.random_201bp_shuf${NUM}_shuffled.fa

ln -s ~/Zhang_Scripts/Negetive_sequence_matched_length_GC.sh .
for file in cut_motif_+-100bp_shuf100000.bed cut_unmotif_+-100bp_shuf100000.bed motif_cut_+-100bp_shuf100000.bed motif_uncut_+-100bp_shuf100000.bed 
do
    split -n l/10 ${file} ${file}_split_
    for j in ${file}_split_*
    do
      nohup bash Negetive_sequence_matched_length_GC.sh ${j} ${GENOME_SIZE} ${FASTA_FILE} 500 &
    done
done

for i in *shuf100000_matched.bed
do
        base=`basename ${i} .bed`
    bedtools getfasta -fi ${FASTA_FILE} -bed  ${i} -fo ${base}.fa
done

#Calculate DNA shape
ln -s ~/Tn5_bias/scripts/generate_shapes.R .
ls *fa | while read i; do nohup Rscript generate_shapes.R -f ${i} & ; done