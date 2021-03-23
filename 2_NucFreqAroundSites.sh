#!/bin/bash
#This script is for Tn5 local pattern analysis

##USAGE EXAMPLE
#bash Tn5_pattern_analysis.sh 
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh
BAM_DIR=mapping && cd ${BAM_DIR}
RUN_THREASHOLD=30

#=========================================================================================
#1. Get random sequence for reference, 
# searching motif using Weblogo and MEME, get Nucleotide composition around
#=========================================================================================
[ ! -d random ] && mkdir random
cd random

for i in ce11 mm10 dm6 hg38 danRer11 tair10 pfa2 zm3 
do

  Configuration_info ${i}

  lenRan=500000
  bpRan=51
  prefixRan=${GENOME_NAME}.random_${bpRan}bp_${lenRan}
  
  if [[ ! -f ${prefixRan}.bed ]]
  then

    #Get random 51 sequence for motif control
    bedtools random -l ${bpRan} -n ${lenRan} -g ${GENOME_SIZE_chrM} > ${prefixRan}.bed
    bedtools getfasta -fi ${FASTA_FILE} -bed ${prefixRan}.bed -fo ${prefixRan}.fa
  
    weblogo --composition ${COMPOSITION} -F pdf -S 0.6 --errorbars NO --color '#109648' A 'A' --color '#255C99' C 'C' --color '#D62839' T 'T' \
    --color '#F7B32B' G 'G' < ${prefixRan}.fa > ${prefixRan}_weblogo.pdf
    
    meme -mod zoops -pal -p ${RUN_THREASHOLD} -nmotifs 3 ${prefixRan}.fa -searchsize 1000000 -oc ${prefixRan}_meme -seed 0616 -dna -revcomp 
    
    #Get random 200bp sequence for Nuc background
    bedtools random -l 200 -n ${lenRan} -g ${GENOME_SIZE_chrM} \
    | bedtools getfasta -fi ${FASTA_FILE} -bed - -fo ${GENOME_NAME}.random_200bp_${lenRan}.fa
    Rscript ~/Zhang_Scripts/Zhang/Fasta_GC_content.R -f ${GENOME_NAME}.random_200bp_${lenRan}.fa -s 1
  fi
done

cd ../

#=========================================================================================
#2. Search motifs using all pairsites
#=========================================================================================
[ ! -d all_dups ] && mkdir all_dups
cd all_dups

for i in ../*forfindcutsites.bed
do
  base=`basename ${i} _forfindcutsites.bed`

  #Get species information for current sample
  id=`echo ${base} | sed "s/_.*//"`
  Configuration_info ${id} 
  echo "Processing ${i} for ${GENOME_NAME} ..."

  grep -vE ${CONSENSUS_REMOVE} ${i} > tmp && mv tmp ${i}

  #Find all 9bp duplicated sites
  python ~/Zhang_Scripts/Zhang/Find_9bp_duplicated_region.py -b ${i}
  mv ${i}_9bp_duplicates.bed ${base}_9bp_duplicates.bed

  prefixDupExt=${base}_9bp_duplicates_ext21bp
  bedtools slop -b 21 -i ${base}_9bp_duplicates.bed -g ${GENOME_SIZE_chrM} > ${prefixDupExt}.bed

  #. Using shuffed 500,000 sequences, three times
  for j in 1 2 3
  do
    LenShufDup=500000
    cat ${prefixDupExt}.bed | shuf | head -${LenShufDup} > ${prefixDupExt}_shuf${LenShufDup}_rep${j}.bed

    #Find motif using +-25bp
    bedtools getfasta -fi ${FASTA_FILE} -bed ${prefixDupExt}_shuf${LenShufDup}_rep${j}.bed -fo ${prefixDupExt}_shuf${LenShufDup}_rep${j}.fa
    weblogo --composition ${COMPOSITION} -F pdf -S 0.6 --errorbars NO --color '#109648' A 'A' --color '#255C99' C 'C' --color '#D62839' T 'T' \
    --color '#F7B32B' G 'G' < ${prefixDupExt}_shuf${LenShufDup}_rep${j}.fa > ${prefixDupExt}_shuf${LenShufDup}_rep${j}_weblogo.pdf

    meme -mod zoops -pal -p ${RUN_THREASHOLD} -nmotifs 3 ${prefixDupExt}_shuf${LenShufDup}_rep${j}.fa -searchsize 1000000 -oc ${prefixDupExt}_shuf${LenShufDup}_rep${j}_meme \
    -seed 0616 -dna -revcomp -w 51 

    bedtools slop -b 175 -i ${base}_9bp_duplicates_ext21bp_shuf500000_rep${j}.bed -g ${GENOME_SIZE_chrM} \
    | bedtools getfasta -fi ${FASTA_FILE} -bed - -fo ${base}_9bp_duplicates_ext200bp_shuf500000_rep${j}.fa
    Rscript ~/Zhang_Scripts/Zhang/Fasta_GC_content.R -f ${base}_9bp_duplicates_ext200bp_shuf500000_rep${j}.fa -s 1
  done
done

[ ! -d shuf ] && mkdir shuf
mv *shuf* shuf
cd ../

#=========================================================================================
#3. Get cut sites
#=========================================================================================
[ ! -d all_sites ] && mkdir all_sites
cd all_sites

for i in ../*filtered.dedup.shifted.insertSites.bed
do
  base=`basename ${i} .filtered.dedup.shifted.insertSites.bed`

  id=`echo ${base} | sed "s/_.*//"`
  Configuration_info ${id} 
  echo "Processing ${i} for ${GENOME_NAME} ..."

  grep -vE ${CONSENSUS_REMOVE} ${i} > tmp && mv tmp ${i}

  for t in `seq 0 1`
    do

    #1. make separte threshold files
    T=`expr ${t} + 2` 
    python ~/Zhang_Scripts/Zhang/Tn5_Cuttingsites_frequency.py -i ${i} -l 0 -r 1 -t ${t} -T ${T} -o ${base}_t${t}T${T}_l0r1.bed
    sort -k1,1 -k2,2n ${base}_t${t}T${T}_l0r1.bed > tmp && mv tmp ${base}_t${t}T${T}_l0r1.bed
    
    #2. make accumulate threshold files
    python ~/Zhang_Scripts/Zhang/Tn5_Cuttingsites_frequency.py -i ${i} -l 0 -r 1 -t ${t} -T 1000 -o ${base}_t${t}Tinf_l0r1.bed
    sort -k1,1 -k2,2n ${base}_t${t}Tinf_l0r1.bed > tmp && mv tmp ${base}_t${t}Tinf_l0r1.bed
    done

    awk '{print $4}' ${base}_t0Tinf_l0r1.bed | sort | uniq -c > ${base}_t0_l0r1.cuttimes
done

echo "All jobs done!!!"
