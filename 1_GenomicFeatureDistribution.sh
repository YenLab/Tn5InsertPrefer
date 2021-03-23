#!/bin/bash
# This script is for calculating Tn5 distribution across genomic features
source /public/home/zhy/Tn5_bias/scripts/0.utilities.sh

for i in *filtered.dedup.shifted.insertSites.bed
do
  base=`basename ${i} .bed`

  #Get species information for current sample
  id=`echo ${base} | sed "s/_.*//"`
  Configuration_info ${id} 
  echo "Processing ${i} for ${GENOME_NAME} ..."

  #Get blacklist removed cut sites
  if [[ ! -f ${base}_rb.bed ]]; then
    sort -k1,1 -k2,2n ${i} | bedtools intersect -nonamecheck -a - -b ${blacklist} -v -nobuf > ${base}_rb.bed
  fi

  #All cut sites for calculation
  python ~/Zhang_Scripts/Zhang/Genomic_feature_occupancy_significance_v2.py -g ${GENOME_SIZE_NUM} -f ${FEATURE_PATH} \
  -s ${base}_rb.bed -o ${base}.feature.out -r 0 

  #Remove outlier cut sites for calculation
  python ~/Zhang_Scripts/Zhang/Genomic_feature_occupancy_significance_v2.py -g ${GENOME_SIZE_NUM} -f ${FEATURE_PATH} \
  -s ${base}_rb.bed -o ${base}_r0.05.feature.out -r 0.05

  #Get unique cut sites
  if [[ ! -f ${base}_rb_uq.bed ]]; then
    sort -k1,1 -k2,2n ${base}_rb.bed | awk '!a[$1$2$3]++' > ${base}_rb_uq.bed
  fi

  #Unique cut sites for calculation
  python ~/Zhang_Scripts/Zhang/Genomic_feature_occupancy_significance_v2.py -g ${GENOME_SIZE_NUM} -f ${FEATURE_PATH} \
  -s ${base}_rb_uq.bed -o ${base}_uq.feature.out -r 0

  #Remove outlier unique cut sites for calculation
  nohup python ~/Zhang_Scripts/Zhang/Genomic_feature_occupancy_significance_v2.py -g ${GENOME_SIZE_NUM} -f ${FEATURE_PATH} \
  -s ${base}_rb_uq.bed -o ${base}_uq_r0.05.feature.out -r 0.05 &
done
