# Genome-wide determinants of Tn5 insertion preference

**Author**: Houyu Zhang

**Email**: Hughiez047@gmail.com

## Introduction

Tn5 transposase has been widely adopted as a molecular tool in next-generation sequencing (NGS) for tagmentation. Although some work has been done to understand Tn5 insertion bias, a comprehensive investigation of Tn5 insertion preference genome-wide is lacking. 

Here we systematically explored Tn5 insertion characteristics across 8 model organisms, in both the naked genomic DNA context to determine its intrinsic preference as well as in a chromatin context, which has relevance for its application in NGS approaches. Evaluation of Tn5 insertion distribution along naked genomic DNA revealed that Tn5 prefers specific genomic regions and the insertion patterns around transcription start sites (TSS) are cell-type specific. Leveraging a machine learning framework on abundant Tn5 insertion sites, we found DNA motif and DNA shape independently or cooperatively contribute to Tn5 preference. The DNA motif and DNA shape effects were also reflected in flexible chromatin contexts. Using a context-dependent approach, we unexpectedly found that DNA methylation, a transcription-associated feature, promotes Tn5 insertion in a highly quantitative manner in naked genomic DNA. 

Our results argue against the random nature proposed for Tn5 insertion and suggest analysis of data generated via Tn5 tagmentation should take this into account.  



## Scripts organization

- Scripts head with `0_` is for general preprocessing NGS data, include

   `ATAC-seq`,  `WGBS`,  `MNase-seq` 

- Script head with series number is for specific analysis follows the order in the Tn5 paper.

  Specifically, the `bash script` and `python script` is for processing data and corresponding `R scripts` is used for downstream analysis, statistics and plotting.
  
- Scripts without number prefix were used as standalones

   
