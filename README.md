# Comprehensive understanding of Tn5 insertion preference recovers expansive transcription regulatory lexicon

**Author**: Houyu Zhang

**Email**: Hughiez047@gmail.com

#Copyright (c) 2021 YenLab@SKLEH. All rights reserved.

## Introduction

Tn5 transposase has been widely adopted as a molecular tool in next-generation sequencing. With more types of Tn5 segmented data subjected to precision genomics, the impact on Tn5 insertion bias comes into focus. Besides, whether explicitly correct the Tn5 bias benefits its general applications is of interest for genomic studies. 

Here we systematically explored Tn5 insertion characteristics across several model organisms to find critical parameters affect its insertion preference. An initial evaluation of Tn5 insertion distribution along naked genomic DNA revealed that Tn5 insertion is not uniformly random. Leveraging a machine learning framework, we found that DNA shape independently or cooperatively works with DNA motif to affect Tn5 insertion. These intrinsic preferences can be faithfully modeled for computational correction using nucleotide dependency information from DNA sequences. 

To promote bias-corrected ATAC-seq analysis, we hence developed a pipeline for this purpose. Using our pipeline, we showed the bias correction would improve the overall performance of ATAC-seq peak detection. Importantly, we showed that with bias correction, many potential false-negative peaks could be recovered. By conducting a genome-wide transcription factors (TFs) enrichment analysis, we found these peaks contain abundant TFs. We propose that a comprehensive understanding and precise correction of Tn5 insertion preference could broaden the view of Tn5-assisted sequencing data. 

<img align="left" width=450 src="https://github.com/YenLab/Tn5InsertPrefer/blob/main/StandaloneScripts/GraphicalAbstract.png">  

## Scripts organization

- Scripts head with `0_` is for general preprocessing NGS data, include

   `ATAC-seq`,  `WGBS`,  `MNase-seq` ,`RNA-seq` ,`ChIP-seq` 

- Script head with series number is pipeline for specific analysis follows the order in the Tn5 paper.

  Specifically, the `bash script` and `python script` is for processing data and corresponding `R scripts` is used for downstream analysis, statistics and plotting.
  
- Scripts in the StandaloneScript folder was called by pipelines.


## BiasFreeATAC: A Tn5-bias corrected ATAC-seq pipeline

BiasFreeATAC pipeline is composed of a bash script and configure `yaml` file, please install following these guidelines. 

### Installation

We recommend using a conda environment to run this pipeline and you can install all dependent packages as following:

```bash
#Download scripts
wget https://github.com/YenLab/Tn5InsertPrefer/blob/main/BiasFreeATAC_env.yaml
wget https://github.com/YenLab/Tn5InsertPrefer/blob/main/BiasFreeATAC
#Install dependencies
conda env create -f BiasFreeATAC_env.yaml
source activate BiasFreeATAC_env
#Install BiasFreeATAC pipeline
chmod +x BiasFreeATAC
```

You can find the install instructions of bias correction package seqOutBias (Martins *et al., NAR*, 2017) from 1) our pre-compiled [version 1.3.0](https://github.com/YenLab/Tn5InsertPrefer/blob/main/seqOutBias); 2). its [github repository](https://github.com/guertinlab/seqOutBias) for specific version (>=1.3.0); 3) compile from source using command below. Please note seqOutBias is based on Rust, the compiler and other dependencies were installed above using conda, so, you need only to compile seqOutBias:

```bash
#Install seqOutBias
wget -c https://github.com/guertinlab/seqOutBias/archive/refs/tags/v1.3.0.tar.gz -O seqOutBias.tar.gz
tar xzf seqOutBias.tar.gz
cd seqOutBias
cargo build --release
#Export seqOutBias executable file into local PATH or you can add it into ~/.bashrc
export PATH=$(realpath ./target/release):$PATH
```

### Usage

BiasFreeATAC accept paired end (PE) `raw fastq` file for *de novo* preprocessing; or a `bam` file that subjected to the bias correction procedure. Provide input files with relative/absolute path as you like, BiasFreeATAC works.

BiasFreeATAC take ~3 hour to run on 65M PE reads for mouse ESC, using 30 CPUs.

```bash
Usage: ./BiasFreeATAC [options]

-r <string>             If start from raw fastq file, please provide Read1 (_R1.fastq.gz) with this parameter.
-b <string>             If start from mapped bam file, please provide with this parameter and you can miss the -r/-m/-i parameter.
-m <string>             Indicate which mapper you wish to use (Bowtie2/BWA)
-i <string>             Location of mapping index corresponding to the mapper you provided (-m).
-g <string>             Genome size for each chromosome.
-f <string>             Genome reference fasta file.
-t <string>             tallymer mappability file directory [optional].
-l <string>             Blacklist regions [optional].
-p <string>             Thresholds for parallelly run this pipeline.
-o <string>             Working directory, all output files will be generated here.
```

`-r` The **Read1** of raw fastq file. In order to unify the pipeline, we suggest you name the file following this format: `Sample_R1.fastq.gz` and `Sample_R2.fastq.gz`.

`-b` You can also provide a bam file to skip all preprocessing steps.

`-m` BiasFreeATAC currently support two mappers (Bowtie2/BWA).

`-i` The location of mapping index for mapper.

`-g` Chromosome size. You can obtain it from:

- UCSC using `fetchChromSizes hg38 > hg38.chrom.sizes`
- index of genome fasta file `samtools faidx input.fa` `cut -f1,2 input.fa.fai > chrom.sizes`  

`-f` Genome sequence file in fasta format.

`-t` [Optional] This directory contain two files for Tn5 bias correction:

- ${GenomeFastaPrefix}.tal_36.gtTxt.gz
- ${GenomeFastaPrefix}_36.18.5.5.tbl

You can find our pre-calculated files for `mm10` `hg38` `dm6` `ce11` `danRer11` `tair10` in [ZENODO](https://zenodo.org/record/5115506#.YRIwpGgzaUk). Note, you need to rename the ${GenomeFastaPrefix} same as your genome fasta file prefix as that in `-f`. 
For `mm10` and `hg38`, please use `cat *_Part0* > ${GenomeFastaPrefix}_36.18.5.5.tbl` to merge our splitted .tbl files.

Otherwise, you can skip this parameter and let BiasFreeATAC create this for you. Hope you will be patient, because for mouse (mm10), this step takes ~5h.

`-l` [Optional] Blacklist that removed from this analysis. You can find this for `mm10` `hg38` `dm6` `ce11`  [here](https://github.com/Boyle-Lab/Blacklist).  

`-p` Thresholds for parallel running the pipeline.

`-o` Working directory. All files will be generated under this directory. Will build one if not exist.

### Example

```bash
#From fastq file
./BiasFreeATAC \
-r ./Mouse_ESC_ATACseq_R1.fastq.gz \
-t ./Tallymer \
-m Bowtie2 \
-i ~/Genomes_info/Mus_musculus/bowtie2_index/mm10_index \
-g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-l ~/Genomes_info/Mus_musculus/mm10-blacklist.v2.bed \
-p 30 -o ./ESC &> ESC_BiasFreeATAC.log

#From bam file
./BiasFreeATAC \
-r ./Mouse_ESC_ATACseq.filtered.dedup.bam \
-t ./Tallymer \
-g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-p 30 -o ./ESC &> ESC_BiasFreeATAC.log
```

### Output files

- `./metrics` : metrics files for reads QC, mapping rate, fragments distribution, deduplication rate.
- `.filtered.dedup.bam` `.filtered.dedup.bam.bai` : The mapped, filtered and deduplicated bam file.
- `.*InsertSites.*` : Single-base Tn5 insertion sites in bed/bam format.
- `_uncorrected.bigWig` `uncorrected.bedGraph` `_corrected.bigWig` `_corrected.bedGraph` : Uncorrected and corrected single-base Tn5 insertion signals by seqOutBias.
- `./_uncorrected_peaks` `./_uncorrected_peaks`: .broadPeak for Uncorrected and corrected Tn5 insertion signals by MACS2.
- `_uncorrected_share_peaks.bed` `_uncorrected_specific_peaks.bed` `_corrected_share_peaks.bed` `_corrected_specific_peaks.bed` : Overlaps of uncorrected/corrected peaks for downstream analysis, for example, the TF-enrichment analysis for each peak set.

### Citation