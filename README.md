# Comprehensive understanding of Tn5 insertion preference recovers expansive transcription regulatory lexicon

**Author**: Houyu Zhang

**Email**: Hughiez047@gmail.com

## Introduction

Tn5 transposase has been widely adopted as a molecular tool in next-generation sequencing. With more types of Tn5 segmented data subjected to precision genomics, the impact on Tn5 insertion bias comes into focus. Besides, whether explicitly correct the Tn5 bias benefits its general applications is of interest for genomic studies. 

Here we systematically explored Tn5 insertion characteristics across several model organisms to find critical parameters affect its insertion preference. An initial evaluation of Tn5 insertion distribution along naked genomic DNA revealed that Tn5 insertion is not uniformly random. Leveraging a machine learning framework, we found that DNA shape independently or cooperatively works with DNA motif to affect Tn5 insertion. These intrinsic preferences can be faithfully modeled for computational correction using nucleotide dependency information from DNA sequences. 

To promote bias-corrected ATAC-seq analysis, we hence developed a pipeline for this purpose. Using our pipeline, we showed the bias correction would improve the overall performance of ATAC-seq peak detection. Importantly, we showed that with bias correction, many potential false-negative peaks could be recovered. By conducting a genome-wide transcription factors (TFs) enrichment analysis, we found these peaks contain abundant TFs. We propose that a comprehensive understanding and precise correction of Tn5 insertion preference could broaden the view of Tn5-assisted sequencing data. 

## Scripts organization

- Scripts head with `0_` is for general preprocessing NGS data, include

   `ATAC-seq`,  `WGBS`,  `MNase-seq` ,`RNA-seq` ,`ChIP-seq` 

- Script head with series number is pipeline for specific analysis follows the order in the Tn5 paper.

  Specifically, the `bash script` and `python script` is for processing data and corresponding `R scripts` is used for downstream analysis, statistics and plotting.
  
- Scripts in the StandaloneScript folder was called by pipelines.

## Tn5-bias corrected ATAC-seq pipeline

This pipeline is composed of a bash script and configure `yaml` file. 

### Prerequisites

We recommend using a conda environment to run this pipeline and you can install all dependent packages as following:

```bash
conda env create -f BiasFreeATAC_env.yaml
source activate BiasFreeATAC_env
chmod +x BiasFreeATAC
```

You can find the install instructions of bias correction package seqOutBias (Martins *et al., NAR*, 2017) from its [github repository](https://github.com/guertinlab/seqOutBias) for specific version; or use followed command. Please note seqOutBias is on Rust, the compiler and other dependencies were also installed above using conda, so, you need only to compile seqOutBias:

```bash
wget -c https://github.com/guertinlab/seqOutBias/archive/refs/tags/v1.3.0.tar.gz -O seqOutBias.tar.gz
tar xzf seqOutBias.tar.gz
cd seqOutBias
cargo build --release
-> /target/release
```

### Help information:

```bash
Usage: BiasFreeATAC.sh [options]

-r <string>             If start from raw fastq file, please provide Read1 (_R1.fastq.gz) with this parameter.
-b <string>             If start from mapped bam file, please provide with this parameter and you can miss the -m/-i parameter
-m <string>             Indicate the mapper you wish to use (Bowtie2/BWA)
-i <string>             Location of mapping index
-g <string>             Genome size for each chromosome
-f <string>             Genome reference fasta file
-t <string>             tallymer mappability file [optional]
-l <string>             Blacklist regions [optional]
-p <string>             Thresholds for parallelly run this pipeline
-o <string>             Working directory, all output files will be generated here
```

`-r` The Read1 of raw fastq file. In order to unify the pipeline, we suggest you name the file following this format: `Sample_R1.fastq.gz` and `Sample_R2.fastq.gz`.

`-b` You can also provide a bam file to skip all preprocessing steps.

``

`-i` The location of bowtie index for mapping

`-g` Chromosome size. You can obtain it from:

- UCSC using `fetchChromSizes hg38 > hg38.chrom.sizes`
- index of genome fasta file `samtools faidx input.fa` `cut -f1,2 input.fa.fai > chrom.sizes`  

`-f` Genome sequence file in fasta format.

`-l` Blacklist that removed from this analysis. You can found the list in [here](https://github.com/Boyle-Lab/Blacklist).  

`-p` Thresholds for parallel running the pipeline.

`-o` Working directory. All files will be generated under this directory. Will build one if not exist.

### Example on mESC:

```bash
bash BiasFreeATAC.sh \
-r ~/Tn5_bias/pre-processing/PeakCalling/ESC/testpipeline/Mouse_ESC_ATACseq_R1.fastq.gz \
-m Bowtie2 \
-i ~/Genomes_info/Mus_musculus/bowtie2_index/mm10_index \
-g ~/Genomes_info/Mus_musculus/mm10.chrom.sizes \
-f ~/Genomes_info/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa \
-p 30 \
-o test
```

### Output files