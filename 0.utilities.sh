#!/bin/bash
##########################################################################
#This script restore all species information
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2020 __YenLab@SKLEH__. All rights reserved.
##########################################################################
#Usage: Configuration_info ce11/ce/Worm 
# ce11 mm10 dm6 hg38 danRer11 tair10 pfa2 zm3 

function Configuration_info {
	
	#Patches and scaffold in the genome
	CONSENSUS_REMOVE="^GL|^JH|^KI|chrPt|KN|KZ|Scaffold|chr2110|chrrDNA|ctg"
	motifPath=/public/home/zhy/Tn5_bias/fimo/motifs

	if [ ${1} = "ce11" ] || [ ${1} = "ce" ] || [ ${1} = "Worm" ]
	then
		sp="Worm" #Alias 
		GENOME_NAME="ce11" #Genome assembly 
		sp_dir="/public/home/zhy/Genomes_info/Caenorhabditis_elegans/" #Dir for each species
		FASTA_FILE="${sp_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel_chrM.fa" #genome sequence 
		GTF_FILE="${sp_dir}/Caenorhabditis_elegans.WBcel235.95_chrM.gtf" #genome anotation 
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/ce11_index" #Bowtie2 mapping index
		GENOME_SIZE="${sp_dir}/ce11.chrom.sizes.clean" #All chrs except Patch/scaffold/chrM
		GENOME_SIZE_chrM="${sp_dir}/ce11.chrom.sizes.clean_chrM" #All chrs except Patch/scaffold
		PATCH_PATTERN="" #Patch/scaffold in genome which need filtering
		FEATURE_PATH="${sp_dir}/basic" # basic annotation of this species
		GENOME_SIZE_NUM=100286401 #Effective genome size 
		COMPOSITION=35.44% #GC content
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist="${sp_dir}/ce11-blacklist.v2.bed" #Blacklist if avaliable
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0" # nucleotide frequency
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1" # Di-nucleotide frequency
		Liftover_file="${sp_dir}/ce10ToCe11.over.chain"
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
	
	elif [ ${1} = "mm10" ] || [ ${1} = "mm" ] || [ ${1} = "Mouse" ]
	then
		sp="Mouse"
		GENOME_NAME="mm10"
		sp_dir="/public/home/zhy/Genomes_info/Mus_musculus/"
		FASTA_FILE="${sp_dir}/Mus_musculus.GRCm38.dna.primary_assembly_chrM.fa"
		GTF_FILE="${sp_dir}/Mus_musculus.GRCm38.95_chrM.gtf"
		GTF_FILE_rs="${sp_dir}/mm10.refGene.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/mm10_index"
		STAR_INDEXES="${sp_dir}/STAR_index_150bp/"
		RSEM_INDEX="${sp_dir}/rsem/mm10"
		GENOME_SIZE="${sp_dir}/mm10.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/mm10.chrom.sizes.clean_chrM"
		PATCH_PATTERN="GL|JH"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=2730871774
		COMPOSITION=40.48%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist="${sp_dir}/mm10-blacklist.v2.bed"
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file="${sp_dir}/mm9ToMm10.over.chain.gz"
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt

	elif [ ${1} = "hg38" ] || [ ${1} = "hg" ] || [ ${1} = "Human" ]
	then
		sp="Human"
		GENOME_NAME="hg38"
		sp_dir="/public/home/zhy/Genomes_info/homo_sapiens"
		FASTA_FILE="${sp_dir}/Homo_sapiens.GRCh38.dna.primary_assembly_chrM.fa"
		GTF_FILE="${sp_dir}/Homo_sapiens.GRCh38.95_chrM.gtf"
		GTF_FILE_rs="${sp_dir}/hg38.refGene.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/hg38_index"
		STAR_INDEXES="${sp_dir}/STAR_index_150bp/"
		RSEM_INDEX="${sp_dir}/rsem/hg38"
		GENOME_SIZE="${sp_dir}/hg38.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/hg38.chrom.sizes.clean_chrM"
		PATCH_PATTERN="KI|GL"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=3099750718
		COMPOSITION=38.84%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist="${sp_dir}/hg38-blacklist.v2.bed"
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file="${sp_dir}/hg19ToHg38.over.chain.gz"
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt

	elif [ ${1} = "tair10" ] || [ ${1} = "tair" ] || [ ${1} = "Plant" ]
	then
		sp="Plant"
		GENOME_NAME="tair10"
		sp_dir="/public/home/zhy/Genomes_info/arabidopsis_thaliana"
		FASTA_FILE="${sp_dir}/Arabidopsis_thaliana.TAIR10.dna.toplevel_chrM.fa"
		GTF_FILE="${sp_dir}/Arabidopsis_thaliana.TAIR10.43_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/tair10_index"
		GENOME_SIZE="${sp_dir}/tair10.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/tair10.chrom.sizes.clean_chrM"
		PATCH_PATTERN="chrPt"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=119667750
		COMPOSITION=36.00%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist=""
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file=""
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
				
	elif [ ${1} = "pfa2" ] || [ ${1} = "pfa" ] || [ ${1} = "huffia" ]
	then
		sp="huffia"
		GENOME_NAME="pfa2"
		sp_dir="/public/home/zhy/Genomes_info/plasmodium_falciparum"
		FASTA_FILE="${sp_dir}/Plasmodium_falciparum.EPr1.dna.toplevel_chrM.fa"
		GTF_FILE="${sp_dir}/Plasmodium_falciparum.EPr1.43_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/pfa2_index"
		GENOME_SIZE="${sp_dir}/pfa2.chrom.sizes"
		GENOME_SIZE_chrM="${sp_dir}/pfa2.chrom.sizes"
		PATCH_PATTERN="na"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=23292622
		COMPOSITION=19.34%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist=""
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file=""
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
				
	elif [ ${1} = "danRer11" ] || [ ${1} = "danRer" ] || [ ${1} = "Fish" ]
	then
		sp="Fish"
		GENOME_NAME="danRer11"
		sp_dir="/public/home/zhy/Genomes_info/danio_rerio"
		FASTA_FILE="${sp_dir}/Danio_rerio.GRCz11.dna.primary_assembly_chrM.fa"
		GTF_FILE="${sp_dir}/Danio_rerio.GRCz11.96_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/danRer11_index"
		GENOME_SIZE="${sp_dir}/danRer11.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/danRer11.chrom.sizes.clean_chrM"
		PATCH_PATTERN="KN|KZ"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=1373471384
		COMPOSITION=36.52%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist=""
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file="${sp_dir}/danRer10ToDanRer11.over.chain.gz"
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
				
	elif [ ${1} = "dm6" ] || [ ${1} = "dm" ] || [ ${1} = "Fly" ]
	then
		sp="Fly"
		GENOME_NAME="dm6"
		sp_dir="/public/home/zhy/Genomes_info/Drosophila_melanogaster"
		FASTA_FILE="${sp_dir}/Drosophila_melanogaster.BDGP6.22.dna.toplevel_chrM.fa"
		GTF_FILE="${sp_dir}/Drosophila_melanogaster.BDGP6.22.96_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/dm6_index"
		GENOME_SIZE="${sp_dir}/dm6.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/dm6.chrom.sizes.clean_chrM"
		PATCH_PATTERN="Scaffold|chr2110|chrrDNA"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=143726002
		COMPOSITION=41.67%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist="${sp_dir}/dm6-blacklist.v2.bed"
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file=""
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
				
	elif [ ${1} = "zm3" ] || [ ${1} = "zm" ] || [ ${1} = "maize" ]
	then
		sp="maize"
		GENOME_NAME="zm3"
		sp_dir="/public/home/zhy/Genomes_info/zea_mays"
		FASTA_FILE="${sp_dir}/Zea_mays.B73_RefGen_v4.dna.toplevel_chrM.fa"
		GTF_FILE="${sp_dir}/Zea_mays.B73_RefGen_v4.43_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/zm3_index"
		GENOME_SIZE="${sp_dir}/zm3.chrom.sizes.clean"
		GENOME_SIZE_chrM="${sp_dir}/zm3.chrom.sizes.clean_chrM"
		PATCH_PATTERN="ctg|chrPt"
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=2135083061
		COMPOSITION=46.19%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist=""
		markov_bg0="${sp_dir}/${GENOME_NAME}_markov_bg_0"
		markov_bg1="${sp_dir}/${GENOME_NAME}_markov_bg_1"
		Liftover_file=""
		Tn5MotifProb=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME_prob.txt
		Tn5MotifMEME=${motifPath}/${GENOME_NAME}_Tn5_motif_MEME.txt
		Tn5MotifMEME19=${motifPath}/${GENOME_NAME}_Tn5_motif19bp_MEME.txt
			
	elif [ ${1} = "sacCer3" ] || [ ${1} = "sacCer" ] || [ ${1} = "yeast" ]
	then
		sp="yeast"
		GENOME_NAME="sacCer3"
		sp_dir="/public/home/zhy/Genomes_info/saccharomyces_cerevisiae"
		FASTA_FILE="${sp_dir}/Saccharomyces_cerevisiae.EF4.67.dna.toplevel_chrM.fa"
		GTF_FILE="${sp_dir}/Saccharomyces_cerevisiae.EF4.67_chrM.gtf"
		BOWTIE2_INDEXES="${sp_dir}/bowtie2_index/sacCer3_index"
		GENOME_SIZE="${sp_dir}/sacCer3.chrom.sizes"
		GENOME_SIZE_chrM="${sp_dir}/sacCer3.chrom.sizes"
		PATCH_PATTERN=""
		FEATURE_PATH="${sp_dir}/basic"
		GENOME_SIZE_NUM=12157105
		COMPOSITION=38.15%
		PREFIX="/public/home/zhy/Tn5_bias/${sp}/"
		blacklist=""
		markov_bg0=""
		Liftover_file=""
	
	else
	  echo "Your genome is not found!"
	  exit 1
	fi
}


