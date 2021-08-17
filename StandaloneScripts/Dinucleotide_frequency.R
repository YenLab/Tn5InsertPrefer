#This script is for calculating all Dinucleotide frequency from fasta file 

setwd("./")
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-f", "--FASTA"), type="character", default=NULL,metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#===========================================================
#1. Get di-nucleotides freqs
#===========================================================
fa <- readDNAStringSet(opt$FASTA)
cat(paste0(length(fa)," sequences for processing..."))

#filter out seqs with abnormal length
seq_len <- max(width(fa))
index <- which(width(fa) != seq_len)
if(length(index)>0){fa <- fa[-index]}

data <- as.matrix(fa)

patts <- c("AA","AT","AC","AG","TA","TT","TC","TG","CA","CT","CC","CG","GA","GT","GC","GG")

result_fre_matrix <- matrix(0, nrow = length(patts), ncol = ncol(data)) %>% as.data.frame()
rownames(result_fre_matrix) <- patts

for (step in 1:(ncol(data)-1)){

	if(step %% 10 == 0){
		cat("Processing col ",step,"\n")
	}

  tmp <- table(paste0(data[,step],data[,step+1]))
  for (patt in patts){
  	result_fre_matrix[patt,step] <- tmp[patt][[1]]
  }
}
result_fre_matrix <- round(result_fre_matrix/nrow(data),3)

write.table(result_fre_matrix,paste0("Dinucleotide_frequency_of_",opt$FASTA,".txt"), col.names = F, quote = F)

