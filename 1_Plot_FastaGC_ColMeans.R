#This script read a fasta file and calculate the nucleotide freq of each column

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-f","--FASTA"), type="character", default=NULL, help="Input a fasta file", metavar="character"),
  make_option(c("-s","--SAMPLE"), type="double", default=NULL, help="Sample some fraction of sequence, to save memory")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

FASTA = opt$FASTA

#Get a fraction of original fasta
if (opt$SAMPLE != 1){
	shell_cmd<-paste0("reformat.sh in=",FASTA," out=",FASTA,"_sample.fa fastawrap=2000 ow=t samplerate=",opt$SAMPLE)
	system(shell_cmd)
	FASTA <- paste0(FASTA,"_sample.fa")
}

fa <- readDNAStringSet(FASTA)
cat(paste0(length(fa)," sequences for processing..."))

#Remove fasta with abnormal length
seq_len <- max(width(fa))
index <- which(width(fa) != seq_len)
if(length(index)>0){fa <- fa[-index]}

fa_mt <- as.matrix(fa)
rm(fa)

num_col <- ncol(fa_mt)
num_row <- nrow(fa_mt)
res_mt <- matrix(NA,nrow = 4,ncol = num_col)

for (i in seq(1:num_col)){
  res_mt[1,i] <- sum(fa_mt[,i] %in% "A")/(num_row - sum(fa_mt[,i] %in% "N"))
  res_mt[2,i] <- sum(fa_mt[,i] %in% "C")/(num_row - sum(fa_mt[,i] %in% "N"))
  res_mt[3,i] <- sum(fa_mt[,i] %in% "T")/(num_row - sum(fa_mt[,i] %in% "N"))
  res_mt[4,i] <- sum(fa_mt[,i] %in% "G")/(num_row - sum(fa_mt[,i] %in% "N"))
}

rm(fa_mt)
rownames(res_mt) <- c("A","C","T","G")
colnames(res_mt) <- paste0(seq(-(seq_len-1)/2,(seq_len-1)/2))
write.table(res_mt,paste0(opt$FASTA,"_nuc.txt"), sep="\t",quote = F)

#system(paste0("rm -f ",FASTA,"_sample.fa"))
