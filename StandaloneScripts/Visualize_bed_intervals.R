#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-b","--bedfile"), type="character", default=NULL, help="Input bed file", metavar="character"),
  make_option(c("-c","--chr"), type="character", default=NULL, help="Choose which chr to plot (chr, chr.*, Separate )")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# ==============================================================================
# Start Plot for bed file
# ==============================================================================
options(readr.num_columns = 0)
data <- read_delim(opt$bedfile, delim = "\t",col_names = F )
data$length <- data$X3 - data$X2

if (opt$chr == "Separate"){
  data_selected <- data
  pdf(paste0(opt$bedfile,"_",opt$chr,".pdf"), width = 12, height = 12)
  p <- ggplot(data_selected) +
    geom_histogram(aes(x=length), binwidth = 10) +
    geom_vline(xintercept = median(data_selected$length), color="red", linetype = "dashed", size = 2) + 
    labs(x="Bed interval length (bp)", y = "Counts", title = paste0(opt$chr," n=",nrow(data_selected))) + 
    facet_wrap(~X1, ncol = 5) + 
    theme_bw() +
    theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=11, face="bold",angle = 70, hjust=1),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_text(color="black", size=16, face="bold"),
      legend.text = element_text(color="black", size=12, face="bold"),
      strip.text.x = element_text(size = 12, face="bold")
    )
  plot(p)
  dev.off()
} else {
  data_selected <- data[grep(opt$chr, data$X1),]
  pdf(paste0(opt$bedfile,"_",opt$chr,".pdf"))
  p <- ggplot(data_selected) +
    geom_histogram(aes(x=length), binwidth = 10) +
    geom_vline(xintercept = median(data_selected$length), color="red", linetype = "dashed", size = 2) + 
    labs(x="Bed interval length (bp)", y = "Counts", title = paste0(opt$chr," n=",nrow(data_selected))) + 
    theme_bw() +
    theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=11, face="bold",angle = 70, hjust=1),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_text(color="black", size=16, face="bold"),
      legend.text = element_text(color="black", size=12, face="bold"),
      strip.text.x = element_text(size = 12, face="bold")
    )
  plot(p)
  dev.off()
}
cat("Processing", opt$bedfile,"...\n")
cat("Totally",nrow(data_selected),"intervals Processed,Median length is [",median(data_selected$length),"bp]\n")


