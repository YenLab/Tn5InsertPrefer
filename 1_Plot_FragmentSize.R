#This script is for plotting fragment size distribution
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

largest_fs <- 2000
files <- list.files(path = "./", pattern = ".*filtered.dedup.insertSizes.txt$")
final_mat <- matrix(0,ncol = length(files), nrow = largest_fs) %>% as.data.frame()
colnames(final_mat) <- gsub("_ATACseq.*","",files)

for (file in files){
  #Read in the Picard result
  insertsize <- readLines(file)
  num <- grep("## HISTOGRAM",insertsize)
  data <- as.data.frame(insertsize[(num+2):(length(insertsize)-1)])
  data$position <- gsub("\t.*","",data[,1])
  data$freqs <- as.integer(sub("\t.*","",sub(".*?\t","",data[,1])))
  
  #For same result contain FR/RF
  if (str_detect(sub(".*?\t","",data[,1])[1],"\t")){
    data$freqs2 <- as.integer(sub(".*\t","",sub(".*?\t","",data[,1])))
    data$freqs <- data$freqs + data$freqs2
  }
  data$freqs2 <- 0
  final_mat[data$position, gsub("_ATACseq.*","",file)] <- data$freqs
}

#Max-min normalization of frequency for comparison among samples with different sequencing depth
normalize <- function(x){ return((x- min(x)) /(max(x)-min(x))) }
final_mat <- as.data.frame(apply(final_mat,2, normalize))

final_mat$position <- 1:largest_fs
final_mat_melt <- reshape2::melt(final_mat,c("position"))
final_mat_melt$species <- gsub("_.*","",final_mat_melt$variable)
final_mat_melt$cond <- "Chromatin"
final_mat_melt$cond[grep("input",final_mat_melt$variable)] <- "Naked_DNA"
final_mat_melt$gg <- paste0(final_mat_melt$species," ",final_mat_melt$cond)
final_mat_melt$gg <- factor(final_mat_melt$gg,levels = unique(final_mat_melt$gg)[c(1,3,2,4)])

pdf("FigS1A.Fragment_size_distribution.pdf", height = 3, width = 6)
p <- ggplot(final_mat_melt, aes(x = position, y = value, color=cond, group=variable)) + 
  geom_line(size=0.5, alpha=0.7) + 
  theme_bw() + xlim(0,1000) + ylim(0,1) +
  labs(x="Fragment length (bp)", y = "Releative frequency") +
  scale_color_manual(values = c("grey","blue")) +
  theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=11, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(color="black", size=11, face="bold"),
    legend.title = element_text(color="black", size=14, face="bold"), 
    legend.text = element_text(color="black", size=12, face="bold"),
    strip.text.x = element_text(size = 12, face="bold"),
    legend.position="none"
  ) +  facet_wrap(~species, ncol = 2)
plot(p)
dev.off()
