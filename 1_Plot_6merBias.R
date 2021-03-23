#This script is for analyze Tn5 insertion preference in 6mers from (Calviello et al., 2019)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(stringr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(Biostrings)
library(ggpubr)

mer6 <- read.table("SeqBias_ATAC.txt")
mer6$V2 <- (mer6$V2 - min(mer6$V2))/(max(mer6$V2) - min(mer6$V2))
mer6$V2 <- mer6$V2[order(mer6$V2, decreasing = T)]

mer6$GC_content <- letterFrequency(DNAStringSet(mer6$V1), "GC", as.prob=T)
mer6$rich <- cut(mer6$GC_content,breaks=c(unique(mer6$GC_content),1.1), include.lowest = T, 
    right = F, labels = c("AT(6)","AT(5)","AT(4)","AT(3)","AT(2)","AT(1)","AT(0)"))
mer6$cutQ <- cut(mer6$V2,breaks=c(0,mer6$V2[3596],mer6$V2[500],1), include.lowest = T, 
                 right = F, labels = c("Bottom 6mers (n = 500)","Middle 6mers (n = 3596)","Top 6mers (n = 500)"))

mer6$cutQ <- factor(mer6$cutQ, 
                    levels = c("Top 6mers (n = 500)","Middle 6mers (n = 3596)","Bottom 6mers (n = 500)"))

pdf("FigS1C.ATAC_6mer_bias.pdf", height = 6, width = 10)
ggplot(mer6,aes(rich,V2)) + 
  geom_boxplot(aes(fill = rich), color = c("#66A61E","#A6761D","#E7298A")[3])+
  ylim(0,1) + 
  stat_compare_means(method = "wilcox.test", paired = F, 
                     comparisons = list(c("AT(6)", "AT(3)"),
                                        c("AT(6)", "AT(1)"))) +
  facet_wrap(~ cutQ, ncol = 3, scales = "free_x") +
  labs(x = "", y = "Relative Tn5 insert. freq.") + 
  scale_fill_manual(values = rev(colorRampPalette(brewer.pal(3,"Greys"))(7))) + 
  theme_bw() + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=11, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(color="black", size=11, face="bold"),
    legend.title = element_text(color="black", size=14, face="bold"), 
    legend.text = element_text(color="black", size=12, face="bold"),
    strip.text.x = element_text(size = 14, face="bold"),
    strip.text.y = element_text(size = 14, face="bold")
  )
dev.off()










