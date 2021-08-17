setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(pracma)
library(ggstatsplot)
options(max.print = 100)

PlotGiggle <- function(prefix){
  giggle1 <- read.table(paste0("SOB_corrected_share_",prefix,".tsv"))
  giggle1$group <- "corrected_share"
  
  giggle2 <- read.table(paste0("SOB_corrected_specific_",prefix,".tsv"))
  giggle2$group <- "corrected_specific"
  
  giggle3 <- read.table(paste0("SOB_uncorrected_shared_",prefix,".tsv"))
  giggle3$group <- "uncorrected_shared"
  
  giggle4 <- read.table(paste0("SOB_uncorrected_specific_",prefix,".tsv"))
  giggle4$group <- "uncorrected_specific"
  
  giggle <- rbind(giggle1,giggle2,giggle3,giggle4)
  colnames(giggle) <- c("file","file_size","overlaps","odds_ratio",
                        "fishers_two_tail","fishers_left_tail","fishers_rigth_tail","combo_score","group")
  giggle$group <- factor(giggle$group, levels = unique(giggle$group)[c(3,1,4,2)])
  giggle$file <- sub(".*/","", giggle$file)
  giggle <- separate(data = giggle, col = file, into = c("Mouse", "ESC","Cellline","TF","GSM","ID"), sep = "_")
  cat(length(unique(giggle$TF)),"types of factors was plotted...\n")
  
  pdf(paste0("SOB_",prefix,"_Peaks.pdf"),height = 8, width = 7)
  pf <- ggplot(giggle, aes(x = group, y = combo_score)) + 
    geom_boxplot(aes(fill=group),notch=F, outlier.shape = NA, size=0.8) +
    geom_jitter(aes(color=group),alpha = 0.1,width = 0.2) + 
    theme_bw() + #ylim(-1000,1000) +
    labs(x="",y="Giggle Score") +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = combn(as.vector(unique(giggle$group)), 2, simplify = FALSE)) +
    scale_color_manual(values = c("grey","grey","#FF7F00","#33A02C","black")) +
    scale_fill_manual(values = c("grey","grey","#FF7F00","#33A02C","black")) +
    theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=14, face="bold",angle = 40, hjust = 1),
      axis.title.y = element_text(color="black", size=18, face="bold"),
      axis.text.y = element_text(color="black", size=14, face="bold"),
      legend.position="none"
    )
  plot(pf)
  dev.off()
  
  for (i in 1:4){
    plan = as.vector(unique(giggle$group))[i]
    color=c("grey","#33A02C","grey","#FF7F00")
    lim=c(2500,1200,2500,1200)
    pdf(paste0("SOB_",prefix,"_giggle_score_rankedByTF_",plan,".pdf"), height = 20, width = 4)
    p <- ggplot(giggle[giggle$group == plan,], aes(x = reorder(TF, combo_score, FUN = median), y=combo_score)) + 
      geom_boxplot(notch=F, outlier.shape = NA, size=0.5, fill=color[i]) +
      geom_jitter(alpha = 0.3, width = 0.1,fill=color[i]) + 
      theme_bw() + coord_flip() +
      ylim(0,lim[i]) +
      labs(x="",y="Giggle Score") +
      theme(
        plot.title = element_text(color="black", size=14, face="bold"),
        axis.title.x = element_text(color="black", size=24, face="bold"),
        axis.text.x = element_text(color="black", size=18, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(color="black", size=10),
        legend.position="none"
      )
    plot(p)
    dev.off()
  }
}

PlotGiggle(prefix = "TFs_giggle")
PlotGiggle(prefix = "equal_giggle")
PlotGiggle(prefix = "matched_equal_giggle")
