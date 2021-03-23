# This script is used for plotting DNA shapes around defined groups

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

#=========================================================================================
#1. Batch run on all samples
#=========================================================================================
Inter_bp_shape <- c("Rise","Shift","Slide","HelT","Roll","Tilt") #long_shape
Intra_bp_shape <- c("ProT","Stretch","Buckle","Shear","Opening","Stagger") #short_shape
groove_shape <- c("MGW","EP")
All_shapes <- c(Inter_bp_shape,Intra_bp_shape,groove_shape)

Plot_DNA_shape <- function(path, flank_window = 50){

   if (! file.exists("mean_DNAshape.txt")){
    files <- list.files(path, ".rds")
    for (file in files){
      cat("Processing",file,"...\n")
      tmp <- readRDS(paste0(path,"/",file))
      shapes <- names(tmp)
      
      NUMC <- ncol(tmp[[shapes[1]]])
      mean_shape <- matrix(0, nrow = length(shapes), ncol = NUMC) %>% as.data.frame()
      rownames(mean_shape) <- shapes
      
      for (shape in shapes){
        #shape <- shapes[1]
        if (shape %in% Inter_bp_shape){
          mean_shape[shape,1:(NUMC-1)] <- colMeans(tmp[[shape]], na.rm = T)
        } else {
          mean_shape[shape,] <- colMeans(tmp[[shape]], na.rm = T)
        }
      }
      rownames(mean_shape) <- paste0(rownames(mean_shape), "_" ,file)
      assign(paste0(file,"_mean_shape"), mean_shape)
    }
    
      v <- paste0(files,"_mean_shape")
      final_matrix <- do.call(rbind,mget(v))
      write.table(final_matrix,"mean_DNAshape.txt",row.names = T, col.names = T,quote = F)
   }
   
  final_matrix <- read.table("mean_DNAshape.txt", header = T, row.names = 1)
  n <- (ncol(final_matrix)-1)/2 + 1
  final_matrix <- t(final_matrix[,(n-flank_window):(n+flank_window)]) %>% as.data.frame()
  final_matrix$position <- -flank_window:flank_window
  
  final_matrix_melt <- reshape2::melt(final_matrix,c("position"))
  
  final_matrix_melt$variable <- gsub(".*_mean_shape.|_ext.*","",final_matrix_melt$variable)
  final_matrix_melt$shape <- gsub("_.*","",final_matrix_melt$variable)
  final_matrix_melt$sample <- sub(".*?_","",final_matrix_melt$variable)
  final_matrix_melt$sample <- gsub("_DNAshape.rds|_\\+\\-100bp_shuf100000|_201bp_shuf100000","",final_matrix_melt$sample)
  final_matrix_melt <- final_matrix_melt[-grep("cut_unmotif_shuffled|mm10.random_shuffled|motif_cut_shuffled|motif_uncut_shuffled|cut_unmotif_matched|motif_cut_matched|motif_uncut_matched",final_matrix_melt$sample),]

  pdf(paste0("FigS2A.DNA shapes around +-", flank_window,"bp.pdf"), height = 5, width = 8)
  p <- ggplot(final_matrix_melt[grep("MGW|HelT|ProT|Roll",final_matrix_melt$shape),], aes(position, value, color = sample)) + 
    geom_smooth(method="loess", se=F, fullrange=T, level=0.95, span=0.5) + 
    facet_wrap(~shape, ncol = 2, scales = "free_y") +
    theme_bw() +
    scale_colour_manual(values = c("#1F78B4","#525252","#1F78B4","#A6CEE3", "grey", "#33A02C", "#B2DF8A")) + 
    labs(x="Relative to motif/insertion sites center (bp)", y = "Mean shape values") +
    theme(
      plot.title = element_text(color="black", size=20, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=11, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_text(color="black", size=14, face="bold"), 
      legend.text = element_text(color="black", size=12, face="bold"),
      strip.text.x = element_text(size = 14, face="bold")
    ) 
  plot(p)
  dev.off()
}

for (i in c(50,98)){
  Plot_DNA_shape("./", flank_window = i)
}

#=========================================================================================
#2. Do statistical test
#=========================================================================================
stat <- final_matrix
colnames(stat) <- gsub(".*_mean_shape.|_ext.*","",colnames(stat))
stat <- stat[,grep("MGW|HelT|ProT|Roll",colnames(stat))]
colnames(stat) <- gsub("\\+\\-100bp_shuf100000_|_DNAshape.rds|DNAshape.rds","",colnames(stat))

sss <- final_matrix_melt[grep("MGW|HelT|ProT|Roll",final_matrix_melt$shape),]
pdf("stat.pdf", height = 6,width = 12)
ggplot(sss[grep("motif_cut$|motif_uncut$|random",sss$sample),], aes(y=value ,x=sample, fill=sample)) +
  geom_boxplot(notch=T, outlier.shape = NA) +
  coord_flip() +
  stat_compare_means(mapping=aes(as_label=..p.adj..), method = "wilcox.test",
                     comparisons = list(c("motif_uncut","motif_cut"))) +
  labs(x="", y = "Tn5 cut relative to mathced median") +
  facet_wrap(~shape, ncol=1, scales = "free") + 
  theme_bw() +
  scale_fill_manual(values = c("grey","#33A02C", "#B2DF8A")) + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_text(color="black", size=11, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(color="black", size=11, face="bold"),
    legend.title = element_text(color="black", size=16, face="bold"),
    legend.text = element_text(color="black", size=12, face="bold"),
    strip.text.x = element_blank(),
    legend.position="bottom")
dev.off()





