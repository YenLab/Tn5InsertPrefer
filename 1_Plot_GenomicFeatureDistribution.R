setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(RColorBrewer)
library(qdap)

#=========================================================================================
#1. Integrate data matrix
#=========================================================================================
allfiles <- list.files(".", pattern = "*feature.out")
for (species in c("Mouse","Human")){
  
  files <- allfiles[c(grep(species,allfiles))]
  
  nucleo <- read.table(paste0(species,".feature_nucleotide_freqs.txt"), row.names = 1,sep = "\t",header = T)
  #Get all unique feature names list
  feature_list <- c()
  for (file in files){
    current_data <- read.table(file, header = T, row.names = 1)
    #Remove genome name mm10/hg38 ...
    rownames(current_data) <- sub(".*?\\.", "", rownames(current_data))
    #Get a consensus feature vector 
    feature_list <- unique(c(feature_list, rownames(current_data)))
    #Assign the matrix to variable for downstream use
    assign(file, current_data)
  }
  
  #Set the fold change/P-value dataframe
  df_value <- as.data.frame(matrix(0, nrow = length(feature_list), ncol = length(files)))
  rownames(df_value) <- feature_list
  colnames(df_value) <- files
  df_p <- df_value
  
  for (i in files){
    #i <- files[1]
    tmp <- eval(parse(text = i))
    sup <- setdiff(feature_list, rownames(tmp))
    sup <- as.data.frame(matrix(NA, nrow = length(sup), ncol = 2), row.names = sup)
    colnames(sup) <- colnames(tmp)
    tmp <- rbind(tmp,sup)
    tmp_ <- tmp[feature_list,]
    df_value[which(colnames(df_value) == i)] <- tmp_$Statistics
    df_p[which(colnames(df_p) == i)] <- tmp_$P.value
  }
  
  p_adjust_dataframe <- function(df){
    dims <- dim(df)
    vec <- c()
    for(ncol in 1:dims[2]){vec <- c(vec, df[,ncol])}
    t <- matrix(p.adjust(vec, method = "fdr"), nrow = dims[1], ncol=dims[2])
    rownames(t) <- rownames(df)
    colnames(t) <- colnames(df)
    return(t)
  }
  df_p <- p_adjust_dataframe(df_p)

  for (i in c("all","naked","chromatin")){

    df_value_selected <- df_value
    df_p_selected <- df_p
    if (i == "naked"){
      id <- grep("input", colnames(df_p))
      df_value_selected <- df_value_selected[,id]
      colnames(df_value_selected) <- sub("_input","",colnames(df_value_selected))
      df_p_selected <- df_p_selected[,id]
    } else if (i == "chromatin"){
      id <- grep("input",colnames(df_p))
      df_value_selected <- df_value_selected[,-id]
      df_p_selected <- df_p_selected[,-id]
    }
    
    #Exclude some features
    features_exclued <- c("repeats","Retroposon","centromeres","gaps","TTS","TSS","repeats_Other",
                          "repeats_Unknown","promoters","all_exons","miscRNA","snoRNA","rna","repeats_scRNA",
                          "repeats_srpRNA","repeats_RNA","miRNA","pseudogen","start_codons","repeats_Low_complexity_polypyrimidin",
                          "repeats_Low_complexity_polypurin","splice5p","splice3p","ncRNA","repeats_snRNA",
                          "stop_codons","repeats_RC","TTS+-0.5k","TTS+-1k","protein_coding","repeats_DNA",
                          "repeats_Low_complexity_GA_rich","repeats_Low_complexity_G_rich"
    )
    df_value_selected <- df_value_selected[setdiff(feature_list,features_exclued),]
    df_p_selected <- df_p_selected[setdiff(feature_list,features_exclued),]
    
    nucleo <- nucleo[rownames(df_value_selected),]
    df_value_selected$GC <- (nucleo$C + nucleo$G)*100
    
    df_value_selected$Anno <- "Basic Annotation"
    df_value_selected$Anno[grep("ENCODE",rownames(df_value_selected))] <- "ENCODE cCREs"
    df_value_selected$Anno[grep("repeats",rownames(df_value_selected))] <- "Repeats Regions"
    
    #replace sample names for better visualization
    meta_patern <-c("all_coding","ENCODE_cCREs_dELS","TSS+-0.5k","ENCODE_cCREs_DNaseH3K4me3_CTCFboun","repeats_Satellit",
                    "repeats_SINE","repeats_rRNA","ENCODE_cCREs_PLS","repeats_Low_complexity_AT_rich","repeats_Low_complexity_GC_rich",
                    "ENCODE_cCREs_DNaseH3K4me3","TTS+-2k","repeats_Low_complexity_A_rich","repeats_Simple_repeat","repeats_Low_complexity",
                    "ENCODE_cCREs_pELS","repeats_Low_complexity_C_rich","repeats_LTR","intergenic","random_500bp_20000",
                    "ENCODE_cCREs_CTCFonly_CTCFboun","introns","utr5","ENCODE_cCREs_PLS_CTCFboun","ENCODE_cCREs_pELS_CTCFboun",
                    "repeats_LINE","utr3","TSS+-1k","ENCODE_cCREs_dELS_CTCFboun","repeats_tRNA",
                    "TSS+-2k","cpgIslan","repeats_Low_complexity_T_rich")
    
    meta_replacement <- c("Exons","dELS","TSS (+-0.5kb)","H3K4me3-DNase (CTCF)","Satellite",
                          "SINE","rRNA","PLS","Low Comp. (AT-rich)","Low Comp. (GC-rich)",
                          "H3K4me3-DNase","TTS (+-2kb)","Low Comp. (A-rich)","Simple Repeat","Low Comp. (ALL)",
                          "pELS","Low Comp. (C-rich)","LTR","Intergenic","Random Regions",
                          "CTCF only","Introns","5' UTR","PLS (CTCF)","pELS (CTCF)",
                          "LINE","3' UTR","TSS (+-1kb)","dELS (CTCF)","tRNA",
                          "TSS (+-2kb)","CpG Island","Low Comp. (T-rich)")
    
    rownames(df_value_selected) <- mgsub(meta_patern, meta_replacement, rownames(df_value_selected))
    colnames(df_value_selected) <- gsub("_ATACseq.*|_THSseq.*|Mouse_|Human_|filtered.*","", colnames(df_value_selected))
    
    colnames(df_p_selected) <- colnames(df_value_selected)[-c((ncol(df_value_selected)-1 ):ncol(df_value_selected))]
    rownames(df_p_selected) <- rownames(df_value_selected)
    
    #Get log2 values for plotting
    scaled <- df_value_selected
    scaled[,-c((ncol(df_value_selected)-1 ):ncol(df_value_selected))] <- 
      round(log2(df_value_selected[,-c((ncol(df_value_selected)-1 ):ncol(df_value_selected))] + 0.0001), 2)
    
  
    if (species == "Human"){
      rank_plot <- meta_replacement[c(30,2,21,4,29,11,16,1,6,27,31,12,7,25,22,18,8,14,
                                      28,3,23,24,32,20,19,26,5,15,13)]
    } else {
      rank_plot <- meta_replacement[c(30,2,21,4,29,11,16,1,6,27,31,12,22,14,5,7,25,
                                      28,3,8,23,24,32,10,20,18,19,26,15,17,9,33,13)]}
    
    scaled <- scaled[rank_plot,]; df_p_selected <- df_p_selected[rank_plot,]
    
    write.table(scaled,paste0(species,"_OE.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
    write.table(df_p_selected,paste0(species,"_FDR.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
      
    pdf(paste0(species,"_feature_distribution_",i,".pdf"), height = 10, width = 8)
    ha = rowAnnotation(GC = scaled$GC,type = scaled$Anno, 
                       gp = gpar(col = "#FFFFFF"), 
                       col = list(
                         type = c("Basic Annotation" = "#EEDC82", "ENCODE cCREs" = "#2E8B57", "Repeats Regions" = "#8B1C62"),
                         GC = colorRamp2(breaks = c(0,80),colors = c("white","black"))
                         )
                       )

   p <- Heatmap(scaled[,-c((ncol(scaled)-1 ):ncol(scaled))],
                 show_row_names = T, show_column_names = T,
                 row_names_gp = gpar(fontsize = 12, fontface = "bold"),
                 column_names_gp = gpar(fontsize = 11, fontface = "bold"),
                 column_names_rot = 70, column_names_centered = F,
                 column_dend_height = unit(3, "cm"),
                 #column_title = "Distribution of Tn5 cut sites\n(*FDR < 0.001, Chi-square test)",
                 clustering_method_columns = "complete", 
                 col = colorRamp2(breaks = c(-1,0,1),colors = c("blue","white","red")),
                 rect_gp = gpar(col = "black", lwd = 1),
                  show_row_dend = F,
                 #column_split = 2, #row_split = 5,
                 cluster_columns = T, cluster_rows = T,
                 name = "Log2(O/E)",
                 right_annotation = ha,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if(!is.na(df_p_selected[i, j])){
                     if(df_p_selected[i, j] > 0.01){
                       grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "black", fill = "Grey"))
                       }}
                 }
    )
    draw(p)
  dev.off()
  }
}
