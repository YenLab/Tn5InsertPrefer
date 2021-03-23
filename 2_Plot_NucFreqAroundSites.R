#This script is for plotting nucleotide frequency and motif around Tn5 insertion sites

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#=========================================================================================
#1. This function is for plotting nucleotide frequency given frequency table
#=========================================================================================
plot_Nuc_freqs_around <- function(file, flank_window = 25, plot_type = "dinuc"){
  
  #get plotting region from given flank_window 
  freq <- read.table(file)
  prefix <- sub(".fa_nuc.txt","",file)
  species <- sub("_.*","",file)
  
  n <- ncol(freq)
  st <- (n+1)/2 - flank_window
  ed <- (n+1)/2 + flank_window
  freq <- freq[,c(st:ed)]
  
  #get plotting region from given flank_window (random file)
  freq_ran <- read.table(paste0(species,".random_200bp_500000.fa_nuc.txt"))
  n_ran <- ncol(freq_ran)
  st_ran <- (n_ran-1)/2 - flank_window
  ed_ran <- (n_ran-1)/2 + flank_window
  freq_ran <- freq_ran[,c(st_ran:ed_ran)]
  n <- ncol(freq)
  rownames(freq) <- c("A","C","T","G")
  colnames(freq) <- paste0(seq(-(n-1)/2,(n-1)/2))
  colnames(freq_ran) <- paste0(seq(-(n-1)/2,(n-1)/2))
  
  tmp1 <- t(freq_ran) %>% as.data.frame()
  tmp1$position <- seq(-(n-1)/2,(n-1)/2)
  tmp1$cond <- "Random"
  
  tmp2 <- t(freq)%>% as.data.frame()
  tmp2$position <- seq(-(n-1)/2,(n-1)/2)
  tmp2$cond <- "true"
  
  combined <- rbind(tmp1, tmp2)
  combined$GC <- combined$G + combined$C
  combined$AT <- combined$A + combined$T
  
  if(plot_type == "dinuc"){
    selected <- combined[,5:8]
  } else if (plot == "mononuc"){
    selected <- combined[,1:6]
  }
  selected_melt <- reshape2::melt(selected,c(c("cond","position")))
  selected_melt$group <- paste0(selected_melt$cond,"_",selected_melt$variable)
  
  
  #pdf(paste0(sub(".tab","",file),"_+-",flank,"bp.pdf"), height = 4, width = 16)
  ggplot(selected_melt, aes(x = position, y = value,group = group, color = group)) + 
    geom_line(aes(linetype=cond)) + 
    ylim(0.2,0.8) + 
    theme_bw() +
    labs(size= "Nitrogen", x="Relative postion to Tn5 cut center (bp)", y = "Nucleotide frequency") +
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
}

#Plot Nuc
files <- list.files(path = "./", pattern = "*_nuc.txt")
for (file in files){
  print(paste0("Processing Nuc ",file))
  plot_Nuc_freqs_around(file,flank_window = 25, plot_type = "dinuc")

}
#=========================================================================================
#3. This function is for plotting MEME result - using exported 'Probability Matrix' file
#=========================================================================================

plot_MEME_PSSM <- function(file, flank_window = 5){
  freqs <- read.table(file)
  prefix <- sub(".txt","",file)
  
  colnames(freqs) <- c("A","C","G","T")
  freqs <- t(freqs)
  
  st <- (ncol(freqs) + 1)/2 - flank_window
  ed <- (ncol(freqs) + 1)/2 + flank_window
  
  freqs <- freqs[,c(st:ed)]
  colnames(freqs) <- seq(-flank_window,flank_window)
  
  #GetConsensusSeq(freqs)
  pdf(paste0(prefix,"_MEME_motif_bits.pdf"),height = 4,width = 8)
  logomaker(freqs, type = "Logo", color_type = "per_row",
            logo_control = list(main_fontsize=20, yscale_change=F,
                                xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
                                xlab = "Position relative to Tn5 cut sites (bp)",ylab = "Bits"),
            colors = c(nuc_color[c(4,2,3,1)]))
  dev.off()
}

#plot_MEME_result_txt(meme_file, species="Mouse")
plot_MEME_result_txt <- function(meme_file, species=NULL){
  #Parse meme.txt result
  prefix <- dirname(meme_file)
  bn <- basename(dirname(meme_file))
  
  motif <- readLines(meme_file)
  
  num <- grep("MEME-1 position-specific probability matrix",motif)
  w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  
  if (w < 19){
    
    num <- grep("MEME-2 position-specific probability matrix",motif)
    w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  }
  
  if (w < 19){
    num <- grep("MEME-3 position-specific probability matrix",motif)
    w <- as.integer(str_extract(str_extract(motif[(num+2)],"w= [0-9]+"),"[0-9]+"))
  }
  
  Tn5_motif <- motif[(num+3):(num+2+w)]
  
  freqs <- matrix(0,nrow = w,ncol = 4)
  for (line in seq(1:length(Tn5_motif))){
    x <- Tn5_motif[line]
    freqs[line,] <- matrix(scan(text = x,quiet=T),nrow = 1,byrow = TRUE)[1,]
  }
  
  write.table(freqs,paste0(bn,"_Tn5_motif_1_freqs.txt"),quote = F,sep = "\t",col.names = F, row.names = F)
  
  colnames(freqs) <- c("A","C","G","T")
  freqs <- t(freqs)
  flank <- w - (ncol(freqs) + 1)/2
  colnames(freqs) <- seq(-flank,flank)
  
  nuc_color <- c('#109648','#255C99', '#D62839','#F7B32B') #ACTG
  
  if (!missing(species)){
    if (species == "Fish"){bg <- c(0.317, 0.183, 0.183, 0.317)}
    else if (species == "Fly"){bg <- c(0.29, 0.21, 0.21, 0.29)}
    else if (species == "huffia"){bg <- c(0.403, 0.097, 0.097, 0.403)}
    else if (species == "Human"){bg <- c(0.296, 0.204, 0.204, 0.296)}
    else if (species == "maize"){bg <- c(0.266, 0.234, 0.234, 0.266)}
    else if (species == "Mouse"){bg <- c(0.292, 0.208, 0.208, 0.292)}
    else if (species == "Plant"){bg <- c(0.32, 0.18, 0.18, 0.32)}
    else if (species == "Worm"){bg <- c(0.323, 0.177, 0.177, 0.323)}
  } else {bg <- c(0.25, 0.25, 0.25, 0.25)}
  names(bg) <- c("A", "C", "G", "T")
  
  #GetConsensusSeq(freqs)
  pdf(paste0(bn,"_MEME_motif_1_bits.pdf"), height = 4,width = 8)
  logomaker(freqs, type = "Logo", color_type = "per_row", bg=bg,
            logo_control = list(main_fontsize=20, yscale_change=F,
                                xaxis_fontsize=10, xlab_fontsize=15, y_fontsize=15,
                                xlab = "Position relative to Tn5 cut sites (bp)", ylab = "Bits"),
            colors = c(nuc_color[c(4,2,3,1)]))
  dev.off()
}

#Plot motif
dirs <- list.dirs(path = ".", full.names = F, recursive = F)
for (dir in dirs[3]){
  print(paste0("Processing motif ",dir))
  meme_files <- dir(path = dir, pattern = "meme.txt", full.names = T, recursive = T)
  for (meme_file in meme_files){
    print(meme_file)
    plot_MEME_result_txt(meme_file, species = dir)
  }
}

