setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(networkD3)

#This script is for analyzing DNAme effect on Tn5 cut
#=========================================================================================
# Integrate data matrix and normalization
#=========================================================================================
#Read four categories and merge into a big file
A_higher <- read.table("E14_higher_unique.bedfa", header = F)
A_higher$condition <- "A_Only"
B_higher <- read.table("Germ_higher_unique.bedfa", header = F)
B_higher$condition <- "B_Only"
All_zeros <- read.table("All_zeros_unique.bedfa", header = F)
All_zeros$condition <- "All_zeros"
All_high <- read.table("All_high_unique.bedfa", header = F)
All_high$condition <- "All_high"

all_df <- rbind(A_higher, B_higher, All_zeros, All_high)
colnames(all_df) <- c("Kmers","DNAme_A","DNAme_B","Tn5_A","Tn5_B","Condition")

#Only select 9mer shared by four categories
common_kmer <- Reduce(intersect, list(A_higher$V1,B_higher$V1,All_zeros$V1,All_high$V1))
all_df <- all_df[all_df$Kmers %in% common_kmer,]
length(unique(all_df$Kmers))
trim_outliers <- function(all_df_scaled, trim = 1){
  #Trim outliers based on percentage
  trim = trim
  Bottom5_A <- quantile(all_df_scaled$Tn5_A, probs = seq(0,1,0.01))[[trim + 1]]
  Top5_A <- quantile(all_df_scaled$Tn5_A, probs = seq(0,1,0.01))[[101 - trim]]
  
  Bottom5_B <- quantile(all_df_scaled$Tn5_B, probs = seq(0,1,0.01))[[trim + 1]]
  Top5_B <- quantile(all_df_scaled$Tn5_B, probs = seq(0,1,0.01))[[101 - trim]]
  
  index <- all_df_scaled$Tn5_A > Top5_A | all_df_scaled$Tn5_B > Top5_B | all_df_scaled$Tn5_A < Bottom5_A | all_df_scaled$Tn5_B < Bottom5_B 
  remove_9mer <- unique(all_df_scaled[index,]$Kmers)
  all_df_scaled <- all_df_scaled[!(all_df_scaled$Kmers %in% remove_9mer),]
  return(all_df_scaled)
}

#Z-scale
all_df[,4:5] <- scale(all_df[,4:5], center = T, scale = T)

#Trim top 1% outliers to relieve data skewness
all_df_scaled <- trim_outliers(all_df, trim = 1)
all_df_selected_1 <- reshape2::melt(all_df_scaled[c("Tn5_A","Tn5_B","Condition")],c("Condition"))

pdf("Fig2A.DNAmeEffect_4Groups_E14GermNakedDNA.pdf", height = 5, width = 8)
p <- ggplot(all_df_selected_1, aes(x = variable, y = value, fill = variable)) + 
  geom_boxplot(notch=F, outlier.shape = NA, size=0.8) +
  stat_compare_means(method = "t.test", paired = F, comparisons = list(c("Tn5_A", "Tn5_B"))) +
  facet_wrap(~Condition, ncol=6, scales = "free_x") + 
  scale_fill_manual(values = c("#66A61E","#D95F02")) + 
  scale_y_continuous(breaks = seq(-2,2,1)) +
  labs(x = paste0(length(unique(all_df_scaled$Kmers))," 9mers"), 
       y = "Relative Tn5 insert. freq.\n(Z-score)") + 
  theme_bw() +
  theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.text.x = element_blank(),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(color="black", size=11, face="bold"),
    legend.title = element_text(color="black", size=14, face="bold"), 
    legend.text = element_text(color="black", size=12, face="bold"),
    strip.text.x = element_text(size = 20, face="bold"),
    strip.text.y = element_text(size = 14, face="bold")
  )
plot(p)
dev.off()

#=========================================================================================
# Titrate DNA methylation level in each 9mer
#=========================================================================================
train <- read.table("refined10_ESC_collapsed")
train <- train[-grep("11|^0",train$V2),]

trim_outliers <- function(all_df_scaled=train, trim = 1){
  #Trim outliers based on percentage
  cat(length(unique(all_df_scaled[,1])),"kmers input !!!\n")
  B1<- quantile(all_df_scaled[,3], probs = seq(0,1,0.01))[[trim + 1]]
  T1 <- quantile(all_df_scaled[,3], probs = seq(0,1,0.01))[[101 - trim]]
  index1 <- all_df_scaled[,3] > T1 | all_df_scaled[,3] < B1 
  
  remove <- index1
  all_df_scaled <- all_df_scaled[!remove,]
  cat(sum(remove),"kmers were removed and",length(unique(all_df_scaled[,1])),"kmers were kept !!!\n" )
  return(all_df_scaled)
}
train <- trim_outliers(train,trim = 1)
tabs <- as.data.frame(table(train$V1))
remove <- tabs$Var1[which(tabs$Freq!= 10)]
train_selected <- train[! train$V1 %in% remove,]
train_selected$V2 <- as.factor(train_selected$V2)

pdf("Fig2B.DNAmeEffectDeciles_ESCnakedDNA.pdf", height = 8, width = 8)
ggplot(train_selected, aes(x = V2, y = V3)) + 
  geom_boxplot(fill = "#66A61E", notch=F, outlier.shape = NA, size=0.8) +
  stat_compare_means(method = "t.test", paired = F, 
                     comparisons = list(c("1", "2"),c("1", "5"),c("5", "10"),c("9", "10"),c("1", "10"))) +
  scale_y_continuous(breaks = seq(-1,3,1)) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 5, name = "Greys"))(13)) + 
  labs(x = "DNA methylation (Deciles)", y = "Relative Tn5 insert. freq.\n(Z-score)") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.text.x  = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y  = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(color="black", size=18, face="bold"), 
    legend.text  = element_text(color="black", size=14, face="bold"),
    strip.text.x = element_text(size = 20, face="bold"),
    strip.text.y = element_text(size = 20, face="bold")
  )
dev.off()

#=========================================================================================
# Titrate DNA methylation level in each 9mer (chromatin)
#=========================================================================================
train <- read.table("refined10_ESC_peaks_collapsed.txt")
train <- train[-grep("11",train$V2),]

trim_outliers <- function(all_df_scaled=train, trim = 1){
  #Trim outliers based on percentage
  cat(length(unique(all_df_scaled[,1])),"kmers input !!!\n")
  B1<- quantile(all_df_scaled[,3], probs = seq(0,1,0.01))[[trim + 1]]
  T1 <- quantile(all_df_scaled[,3], probs = seq(0,1,0.01))[[101 - trim]]
  index1 <- all_df_scaled[,3] > T1 | all_df_scaled[,3] < B1 
  
  remove <- index1
  all_df_scaled <- all_df_scaled[!remove,]
  cat(sum(remove),"kmers were removed and",length(unique(all_df_scaled[,1])),"kmers were kept !!!\n" )
  return(all_df_scaled)
}

train <- trim_outliers(train,trim = 1)
tabs <- as.data.frame(table(train$V1))
remove <- tabs$Var1[which(tabs$Freq!=11)]
train_selected <- train[! train$V1 %in% remove,]
train_selected$V2 <- as.factor(train_selected$V2)

pdf("DNAmeEffectDeciles_ESCChromatin.pdf", height = 8, width = 8)
ggplot(train_selected, aes(x = V2, y = V3)) + 
  geom_boxplot(fill = "#66A61E", notch=F, outlier.shape = NA, size=0.8) +
  stat_compare_means(method = "t.test", paired = F, 
                     comparisons = list(c("1", "2"),c("1", "5"),c("5", "10"),c("9", "10"),c("1", "10"))) +
  scale_y_continuous(breaks = seq(-1,3,1)) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 5, name = "Greys"))(13)) + 
  labs(x = "DNA methylation (Deciles)", y = "Relative Tn5 insert. freq.\n(Z-score)") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.text.x  = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y  = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(color="black", size=18, face="bold"), 
    legend.text  = element_text(color="black", size=14, face="bold"),
    strip.text.x = element_text(size = 20, face="bold"),
    strip.text.y = element_text(size = 20, face="bold")
  )
dev.off()

#Nucleosome
train <- read.table("refined10_ESC_peaksMNase_collapsed.txt")
train <- train[-grep("11",train$V2),]

trim_outliers <- function(all_df_scaled=train, trim = 1){
  #Trim outliers based on percentage
  cat(length(unique(all_df_scaled[,1])),"kmers input !!!\n")
  B1<- quantile(all_df_scaled[,5], probs = seq(0,1,0.01))[[trim + 1]]
  T1 <- quantile(all_df_scaled[,5], probs = seq(0,1,0.01))[[101 - trim]]
  index1 <- all_df_scaled[,5] > T1 | all_df_scaled[,5] < B1 
  
  remove <- index1
  all_df_scaled <- all_df_scaled[!remove,]
  cat(sum(remove),"kmers were removed and",length(unique(all_df_scaled[,1])),"kmers were kept !!!\n" )
  return(all_df_scaled)
}

train <- trim_outliers(train,trim = 1)
tabs <- as.data.frame(table(train$V1))
#Use shared 9mers
remove <- tabs$Var1[which(tabs$Freq!=11)]
train_selected <- train[! train$V1 %in% remove,]
train_selected$V2 <- as.factor(train_selected$V2)

pdf("NucleosomeDeciles_ESCChromatin.pdf", height = 8, width = 8)
ggplot(train_selected, aes(x = V2, y = V5)) + 
  geom_boxplot(fill = "#D95F02", notch=F, outlier.shape = NA, size=0.8) +
  stat_compare_means(method = "t.test", paired = F, 
                     comparisons = list(c("1", "2"),c("1", "5"),c("5", "10"),c("9", "10"),c("1", "10"))) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 5, name = "Greys"))(13)) + 
  labs(x = "DNA methylation (Deciles)", y = "Relative Tn5 insert. freq.\n(Z-score)") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.text.x  = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y  = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(color="black", size=18, face="bold"), 
    legend.text  = element_text(color="black", size=14, face="bold"),
    strip.text.x = element_text(size = 20, face="bold"),
    strip.text.y = element_text(size = 20, face="bold")
  )
dev.off()

#
train <- read.table("refined10_ESC_peaksdeMNase_collapsed.txt")
train <- train[-grep("11",train$V2),]

trim_outliers <- function(all_df_scaled=train, trim = 1){
  #Trim outliers based on percentage
  cat(length(unique(all_df_scaled[,1])),"kmers input !!!\n")
  B1<- quantile(all_df_scaled[,4], probs = seq(0,1,0.01))[[trim + 1]]
  T1 <- quantile(all_df_scaled[,4], probs = seq(0,1,0.01))[[101 - trim]]
  index1 <- all_df_scaled[,4] > T1 | all_df_scaled[,4] < B1 
  
  remove <- index1
  all_df_scaled <- all_df_scaled[!remove,]
  cat(sum(remove),"kmers were removed and",length(unique(all_df_scaled[,1])),"kmers were kept !!!\n" )
  return(all_df_scaled)
}

train <- trim_outliers(train,trim = 1)
tabs <- as.data.frame(table(train$V1))
remove <- tabs$Var1[which(tabs$Freq!=11)]
train_selected <- train[! train$V1 %in% remove,]
train_selected[,4] <- scale(train_selected[,4] )
train_selected$V2 <- as.factor(train_selected$V2)

pdf("DNAmeEffectDeciles_ESCChromatin.pdf", height = 8, width = 8)
ggplot(train_selected, aes(x = V2, y = V4)) + 
  geom_boxplot(fill = "#66A61E", notch=F, outlier.shape = NA, size=0.8) +
  stat_compare_means(method = "t.test", paired = F, 
                     comparisons = list(c("1", "2"),c("1", "5"),c("5", "10"),c("9", "10"),c("1", "10"))) +
  scale_y_continuous(breaks = seq(-1,3,1)) + 
  scale_fill_manual(values = colorRampPalette(brewer.pal(n = 5, name = "Greys"))(13)) + 
  labs(x = "DNA methylation (Deciles)", y = "Relative Tn5 insert. freq.\n(Z-score)") + 
  theme_bw() +
  theme(
    axis.title.x = element_text(color="black", size=18, face="bold"),
    axis.text.x  = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=18, face="bold"),
    axis.text.y  = element_text(color="black", size=14, face="bold"),
    legend.title = element_text(color="black", size=18, face="bold"), 
    legend.text  = element_text(color="black", size=14, face="bold"),
    strip.text.x = element_text(size = 20, face="bold"),
    strip.text.y = element_text(size = 20, face="bold")
  )
dev.off()