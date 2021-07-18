# This script is used for running GLM model on Tn5

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./")
suppressPackageStartupMessages(library(DNAshapeR))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qdap))

cl <- makePSOCKcluster(10)
registerDoParallel(cl)

TrainModelsTn5 <- function(trainning_final=trainning,sample_prefix=sample_prefix,
                           trainning_plan=trainning_plan, featureName=featureName,
                           out_file=out_file,varImp_list=varImp_list){
  
  #Cut sites as 1, uncut sites as 0
  trainning_final$Tn5_cut <- factor(trainning_final$Tn5_cut, levels=c(0,1),labels = c("uncut","cut"))
  result_prefix <- paste0(sample_prefix, "_", trainning_plan,"_",paste(featureName, collapse = ","))
  result_prefix <- sub("1-MGW,1-HelT,1-ProT,1-Roll,1-EP,1-Stretch,1-Tilt,1-Buckle,1-Shear,1-Opening,1-Rise,1-Shift,1-Stagger,1-Slide","14shapes",result_prefix)
  
  cat("-> Start trainning on", result_prefix, "<- ")

  train_grid <- expand.grid(.alpha = seq(0, 0.5, length = 5), .lambda = c((1:3)/10))
  set.seed(0218)
  trcl <- trainControl(method = "cv", number = 10, search = "grid", returnResamp = "final", 
                       savePredictions = "final", classProbs = TRUE, summaryFunction = prSummary, allowParallel = T)
  
  fitted_model <- train(Tn5_cut ~ ., data = trainning_final,
                        method = "glmnet", preProc = c("range"), trControl = trcl,
                        metric = "AUC", tuneGrid = train_grid)
  
  cat("Trainning finished\n")
  
  res_table <- fitted_model$pred[,c("pred","obs","uncut","cut")]
  AUROC <- twoClassSummary(res_table, lev = levels(res_table$obs))
  AUPRC <- prSummary(res_table, lev = rev(levels(res_table$obs)))
  acc <- sum(res_table$pred == res_table$obs)/nrow(res_table)
  
  result <- paste(AUROC[[1]],AUROC[[2]],AUROC[[3]],AUPRC[[1]], AUPRC[[2]],AUPRC[[3]],AUPRC[[4]],acc, result_prefix)
  write(result, file = out_file, append = TRUE)
  
  #Store this model
  saveRDS(fitted_model, paste0(result_prefix, "_GLM_model.rds"))

  #Add predictor importance (I have modified the source code in caret to report pos/neg values)
  varImp_list[[result_prefix]] <- varImp(fitted_model, scale = F)$importance
}

# ==============================================================================
# 1. Main trainning codes
# ==============================================================================
files <- list.files("./",pattern = ".*.fa")
prefixs <- sub("\\.fa","",files)

prefixs <- prefixs[grep("shuffled",prefixs)]

for (sample_prefix in prefixs){
  
  cat("--->>> Start Processing", sample_prefix, "...\n")
  cat("Calculating Motif affinity ...\n")

  motifPrefix="/public/home/zhy/Tn5_bias/fimo/motifs/"
  motif <- read.table(paste0(motifPrefix,str_extract(sample_prefix,".*?_"),"Tn5_motif_MEME_prob.txt"))
  colnames(motif) <- c("A","C","G","T")

  fa <- paste0(sample_prefix,".fa")
  fastr_mt <- as.matrix(readDNAStringSet(fa))
  res_matrix <- matrix(NA, nrow = nrow(fastr_mt), ncol = 51)
  
  for (j in 1:ncol(fastr_mt)){
    res_matrix[,j] <- mgsub(c("A","C","G","T","N"), c(motif[j,1],motif[j,2],motif[j,3],motif[j,4],0.25), as.vector(fastr_mt[,j]))
  }
  mode(res_matrix) <- "numeric"

  bed_file <- paste0(sub("_shuffled","",sample_prefix),".bed")
  bed <- read.table(bed_file)

  #Set cut sites into 1 for a cut event (classification)
  bed[which(bed$V4 > 0), 4] <- 1
  
  #Merge motif data into bed backbone file
  bed <- cbind(bed,as.data.frame(res_matrix))
  colnames(bed) <- c("chr","st","ed","cut","strand", paste0("motif",1:51))
  colna <- paste0(bed$chr,bed$st,bed$ed)
  
  # 2. Get DNA shape ----
  shape <- readRDS(paste0(sample_prefix, "_DNAshape.rds"))
  all_shape_mer <- c("1-MGW", "1-HelT","1-ProT","1-Roll","1-EP",  
                     "1-Stretch","1-Tilt","1-Buckle","1-Shear","1-Opening","1-Rise","1-Shift","1-Stagger","1-Slide")
  
  featureNames <- list(
    c("1-MGW"), c("1-HelT"), c("1-ProT"), c("1-Roll"), c("1-EP"), 
    c("1-Stretch"), c("1-Tilt"), c("1-Buckle"), c("1-Shear"), c("1-Opening"),
    c("1-Rise"), c("1-Shift"), c("1-Stagger"), c("1-Slide"), 
    c("1-mer"),c("2-mer"),c("3-mer"),
    c("1-mer","2-mer"),c("2-mer","3-mer"),c("1-mer","3-mer"),
    c("1-mer","2-mer","3-mer"),
    c("1-MGW","1-ProT","1-Roll","1-HelT"), #4shapes
    c(all_shape_mer)  #14shapes
  )
  
  # 3. Start trainning ----
  #For store feature importance
  varImp_list <- list()
  #For store result ML performance matrices
  out_file <- paste0(sample_prefix,"_matrice.txt")
  
  #Only for motif
  trainning <- data.frame(Tn5_cut = bed$cut, bed[,6:56])
  TrainModelsTn5(trainning_final=trainning,sample_prefix=sample_prefix,
                 trainning_plan="motif", featureName="", out_file= out_file, varImp_list=varImp_list)
  
  for(featureName in featureNames){
    #featureName <- featureNames[[1]]
    featureVector <- encodeSeqShape(fa, shape, featureName, normalize = F)
    
    trainning_plans <- c("shape", "motif_shape")
    for (trainning_plan in trainning_plans){
      #trainning_plan <- trainning_plans[6]
      if (trainning_plan == "shape"){ 
        trainning <- data.frame(Tn5_cut = bed$cut, featureVector)
      } else if (trainning_plan == "motif_shape"){
          trainning <- data.frame(Tn5_cut = bed$cut, bed[,6:56], featureVector)
        } 
      TrainModelsTn5(trainning_final=trainning,sample_prefix=sample_prefix,
                     trainning_plan=trainning_plan, featureName=featureName,
                     out_file= out_file,varImp_list=varImp_list)
    }
  }
  #save varImp_list
  saveRDS(varImp_list,paste0(sample_prefix,"_varImp_list.rds"))
}

stopCluster(cl)

#tail -n 1  *txt | grep -v "=" | grep -vE "^$" | cut -d " " -f 8 | sort -k1nr