# This script is used for calculate DNA shapes from given fasta file and 
# return a rds file contain all DNA shapes (14 types)

setwd("./")
suppressPackageStartupMessages(library(DNAshapeR))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-f", "--FASTA"), type="character", default=NULL,metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Inter_bp_shape <- c("Rise","Shift","Slide","HelT","Roll","Tilt")
Intra_bp_shape <- c("ProT","Stretch","Buckle","Shear","Opening","Stagger")
groove_shape <- c("MGW","EP")

shapes <- c(Inter_bp_shape,Intra_bp_shape,groove_shape)
long_shape <- c("HelT","Roll","Tilt","Rise","Shift","Slide")
short_shape <- setdiff(shapes, long_shape)

#===========================================================
#1. Calculate given fasta file
#===========================================================
prefix <- sub("\\.fa","",opt$FASTA)
fa <- paste0(prefix,".fa")
rds <- paste0(prefix, "_DNAshape.rds")

if (file.exists(rds)){
  cat(rds, "existed!")
} else {
  cat("DNA shapes are calculated for ",opt$FASTA,"...\n")
  shape <- getShape(fa, parse = TRUE)
  saveRDS(shape, file = rds)
}

#===========================================================
#2. Calculate all fasta file listed in current directory
#===========================================================
files <- rev(list.files("./",pattern = ".*fa$"))
prefixs <- sub("\\.fa","",files)

for (sample_prefix in prefixs){
  fa <- paste0(sample_prefix,".fa")
  rds <- paste0(sample_prefix, "_DNAshape.rds")
  shape <- getShape(fa, parse = TRUE)
  saveRDS(shape, file = rds)
}

#===========================================================
#3. Generate mean shape value of columns
#===========================================================
tmp <- readRDS(rds)
shapes <- names(tmp)
long_shape <- c("HelT","Roll","Tilt","Rise","Shift","Slide")
short_shape <- setdiff(shapes, long_shape)

NUMC = ncol(tmp[[shapes[1]]])
NUMR = nrow(tmp[[shapes[1]]])

mean_shape <- matrix(0, nrow = length(shapes), ncol = NUMC ) %>% as.data.frame()
rownames(mean_shape) <- shapes

for (shape in shapes){
  if (shape %in% short_shape){
    mean_shape[shape,] <- colMeans(tmp[[shape]], na.rm = T)
  } else {
    mean_shape[shape,1:(NUMC-1)] <- colMeans(tmp[[shape]], na.rm = T)
  }
}

write.table(mean_shape, paste0(rds,"_mean_shape_col.txt"), col.names = F, quote = F)

#===========================================================
#4. Generate mean shape value of columns and rows
#===========================================================
files <- list.files("./",pattern = ".*DNAshape.rds$")
prefixs <- sub("_DNAshape\\.rds","",files)

shapes <- c("Rise","Shift","Slide","HelT","Roll","Tilt",
            "ProT","Stretch","Buckle","Shear","Opening","Stagger","MGW","EP")

res_mat <- as.data.frame(matrix(0, ncol = length(shapes), nrow = length(prefixs)))
rownames(res_mat) <- prefixs
colnames(res_mat) <- shapes

for (prefix in prefixs){
  cat("Processing ",prefix,"...\n")
  tmp <- readRDS(paste0(prefix,"_DNAshape.rds"))
  
  for(shape in shapes){
    res_mat[prefix,shape] <- mean(tmp[[shape]],na.rm = T)
  }
}
write.table(round(res_mat,4), "mean_shape.txt", quote = F, sep = "\t")
