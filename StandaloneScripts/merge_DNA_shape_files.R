#This script is for merging DNA shapes from DiProDB 
shape_files <- list.files(path = ".", pattern = "fa.",full.names = F, recursive = F)
samples <- unique(gsub(".fa.*","",shape_files))
full_sample <- paste0(samples,".fa.")
shapes <- unique(gsub(full_sample[1],"",shape_files[grep(samples[1],shape_files)]))

for (ss in full_sample){
  tmp <- c()
  for (shape in shapes){
    a <- readLines(paste0(ss,shape))
    value <- numeric(length=length(a)/2)
    for (i in 1:length(a)){
      if (i%%2 != 1){value[i/2] <- round(mean(as.numeric(strsplit(a[i],",")[[1]][-1])),2)}
    }
    value <- as.data.frame(value)
    colnames(value)[1] <- sub("\\(DNA-protein_complex\\)","DBP",gsub(" ","_",shape))
    if(!is.data.frame(tmp)){tmp <- rbind(tmp,value)
    } else {tmp <- cbind(tmp,value)}
  }
  names <- character(length=length(a)/2)
  for (i in 1:length(a)){
    if (i%%2 == 1){names[(i+1)/2] <- a[i]}
  }
  rownames(tmp) <- names
  assign(ss,tmp)
  out <- eval(parse(text = ss))
  write.table(out,paste0(ss,"full_shapes.txt"),sep = "\t",quote = F)
}





