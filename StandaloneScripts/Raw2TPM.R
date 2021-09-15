#This script is for converting Gene IDs from Ensembl to RefSeq

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(max.print = 100)
#=========================================================================================
#1. Convert rawcout to TPM
#=========================================================================================
RawCount2TPM <- function(Rawcount=Rawcount, Reference=Reference, result="TPM.txt") {
  Raw <- read.table(Rawcount, row.names = 1)
  Ref <- read.table(Reference)
  Ref_sorted <- Ref[match(rownames(Raw),Ref$V4),]
  combined <- cbind(Raw,Ref_sorted)
  colnames(combined) <- c("Rawcount","Chr","Start","End","ENSEM","NONE","Strand","Name","GeneType")
  combined_filtered <- combined[!is.na(combined$ENSEM),]
  filterNUM <- nrow(combined) - nrow(combined_filtered)
  cat(filterNUM,"genes are filtered because not mapped in your reference!!!\n")
  
  combined_filtered$GeneLength <- combined_filtered$End - combined_filtered$Start
  combined_filtered$ss <- combined_filtered$Rawcount / combined_filtered$GeneLength
  combined_filtered$TPM <- combined_filtered$ss * 1e6 / sum(combined_filtered$ss, na.rm = T)
  #return(Rawcount)
  write.table(combined_filtered[,"TPM", drop=F], result, col.names = F, row.names = T, sep = "\t", quote = F)

}

Rawcount="your_Rawcount.txt" #ENSGXXX\tnnnn
Reference="mm10.EMSEMBLE_54838.bed"
RawCount2TPM(Rawcount, Reference, paste0(sub("_Rawcount.txt","",Rawcount),"_TPM_ENSEMBL.txt"))

#=========================================================================================
#2. Convert ENSEMBL gene expression to Refseq (to TSS lists)
#=========================================================================================
ENSG_TPM <- "your_TPM_ENSEMBL.txt"
TSS_list <- "mm10.Refseq_TSS.bed"
martfiles <- "Mouse_mart_export.txt"
species <- c("Mouse","Human")[1]

if (file.exists(martfiles)){
  ENSG_REF <- read.table(martfiles, sep = "\t")
} else {
  library(biomaRt)
  mymart <- useEnsembl(biomart = "ensembl")
  #searchDatasets(mart = mymart, pattern = "mmu*")
  if (species == "Mouse"){dataset = "mmusculus_gene_ensembl"
  } else {dataset = "hsapiens_gene_ensembl"}
  Myensembl = useMart(biomart="ensembl", dataset = dataset)
  #listDatasets(Myensembl)
  #searchAttributes(mart = Myensembl, pattern = "refseq")
  ENSG_REF <- getBM(attributes = c("ensembl_gene_id","refseq_mrna","refseq_ncrna"), mart = Myensembl)
  write.table(ENSG_REF,martfiles,quote = F, col.names = F, row.names = F)
}

TSS <- read.table(TSS_list)
TSS_mapped <- ENSG_REF[(ENSG_REF$V2 %in% TSS$V4) | (ENSG_REF$V3 %in% TSS$V4),]
TSS_mapped$refseq <- TSS_mapped$V2
TSS_mapped$refseq[TSS_mapped$refseq==""] <- TSS_mapped$V3[TSS_mapped$refseq==""]

expr <- read.table(ENSG_TPM, header = F,sep = "\t", row.names = 1)
sum(!rownames(expr) %in% ENSG_REF$ensembl_gene_id)

Ref_expr <- data.frame(RefseqID = TSS$V4, geneExpr=NA)

for (i in 1:nrow(Ref_expr)){
  rID <- Ref_expr$RefseqID[i]
  if (rID %in% TSS_mapped$refseq){
    eID <- TSS_mapped$V1[grep(rID,TSS_mapped$refseq)]
    if(eID %in% rownames(expr)){
      value <- expr[eID,]
    }
  } else {value = NA}
  
  Ref_expr$geneExpr[i] <- value
}

Refseq_TPM <- "your_TPM_Refseq.txt"
write.table(Ref_expr, file = Refseq_TPM, row.names = F, col.names = F, quote = F, sep = "\t")
