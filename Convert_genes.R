library(org.Hs.eg.db)
library(SummarizedExperiment)

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }

  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

load("RNATableL.Rdata")
#load("RNAByAliquotL.Rdata")
ensembles <- rownames(assay(RNATable,1))

entrez <- convertIDs(ensembles, "ENSEMBL", "ENTREZID", org.Hs.eg.db, ifMultiple = "useFirst")

## Number of genes
print(length(unique(entrez)))

## Only entrez genes with value
filtered <- assay(RNATable,1)[- which(is.na(entrez)),]

## Change rownames ensemble to entrezID
rownames(filtered) <- entrez[- which(is.na(entrez))]

## Aggregate by entrezID
RNAByAliquot <- as.data.frame(filtered)
##----------------------------------------------
RNAByAliquot$genes <- row.names(RNAByAliquot)
RNAByAliquot <- aggregate(. ~ genes,FUN = median, data=RNAByAliquot)

#save(RNAByAliquot,file = "RNAByAliquotL.Rdata")
save(RNAByAliquot,file = "RNAByAliquot_510.Rdata")
load("RNAByAliquot_510.Rdata")

## Changed Barcode with patient
RNASeqData <- colData(RNATable)

### sp
RNASeqData <- RNASeqData[which(rownames(RNASeqData) %in% colnames(RNAByAliquot)),]

names(RNAByAliquot) <- c("genes",RNASeqData$patient)
row.names(RNAByAliquot) <- RNAByAliquot$genes
RNAByAliquot <- RNAByAliquot[,-1]

### sp
tRNA_final$V14 <- colnames(RNAByAliquot)

## Filtering out patients with not reported tumor stage 
RNAByAliquot <- RNAByAliquot[,-which(RNASeqData$tumor_stage == "not reported")]

## sp
tRNA_final <- tRNA_final[-which(RNASeqData$tumor_stage == "not reported"),]

stages <- RNASeqData$tumor_stage[-which(RNASeqData$tumor_stage == "not reported")]

## Calculating index for tumor stages and removing patients at stage 4
stage4 <- grep("stage iv$", stages, perl = TRUE, value = FALSE)
RNAByAliquot <- RNAByAliquot[,-stage4]
stages <- stages[-stage4]

stage1 <- grep("stage i[^iv]*$", stages, perl = TRUE, value = FALSE)
stage2 <- grep("stage ii[^i]*$", stages, perl = TRUE, value = FALSE)
stage3 <- grep("stage iii[ab]*$", stages, perl = TRUE, value = FALSE)

## Grouping subclasses for tumor stages
stages[stage1] <- "stage i"
stages[stage2] <- "stage ii"
stages[stage3] <- "stage iii"

## sp
tRNA_final <- tRNA_final[-stage4,]

save(RNAByAliquot,file = "RNAByAliquotL.Rdata")

tRNA_final$cases_0_diagnoses_0_tumor_stage <- stages

save(tRNA_final,file = "tRNA_final.Rdata")

length(which(nchar(colnames(RNAByAliquot)) ==12))
length(unique(colnames(RNAByAliquot)))