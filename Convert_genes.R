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
RNAByAliquot$genes <- row.names(RNAByAliquot)
RNAByAliquot <- aggregate(. ~ genes,FUN = median, data=RNAByAliquot)

save(RNAByAliquot,file = "RNAByAliquotL.Rdata")

## Changed Barcode with patient
RNASeqData <- colData(RNATable)
names(RNAByAliquot) <- c("genes",RNASeqData$patient)
row.names(RNAByAliquot) <- RNAByAliquot$genes
RNAByAliquot <- RNAByAliquot[,-1]

## Filtering out patients with not reported tumor stage 
