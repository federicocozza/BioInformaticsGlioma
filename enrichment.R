library(GSA)
#gmtfile <- system.file("extdata", "msigdb.v5.2.entrez.gmt",package="clusterProfiler")
gmtfile <- system.file("extdata", "msigdb.v5.2.symbols.gmt",package="clusterProfiler")
c2kegg <- read.gmt(gmtfile)

#sum(gene %in% c2kegg$gene)

### Glioma - Gene
load("D:/R_workspace/BioInformaticsGlioma/paramList10.Rdata")
gene <- rownames(paramList$dfps)

egmt <- enricher(gene, TERM2GENE=c2kegg,pvalueCutoff = 0.05)
#head(egmt)

#View(egmt@result[grep(pattern = "KEGG",x = rownames(egmt@result),ignore.case = F),])

reactome <- egmt@result[grep(pattern = "reactome",x = rownames(egmt@result),ignore.case = T),]
kegg <- egmt@result[grep(pattern = "kegg",x = rownames(egmt@result),ignore.case = T),]
biocarta <- egmt@result[grep(pattern = "biocarta",x = rownames(egmt@result),ignore.case = T),]

# in c2
# estrarre i primi 30 pathway con p-value più alto
# creare file per ogni pathway con i soli geni dei pazienti

dim(egmt@result)
View(egmt@result)
plot(egmt@result$pvalue)
c2 <- read.gmt('D:\\Download\\c2.all.v5.2.symbols.gmt')
"VERHAAK_GLIOBLASTOMA_MESENCHYMAL" %in% c2$ont

egmt_result <- egmt@result
save(egmt_result, file="egmt_result.Rdata")

# reactomeFinal <- reactome[which(reactome$Count >=10),]
# keggFinal <- kegg[which(kegg$Count >=10),]
# c6Final <- egmt@result[egmt@result$ID %in% c6 & egmt@result$Count >= 10,]

save(reactomeFinal,file="reactomeFinal.Rdata")
save(keggFinal,file="keggFinal.Rdata")
save(c6Final,file="c6Final.Rdata")