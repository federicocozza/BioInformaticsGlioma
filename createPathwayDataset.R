# trainAliquot : genes X patients [12042, 124]
# tRNA : genes X patients [12042, 167]
# egmt_top30 : 30 pathway with highest p-value [30, 9]

load("indici_training.Rdata")
load("egmt_result.Rdata")
#pathways <- egmt_top30

#pathways <- egmt_result[which(egmt_result$pvalue < 0.05),]

pathways <- egmt_result

pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
    cur <- pathways[i,]
    genes <- strsplit(cur$geneID,"/")[[1]]
    if(length(genes) < 10){
       break()
    }
    data = t(tRNA[genes,indici_training])
    write.table(data,file=paste("Pathways_more_than_10/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(patient_classes[indici_training,]$x))
names(class) <- "x"
row.names(class) <- rownames(patient_classes[indici_training,])
write.table(class,file="Pathways_more_than_10/labels",quote=F)

pathway_accuracy <- read.table("ourPathways_res/test_avg_accuracy.txt", header = TRUE, sep = ",")

pathway_accuracy <- pathway_accuracy[order(pathway_accuracy$X),]

egmt_result_sorted <- egmt_result[order(egmt_result$ID),]
egmt_result_sorted2 <- egmt_result_sorted[-which(!(egmt_result_sorted$ID %in% pathway_accuracy$X)),]

pathway_accuracy$count <- egmt_result_sorted2$Count

save(pathway_accuracy,file = "pathway_accuracy.Rdata")

