# trainAliquot : genes X patients [12042, 124]
# tRNA : genes X patients [12042, 167]
# egmt_top30 : 30 pathway with highest p-value [30, 9]

load("indici_training.Rdata")
pathways <- egmt_top30

pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
    cur <- pathways[i,]
    genes <- strsplit(cur$geneID,"/")[[1]]
    data = t(tRNA[genes,indici_training])
    write.table(data,file=paste("ourPathways/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(patient_classes[indici_training,]$x))
names(class) <- "x"
row.names(class) <- rownames(patient_classes[indici_training,])
write.table(class,file="ourPathways/labels",quote=F)