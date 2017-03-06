library(ranger)

varGenes <- apply(log(RNAByAliquot + 8),1,var)
varGenes <- as.data.frame(varGenes)

varGenes$genes <- row.names(varGenes)
varGenes <- varGenes[order(varGenes$varGenes, decreasing=TRUE),]
plot(varGenes$varGenes, type='l')

selectedGenes <- varGenes[which(varGenes$varGenes > 0),]
rfTrain <- RNAByAliquot_ordered[,which(selectedGenes$genes %in% colnames(RNAByAliquot_ordered))]
rfTrain$stage <- RNAByAliquot_ordered$stage

rf <- ranger(stage ~ ., data = rfTrain, num.trees=5001, min.node.size = 1, importance = "impurity")
