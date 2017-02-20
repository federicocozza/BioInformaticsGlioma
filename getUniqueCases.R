tRNA <- read.table("all_rna.txt", header = TRUE, sep = "\t")
tmiRNA <- read.table("all_mirna.txt", header = TRUE, sep = "\t")

##tmiRNA with only mirnas.quantification.txt
tmiRNA <- tmiRNA[which(tmiRNA$file_name == "mirnas.quantification.txt"),]

##length(which(tmiRNA$cases_0_samples_0_sample_type == "Primary Tumor"))

#removed Metastatic
#tmiRNA <- tmiRNA[-which(tmiRNA$cases_0_samples_0_sample_type == "Metastatic"),]

#tRNA filtered on FPKM-UQ.txt.gz
tRNA <- tRNA[which(endsWith(as.character(tRNA$file_name) ,".FPKM-UQ.txt.gz") == TRUE),]

# Now we have 594 obs and we add a column = col names of RNAByAliquot

tRNA[14] <- colnames(RNAByAliquot)

#removed Metastatic
#tRNA <- tRNA[-which(tRNA$cases_0_samples_0_sample_type == "Metastatic"),]

#get common case_id
tmerged <- intersect(tRNA$cases_0_case_id,tmiRNA$cases_0_case_id)

# get just the miRNA with common case_id
l <- c()
tmiRNA2 = tmiRNA
tmiRNA_final <- tmiRNA2[0,]
z <- length(tmiRNA2$cases_0_case_id)
tmerged2 <- tmerged
for(i in 1:z){
    print(i)
    if((tmiRNA2$cases_0_case_id[i] %in% tmerged2) == FALSE){
        #tmiRNA_final <- tmiRNa
        l<-c(l,i)
        #tmiRNA_final<- rbind(tmiRNA_final, tmiRNA_tw[i,])
    }
    else{
        toBeRemoved<-which(tmerged2 == toString(tmiRNA2$cases_0_case_id[i]))
        tmerged2<-tmerged2[-toBeRemoved]
    }
}

tmiRNA_final <- tmiRNA2[-l,]

# get just the RNA with common case_id
l <- c()
tRNA2 = tRNA
tRNA_final <- tRNA2[0,]
z <- length(tRNA2$cases_0_case_id)
tmerged2 <- tmerged
for(i in 1:z){
    print(i)
    if((tRNA2$cases_0_case_id[i] %in% tmerged2) == FALSE){
        #tRNA_final <- tRNA
        l<-c(l,i)
        #tRNA_final<- rbind(tRNA_final, tRNA_tw[i,])
    }
    else{
        toBeRemoved<-which(tmerged2 == toString(tRNA2$cases_0_case_id[i]))
        tmerged2<-tmerged2[-toBeRemoved]
    }
}

tRNA_final <- tRNA2[-l,]

save(tRNA_final, file = "tRNA_final.RData")
save(tmiRNA_final, file = "tmiRNA_final.RData")