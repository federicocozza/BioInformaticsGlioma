tRNA <- read.table("all_rna.txt", header = TRUE, sep = "\t")
tmiRNA <- read.table("all_mirna.txt", header = TRUE, sep = "\t")


tmerged <- intersect(unique(tRNA$cases_0_case_id), unique(tmiRNA$cases_0_case_id))
tmerged

##tmiRNA with only mirnas.quantification.txt
tmiRNA <- tmiRNA[which(tmiRNA$file_name == "mirnas.quantification.txt"),]

##length(which(tmiRNA$cases_0_samples_0_sample_type == "Primary Tumor"))

#removed Metastatic
tmiRNA <- tmiRNA[-which(tmiRNA$cases_0_samples_0_sample_type == "Metastatic"),]

#tRNA filtered on FPKM-UQ.txt.gz
tRNA <- tRNA[which(endsWith(as.character(tRNA$file_name) ,".FPKM-UQ.txt.gz") == TRUE),]

#removed Metastatic
tRNA <- tRNA[-which(tRNA$cases_0_samples_0_sample_type == "Metastatic"),]

tmerged <- intersect(tRNA$cases_0_case_id,tmiRNA$cases_0_case_id)
tmerged <- as.data.frame(tmerged)

tRNA2 <- which(tRNA$cases_0_case_id %in% tmiRNA$cases_0_case_id)
