tRNA <- read.table("Dataset//all_genes.txt", header = TRUE, sep = "\t")
tmiRNA <- read.table("Dataset/all_mirnas.txt", header = TRUE, sep = "\t")
patientsID <- read.table("Dataset/patientsID.txt", header = TRUE, sep = "\t")

ORIG_tRNA <- tRNA

tRNA <- tRNA[,which(names(tRNA) %in% patientsID$x)]

patient_classes <- read.table("Dataset/patient_classes.txt", header = TRUE, sep = "\t")

patientsID[,2] <- patient_classes[,1]

colnames(patientsID) <- c('ID', 'class')
names(patientsID) <- patientsID$ID
rownames(patientsID) <- patientsID$ID