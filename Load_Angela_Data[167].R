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

### DFP

class1PatientsRNA <- patientsID[patientsID$class == "Classical",]
class2PatientsRNA <- patientsID[patientsID$class == "Mesenchymal",]
class3PatientsRNA <- patientsID[patientsID$class == "Neural",]
class4PatientsRNA <- patientsID[patientsID$class == "Proneural",]

## 75% of the sample size

smp_size_1 <- floor(0.75 * nrow(class1PatientsRNA))
smp_size_2 <- floor(0.75 * nrow(class2PatientsRNA))
smp_size_3 <- floor(0.75 * nrow(class3PatientsRNA))
smp_size_4 <- floor(0.75 * nrow(class4PatientsRNA))

## set the seed to make your partition reproductible

set.seed(123)
train_ind_1 <- sample(seq_len(nrow(class1PatientsRNA)), size = smp_size_1)
train_ind_2 <- sample(seq_len(nrow(class2PatientsRNA)), size = smp_size_2)
train_ind_3 <- sample(seq_len(nrow(class3PatientsRNA)), size = smp_size_3)
train_ind_4 <- sample(seq_len(nrow(class4PatientsRNA)), size = smp_size_4)

train_1 <- class1PatientsRNA[train_ind_1, ]
test_1 <- class1PatientsRNA[-train_ind_1, ]

train_2 <- class2PatientsRNA[train_ind_2, ]
test_2 <- class2PatientsRNA[-train_ind_2, ]

train_3 <- class3PatientsRNA[train_ind_3, ]
test_3 <- class3PatientsRNA[-train_ind_3, ]

train_4 <- class4PatientsRNA[train_ind_4, ]
test_4 <- class4PatientsRNA[-train_ind_4, ]

trainFinalRNA <- rbind(train_1, train_2, train_3, train_4)
testFinalRNA <- rbind(test_1, test_2, test_3, test_4)

trainAliquot <- tRNA[,trainFinalRNA$ID]
testAliquot <- tRNA[,testFinalRNA$ID]

row.names(trainFinalRNA) <- trainFinalRNA$ID
row.names(testFinalRNA) <- testFinalRNA$ID

save(trainAliquot,file = "trainAliquot_GLIOMA.Rdata")
save(testAliquot,file = "testAliquot_GLIOMA.Rdata")
save(trainFinalRNA,file = "trainFinalRNA_GLIOMA.Rdata")
save(testFinalRNA,file = "testFinalRNA_GLIOMA.Rdata")

load("trainAliquot_GLIOMA.Rdata")
load("testAliquot_GLIOMA.Rdata")
load("trainFinalRNA_GLIOMA.Rdata")
load("testFinalRNA_GLIOMA.Rdata")

multiDFP(trainAliquot, trainFinalRNA, "Glioma", core = 4, overlapping = c(1, 2), piVal = 0.9, z = c(0.35, 0.4, 0.45, 0.5), skipFactor = 2)
#################
multiDFP(trainAliquot, trainFinalRNA, "Glioma", core = 4, overlapping = c(1, 2), piVal = c(0.5, 0.6, 0.7, 0.8),  skipFactor = 1)


#### da rivedere

indici_training <- which(names(tRNA) %in% names(trainAliquot))

save(indici_training,file = "indici_training.Rdata")
