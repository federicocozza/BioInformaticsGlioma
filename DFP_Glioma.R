load("tRNA_final.RData")
load("tmiRNA_final.RData")
load("RNAByAliquotL.RData")

movetolast <- function(data, move) {
  data[c(setdiff(names(data), move), move)]
}

movetolast(tRNA_final, 'cases_0_diagnoses_0_tumor_stage')

stage1PatientsRNA <- tRNA_final[tRNA_final$cases_0_diagnoses_0_tumor_stage == "stage i",]
stage2PatientsRNA <- tRNA_final[tRNA_final$cases_0_diagnoses_0_tumor_stage == "stage ii",]
stage3PatientsRNA <- tRNA_final[tRNA_final$cases_0_diagnoses_0_tumor_stage == "stage iii",]

## 75% of the sample size
smp_size_1 <- floor(0.75 * nrow(stage1PatientsRNA))
smp_size_2 <- floor(0.75 * nrow(stage2PatientsRNA))
smp_size_3 <- floor(0.75 * nrow(stage3PatientsRNA))

## set the seed to make your partition reproductible
set.seed(123)
train_ind_1 <- sample(seq_len(nrow(stage1PatientsRNA)), size = smp_size_1)
train_ind_2 <- sample(seq_len(nrow(stage2PatientsRNA)), size = smp_size_2)
train_ind_3 <- sample(seq_len(nrow(stage3PatientsRNA)), size = smp_size_3)

train_1 <- stage1PatientsRNA[train_ind_1, ]
test_1 <- stage1PatientsRNA[-train_ind_1, ]

train_2 <- stage2PatientsRNA[train_ind_2, ]
test_2 <- stage2PatientsRNA[-train_ind_2, ]

train_3 <- stage3PatientsRNA[train_ind_3, ]
test_3 <- stage3PatientsRNA[-train_ind_3, ]

trainFinalRNA <- rbind(train_1, train_2, train_3)
testFinalRNA <- rbind(test_1, test_2, test_3)

trainAliquot <- RNAByAliquot[,trainFinalRNA$V14]
testAliquot <- RNAByAliquot[,testFinalRNA$V14]

row.names(trainFinalRNA) <- trainFinalRNA$V14
row.names(testFinalRNA) <- testFinalRNA$V14

multiDFP(trainAliquot, trainFinalRNA, "LungDFP", core = 4, overlapping = c(1, 2), piVal = c(0.4, 0.5, 0.6, 0.7), z = c(0.35, 0.4, 0.45, 0.5, 0.55), skipFactor = c(0, 1, 2, 3))
