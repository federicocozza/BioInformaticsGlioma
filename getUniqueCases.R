tRNA <- read.table("all_rna.txt", header = TRUE, sep = "\t")
tmiRNA <- read.table("all_mirna.txt", header = TRUE, sep = "\t")

tmerged <- intersect(unique(tRNA$cases_0_case_id), unique(tmiRNA$cases_0_case_id))
tmerged

miRNA 34464
RNA 46329

totale 80793

intersezione 10140