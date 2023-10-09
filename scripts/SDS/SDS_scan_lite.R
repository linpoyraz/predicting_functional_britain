# Lin Poyraz 
# Usage SDS_scan_lite.R <chr> <ensembl name> 

library(data.table)
library(RSQLite)

args = commandArgs(trailingOnly = T)

# read in weights database and extract
con =  dbConnect(RSQLite::SQLite(), "/project/mathilab/colbranl/data/zhou2020_JTI_models/Best_Models_ProteinCoding.db")
gene = args[2]
chr = args[1]

weights = dbGetQuery(con, "SELECT * FROM weights")
i = which(weights$gene == gene)
weights = weights[i,]

# SDS scores
SDS_f = paste("/project/mathilab/ppoyraz/scan/SDS/SDS_", chr, ".tab.gz", sep = "") 
SDS = read.table(SDS_f, header = F, sep = "\t")
# extract SNPs relevant to this gene 
i = which(SDS[,3] %in% weights$rsid)
extracted_SDS = SDS[i,]
in_db = length(i)
snp_num = length(weights$rsid)
# generate score 
score = c()
normalization = c()
sub = c() 
if (nrow(extracted_SDS) == 0){
	z = NA
} else {
	for (i in 1:nrow(extracted_SDS)){
		j = which(weights$rsid == extracted_SDS[i,3])
		sub = c(sub, j) 
		# match effect alleles between weight and SDS
		if (extracted_SDS[i,5] != weights$eff_allele[j]){
			weights$weight[j] = -weights$weight[j]
		}
		score = c(score, weights$weight[j] * extracted_SDS[i,7]) 
		normalization = c(normalization, weights$weight[j]^2)
	}
	z = sum(score) / sqrt(sum(normalization))
}
results = paste(z, snp_num, in_db)
print(results)

