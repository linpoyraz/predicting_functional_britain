# Lin Poyraz 
# Usage: Rscript iHS_scan_lite.R <chr> <ensembl name> 

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

# iHS scores
iHS_f = paste("/project/mathilab/ppoyraz/scan/iHS_selscan/results_AA/GBR_chr", chr, "_ihs.txt", sep = "") 
iHS = read.table(iHS_f, header = T, sep = "\t")
# extract SNPs relevant to this gene 
i = which(iHS$snp_id  %in% weights$rsid)
in_db = length(i)
snp_num = length(weights$rsid)
extracted_iHS = iHS[i,]
extracted_iHS$iHS = as.numeric(extracted_iHS$iHS)
# generate score 
score = c()
normalization = c()
sub = c() 
if (nrow(extracted_iHS) == 0){
	z = NA
} else {
	for (i in 1:nrow(extracted_iHS)){
		j = which(weights$rsid == extracted_iHS$snp_id[i])
		sub = c(sub,j) 
		# match effect alleles between weight and SDS
		if (extracted_iHS$a2[i] != weights$eff_allele[j]){
			weights$weight[j] = -weights$weight[j]
			weights$eff_allele[j] = extracted_iHS$a2[i] 
			weights$ref_allele[j] = extracted_iHS$a1[i]
		}
		score = c(score, weights$weight[j] * extracted_iHS$iHS[i]) 
		normalization = c(normalization, weights$weight[j]^2) 
	}
	z = sum(score) / sqrt(sum(normalization)) 
}
results = paste(z, snp_num, in_db)
print(results)

