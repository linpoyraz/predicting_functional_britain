# Lin Poyraz 
# Usage imputation_quality_by_gene.R <chr> <ensembl name> 

library(data.table)
library(RSQLite)

args = commandArgs(trailingOnly = T)
gene = args[2]
chr = args[1]

# read in weights database and extract
con =  dbConnect(RSQLite::SQLite(), "./data/zhou2020_JTI_models/Best_Models_ProteinCoding.db")
weights = dbGetQuery(con, "SELECT * FROM weights")
i = which(weights$gene == gene)
weights = weights[i,]
total = nrow(weights) 

# imputation quality: extract quality pertaining to gene 
imp = read.table(paste("./data/imputation_quality/chr", chr, ".txt.gz", sep = ""), header = T, sep = "\t")
idx = which(imp$RSID %in% weights$rsid)
imp = imp[idx,]
with_id = nrow(imp)
genotyped = sum(imp$Genotyped == "Genotyped")

# parse through imp file 
score = 0 
sum = 0 
for (j in 1:nrow(imp)){
	# find on table
	weight = abs(weights$weight[which(weights$rsid == imp$RSID[j])])
	score = score + weight * imp$Rsq[j]
	sum = sum + weight	
}
# weighted average imputation  quality 
mean = sum(score)/sum(sum)

results = paste(mean, total, with_id, genotyped)
print(results)

