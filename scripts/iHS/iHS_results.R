# Lin Poyraz 
# Usage: Rscript iHS_results.R <regression results> 

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
regression_results = read.table(args[1], sep = "\t", header = F)
names(regression_results) = c("chr", "start","end","Gene","Ensembl","P.value","FDR","Bonferroni","GC","Beta") 
genes = regression_results$Ensembl
chrs = gsub("chr", "", regression_results$chr)

regression_results$iHS = rep(NA, nrow(regression_results))
regression_results$iHS.p = rep(NA, nrow(regression_results))
regression_results$SNP.total = rep(NA, nrow(regression_results))
regression_results$SNP.iHS = rep(NA, nrow(regression_results))

for (i in 1:length(genes)){
	com = paste("Rscript iHS_scan_lite.R", chrs[i], genes[i], sep = " ")
	score = system(command = com, intern = T)
	print(score) 
	score = gsub("\\[1\\] ", "", score)
	score = gsub('"', "", score)
	print(score) 
	score = strsplit(score, " ")[[1]]
	score = as.numeric(score)
	print(score)
	regression_results$iHS[i] = score[1]
	# convert to p value 
	p = 2*pnorm(q = abs(score[1]), lower.tail=FALSE)
	regression_results$iHS.p[i] = p
	# note snps 
	regression_results$SNP.total[i] = score[2]
	regression_results$SNP.iHS[i] = score[3]
}
print("saving") 
f_results = paste("./data/iHS_split/", chrs[1], "_results.txt", sep = "")
fwrite(regression_results, f_results, quote = F, row.names = F, col.names = T, sep = "\t", na = NA)


