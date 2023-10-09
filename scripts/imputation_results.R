# Lin Poyraz 
# Usage: Rscript imputation_results.R <regression results> 
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
options(stringsAsFactors = FALSE)

regression_results = read.table(args[1], sep = "\t", header = F)
names(regression_results) = c("chr", "start","end","Gene","Ensembl","P.value","FDR","Bonferroni","GC","Beta") 
genes = regression_results$Ensembl
chrs = gsub("chr", "", regression_results$chr)

regression_results$Mean.Imputation.Rsq = rep(NA, nrow(regression_results))
regression_results$SNP.total = rep(NA, nrow(regression_results))
regression_results$SNP.quality = rep(NA, nrow(regression_results))
regression_results$SNP.Genotyped = rep(NA, nrow(regression_results))

for (i in 1:length(genes)){
	print(genes[i]) 
	com = paste("Rscript imputation_quality_by_gene.R", chrs[i], genes[i], sep = " ")
        score = system(command = com, intern = T)
	score = gsub("\\[1\\] ", "", score)
	score = gsub('"', "", score)
	score = strsplit(score, " ")[[1]]
	score = as.numeric(score)
	regression_results$Mean.Imputation.Rsq[i] = score[1]
	# add snp information 
	regression_results$SNP.total[i] = score[2]
	regression_results$SNP.quality[i] = score[3]
	regression_results$SNP.Genotyped[i] = score[4]

}
f_results = paste("./output/imputation_split/", chrs[1], "_results.txt", sep = "")
fwrite(regression_results, f_results, quote = F, row.names = F, col.names = T, sep = "\t", na = NA)


