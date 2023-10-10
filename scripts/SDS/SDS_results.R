# Lin Poyraz 
# Usage: Rscript SDS_results.R <regression results> 

library(qqman)
library(data.table)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)
regression_results = read.table(args[1], sep = "\t", header = F)
names(regression_results) = c("chr", "start","end","Gene","Ensembl","P.value","FDR","Bonferroni","GC","Beta") 
genes = regression_results$Ensembl
chrs = gsub("chr", "", regression_results$chr)

regression_results$SDS = rep(1000, nrow(regression_results))
regression_results$SDS.p = rep(1000, nrow(regression_results))
regression_results$SNP.total = rep(1000, nrow(regression_results))
regression_results$SNP.SDS = rep(1000, nrow(regression_results))

for (i in 1:length(genes)){
	print(genes[i]) 
	com = paste("Rscript SDS_scan_lite.R", chrs[i], genes[i], sep = " ")
        score = system(command = com, intern = T)
	score = gsub("\\[1\\] ", "", score)
	score = gsub('"', "", score)
	score = strsplit(score, " ")[[1]]
	score = as.numeric(score)
	print(score)
	regression_results$SDS[i] = score[1]
	# convert to p value 
	p = 2*pnorm(q = abs(score[1]), lower.tail=FALSE)
	regression_results$SDS.p[i] = p
	# add snp information 
	regression_results$SNP.total[i] = score[2]
	regression_results$SNP.SDS[i] = score[3]

}
f_results = paste("/project/mathilab/ppoyraz/time_series_predixcan/results/data/scan/SDS_split/", chrs[1], "_results.txt", sep = "")
fwrite(regression_results, f_results, quote = F, row.names = F, col.names = T, sep = "\t", na = NA)


