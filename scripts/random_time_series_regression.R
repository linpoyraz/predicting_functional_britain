# Lin Poyraz 
# Usage: random_time_series_regression <expression file> 


library(data.table)
options(stringsAsFactors = FALSE)
set.seed(444)

args = commandArgs(trailingOnly=TRUE)

print("reading in data...") 

meta = read.table("./data/brit_meta.txt", header = F, sep = "\t")[,c(1,2,4,5)]
names(meta) = c("iid", "time", "lat", "long")
meta$time = -as.numeric(meta$time)
print(head(meta)) 

# change gene names to human-identifiable names 
names = read.table("./data/gencode.v26.GRCh38.all_genes.bed.gz", header = F, sep = "\t")[, c(4,5, 1, 2, 3)]
names[,2] = as.character(gsub("\\..*", "", names[,2]))

exp = read.table(paste(args[1], "/predicted_expression.txt.gz", sep = ""), header = T, sep = "\t")
exp = exp[,2:ncol(exp)]
names(exp)[1] = "iid"
exp$time = rep(0, nrow(exp))

print ("adding time information...") 

# this sample is a problem we need to fix: 
problem = which(exp$iid == "I5473_published")
exp$iid[problem] = "I5473"

#fwrite(list(exp$iid[which(!grepl(exp$iid, pattern = "^HG0|_published"))]), "samples_all", quote = F, row.names = F, col.names = F, sep = "\t", na = NA)

# read in samples to subset to 
subset = read.table("./data/samples_01", header = F)
print(head(subset))
print("subset sample number")
print(nrow(subset))
# subset to samples in subset file 
i = which(!exp$iid %in% subset[,1])
print("not in subset") 
k = subset[which(!subset[,1] %in% exp$iid),1]
print(k)
print("subsetted to sample number")
exp = exp[i,]
print("should me this many modern individuals") 
print(sum(grepl(exp$iid, pattern = "^HG0")))
remove = c()
# add time information for each individual to exp 
for (i in 1:nrow(exp)){
	if (exp$iid[i] %in% meta$iid){	
		idx = which(meta$iid == exp$iid[i])
		exp$time[i] = meta$time[idx]
	} else if (!grepl(exp$iid[i], pattern = "^HG0")) {
		remove = c(remove, exp$iid[i])
		print(exp$iid[i])
	}  
}

exp = exp[which(!exp$iid %in% remove),]
print("number of ancient individuals")
print(sum(exp$time != 0))
print("number of modern individuals")
print(sum(exp$time == 0))
print("running linear regression, updating gene names...") 


print(sum(exp$time != 0))

# randomize time
exp$time = sample(exp$time)

print("running linear regression, updating gene names...") 

# run linear regression for each gene vs. time, store all p values and betas 
genes = names(exp)[which(!names(exp) %in% c("iid", "time"))]
results = data.frame("chr" = genes, "start" = rep(1000, length(genes)), "end" = rep(1000, length(genes)), 
	"Gene" = genes, "Ensembl" = genes, "P.value" = rep(1000, length(genes)), "FDR" = rep(1000, length(genes)), 
	"Bonferroni" = rep(1000, length(genes)), "GC" = rep(1000, length(genes)), "Beta" = rep(1000, length(genes)))

results$chr = as.character(results$chr)
results$Gene = as.character(results$Gene)
results$Ensembl = as.character(results$Ensembl)

for (i in 1:length(genes)){
	gene = genes[i]
	names_row = names[which(names[,2] == gene),]
	h = names_row[,1]
	results[i, "Gene"] = h
	dat = exp[, c(gene, "time")]
	colnames(exp)[i+1] = h
	names(dat) = c("gene", "time")
	model = lm(gene ~ time, data = dat)
	results$P.value[i] = summary(model)$coefficients["time", 4]
	results$Beta[i] = summary(model)$coefficients["time", 1]
	results$chr[i] = names_row[,3]
	results$start[i] = names_row[,4]
	results$end[i] = names_row[,5]
}

results = na.omit(results) ##NA omitted 
print(head(results))
print("adjusting p values and preparing output files...") 
p = results$P.value


# save results 
f = paste("./output/samples01_random_time_series_regression_results.txt", sep = "")
fwrite(results, f, quote = F, row.names = F, col.names = T, sep = "\t", na = NA)

