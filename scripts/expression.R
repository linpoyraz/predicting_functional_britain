# Lin Poyraz 
# Usage: expression.R 

library(data.table)
args = commandArgs(trailingOnly = T) 

print("reading in data...") 

meta = read.table("./data/brit_meta.txt", header = F, sep = "\t")[,c(1,2,4,5)]
names(meta) = c("iid", "time", "lat", "long")
meta$time = -as.numeric(meta$time)
print(head(meta)) 

exp = read.table(paste("./data/predicted_expression.txt.gz", sep = ""), header = T, sep = "\t")
exp = exp[,2:ncol(exp)]
names(exp)[1] = "iid"
exp$time = rep(0, nrow(exp))

print ("adding time information...") 

# this sample is a problem we need to fix: 
problem = which(exp$iid == "I5473_published")
exp$iid[problem] = "I5473"

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

dir = "./output/expression/"


genes = read.table("./data/fdr_sig.txt", header = T, sep = "\t")

for (i in 1:nrow(genes)){
	gene = genes$Ensembl[i]
	dat = exp[, c(gene, "time", "iid")]
	names(dat) = c("gene", "time", "iid")
	gene = genes$Gene[i]
	f = paste(dir, gene, ".txt", sep = "") 
	print(head(dat))
	fwrite(dat, f, row.names = F, quote = F, col.names = T, na = NA, sep = "\t")
}

