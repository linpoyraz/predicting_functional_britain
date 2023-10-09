# Lin Poyraz 
# Usage: Rscript allele_frequencies.R <FLAG:"ensembl" or "gene" name> <name>

library(data.table)
library(RSQLite)

options(stringsAsFactors = FALSE) 

args = commandArgs(trailingOnly = T)

# read in gene name and chromosome information  
matched = read.table("/project/mathilab/colbranl/data/gene_lists/gencode.v26.GRCh38.all_genes.bed.gz", header = F, sep = "\t")[, c(4,5, 1)]
names(matched) = c("gene", "ensembl", "chr") 
matched[,2] = gsub("\\..*", "", matched[,2])

# read in weights database and extract 
con =  dbConnect(RSQLite::SQLite(), "/project/mathilab/colbranl/data/zhou2020_JTI_models/Best_Models_ProteinCoding.db")
#con =  dbConnect(RSQLite::SQLite(), "/project/mathilab/colbranl/data/zhou2020_JTI_models/JTI_Brain_Caudate_basal_ganglia.db")
# get ensembl gene name 
if (args[1] == "ensembl"){
	gene = args[2]
	i = which(matched$ensembl == args[2])
	human = matched$gene[i]
	chr = matched$chr[i]

} else{
	i = which(matched$gene == args[2])
	gene = matched$ensembl[i]
	human = args[2]
	chr = matched$chr[i]
}

# extract weights pertaining to specified gene from db
weights = dbGetQuery(con, "SELECT * FROM weights")
i = which(weights$gene == gene)
weights = weights[i,]

# print tissue 
extra = dbGetQuery(con, "SELECT * FROM extra")
i = which(extra$gene == gene)
extra = extra[i,]

# read in dosage information 
exp = read.table(paste("/project/mathilab/ppoyraz/dosage_updated/dosage_files/", chr, ".rs_updated.dos.txt.gz", sep = ""), header = T, sep = "\t", comment.char = "")
names(exp)[c(1,2,4,5)] = c("chr", "rsid", "a1", "a2") 

# extract snps pertaining to specified gene from dosage 
i = which(exp$rsid %in% weights$rsid)
print(weights$rsid)
exp = exp[i,]
row.names(exp) = exp$rsid

#print(length(exp$rsid))
#print(exp$rsid)
#print(weights) 
print("weights length") 
print(nrow(weights))
print("exp length") 
print(nrow(exp))

sub = c() 
# make sure a2 and effect allele match 
for (i in 1:nrow(exp)){
	j = which(weights$rsid == exp$rsid[i])
	sub = c(sub, j)
	if (weights$eff_allele[j] != exp$a2[i]){
		print(exp[i,1:5])
		print(weights[j,])
		weights$weight[j] = -weights$weight[j]
		weights$eff_allele[j] = exp$a2[i] 
		weights$ref_allele[j] = exp$a1[i]
	}
}
weights = weights[sub,]
#print(weights$weight)
iid = read.table("/project/mathilab/ppoyraz/dosage_updated/dosage_files/sample.txt", header = F, sep = "\t")
exp = exp[,c(7:ncol(exp))]
gene_dosages = as.data.frame(t(exp))
gene_dosages$iid = iid[,1]
gene_dosages$time = rep(0, nrow(gene_dosages))

#print(nrow(gene_dosages))
#print(ncol(gene_dosages))

# read in individual and time information 
meta = read.table("/project/mathilab/data/brit_data2/brit_meta.txt", header = F, sep = "\t")[,c(1,2)]
names(meta) = c("iid", "time")
meta$time = -as.numeric(meta$time)

# read in samples to subset to 
subset = read.table("/project/mathilab/ppoyraz/time_series_predixcan/analysis/expression_time_series/samples_01", header = F)
print(head(subset))
print("subset sample number")
print(nrow(subset))
# subset to samples in subset file 
i = which(!gene_dosages$iid %in% subset[,1])
print("not in subset") 
k = subset[which(!subset[,1] %in% gene_dosages$iid),1]
print(k)
print("subsetted to sample number")
gene_dosages = gene_dosages[i,]
print("should me this many modern individuals") 
print(sum(grepl(gene_dosages$iid, pattern = "^HG0")))
remove = c()
# add time information for each individual 
for (i in 1:nrow(gene_dosages)){
	if (gene_dosages$iid[i] %in% meta$iid){	
		idx = which(meta$iid == gene_dosages$iid[i])
		gene_dosages$time[i] = meta$time[idx]
	} else if (!grepl(gene_dosages$iid[i], pattern = "^HG0")) {
		remove = c(remove, gene_dosages$iid[i])
		print(gene_dosages$iid[i])
	}  
}

gene_dosages = gene_dosages[which(!gene_dosages$iid %in% remove),]
print("number of ancient individuals")
print(sum(gene_dosages$time != 0))
print("number of modern individuals")
print(sum(gene_dosages$time == 0))

print(nrow(meta))
#print(weights) 
f = paste("/project/mathilab/ppoyraz/time_series_predixcan/results_updated/data/allele_frequency/", human, ".txt", sep = "") 
fwrite(gene_dosages, f, row.names = F, col.names = T, quote = F, sep = "\t")
f = paste("/project/mathilab/ppoyraz/time_series_predixcan/results_updated/data/allele_frequency/", human, "_weights.txt", sep = "") 
fwrite(weights, f, row.names = F, col.names = T, quote = F, sep = "\t")




