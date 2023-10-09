# Lin Poyraz
# Usage: Rscript iHS_split.R
# Adds rsid field to iHS output files. Formats iHS output files 

library(data.table) 

f = paste("/project/mathilab/ppoyraz/scan/iHS_selscan/results_AA/GBR_chr", 1, ".ihs.out.100bins.norm", sep = "") 
file = read.table(f, header = F, sep = "\t", comment.char = "")
print(head(file))
file = file[,c(2,7)]
print(head(file))
names(file) = c("pos", "iHS") 
file$chr = rep(1, nrow(file))
all = file

for (i in 1:22){
	print(i)
        # read iHS results 
	f = paste("/project/mathilab/ppoyraz/scan/iHS_selscan/results_AA/GBR_chr", i, ".ihs.out.100bins.norm", sep = "") 
	file = read.table(f, header = F, sep = "\t", comment.char = "")
	file = file[,c(2,7)]
	names(file) = c("pos", "iHS") 
	print(head(file))
	file$chr = rep(i, nrow(file))	
	# read jti dosages to extract ancestral/derived information. assume ref (a1) is ancestral and alt (a2) is derived. 
	jti_f = paste("/project/mathilab/ppoyraz/dosage/dosage_files/chr", i, ".rs_updated.dos.txt.gz", sep = "") 
	jti = read.table(jti_f, header = T, sep = "\t", comment.char = "")
        jti = subset(jti, select = c("snp_id", "pos", "a1", "a2")) 	
	if (i != 1){
		all = rbind(all, file)
	}
	# match on pos. new file will have a1 and a2 information, as well as rsids. 
	file = merge(file, jti, by = "pos") 
        f_new = paste("/project/mathilab/ppoyraz/scan/iHS_selscan/results_AA/GBR_chr", i, "_ihs.txt", sep = "")
	fwrite(file, f_new, col.names = T, row.names = F, na = NA, quote = F, sep = "\t") 
	print("file")
	print(nrow(file))
	print("in databases") 
	print(nrow(jti))
}
f_all = "/project/mathilab/ppoyraz/scan/iHS_AA_GBR_all.txt"
fwrite(all, f_all, col.names = T, row.names = F, na = NA, quote = F, sep = "\t") 

