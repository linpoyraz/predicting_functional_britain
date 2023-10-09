# Lin Poyraz 
# Usage: Rscript ancestral_match.R <lifted file> <unlifted file> <vcf to match>  

library(seqinr)
library(data.table)

args = commandArgs(trailingOnly = T) 

# read files: only first 5 columns of the vcf file is read 
vcf = fread(cmd = paste("gunzip -c", args[3], "| grep -v '^#'"), sep = "\t", select = 1:5)
names(vcf) = c("chr", "POS", "ID", "REF", "ALT") 
lifted = read.table(args[1], header = F, sep = "\t") 
names(lifted) = c("chr", "pos", "pos1", "id") 
unlifted = read.table(args[2], header = F, sep = "\t") 
names(unlifted) = c("chr", "pos", "pos1", "id") 

# initialize output ancestral allele file 
output = data.frame(CHR = rep(vcf$chr[1], nrow(vcf)), POS = vcf$POS, REF = vcf$REF, ALT = vcf$ALT, "RSID" = vcf$ID, "INFO/AA" = rep("F", nrow(vcf)))

# read chimpanzee fasta sequences: initially only read in the chromosome the first variant was lifted to 
read = read.fasta(paste("/project/mathilab/ppoyraz/lift/pantro6/", lifted$chr[1], ".fasta-1.gz", sep = ""))

# iterate over vcf 
lifted_i = 1 
unlifted_i = 1
unmatched = c()
for (i in 1:nrow(vcf)){
	if (unlifted_i > nrow(unlifted)){ 
		unlifted_i = nrow(unlifted)
	}
	if (vcf$POS[i] != unlifted$pos[unlifted_i]){ # check if the variant was lifted 
		pt_chr = lifted$chr[lifted_i]
		if (pt_chr %in% names(read)){ # check if the chromosome the variant was lifted to has been read in 
			aa = toupper(read[[pt_chr]][lifted$pos[lifted_i]])
			if(!(aa == vcf$REF[i] | aa == vcf$ALT[i])){ # ensure that the ancestral allele matches either ref or alt in the vcf 
				unmatched = c(unmatched, vcf$POS[i])
				aa = "."
			} else {
				output[i,"INFO/AA"] = aa 
			}
		} else {
			read = append(read, read.fasta(paste("/project/mathilab/ppoyraz/lift/pantro6/", pt_chr, ".fasta-1.gz", sep = ""))) # if the chromosome hasn't been read in yet 
			aa = toupper(read[[pt_chr]][lifted$pos[lifted_i]])
			if(!(aa == vcf$REF[i] | aa == vcf$ALT[i])){
				unmatched = c(unmatched, vcf$POS[i])
				aa = "."
			} else {
				output[i,"INFO/AA"] = aa 
			}

		}
		lifted_i = lifted_i + 1
	} else {
		unlifted_i = unlifted_i + 1 
		aa = "." 
	}
}

print(head(output))
f_out = paste("/project/mathilab/ppoyraz/scan/iHS_selscan/temp/", vcf$chr[1], "_AA.txt", sep = "") # ancestral allele file 
output_fil = output[which(output[,"INFO/AA"] %in% c("A", "C", "T", "G") & output[,"RSID"] != "."),] # remove unmatched variants and variants without rsids 
fwrite(output_fil[,c("RSID", "INFO/AA")], f_out, sep = " ", quote = F, col.names = F, row.names = F) 
f_unmatched =  paste("/project/mathilab/ppoyraz/scan/iHS_selscan/temp/", vcf$chr[1], "_unmatched.tab", sep = "") # unmatched file 
# add variants without rsids to unmatched file 
unmatched = c(unmatched, subset(output, RSID == ".")$POS)
unmatched_df = data.frame(chr = rep(vcf$chr[1], length(unmatched)), pos = unmatched) 
fwrite(unmatched_df, f_unmatched, sep = "\t", quote = F, col.names = T, row.names = F)
sink(paste("/project/mathilab/ppoyraz/scan/iHS_selscan/temp/", vcf$chr[1], "_AA.log", sep = "")) # log file 
print("length of output") 
print(nrow(output_fil))
print("length of unlifted")
print(nrow(unlifted))
print("length of unmatched") 
print(length(unmatched))
print("length of vcf") 
print(nrow(vcf))
print("vcf - output = unlifted + unmatched") 
print(nrow(vcf) - nrow(output) == nrow(unlifted) + length(unmatched))
print("num lost due to missing RSID") 
print(sum(output$RSID == ".")) 
sink()
