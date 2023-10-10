# Lin Poyraz 
# Usage: Rscript prepare_predixcan_dir.R <dosage dir> 
# Creates sample.txt file from a given directory 
library(data.table)

args = commandArgs(trailingOnly = TRUE)
dir = args[1]

#print(dir)
# get all files 
files = list.files(path = dir, full.names = TRUE)
print(files[1])

# read first line of file 
line = readLines(files[1], n = 1)
print(line)

# separate at \t
split = unlist(strsplit(line, split  = "\t"))
#print(split)
split = split[7:length(split)]
print(split) 

# write samples file 
f = paste(dir, "/sample.txt", sep = "")
res = as.data.frame(matrix(c(split, split), ncol = 2))
fwrite(res, f, row.names = F, col.names = F, quote = F, sep = "\t")

