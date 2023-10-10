#!/usr/bin/env Rscript
#tweaked from iain's version to take arguments.

args = commandArgs(trailingOnly=TRUE)

chr=args[1]
map<-read.table(args[2], as.is=T, header=F)
new.pos<-scan(args[3])
mapfun<-approxfun(map[,4], map[,3], rule=2)
new.map<-mapfun(new.pos)
new.map.file<-data.frame(chr=chr, id=".", map=new.map, pos=new.pos)

new.map.file$pos <- format(new.map.file$pos, scientific = FALSE)
write.table(new.map.file, args[4], row.name=F, col.name=F, sep="\t", quote=F)
