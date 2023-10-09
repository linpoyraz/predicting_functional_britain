Description
-----------

This file produces all main figures included in “Predicting functional
consequences of recent natural selection in Britain”

Packages
--------

    library(ggplot2)
    library(dplyr)
    library(data.table)
    library(ggpubr)
    library(calibrate)
    library(scales)
    library(ggrepel)
    library(gtable)
    library(gt)

Functions to make figures and tables
------------------------------------

### Functions to make manhattan plots

This is the general manhattan plot function. Edited minimally and
adapted from R-package qqman
<a href="https://rdrr.io/cran/qqman/src/R/manhattan.R" class="uri">https://rdrr.io/cran/qqman/src/R/manhattan.R</a>

    # PACKAGE: qqman https://rdrr.io/cran/qqman/src/R/manhattan.R 
    manhattan <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                          col=c("gray10", "gray60"), chrlabs=NULL,
                          suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                          highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
      
      par(mar = c(5, 5, 4, 2) + 0.1)
      # Not sure why, but package check will warn without this.
      CHR=BP=P=index=NULL
      
      # Check for sensible dataset
      ## Make sure you have chr, bp and p columns.
      if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
      if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
      if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
      ## warn if you don't have a snp column
      if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
      ## make sure chr, bp, and p columns are numeric.
      if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
      if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
      if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
      
      # Create a new data.frame with columns called CHR, BP, and P.
      # d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA) # with millions of SNPs, create dataframe at once 
      #  rather than dynamically allocated(see line 72-73, and remove line 87 and line 91 )
      
      # If the input data frame has a SNP column, add it to the new data frame you're creating.
      if (!is.null(x[[snp]])) d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA ,SNP=x[[snp]], stringsAsFactors = FALSE) else 
        d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA)
      
      
      # Set positions, ticks, and labels for plotting
      ## Sort and keep only values where is numeric.
      #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
      #  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))       ## unused, all three variables are numeric, line:63-65 
      d <- d[order(d$CHR, d$BP), ]
      #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
      if (logp) {
        d$logp <- -log10(d$P)
      } else {
        d$logp <- d$P
      }
      # d$pos=NA
      
      
      # Fixes the bug where one chromosome is missing by adding a sequential index column.
      # d$index=NA
      # ind = 0
      # for (i in unique(d$CHR)){
      #     ind = ind + 1
      #     d[d$CHR==i,]$index = ind
      # }
      d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency
      
      # This section sets up positions and ticks. Ticks should be placed in the
      # middle of a chromosome. The a new pos column is added that keeps a running
      # sum of the positions of each successive chromsome. For example:
      # chr bp pos
      # 1   1  1
      # 1   2  2
      # 2   1  3
      # 2   2  4
      # 3   1  5
      nchr = length(unique(d$CHR))
      if (nchr==1) { ## For a single chromosome
        ## Uncomment the next two linex to plot single chr results in Mb
        #options(scipen=999)
        #d$pos=d$BP/1e6
        d$pos=d$BP
        #ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
        xlabel = paste('Chromosome',unique(d$CHR),'position')
        #labs = ticks          ## unused, from code line: 169
      } else { ## For multiple chromosomes
        lastbase=0
        ticks=NULL
        for (i in unique(d$index)) {
          if (i==1) {
            d[d$index==i, ]$pos=d[d$index==i, ]$BP
          } else {
            ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
            lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
            d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
            d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
            # lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
            # d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
            
          }
          # Old way: assumes SNPs evenly distributed
          #ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
          # New way: doesn't make that assumption
          ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
        }
        #ticks <-tapply(d$pos,d$index,mean)   # replace line 135
        xlabel = 'Chromosome'
        #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
        labs <- unique(d$CHR)
      }
      
      # Initialize plot
      xmax = ceiling(max(d$pos) * 1.03)
      xmin = floor(max(d$pos) * -0.03)
      
      # The old way to initialize the plot
      # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
      #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
      
      
      # The new way to initialize the plot.
      ## See http://stackoverflow.com/q/23922130/654296
      ## First, define your default arguments
      def_args <- list(xaxt='n', bty='l', xaxs='i', yaxs='i', las=1, pch=20,
                       xlim=c(xmin,xmax), cex.lab=1.3, ylim=c(0,ceiling(max(d$logp))),
                       xlab=xlabel, ylab=expression(-log[10](italic(p))))
      ## Next, get a list of ... arguments
      #dotargs <- as.list(match.call())[-1L]
      dotargs <- list(...)
      ## And call the plot function passing NA, your ... arguments, and the default
      ## arguments that were not defined in the ... arguments.
      do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

      # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
      if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
          if (length(chrlabs)==length(labs)) {
            labs <- chrlabs
          } else {
            warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
          }
        } else {
          warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
      }
      
      # Add an axis. 
      if (nchr==1) { #If single chromosome, ticks and labels automatic.
        axis(1, ...)
      } else { # if multiple chrs, use the ticks and labels you created above.
        axis(1, at=ticks, labels=labs, cex.axis = 0.75, ...)
      }
      
      # Create a vector of alternatiting colors
      #col=rep(col, max(d$CHR))  # replaced by line 187
      col = rep_len(col, max(d$index))  ## mean this one?  the results are same
      
      # Add points to the plot
      if (nchr==1) {
        with(d, points(pos, logp, pch=20, col=col[1], ...))
      } else {
        # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
        icol=1
        for (i in unique(d$index)) {
          #with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
          points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=20, ...)
          icol=icol+1
        }
      }
      
      # Add suggestive and genomewide lines
      if (suggestiveline) abline(h=suggestiveline, col="blue", lty = 2)
      if (genomewideline) abline(h=genomewideline, col="red", lty = 2)
      
      # Highlight snps from a character vector
      if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight=d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col="green3", pch=20, ...)) 
      }
      
      # Highlight top SNPs
      if (!is.null(annotatePval)) {
        # extract top SNPs at given p-val
        if (logp) {
          topHits = subset(d, P <= annotatePval)
        } else
          topHits = subset(d, P >= annotatePval)
        par(xpd = TRUE)
        # annotate these SNPs
        if (annotateTop == FALSE) {
          if (logp) {
            with(subset(d, P <= annotatePval), 
                 textxy(pos, -log10(P), labs = topHits$SNP, cex = 0.45), ...)
          } else
            with(subset(d, P >= annotatePval), 
                 textxy(pos, P, labs = topHits$SNP, cex = 0.45), ...)
        }
        else {
          # could try alternative, annotate top SNP of each sig chr
          topHits <- topHits[order(topHits$P),]
          topSNPs <- NULL
          
          for (i in unique(topHits$CHR)) {
            
            chrSNPs <- topHits[topHits$CHR == i,]
            topSNPs <- rbind(topSNPs, chrSNPs[1,])
            # for selection scan
            #if (i == 5){
            #  topSNPs <- rbind(topSNPs, chrSNPs[4,])
            #}
            
          #}
          if (logp ){
            # pretty labels 
            topSNPs = topSNPs[order(topSNPs$pos, decreasing = FALSE), ]
            #for(r in 1:(nrow(topSNPs)-1)){
             # if(abs(topSNPs$pos[r] - topSNPs$pos[r+1]) < 5e7){
              #  topSNPs$pos[r] = topSNPs$pos[r] - 1e7
               # topSNPs$pos[r+1] = topSNPs$pos[r+1] + 1e7
              #}
            #}
            # for selection scan manhattan 
            #topSNPs$SNP = c("LCT", "MANF", "SLC45A2", "PDLIM4", "HLA", "SLC26A5", "SHOC2", "DHCR7", "OAS1", "HERC2", "PAPD5", "NF1", "chr18:41.5Mb", "VAV1")
            y = -log10(topSNPs$P) + 3
            long = which(nchar(topSNPs$SNP) >= 5)
            y[long] = y[long] + 1
            longer = which(nchar(topSNPs$SNP) >= 9)
            y[longer] = y[longer] + 1.3

            textxy(topSNPs$pos, y , labs = topSNPs$SNP, cex = 1, font = 3,srt=90,offset = 0, ...)
          } else
            text(topSNPs$pos, -log10(topSNPs$P) + 0.5,labs = topSNPs$SNP, cex = 0.9, font = 3, srt=90, ...)
        }
      }  
      par(xpd = FALSE)
      }
    }

This is the manhattan plot function edited to annotate gene peaks in the
genome-wide selection scan.

    # PACKAGE: qqman https://rdrr.io/cran/qqman/src/R/manhattan.R 
    manhattan_ss <- function(x, chr="CHR", bp="BP", p="P", snp="SNP", 
                          col=c("gray10", "gray60"), chrlabs=NULL,
                          suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), 
                          highlight=NULL, logp=TRUE, annotatePval = NULL, annotateTop = TRUE, ...) {
      
      par(mar = c(5, 5, 4, 2) + 0.1)
      # Not sure why, but package check will warn without this.
      CHR=BP=P=index=NULL
      
      # Check for sensible dataset
      ## Make sure you have chr, bp and p columns.
      if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
      if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
      if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
      ## warn if you don't have a snp column
      if (!(snp %in% names(x))) warning(paste("No SNP column found. OK unless you're trying to highlight."))
      ## make sure chr, bp, and p columns are numeric.
      if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
      if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
      if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
      
      # Create a new data.frame with columns called CHR, BP, and P.
      # d=data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA) # with millions of SNPs, create dataframe at once 
      #  rather than dynamically allocated(see line 72-73, and remove line 87 and line 91 )
      
      # If the input data frame has a SNP column, add it to the new data frame you're creating.
      if (!is.null(x[[snp]])) d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA ,SNP=x[[snp]], stringsAsFactors = FALSE) else 
        d = data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos = NA, index = NA)
      
      
      # Set positions, ticks, and labels for plotting
      ## Sort and keep only values where is numeric.
      #d <- subset(d[order(d$CHR, d$BP), ], (P>0 & P<=1 & is.numeric(P)))
      #  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))       ## unused, all three variables are numeric, line:63-65 
      d <- d[order(d$CHR, d$BP), ]
      #d$logp <- ifelse(logp, yes=-log10(d$P), no=d$P)
      if (logp) {
        d$logp <- -log10(d$P)
      } else {
        d$logp <- d$P
      }
      # d$pos=NA
      
      
      # Fixes the bug where one chromosome is missing by adding a sequential index column.
      # d$index=NA
      # ind = 0
      # for (i in unique(d$CHR)){
      #     ind = ind + 1
      #     d[d$CHR==i,]$index = ind
      # }
      d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency
      
      # This section sets up positions and ticks. Ticks should be placed in the
      # middle of a chromosome. The a new pos column is added that keeps a running
      # sum of the positions of each successive chromsome. For example:
      # chr bp pos
      # 1   1  1
      # 1   2  2
      # 2   1  3
      # 2   2  4
      # 3   1  5
      nchr = length(unique(d$CHR))
      if (nchr==1) { ## For a single chromosome
        ## Uncomment the next two linex to plot single chr results in Mb
        #options(scipen=999)
        #d$pos=d$BP/1e6
        d$pos=d$BP
        #ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
        xlabel = paste('Chromosome',unique(d$CHR),'position')
        #labs = ticks          ## unused, from code line: 169
      } else { ## For multiple chromosomes
        lastbase=0
        ticks=NULL
        for (i in unique(d$index)) {
          if (i==1) {
            d[d$index==i, ]$pos=d[d$index==i, ]$BP
          } else {
            ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
            lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
            d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
            d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
            # lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
            # d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
            
          }
          # Old way: assumes SNPs evenly distributed
          #ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
          # New way: doesn't make that assumption
          ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
        }
        #ticks <-tapply(d$pos,d$index,mean)   # replace line 135
        xlabel = 'Chromosome'
        #labs = append(unique(d$CHR),'') ## I forgot what this was here for... if seems to work, remove.
        labs <- unique(d$CHR)
      }
      
      # Initialize plot
      xmax = ceiling(max(d$pos) * 1.03)
      xmin = floor(max(d$pos) * -0.03)
      
      # The old way to initialize the plot
      # plot(NULL, xaxt='n', bty='n', xaxs='i', yaxs='i', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
      #      xlab=xlabel, ylab=expression(-log[10](italic(p))), las=1, pch=20, ...)
      
      
      # The new way to initialize the plot.
      ## See http://stackoverflow.com/q/23922130/654296
      ## First, define your default arguments
      def_args <- list(xaxt='n', bty='l', xaxs='i', yaxs='i', las=1, pch=20,
                       xlim=c(xmin,xmax), cex.lab=1.3, ylim=c(0,ceiling(max(d$logp))),
                       xlab=xlabel, ylab=expression(-log[10](italic(p))))
      ## Next, get a list of ... arguments
      #dotargs <- as.list(match.call())[-1L]
      dotargs <- list(...)
      ## And call the plot function passing NA, your ... arguments, and the default
      ## arguments that were not defined in the ... arguments.
      do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

      # If manually specifying chromosome labels, ensure a character vector and number of labels matches number chrs.
      if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
          if (length(chrlabs)==length(labs)) {
            labs <- chrlabs
          } else {
            warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
          }
        } else {
          warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
      }
      
      # Add an axis. 
      if (nchr==1) { #If single chromosome, ticks and labels automatic.
        axis(1, ...)
      } else { # if multiple chrs, use the ticks and labels you created above.
        axis(1, at=ticks, labels=labs, cex.axis = 0.75, ...)
      }
      
      # Create a vector of alternatiting colors
      #col=rep(col, max(d$CHR))  # replaced by line 187
      col = rep_len(col, max(d$index))  ## mean this one?  the results are same
      
      # Add points to the plot
      if (nchr==1) {
        with(d, points(pos, logp, pch=20, col=col[1], ...))
      } else {
        # if multiple chromosomes, need to alternate colors and increase the color index (icol) each chr.
        icol=1
        for (i in unique(d$index)) {
          #with(d[d$index==unique(d$index)[i], ], points(pos, logp, col=col[icol], pch=20, ...))
          points(d[d$index==i,"pos"], d[d$index==i,"logp"], col=col[icol], pch=20, ...)
          icol=icol+1
        }
      }
      
      # Add suggestive and genomewide lines
      if (suggestiveline) abline(h=suggestiveline, col="blue", lty = 2)
      if (genomewideline) abline(h=genomewideline, col="red", lty = 2)
      
      # Highlight snps from a character vector
      if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight=d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col="green3", pch=20, ...)) 
      }
      
      # Highlight top SNPs
      if (!is.null(annotatePval)) {
        # extract top SNPs at given p-val
        if (logp) {
          topHits = subset(d, P <= annotatePval)
        } else
          topHits = subset(d, P >= annotatePval)
        par(xpd = TRUE)
        # annotate these SNPs
        if (annotateTop == FALSE) {
          if (logp) {
            with(subset(d, P <= annotatePval), 
                 textxy(pos, -log10(P), labs = topHits$SNP, cex = 0.45), ...)
          } else
            with(subset(d, P >= annotatePval), 
                 textxy(pos, P, labs = topHits$SNP, cex = 0.45), ...)
        }
        else {
          # could try alternative, annotate top SNP of each sig chr
          topHits <- topHits[order(topHits$P),]
          topSNPs <- NULL
          
          for (i in unique(topHits$CHR)) {
            
            chrSNPs <- topHits[topHits$CHR == i,]
            topSNPs <- rbind(topSNPs, chrSNPs[1,])
            # for selection scan
            if (i == 5){
              topSNPs <- rbind(topSNPs, chrSNPs[3,])
            }
            
          }
          if (logp ){
            # pretty labels 
            topSNPs = topSNPs[order(topSNPs$pos, decreasing = FALSE), ]
            #for(r in 1:(nrow(topSNPs)-1)){
             # if(abs(topSNPs$pos[r] - topSNPs$pos[r+1]) < 5e7){
              #  topSNPs$pos[r] = topSNPs$pos[r] - 1e7
               # topSNPs$pos[r+1] = topSNPs$pos[r+1] + 1e7
              #}
            #}
            # for selection scan manhattan 
            topSNPs$SNP = c("LCT", "MANF", "SLC45A2", "PDLIM4", "HLA", "SLC26A5", "SHOC2", "DHCR7", "OAS1", "HERC2", "PAPD5", "NF1", "chr18:41.3Mb", "VAV1", "chr21:30.2Mb")
            y = -log10(topSNPs$P) + 3
            long = which(nchar(topSNPs$SNP) >= 5)
            y[long] = y[long] + 1
            longer = which(nchar(topSNPs$SNP) >= 9)
            y[longer] = y[longer] + 1.3

            textxy(topSNPs$pos, y , labs = topSNPs$SNP, cex = 1, font = 3,srt=90,offset = 0, ...)
          } else
            text(topSNPs$pos, -log10(topSNPs$P) + 0.5,labs = topSNPs$SNP, cex = 0.9, font = 3, srt=90, ...)
        }
      }  
      par(xpd = FALSE)
      }

This function takes in outputs from time\_series\_regression.R annotated
with imputation qualities with imputation\_results.R and turns them into
MH\_df format, which is the input format for our Manhattan plot
functions. This function filters out the genes with the lowest 20%
imputation scores, performs GC, FDR and Bonferroni correction.

    # Returns regression results file in MH_df format 
    to_MH_df = function(results){
      # filter for imputation quality, first remove those without quality values 
      results = results[which(!is.na(results$Mean.Imputation.Rsq)),]
      # then remove lowest 20% for imputation quality 
      results = results[which(results$Mean.Imputation.Rsq > quantile(results$Mean.Imputation.Rsq, 0.2)),]
      # now, add GC, FDR and Bonferroni 
      p = results$P.value
      chi_stat = qchisq(p, df = 1, lower.tail=F)
      lambda = median(chi_stat)/qchisq(0.5, df=1) #inflation factor
      corrected_p = pchisq(chi_stat/lambda, df=1, lower.tail=F)
      results$GC = corrected_p
      results$FDR = p.adjust(corrected_p, method = "fdr")
      results$Bonferroni = p.adjust(corrected_p, method = "bonferroni")
      # convert to MH format 
      df = data.frame("CHR" = rep(NA, nrow(results)), "POS" = rep(NA, nrow(results)), "PVAL" = rep(NA, nrow(results)))
      # Set values for df 
      df$CHR = gsub(results$chr, pattern = "chr", replacement = "")
      df$POS = round((results$start + results$end) / 2)
      df$START = results$start
      df$END = results$end
      df$PVAL = results$GC
      df$Gene = results$Gene
      df$CHR = as.numeric(df$CHR)
      df$BETA = results$Beta
      df$FDR = results$FDR
      df$RAW = p
      df$Ensembl = results$Ensembl
      return(df)
    }

This function returns FDR and Bonferroni cutoffs from MH\_df files to
input into the manhattan plot functions.

    # Returns FDR and Bonferroni cutoffs from MH_df file 
    cutoffs = function(df){
      # find sig FDR cutoff 
      FDR = p.adjust(df$PVAL, "fdr")
      i = which(FDR <= 0.05)
      fdr_cut = max(df$PVAL[i])
      # find sig bonferroni cutoff 
      b = p.adjust(df$PVAL, "bonferroni")
      i = which(b <= 0.05)
      b_cut = max(df$PVAL[i])
     return(c(fdr_cut, b_cut)) 
    }

### Functions to combine GWSS and TWSS results

This function generates a list of selected regions from selection scan
by combining consecutive (&lt;5e6 bp apart) FDR-significant windows. It
adds 1e5 bp buffers to each “consecutive region under selection.”

    find_under_selection = function(selection, fdr_cutoff){
      i = 1 
      under_selection = data.frame("chr" = "a", "start" = 0, "end" = 0)
      # scan across selection scan file 
      while (i < nrow(selection)){
        # if a region is under selection (p value is < fdr cutoff)
        if (selection$P1[i] <= fdr_cutoff){
          # start 
          start = selection$start[i]
          j = i + 1
          # if the consecutive regions are also significant 
          while (j < nrow(selection) & selection$P1[j] <= fdr_cutoff) {
            j = j + 1
          }
          end = selection$start[j]
          # add to under selection 
          if (selection$CHR[i] != selection$CHR[j]){
            print("error")
          }
          under_selection = rbind(under_selection, data.frame("chr" = selection$CHR[i], "start" = start, "end" = end))
          i = j 
        } else {
          i = i + 1
        }
      }
      under_selection = under_selection[2:nrow(under_selection),]
      consecutive_under_selection = data.frame("chr" = "a", "start" = 0, "end" = 0)

      # extract consecutive regions under selection
      for (i in unique(under_selection$chr)){
        pos = sort(c(under_selection$start[which(under_selection$chr == i)], under_selection$end[which(under_selection$chr == i)]))
        k = 1
        while (k < length(pos)){
          j = length(pos)
          while (pos[j] - pos[k] > 5e6 & j > k){
            j = j - 1
          }
          consecutive_under_selection = rbind(consecutive_under_selection, data.frame("chr" = i, "start" = pos[k], "end" = pos[j]))
          k = j + 1
        }
      }
      consecutive_under_selection = consecutive_under_selection[2:nrow(consecutive_under_selection),]

        # add 1e5 buffers 
      consecutive_under_selection$start = consecutive_under_selection$start - 1e5
      consecutive_under_selection$end = consecutive_under_selection$end + 1e5
        
      return(consecutive_under_selection)
    }

Functions to create combined selection scan and expression time-series
plots.

    # This function creates manhattan plot backdrop for the segment plots 
    segment_manhattan_plot = function(selection_scan, chr, start_end, p, main, fdr_cutoff, b_cutoff){
      manhattan(subset(selection_scan, CHR == chr), ylim = p, xlim = start_end, chr = "CHR", bp = "POS", snp = "lead.snp", p = "P1", suggestiveline = -log10(fdr_cutoff), genomewideline = -log10(b_cutoff), cex = 0.001, main = main) 
      subset_segments = subset(selection_scan, CHR == chr)
      segments(x0 = subset_segments$start, y0 = -log10(subset_segments$P1), x1 = subset_segments$end, lwd = 4, col = "black")
    }
      
    # selection_scan: hg38 selection scan results file 
    # df: hg38 time series regression results file in MH_df format 
    # start_end = c(start position for manhattan plot, end position for manhattan plot)
    # fdr_cutoff, b_cutoffs are horizontal lines
    # main is the title of the plot
    scan_regression_overlap = function(selection_scan, df, chr, start_end, fdr_cutoff, b_cutoff, main){
      subset_scan = subset(selection_scan, CHR == chr & POS >= start_end[1] & POS <= start_end[2])
      subset = subset(df, CHR == chr & START >= start_end[1] & END <= start_end[2]) 
      p_cut = max(c(-log10(b_cutoff), max(-log10(subset_scan$P1)), max(-log10(subset$PVAL))), na.rm = T) + 5
      segment_manhattan_plot(selection_scan, chr, start_end, c(0,p_cut), main, fdr_cutoff, b_cutoff) 
      if (nrow(subset) > 0){
        color = ifelse(subset$BETA > 0, yes = ifelse(subset$PVAL <= fdr_cutoff, yes = "red", no = alpha("gray", 0.6)), no = ifelse(subset$PVAL <= fdr_cutoff, yes = "blue", no = alpha("gray", 0.6)))
        color_text = ifelse(subset$BETA > 0, yes = ifelse(subset$PVAL <= fdr_cutoff, yes = "red", no = "gray27"), no = ifelse(subset$PVAL <= fdr_cutoff, yes = "blue", no = "gray27"))
        shape = ifelse(subset$BETA > 0, yes = 24, no = 25)
        points(x = subset$POS, y = -log10(subset$PVAL), pch = shape, col = color, bg = color, cex = 4)
        text(x = subset$POS, y = -log10(subset$PVAL), labels = subset$Gene, cex = 1.2, offset = 1.45, pos = 3, col = color_text, font = 3)
      }
    }

### Functions to visualize expression changes across time

Function to plot expression against time.

    # Reads in the expression file of given gene and plots expression against time. 
    # Note that expression files must be generated beforehand using ./analysis/expression_time_series/sig_genes.R 
    expression_plot = function(gene){
      gene_df = read.table(paste("./data/expression/", gene, ".txt.gz", sep = ""), header = T, sep = "\t")
      gene_df = na.omit(gene_df)
      p = ggplot(data = gene_df, aes(x = time, y = gene)) + geom_point(size = 1.3) + theme_bw() + theme(text = element_text(size = 15))+ labs(y = "Predicted expression", x = "Time", title = gene) + 
        scale_color_gradient(low = "black", high = "black") + geom_smooth(method = "lm", color = "blue") + lims(y = c(-2.5, 2.5)) + theme(legend.position = "none",  plot.title=element_text(face="italic"))
      return(p)
    }

Functions to create allele frequency against time plots.

    # Calculates change in time of each snp (beta in snp ~ time), Ensures that all alleles are trait increasing 
    allele_freq_paper = function(gene){
        gene_df = read.table(paste("./data/allele_frequency/", gene, ".txt.gz", sep = ""), header = T, sep = "\t")
        weights = read.table(paste("./data/allele_frequency/", gene, "_weights.txt.gz", sep = ""), header = T, sep = "\t")
        expression = read.table(paste("./data/expression/", gene, ".txt.gz", sep = ""), header = T, sep = "\t")
      joint = merge(gene_df, expression, by = "iid")
      gene_df$time = as.numeric(gene_df$time)
      num_snp = ncol(gene_df) - 2
      af_df = data.frame(SNP = colnames(gene_df)[1:num_snp], time_beta = rep(NA, num_snp))
      for (i in 1:num_snp){
        # make all trait increasing 
        snp = colnames(gene_df)[i] 
        weight = weights$weight[which(weights$rsid == snp)]
        m_df = gene_df[,c(i, which(colnames(gene_df) == "time"))]
        names(m_df) = c("snp", "time")
         if (weight < 0){
          m_df$snp = 2 - m_df$snp
        }
        m = summary(lm(snp ~ time, data = m_df))$coefficients["time", "Estimate"]
        af_df$time_beta[i] = m
      }
      return(af_df)
    }

    af_plot_paper = function(gene){
      weights = read.table(paste("./data/allele_frequency/", gene, "_weights.txt.gz", sep = ""), header = T, sep = "\t")
      af_df = allele_freq_paper(gene)
      af_df$beta - rep(1000, nrow(af_df))
      af_df$dire = rep(NA, nrow(af_df))
      for (i in 1:nrow(af_df)){
        j = which(weights$rsid == af_df$SNP[i])
        af_df$beta[i] = abs(as.numeric(weights$weight[j]))
        if ((af_df$beta[i] * af_df$time_beta[i]) > 0){
          af_df$dire[i] = "up"
        } else {
          af_df$dire[i] = "down"
        }
      }
      p = ggplot(data = af_df, aes(x = beta, y = time_beta, fill = dire, color = dire, shape = dire, label = SNP)) + geom_point(size = 2.5) + theme_bw() +
        scale_color_manual(breaks = c("down", "up"), values = c("blue", "red")) + labs(x = "Effect on normalized expression", y = "Average change in frequency per-year", fill = "SNP", title = gene) + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + coord_cartesian(xlim = c(0, 1.7), ylim = c(-3e-4, 3e-4)) +
        theme(legend.position = "none", plot.title=element_text(face="italic")) + scale_shape_manual(values = c("up" = 24, "down" = 25)) + scale_fill_manual(breaks = c("down", "up"), values = c("blue", "red")) 
      return(p)
    }

Figure 1
--------

We begin by preparing the data:

    # Read in imputation qualities: 
    imp_qualities = read.table("./data/imputation_results.txt.gz", header = T, sep = "\t") 
    imp_scores = imp_qualities[,c("Gene", "Mean.Imputation.Rsq")]
    # Read in output from time_series_regression limited to samples with >0.01 coverage
    samples_01 = read.table("./data/samples01_time_series_regression_results.txt.gz", header = T) 
    # Populate output with imputation quality scores 
    samples_01$Mean.Imputation.Rsq = NA
    for (i in 1:nrow(samples_01)){
      g = which(imp_scores$Gene == samples_01$Gene[i])[1]
      samples_01$Mean.Imputation.Rsq[i] = imp_scores$Mean.Imputation.Rsq[g]
    }
    MH_df_whole = to_MH_df(samples_01)
    r_fdr_cut = cutoffs(MH_df_whole)[1]
    r_b_cut = cutoffs(MH_df_whole)[2]

    # Read in randomized p-values 
    random  = read.table("./data/samples01_random_time_series_regression_results.txt.gz", header = T, sep = "\t") 
    random = random[,c("Ensembl", "P.value")]
    names(random) = c("Ensembl", "Random.P.value")
    # merge randomized and non-randomized data frames on Ensembl gene names 
    qq_regression = MH_df_whole[c("Ensembl", "RAW", "PVAL")]
    merged_twss = merge(random, qq_regression, by = "Ensembl")

    # Read in selection scan and identify regions under selection. 
    selection_scan = read.table("./data/s_scan_all_brit2_01_window.txt.gz", header = T, sep = "\t")

    # Read in randomized p-values from selection scan 
    random_selection_scan = read.table("./data/s_scan_all_brit2_01_random_window.txt.gz", header = T, sep = "\t")

We plot the Manhattan plot for the transciptome-wide selection scan
(TWSS).

    manhattan(MH_df_whole, ylim = c(0,20), chrlabs = as.character(1:22), chr = "CHR", bp = "POS", snp = "Gene", p = "PVAL", annotatePval = r_fdr_cut, suggestiveline = -log10(r_fdr_cut), genomewideline = -log10(r_b_cut)) 

![](main_figures_files/figure-markdown_strict/unnamed-chunk-11-1.png)

We generate a qq-plot for the TWSS.

    p_results = as.numeric(merged_twss$RAW)
    p_random = as.numeric(merged_twss$Random.P.value)
    p_GC = as.numeric(merged_twss$PVAL)
    observed_results = -log10(sort(p_results))
    observed_random = -log10(sort(p_random))
    observed_GC = -log10(sort(p_GC))
    expected = -log10(ppoints(length(p_results)))

    df = data.frame("GC.Observed" = observed_GC, "Randomized" = observed_random, "Expected" = expected)
    qqlot = ggplot(data = df, aes(x = Expected, y = GC.Observed)) + 
      geom_point(color = "darkblue") + geom_abline(intercept = 0, slope = 1) + labs(cex = 1.4, y = expression("Observed"-log[10](italic(p))), x = expression("Expected"-log[10](italic(p))))+ theme_bw() +  
      geom_point(data = df, aes(x = Expected, y = Randomized),color = "darkgray") + theme(axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
    qqlot

![](main_figures_files/figure-markdown_strict/unnamed-chunk-12-1.png)

We plot the Manhattan plot for the genome-wide selection scan (GWSS).

    # find sig FDR cutoff 
    FDR = p.adjust(selection_scan$P1, "fdr")
    i = which(FDR <= 0.05)
    scan_fdr_cut = max(selection_scan$P1[i])
    # find sig bonferroni cutoff 
    b = p.adjust(selection_scan$P1, "bonferroni")
    i = which(b <= 0.05)
    scan_b_cut = max(selection_scan$P1[i])
    manhattan_ss(selection_scan, ylim = c(0,20), chrlabs = as.character(1:22), chr = "CHR", bp = "POS", snp = "lead.snp", p = "P1", annotatePval = scan_fdr_cut, logp = T, suggestiveline = -log10(scan_fdr_cut), genomewideline = -log10(scan_b_cut)) 

![](main_figures_files/figure-markdown_strict/unnamed-chunk-13-1.png)

We generate a qq-plot for the GWSS.

    p_random = as.numeric(random_selection_scan$P2)
    p_GC = as.numeric(selection_scan$P1)
    observed_random = -log10(sort(p_random))
    observed_GC = -log10(sort(p_GC))
    expected = -log10(ppoints(length(p_random)))
    df = data.frame("GC.Observed" = observed_GC, "Randomized" = observed_random, "Expected" = expected)

    qqlot = ggplot(data = df, aes(x = Expected, y = GC.Observed)) + 
      geom_point(color = "darkblue") + geom_abline(intercept = 0, slope = 1) + labs( y = expression("Observed"-log[10](italic(p))), x = expression("Expected"-log[10](italic(p))))+ theme_bw() +
      geom_point(data = df, aes(x = Expected, y = Randomized),color = "darkgray") + theme(axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
    qqlot

![](main_figures_files/figure-markdown_strict/unnamed-chunk-14-1.png)

Tables 1 & 2
------------

We begin by Table 2, as Table 1 depends on the outputs of this process.

We identify consecutive regions under selection, and add gene
information to each selected region. For each selected region, we
include genes that start and have a midpoint within the bounds of the
selected region (including buffers).

    results_whole = samples_01 # results before filtering for imputation quality 
    results_whole$pos = round((results_whole$start + results_whole$end)/2)
    results_whole$chr = gsub(pattern = "chr", replacement ="", results_whole$chr)
    # We identify consecutive genes under selection 
    selection_to_gene = find_under_selection(selection_scan, scan_fdr_cut) 
    # This is for the purpose of counting all protein coding genes in the region, and the number for which we have models (meaning that they were not filtered out due to low imputation quality)
    selection_to_gene$sig_genes = NA
    selection_to_gene$all_genes_modeled_num = NaN
    selection_to_gene$all_genes_num = NaN
    selection_to_gene$sig_genes_num = NaN
    for (i in 1:nrow(selection_to_gene)){
      # select genes that start and have a midpoint within the bounds of the selected region
      all_genes = results_whole[which(results_whole$chr == selection_to_gene$chr[i] & results_whole$start >= selection_to_gene$start[i] & results_whole$pos <= selection_to_gene$end[i]),]
      all_genes_modeled = MH_df_whole[which(MH_df_whole$CHR == selection_to_gene$chr[i] & MH_df_whole$START >= selection_to_gene$start[i] & MH_df_whole$POS <= selection_to_gene$end[i]),]
      selection_to_gene$all_genes_modeled_num[i] = nrow(all_genes_modeled)
      selection_to_gene$all_genes_num[i] = nrow(all_genes)
      sig_genes = all_genes_modeled[which(all_genes_modeled$FDR < 0.05),]
      selection_to_gene$sig_genes[i] = paste(sig_genes$Gene, collapse = ", ")
      selection_to_gene$sig_genes_num[i] = nrow(sig_genes)
    }

We then generate Table 2:

    selection_to_gene$chr = as.numeric(selection_to_gene$chr)
    selection_to_gene$start = as.numeric(selection_to_gene$start)
    selection_to_gene = selection_to_gene[order(selection_to_gene$chr, selection_to_gene$start),]
    selection_to_gene_tbl = gt(data = selection_to_gene) %>% cols_label(
        chr = html("Chr"),
        start = html("GWSS Start"),
        end = html("GWSS End"),
        sig_genes = html("TWSS Significant Genes"),
        all_genes_num = "# All Genes",
        all_genes_modeled_num = "# All Modeled Genes",
        sig_genes_num = "# Significant Genes")

    selection_to_gene_tbl

<div id="hjvgojgfjt" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#hjvgojgfjt .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#hjvgojgfjt .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#hjvgojgfjt .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#hjvgojgfjt .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#hjvgojgfjt .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#hjvgojgfjt .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hjvgojgfjt .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#hjvgojgfjt .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#hjvgojgfjt .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#hjvgojgfjt .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#hjvgojgfjt .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#hjvgojgfjt .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#hjvgojgfjt .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#hjvgojgfjt .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#hjvgojgfjt .gt_from_md > :first-child {
  margin-top: 0;
}

#hjvgojgfjt .gt_from_md > :last-child {
  margin-bottom: 0;
}

#hjvgojgfjt .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#hjvgojgfjt .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#hjvgojgfjt .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#hjvgojgfjt .gt_row_group_first td {
  border-top-width: 2px;
}

#hjvgojgfjt .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hjvgojgfjt .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#hjvgojgfjt .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#hjvgojgfjt .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hjvgojgfjt .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#hjvgojgfjt .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#hjvgojgfjt .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#hjvgojgfjt .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#hjvgojgfjt .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#hjvgojgfjt .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#hjvgojgfjt .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#hjvgojgfjt .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#hjvgojgfjt .gt_left {
  text-align: left;
}

#hjvgojgfjt .gt_center {
  text-align: center;
}

#hjvgojgfjt .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#hjvgojgfjt .gt_font_normal {
  font-weight: normal;
}

#hjvgojgfjt .gt_font_bold {
  font-weight: bold;
}

#hjvgojgfjt .gt_font_italic {
  font-style: italic;
}

#hjvgojgfjt .gt_super {
  font-size: 65%;
}

#hjvgojgfjt .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#hjvgojgfjt .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#hjvgojgfjt .gt_indent_1 {
  text-indent: 5px;
}

#hjvgojgfjt .gt_indent_2 {
  text-indent: 10px;
}

#hjvgojgfjt .gt_indent_3 {
  text-indent: 15px;
}

#hjvgojgfjt .gt_indent_4 {
  text-indent: 20px;
}

#hjvgojgfjt .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Chr">Chr</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="GWSS Start">GWSS Start</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="GWSS End">GWSS End</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="TWSS Significant Genes">TWSS Significant Genes</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="# All Modeled Genes"># All Modeled Genes</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="# All Genes"># All Genes</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="# Significant Genes"># Significant Genes</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="chr" class="gt_row gt_right">2</td>
<td headers="start" class="gt_row gt_right">134292811</td>
<td headers="end" class="gt_row gt_right">136385296</td>
<td headers="sig_genes" class="gt_row gt_left">MCM6, LCT, DARS, CXCR4, TMEM163, MAP3K19</td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">11</td>
<td headers="all_genes_num" class="gt_row gt_right">12</td>
<td headers="sig_genes_num" class="gt_row gt_right">6</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">3</td>
<td headers="start" class="gt_row gt_right">50849763</td>
<td headers="end" class="gt_row gt_right">51387494</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">1</td>
<td headers="all_genes_num" class="gt_row gt_right">1</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">5</td>
<td headers="start" class="gt_row gt_right">33777346</td>
<td headers="end" class="gt_row gt_right">34069589</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">2</td>
<td headers="all_genes_num" class="gt_row gt_right">4</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">5</td>
<td headers="start" class="gt_row gt_right">132101060</td>
<td headers="end" class="gt_row gt_right">132454548</td>
<td headers="sig_genes" class="gt_row gt_left">P4HA2, PDLIM4, SLC22A5</td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">5</td>
<td headers="all_genes_num" class="gt_row gt_right">5</td>
<td headers="sig_genes_num" class="gt_row gt_right">3</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">6</td>
<td headers="start" class="gt_row gt_right">28140757</td>
<td headers="end" class="gt_row gt_right">33178885</td>
<td headers="sig_genes" class="gt_row gt_left">PPP1R18, TUBB, HLA-DMA, PBX2, RNF5, APOM, CCHCR1, CDSN, PSORS1C1, ATF6B, HLA-DPA1, C4A</td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">148</td>
<td headers="all_genes_num" class="gt_row gt_right">148</td>
<td headers="sig_genes_num" class="gt_row gt_right">12</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">6</td>
<td headers="start" class="gt_row gt_right">128796160</td>
<td headers="end" class="gt_row gt_right">129062458</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">0</td>
<td headers="all_genes_num" class="gt_row gt_right">0</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">8</td>
<td headers="start" class="gt_row gt_right">33218070</td>
<td headers="end" class="gt_row gt_right">33438568</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">1</td>
<td headers="all_genes_num" class="gt_row gt_right">1</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">10</td>
<td headers="start" class="gt_row gt_right">49643164</td>
<td headers="end" class="gt_row gt_right">50348209</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">5</td>
<td headers="all_genes_num" class="gt_row gt_right">7</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">10</td>
<td headers="start" class="gt_row gt_right">110778270</td>
<td headers="end" class="gt_row gt_right">111086977</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">4</td>
<td headers="all_genes_num" class="gt_row gt_right">4</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">11</td>
<td headers="start" class="gt_row gt_right">61684455</td>
<td headers="end" class="gt_row gt_right">61903876</td>
<td headers="sig_genes" class="gt_row gt_left">FADS1</td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">6</td>
<td headers="all_genes_num" class="gt_row gt_right">6</td>
<td headers="sig_genes_num" class="gt_row gt_right">1</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">11</td>
<td headers="start" class="gt_row gt_right">71302258</td>
<td headers="end" class="gt_row gt_right">71592390</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">4</td>
<td headers="all_genes_num" class="gt_row gt_right">7</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">12</td>
<td headers="start" class="gt_row gt_right">110796571</td>
<td headers="end" class="gt_row gt_right">113044151</td>
<td headers="sig_genes" class="gt_row gt_left">OAS3, FAM109A</td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">17</td>
<td headers="all_genes_num" class="gt_row gt_right">21</td>
<td headers="sig_genes_num" class="gt_row gt_right">2</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">13</td>
<td headers="start" class="gt_row gt_right">111558732</td>
<td headers="end" class="gt_row gt_right">111786334</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">0</td>
<td headers="all_genes_num" class="gt_row gt_right">0</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">15</td>
<td headers="start" class="gt_row gt_right">27951279</td>
<td headers="end" class="gt_row gt_right">29045218</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">1</td>
<td headers="all_genes_num" class="gt_row gt_right">5</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">16</td>
<td headers="start" class="gt_row gt_right">49972683</td>
<td headers="end" class="gt_row gt_right">50325572</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">2</td>
<td headers="all_genes_num" class="gt_row gt_right">4</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">16</td>
<td headers="start" class="gt_row gt_right">82983756</td>
<td headers="end" class="gt_row gt_right">83190020</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">0</td>
<td headers="all_genes_num" class="gt_row gt_right">0</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">17</td>
<td headers="start" class="gt_row gt_right">30984555</td>
<td headers="end" class="gt_row gt_right">31423638</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">4</td>
<td headers="all_genes_num" class="gt_row gt_right">4</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">18</td>
<td headers="start" class="gt_row gt_right">41362995</td>
<td headers="end" class="gt_row gt_right">41666874</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">0</td>
<td headers="all_genes_num" class="gt_row gt_right">0</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
    <tr><td headers="chr" class="gt_row gt_right">21</td>
<td headers="start" class="gt_row gt_right">43286434</td>
<td headers="end" class="gt_row gt_right">44295140</td>
<td headers="sig_genes" class="gt_row gt_left"></td>
<td headers="all_genes_modeled_num" class="gt_row gt_right">11</td>
<td headers="all_genes_num" class="gt_row gt_right">12</td>
<td headers="sig_genes_num" class="gt_row gt_right">0</td></tr>
  </tbody>
  
  
</table>
</div>

We identify consecutive selected regions for each FDR significant gene.
Genes start and have a midpoint within the bounds of the selected region
(including buffers).

    fdr_sig = MH_df_whole[which(MH_df_whole$FDR < 0.05), ]
    gene_to_selection_scan = data.frame("Chr" = fdr_sig$CHR, "Start" = fdr_sig$START, "End" = fdr_sig$END,
                                        "Position" = fdr_sig$POS, "Gene" = fdr_sig$Gene,
                                        "GC_P" = fdr_sig$PVAL, 
                                        "Beta" = fdr_sig$BETA, "GWSS_Start" = NA, "GWSS_End" = NA)
    for (i in 1:nrow(gene_to_selection_scan)){
      j = which(gene_to_selection_scan$Chr[i] == selection_to_gene$chr & gene_to_selection_scan$Start[i] >= selection_to_gene$start & gene_to_selection_scan$Position[i] <= selection_to_gene$end)
      if (length(j) > 1){
        print("ERROR, two matches")
      } else if (length(j) == 0) { 
        gene_to_selection_scan$GWSS_Start[i] = 0
        gene_to_selection_scan$GWSS_End[i] = 0
      } else {
        gene_to_selection_scan$GWSS_Start[i] = selection_to_gene$start[j]
        gene_to_selection_scan$GWSS_End[i] = selection_to_gene$end[j]
      }
    }

Next, we add tissue information for each gene.

    # Add tissue information to gene_to_selection_scan 
    # read in tissue information 
    model_info = read.table("./data/model_info.txt.gz", header = T)
    # make tissue pretty
    model_info$tissue = gsub(pattern = "_", replace = " ", x = model_info$tissue)
    gene_to_selection_scan$Tissue = NA
    # then, add in tissue information to gene_to_selection_scan 
    for (i in 1:nrow(gene_to_selection_scan)){
      j = which(model_info$genename == gene_to_selection_scan$Gene[i])
      if(gene_to_selection_scan$Gene[i] == "FAM109A"){
          j = which(model_info$genename == "PHETA1") #these genes are the same, but have different aliases across these two datasets 
      } 
       gene_to_selection_scan$Tissue[i] = model_info$tissue[j]
    }

We then add loci names to the corresponding selected regions for each
gene.

    # list of loci 
    gwss_loci = selection_to_gene[which(selection_to_gene$sig_genes_num != 0),]
    gwss_loci$locus = c("LCT", "PDLIM4", "HLA", "FADS1", "OAS")
    # Add GWSS Peak to gene_to_selection_scan 
    gene_to_selection_scan$GWSS_Peak = NA
    for (i in 1:nrow(gene_to_selection_scan)){
      chr = gene_to_selection_scan$Chr[i]
      start = gene_to_selection_scan$GWSS_Start[i]
      end = gene_to_selection_scan$GWSS_End[i]
      pos = (start+end)/2
      j = which(gwss_loci$chr == chr & gwss_loci$start <= start & gwss_loci$end >= pos)
      if (length(j) != 0){
        gene_to_selection_scan$GWSS_Peak[i] = gwss_loci$locus[j]
      } else { 
          gene_to_selection_scan$GWSS_Peak[i] = "Novel"
        }
    }

We then generate Table 1:

    gene_to_selection_scan = gene_to_selection_scan[order(gene_to_selection_scan$Chr, gene_to_selection_scan$Start),]
    gene_to_selection_scan$GC_P = signif(gene_to_selection_scan$GC_P, 4)
    gene_to_selection_scan$Beta = signif(gene_to_selection_scan$Beta, 4)
    # which to italicize in the GWSS signal row
    ital_gwss = 
    # Create .tex tables for gene_to_selection_scan and selection_scan_to_gene
    gene_to_selection_scan_tbl = gt(data = gene_to_selection_scan[,c("Chr", "Gene", "Tissue","GC_P", "Beta", "GWSS_Peak")]) %>% cols_label(
        GC_P = "P-value",
        GWSS_Peak = "GWSS Peak"
      ) %>% tab_style(
        style = cell_text(style = "italic"),
        locations = cells_body(
          columns = "Gene"
        )
      )

    gene_to_selection_scan_tbl

<div id="uomxyklmfg" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#uomxyklmfg .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#uomxyklmfg .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#uomxyklmfg .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#uomxyklmfg .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#uomxyklmfg .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#uomxyklmfg .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uomxyklmfg .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#uomxyklmfg .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#uomxyklmfg .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#uomxyklmfg .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#uomxyklmfg .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#uomxyklmfg .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#uomxyklmfg .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#uomxyklmfg .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#uomxyklmfg .gt_from_md > :first-child {
  margin-top: 0;
}

#uomxyklmfg .gt_from_md > :last-child {
  margin-bottom: 0;
}

#uomxyklmfg .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#uomxyklmfg .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#uomxyklmfg .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#uomxyklmfg .gt_row_group_first td {
  border-top-width: 2px;
}

#uomxyklmfg .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uomxyklmfg .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#uomxyklmfg .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#uomxyklmfg .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uomxyklmfg .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#uomxyklmfg .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#uomxyklmfg .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#uomxyklmfg .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#uomxyklmfg .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#uomxyklmfg .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#uomxyklmfg .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#uomxyklmfg .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#uomxyklmfg .gt_left {
  text-align: left;
}

#uomxyklmfg .gt_center {
  text-align: center;
}

#uomxyklmfg .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#uomxyklmfg .gt_font_normal {
  font-weight: normal;
}

#uomxyklmfg .gt_font_bold {
  font-weight: bold;
}

#uomxyklmfg .gt_font_italic {
  font-style: italic;
}

#uomxyklmfg .gt_super {
  font-size: 65%;
}

#uomxyklmfg .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#uomxyklmfg .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#uomxyklmfg .gt_indent_1 {
  text-indent: 5px;
}

#uomxyklmfg .gt_indent_2 {
  text-indent: 10px;
}

#uomxyklmfg .gt_indent_3 {
  text-indent: 15px;
}

#uomxyklmfg .gt_indent_4 {
  text-indent: 20px;
}

#uomxyklmfg .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Chr">Chr</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Gene">Gene</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="Tissue">Tissue</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="P-value">P-value</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Beta">Beta</th>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="GWSS Peak">GWSS Peak</th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Chr" class="gt_row gt_right">1</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">SLC44A5</td>
<td headers="Tissue" class="gt_row gt_left">Skin Sun Exposed Lower leg</td>
<td headers="GC_P" class="gt_row gt_right">8.827e-06</td>
<td headers="Beta" class="gt_row gt_right">-1.178e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">Novel</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">TMEM163</td>
<td headers="Tissue" class="gt_row gt_left">Kidney Cortex</td>
<td headers="GC_P" class="gt_row gt_right">5.012e-13</td>
<td headers="Beta" class="gt_row gt_right">1.400e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">MAP3K19</td>
<td headers="Tissue" class="gt_row gt_left">Testis</td>
<td headers="GC_P" class="gt_row gt_right">1.490e-07</td>
<td headers="Beta" class="gt_row gt_right">1.127e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">LCT</td>
<td headers="Tissue" class="gt_row gt_left">Cells EBV-transformed lymphocytes</td>
<td headers="GC_P" class="gt_row gt_right">1.078e-09</td>
<td headers="Beta" class="gt_row gt_right">8.239e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">MCM6</td>
<td headers="Tissue" class="gt_row gt_left">Esophagus Muscularis</td>
<td headers="GC_P" class="gt_row gt_right">1.736e-18</td>
<td headers="Beta" class="gt_row gt_right">-1.795e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">DARS</td>
<td headers="Tissue" class="gt_row gt_left">Whole Blood</td>
<td headers="GC_P" class="gt_row gt_right">2.137e-07</td>
<td headers="Beta" class="gt_row gt_right">-9.886e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">2</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">CXCR4</td>
<td headers="Tissue" class="gt_row gt_left">Adrenal Gland</td>
<td headers="GC_P" class="gt_row gt_right">1.457e-05</td>
<td headers="Beta" class="gt_row gt_right">-4.051e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">LCT</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">5</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">P4HA2</td>
<td headers="Tissue" class="gt_row gt_left">Thyroid</td>
<td headers="GC_P" class="gt_row gt_right">9.358e-05</td>
<td headers="Beta" class="gt_row gt_right">1.027e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">PDLIM4</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">5</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">PDLIM4</td>
<td headers="Tissue" class="gt_row gt_left">Brain Cortex</td>
<td headers="GC_P" class="gt_row gt_right">1.098e-06</td>
<td headers="Beta" class="gt_row gt_right">-4.210e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">PDLIM4</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">5</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">SLC22A5</td>
<td headers="Tissue" class="gt_row gt_left">Cells Cultured fibroblasts</td>
<td headers="GC_P" class="gt_row gt_right">1.605e-05</td>
<td headers="Beta" class="gt_row gt_right">-1.283e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">PDLIM4</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">PPP1R18</td>
<td headers="Tissue" class="gt_row gt_left">Artery Aorta</td>
<td headers="GC_P" class="gt_row gt_right">4.166e-05</td>
<td headers="Beta" class="gt_row gt_right">-2.423e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">TUBB</td>
<td headers="Tissue" class="gt_row gt_left">Cells EBV-transformed lymphocytes</td>
<td headers="GC_P" class="gt_row gt_right">1.635e-06</td>
<td headers="Beta" class="gt_row gt_right">4.196e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">PSORS1C1</td>
<td headers="Tissue" class="gt_row gt_left">Thyroid</td>
<td headers="GC_P" class="gt_row gt_right">7.309e-05</td>
<td headers="Beta" class="gt_row gt_right">-7.907e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">CDSN</td>
<td headers="Tissue" class="gt_row gt_left">Skin Sun Exposed Lower leg</td>
<td headers="GC_P" class="gt_row gt_right">4.200e-05</td>
<td headers="Beta" class="gt_row gt_right">4.507e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">CCHCR1</td>
<td headers="Tissue" class="gt_row gt_left">Spleen</td>
<td headers="GC_P" class="gt_row gt_right">9.258e-07</td>
<td headers="Beta" class="gt_row gt_right">8.014e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">APOM</td>
<td headers="Tissue" class="gt_row gt_left">Testis</td>
<td headers="GC_P" class="gt_row gt_right">1.935e-06</td>
<td headers="Beta" class="gt_row gt_right">-2.599e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">C4A</td>
<td headers="Tissue" class="gt_row gt_left">Brain Cerebellum</td>
<td headers="GC_P" class="gt_row gt_right">2.017e-06</td>
<td headers="Beta" class="gt_row gt_right">-8.152e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">ATF6B</td>
<td headers="Tissue" class="gt_row gt_left">Brain Spinal cord cervical c-1</td>
<td headers="GC_P" class="gt_row gt_right">6.076e-06</td>
<td headers="Beta" class="gt_row gt_right">-5.952e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">RNF5</td>
<td headers="Tissue" class="gt_row gt_left">Colon Transverse</td>
<td headers="GC_P" class="gt_row gt_right">1.145e-09</td>
<td headers="Beta" class="gt_row gt_right">-6.436e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">PBX2</td>
<td headers="Tissue" class="gt_row gt_left">Whole Blood</td>
<td headers="GC_P" class="gt_row gt_right">1.059e-06</td>
<td headers="Beta" class="gt_row gt_right">-3.343e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">HLA-DMA</td>
<td headers="Tissue" class="gt_row gt_left">Cells Cultured fibroblasts</td>
<td headers="GC_P" class="gt_row gt_right">8.493e-05</td>
<td headers="Beta" class="gt_row gt_right">5.442e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">6</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">HLA-DPA1</td>
<td headers="Tissue" class="gt_row gt_left">Cells Cultured fibroblasts</td>
<td headers="GC_P" class="gt_row gt_right">2.249e-05</td>
<td headers="Beta" class="gt_row gt_right">-1.380e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">HLA</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">11</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">FADS1</td>
<td headers="Tissue" class="gt_row gt_left">Brain Cerebellum</td>
<td headers="GC_P" class="gt_row gt_right">3.310e-05</td>
<td headers="Beta" class="gt_row gt_right">6.236e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">FADS1</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">12</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">FAM109A</td>
<td headers="Tissue" class="gt_row gt_left">Kidney Cortex</td>
<td headers="GC_P" class="gt_row gt_right">6.886e-05</td>
<td headers="Beta" class="gt_row gt_right">-1.752e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">OAS</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">12</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">OAS3</td>
<td headers="Tissue" class="gt_row gt_left">Cells Cultured fibroblasts</td>
<td headers="GC_P" class="gt_row gt_right">8.779e-06</td>
<td headers="Beta" class="gt_row gt_right">-1.264e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">OAS</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">17</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">NUP85</td>
<td headers="Tissue" class="gt_row gt_left">Brain Cerebellar Hemisphere</td>
<td headers="GC_P" class="gt_row gt_right">7.025e-05</td>
<td headers="Beta" class="gt_row gt_right">1.011e-04</td>
<td headers="GWSS_Peak" class="gt_row gt_left">Novel</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">17</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">GGA3</td>
<td headers="Tissue" class="gt_row gt_left">Breast Mammary Tissue</td>
<td headers="GC_P" class="gt_row gt_right">4.907e-05</td>
<td headers="Beta" class="gt_row gt_right">3.052e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">Novel</td></tr>
    <tr><td headers="Chr" class="gt_row gt_right">17</td>
<td headers="Gene" class="gt_row gt_left" style="font-style: italic;">MRPS7</td>
<td headers="Tissue" class="gt_row gt_left">Brain Cerebellar Hemisphere</td>
<td headers="GC_P" class="gt_row gt_right">2.420e-05</td>
<td headers="Beta" class="gt_row gt_right">-3.548e-05</td>
<td headers="GWSS_Peak" class="gt_row gt_left">Novel</td></tr>
  </tbody>
  
  
</table>
</div>

Figure 2
--------

We plot TWSS-GWSS overlap plots for each significant locus.

FADS region

    scan_regression_overlap(selection_scan, MH_df_whole, 11, c(61700000, 61950000), r_fdr_cut, r_b_cut, "FADS region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-21-1.png)

HLA region

    scan_regression_overlap(selection_scan, MH_df_whole, 6, c(31850000, 32300000), r_fdr_cut, r_b_cut, "HLA region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-22-1.png)

LCT region

    scan_regression_overlap(selection_scan, MH_df_whole, 2, c(134300000, 136400000), r_fdr_cut, r_b_cut, "LCT region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-23-1.png)

OAS region

    scan_regression_overlap(selection_scan, MH_df_whole, 12, c(112600000, 113050000), r_fdr_cut, r_b_cut, "OAS region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-24-1.png)

PDLIM4 region !!!!!!

    i = which(fdr_sig$Gene == "PDLIM4")
    scan_regression_overlap(selection_scan, MH_df_whole, fdr_sig$CHR[i], c(fdr_sig$START[i] - 1e5, fdr_sig$END[i] + 1e5), r_fdr_cut, r_b_cut, "PDLIM4 region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-25-1.png)

Figure 3
--------

We plot TWSS-GWSS overlap, expression against time, and allele frequency
against time for out two novel loci.

SLC44A5 region

    i = which(fdr_sig$Gene == "SLC44A5")
    scan_regression_overlap(selection_scan, MH_df_whole, fdr_sig$CHR[i], c(fdr_sig$START[i] - 1e5, fdr_sig$END[i] + 1e5), r_fdr_cut, r_b_cut, "SLC44A5 region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-26-1.png)

    expression_plot("SLC44A5")

    ## `geom_smooth()` using formula 'y ~ x'

![](main_figures_files/figure-markdown_strict/unnamed-chunk-27-1.png)

    af_plot_paper("SLC44A5")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-28-1.png)

NUP85 region

    i = which(fdr_sig$Gene == "NUP85")
    scan_regression_overlap(selection_scan, MH_df_whole, fdr_sig$CHR[i], c(fdr_sig$START[i] - 1e5, fdr_sig$END[i] + 1e5), r_fdr_cut, r_b_cut, "NUP85 region")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-29-1.png)

    expression_plot("NUP85")

    ## `geom_smooth()` using formula 'y ~ x'

![](main_figures_files/figure-markdown_strict/unnamed-chunk-30-1.png)

    af_plot_paper("NUP85")

![](main_figures_files/figure-markdown_strict/unnamed-chunk-31-1.png)

Figure 4
--------

We begin by plotting gene-level SDS results.

We read in SDS results, filter genes with less than half of the SNPs
included, and impose GC.

    # read file 
    sds = read.table("./data/SDS_results.txt.gz", header = T, sep = "\t")

    # filter genes with less than half of snps included 
    i = which(sds$SNP.SDS/sds$SNP.total > 0.5)
    sds = sds[i,]

    # impose GC on SDS
    p = as.numeric(sds$SDS.p)
    chi_stat = qchisq(p, df = 1, lower.tail=F)
    lambda = median(chi_stat)/qchisq(0.5, df=1) #inflation factor
    corrected_p = pchisq(chi_stat/lambda, df=1, lower.tail=F)
    sds$SDS.GC = corrected_p

    sds$SDS.FDR = p.adjust(sds$SDS.GC, "fdr")

We generate a scatter plot of gene-level SDS and TWSS Betas.

    df = data.frame("Gene" = sds$Gene, "SDS.GC" = -log10(sds$SDS.GC), "GC.p" = -log10(sds$GC), "GC.Beta" = sds$Beta, "SDS" = sds$SDS, "SDS.FDR" = sds$SDS.FDR)
    df = merge(df, MH_df_whole[,c("Gene", "BETA", "FDR")], by = "Gene")

    df$Significant = ifelse(df$SDS.FDR <= 0.05, yes = ifelse(df$FDR <= 0.05, yes = "Both", no = "SDS"), no = ifelse(df$FDR <= 0.05, yes = "TWSS", no = "Neither"))
    df$o = ifelse(df$Significant == "Neither", 1, ifelse(df$Significant == "Both", 3, 2))
    df$t = df$Significant == "Both"

    scatter = ggplot(data = df %>% arrange(o), aes(x = BETA, y = SDS, color = Significant)) + geom_point(size=2.5) + theme_bw() + labs(x = "TWSS Beta", y = "Gene-level SDS", color = "FDR Significant") + scale_color_manual(values = c("Both" = "purple", "SDS" = "red", "TWSS" = "blue", "Neither" = alpha("gray", 0.3))) +theme(text = element_text(size = 20)) + geom_text_repel(data=subset(df, t), fontface = "italic", aes(x = BETA, y = SDS, label = Gene), color = "black", cex = 4.5) + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

    scatter 

![](main_figures_files/figure-markdown_strict/unnamed-chunk-33-1.png)

We generate a qq-plot for the gene-level SDS results.

    p_results = as.numeric(sds$SDS.p)
    p_GC = as.numeric(sds$SDS.GC)
    observed_results = -log10(sort(p_results))
    observed_GC = -log10(sort(p_GC))
    expected = -log10(ppoints(length(p_results)))

    df = data.frame("GC.Observed" = observed_GC,  "Expected" = expected)
    qqlot = ggplot(data = df, aes(x = Expected, y = GC.Observed)) + 
      geom_point(color = "darkblue") + geom_abline(intercept = 0, slope = 1) + labs(cex = 1.4, y = expression("Observed"-log[10](italic(p))), x = expression("Expected"-log[10](italic(p))))+ theme_bw() +  theme(axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

    qqlot

![](main_figures_files/figure-markdown_strict/unnamed-chunk-34-1.png)

We then plot gene-level iHS results.

We read in iHS results, filter genes with less than half of the SNPs
included, and impose GC.

    # read file 
    iHS = read.table("./data/iHS_results.txt.gz", header = T, sep = "\t")
    iHS = na.omit(iHS)

    # filter genes with less than half of snps included 
    i = which(iHS$SNP.iHS/iHS$SNP.total > 0.5)
    iHS = iHS[i,]

    # impose GC on SDS
    p = as.numeric(iHS$iHS.p)
    chi_stat = qchisq(p, df = 1, lower.tail=F)
    lambda = median(chi_stat)/qchisq(0.5, df=1) #inflation factor
    corrected_p = pchisq(chi_stat/lambda, df=1, lower.tail=F)
    iHS$iHS.GC = corrected_p
    iHS$iHS.FDR = p.adjust(iHS$iHS.GC, "fdr")

We generate a scatter plot of gene-level SDS and TWSS Betas.

    df = data.frame("Gene" = iHS$Gene, "iHS.GC" = -log10(iHS$iHS.GC), "GC.p" = -log10(iHS$GC), "GC.Beta" = iHS$Beta, "iHS" = iHS$iHS, "iHS.FDR" = iHS$iHS.FDR)
    df = merge(df, MH_df_whole[,c("Gene", "BETA", "FDR")], by = "Gene")

    df$Significant = ifelse(df$iHS.FDR <= 0.05, yes = ifelse(df$FDR <= 0.05, yes = "Both", no = "iHS"), no = ifelse(df$FDR <= 0.05, yes = "TWSS", no = "Neither"))
    df$o = ifelse(df$Significant == "Neither", 1, ifelse(df$Significant == "Both", 3, 2))
    df$t = df$Significant == "Both"

    scatter = ggplot(data = df %>% arrange(o), aes(x = BETA, y = iHS, color = Significant)) + geom_point(size=2.5) + theme_bw() + labs(x = "TWSS Beta", y = "Gene-level iHS", color = "FDR Significant") + scale_color_manual(values = c("Both" = "purple", "iHS" = "red", "TWSS" = "blue", "Neither" = alpha("gray", 0.3))) +theme(text = element_text(size = 20)) + geom_text_repel(data=subset(df, t), aes(y = iHS, x = BETA, label = Gene, fontface = "italic"), color = "black", cex = 4.5) + coord_cartesian(ylim = c(-11,10), xlim = c(-1.75e-4, 1.75e-4)) + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") + theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

    scatter 

![](main_figures_files/figure-markdown_strict/unnamed-chunk-36-1.png)

We generate a qq-plot for the gene-level iHS results.

    p_results = as.numeric(iHS$iHS.p)
    p_GC = as.numeric(iHS$iHS.GC)
    observed_results = -log10(sort(p_results))
    observed_GC = -log10(sort(p_GC))
    expected = -log10(ppoints(length(p_results)))

    df = data.frame("GC.Observed" = observed_GC,  "Expected" = expected)
    qqlot = ggplot(data = df, aes(x = Expected, y = GC.Observed)) + 
      geom_point(color = "darkblue") + geom_abline(intercept = 0, slope = 1) + labs(cex = 1.4, y = expression("Observed"-log[10](italic(p))), x = expression("Expected"-log[10](italic(p))))+ theme_bw() +  theme(axis.line = element_line(colour = "black"),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

    qqlot

![](main_figures_files/figure-markdown_strict/unnamed-chunk-37-1.png)
