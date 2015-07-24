#! /usr/bin/Rscript

## Plot a manhattan plot from a file with fsts out of VCFtools

##------------------------------------------------------------------------------
getFsts <- function(fst.tab){
    fsts = read.table(fst.tab, header=T)
    # Remove all lines with "NAs"
    fsts = fsts[complete.cases(fsts),]
    
    last.max = 0
    
    for (chr in unique(fsts$CHROM)){
        
        # Increment with the maximum of the last chromosome (just for plotting)
        fsts[fsts$CHROM==chr,2] = fsts[fsts$CHROM==chr,2] + last.max
        last.max = max(fsts[fsts$CHROM==chr,2])
    }
    
    # Plot
    pdf("test.pdf")    
    plot(fsts$POS, fsts$WEIR_AND_COCKERHAM_FST, ylim=c(-0.1,1), type="n", ylab="FST", xlab="Chromosomes", xaxt="n")

    # Only for colors
    colors=rep(c("grey","lightgrey"), length(unique(fsts$CHROM)))
    n=1
    for (chr in unique(fsts$CHROM)){ 
        rect(min(fsts[fsts$CHROM==chr,2]),
             -0.2,
             max(fsts[fsts$CHROM==chr,2]),
             1.2,
             col=colors[n], border=NA)
        # only for colors
        n=n+1
    }
    points(fsts$POS, fsts$WEIR_AND_COCKERHAM_FST, pch=20, cex=0.5)
    axis(side=1, at = by (fsts[,2], fsts$CHROM,
                     function(x) {mean(c(max(x),min(x)))})[unique(fsts$CHROM)],
         labels=unique(fsts$CHROM), cex.axis=0.5, las=3)

    dev.off()
    print(length(unique(fsts$CHROM)))

}



##------------------------------------------------------------------------------
## MAIN
##------------------------------------------------------------------------------

## USAGE:
## $1: the fst file from VCFtools 

args = commandArgs(trailingOnly = TRUE)
getFsts(args[1])

