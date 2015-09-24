#! /usr/bin/Rscript

##------------------------------------------------------------------------------
plotDisp <- function(cds, outPrefix){
    ## from http://dwheelerau.com/2013/04/15/how-to-use-deseq-to-analyse-rnaseq-data/

    jpeg(paste(c(outPrefix, '_plotDisp.jpg'), collapse=''))
    
    plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)$perGeneDispEsts,
         pch = '.', log='xy', ylab='dispersion', xlab='mean of normalized counts')
    xg = 10^seq( -.5, 5, length.out=300 )    
    lines( xg, fitInfo(cds)$dispFun( xg ), col='red' )

    dev.off()
}

##------------------------------------------------------------------------------
saveRes <- function(res, FDR, outPrefix){
    resSig = res[res$padj < FDR, ]

    write.csv(res, file=paste(c(outPrefix, '_total.csv'), collapse=''))
    write.csv(resSig, file=paste(c(outPrefix, '_FDR_', FDR, '.csv'), collapse=''))
}

##------------------------------------------------------------------------------
DESEQ.fun <- function(countFile, groupFile, groupChoice, outPrefix){
    library('DESeq')
    countsTab <- read.delim(countFile, header=T, check.names=F, row.names=1)
    groups <- read.delim(groupFile, header=F, row.names=1)[colnames(countsTab), groupChoice]
    groups <- factor(groups)

    cds <- newCountDataSet( countsTab, groups)
    cds <- estimateSizeFactors( cds )
    cds <- estimateDispersions( cds ) 

    # Plot dispersion
    plotDisp(cds, outPrefix)

    # Do DE
    res = nbinomTest( cds, '25B', '25V' )

    # Plot log2 fold change
    jpeg(paste(c(outPrefix, '_plotMA.jpg'), collapse=''))
    plotMA(res)
    dev.off()

    # Save data (one complete res, and one only with FDR <= 0.05):
    saveRes(res, 0.05, outPrefix)
}

##------------------------------------------------------------------------------
## MAIN
##------------------------------------------------------------------------------

## WARNING: this gave me last time FDR all equal to one, so switched to DESeq2

## Performing DE analysis with DESeq

## $1, countFile: the files with counts in, tab delimited
## $2, groupFile: tab delim file with C1: sample name; C2...CN: factor for grouping 
## $3, groupChoice: the number of the group column to use from groupsFile
## $4, outPrefix: the prefix for the outputs created

args = commandArgs(trailingOnly = TRUE)
DESEQ.fun(args[1], args[2], as.numeric(args[3]), args[4])
