#! /usr/bin/Rscript

DESeq2.fun <- function(countFile, groupFile, groupChoice, outPrefix){
    library('DESeq2')
    countsData <- read.table(countFile, header=T, check.names=F)
    groups <- read.table(groupFile, header=T, row.names=1)

    dds <- DESeqDataSetFromMatrix(countData = countsData, colData = groups, design = ~ morphtissue)
    dds <- DESeq(dds)

    gpl=levels(groups[,groupChoice])
    
    for (i in gpl){
        gpl = gpl[-1]

        for (j in gpl){
            res <- results(dds, contrast=c("morphtissue",i,j))
            final.res <- res[]
        }
    }

}

##------------------------------------------------------------------------------
## MAIN
##------------------------------------------------------------------------------

## Performing DE analysis with DESeq2

## $1, countFile: the files with counts in, tab delimited
## $2, groupFile: tab delim file with C1: sample name; C2...CN: factor for grouping
## WARNING $2 is now using a file with HEADER !!!! i.e junco_sample_grps.tab2
## $3 groupChoice
## $4, outPrefix: the prefix for the outputs created

args = commandArgs(trailingOnly = TRUE)
DESeq2.fun(args[1], args[2], args[3], args[4])
