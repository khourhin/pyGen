## INTRO:

## GEO Platform (GPL)
## These files describe a particular type of microarray. They are
## annotation files.

## GEO Sample (GSM)
## Files that contain all the data from the use of a single chip. For
## each gene there will be multiple scores including the main one,
## held in the VALUE column.

## GEO Series (GSE)
## Lists of GSM files that together form a single experiment.

## GEO Dataset (GDS)
## These are curated files that hold a summarised
## combination of a GSE file and its GSM files. They contain
## normalised expression levels for each gene from each sample
## (i.e. just the VALUE field from the GSM file).

library(Biobase)
library(GEOquery)

## WARNING !! SEEMS MUCH EASIER WITH GSM accession
## http://bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html

gse <- getGEO("GSE25282", GSEMatrix=TRUE, AnnotGPL=FALSE,)

gpl5175 <- getGEO('GPL5175')

## Does not seem to work
#head(Meta(gset))
#head(Table(gset))

## Get the expression table
ex <- exprs(gse$`GSE25282-GPL5175_series_matrix.txt.gz`)
an <- Table(gpl5175)

## Use ID column as rownames and remove it from table
rownames(an) <- an$ID
an$ID <- NULL

### DOUBLE CHECK THIS !!!!!!!!!!!!!!!!!
choice <- rownames(ex)
total.table <- cbind(an[choice,], ex)
