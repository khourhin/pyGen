                                        # REQUIREMENTS:
## meta_sample.tab:
## col1: sample names, col2 NAMED "condition" ("condition" string is used in model parameters later on)
## SEEMS LIKE SAMPLE NAMES SHOULD BE IN THE SAME ORDER AS countFiles !!!

suppressPackageStartupMessages(library( "DEXSeq" ))
pdf()

# Threads
NTHREADS = MulticoreParam(workers=8)

## TO EDIT
countFiles <- list.files("../dexseq/counts_considered_stranded", full.names=TRUE)
flattenedFile <- "../dexseq/hg19_ensembl_4DEXSeq.gff"
metaSample <- read.table("../dexseq/meta_sample.tab")
##

message("Count files used:")
print(countFiles)
message("MetaData:")
print(metaSample)
message("GFF flattened file used:")
print(flattenedFile)

message("Importing data ...")
dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=metaSample,
    design= ~ sample + exon + condition:exon,
    flattenedfile=flattenedFile
)

head(colData(dxd))
head(counts(dxd))

message("Show how preceding columns are separated according
 to counts in 'this' exons and counts in 'other' exons:")
split(seq_len(ncol(dxd)), colData(dxd)$exon)

message("See the vignette for more exploring commands ...")

## Normalization (same as DESeq, DESeq2)
message("Estimating size factors ...")
dxd <- estimateSizeFactors(dxd)
message("Estimating dispersions ...")
dxd <- estimateDispersions(dxd, BPPARAM=NTHREADS)

## per-exon dispersion estimates versus the mean normalised count
plotDispEsts(dxd)

## Perform test if condition effect
message("Testing DEU ...")
dxd <- testForDEU(dxd, BPPARAM=NTHREADS)

message("Calculate fold change ...")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=NTHREADS)

dxr1 <- DEXSeqResults(dxd)
print(dxr1)

## Number of significant exons? (called genes in vignette) with false discovery rate of 5%
table(dxr1$padj < 0.05)

## How many genes are affected
table(tapply(dxr1$padj<0.05, dxr1$groupID, any))

## Get genes with DR exons
which(tapply(dxr1$padj<0.05, dxr1$groupID, any))

dev.off()
