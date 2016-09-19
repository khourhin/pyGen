suppressPackageStartupMessages(library("JunctionSeq"))

juncSeqPipe <- function(sampleMeta, gffFile, countFolder, outPrefix, Nthreads){

    sampleMeta <- read.table(sampleMeta, header=TRUE, stringsAsFactors=FALSE)
    message("Using this metadata")
    message(sampleMeta)
    countFiles <- list.files(countFolder, full.names=TRUE)
    message("Using this CountFiles")
    message(countFiles)
    Nthreads <- as.integer(Nthreads)

    jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                   sample.names = sampleMeta[,1],
                                   condition=factor(sampleMeta[,2]),
                                   flat.gff.file = gffFile,
                                   nCores = Nthreads,
                                   analysis.type = "junctionsAndExons" )

    message("Writing results...")
    writeCompleteResults(jscs,outfile.prefix = outPrefix, save.jscs = TRUE)

    message("Plotting results...")
    buildAllPlots(jscs=jscs, outfile.prefix = paste(outPrefix, "plots", sep=""),
                  use.plotting.device = "png", FDR.threshold = 0.01 )
}

##-------------------------------------------------------------------------------
                                        # MAIN
##-------------------------------------------------------------------------------
message("Starting JunctionSeq analysis...")

args = commandArgs(trailingOnly = TRUE)
juncSeqPipe(args[1], args[2], args[3], args[4], args[5])

save.image()

##-------------------------------------------------------------------------------
                                        # USAGE
##-------------------------------------------------------------------------------

## $1: Table with headers:
## col1: sample names(same order as in countFolder)
## col2: conditions
## $2: Flattened (see JunctionSeq manual) GFF.gz file
## $3: Path to folder with counts (can be gz)
## (see JunctionSeq for counting instructions) 
## The path SHOULD NOT end up with "/"
## $4: Prefix for the results
## $5: Number of threads

##-------------------------------------------------------------------------------
                                        # EXAMPLE
##-------------------------------------------------------------------------------
## time Rscript ~/Programming/pyGen/R/JunctionSeq.R sample_meta_2_3.tab JunctionSeq_flat_fr_stranded_copiedFromMuchardtest.gff.gz replicates_2_3 juncSeq_rep23 10
