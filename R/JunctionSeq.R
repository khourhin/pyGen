suppressPackageStartupMessages(library("JunctionSeq"))

juncSeqPipe <- function(sampleMeta, gffFile, countFolder, outPrefix){

    sampleMeta <- read.table(sampleMeta, header=TRUE, stringsAsFactors=FALSE)
    message("Using this metadata")
    message(sampleMeta)
    countFiles <- list.files(countFolder)
    message("Using this CountFiles")
    message(countFiles)


    jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                                   sample.names = sampleMeta[,1],
                                   condition=factor(sampleMeta[,2]),
                                   flat.gff.file = gffFile,
                                   nCores = 1,
                                   analysis.type = "junctionsAndExons" )

    message("Writing results...")
    writeCompleteResults(jscs,outfile.prefix = outPrefix, save.jscs = TRUE)

    message("Plotting results...")
    buildAllPlots(jscs=jscs, outfile.prefix = paste(outPrefix, "plots", sep="/"),
                  use.plotting.device = "png", FDR.threshold = 0.01 )
}

##-------------------------------------------------------------------------------
                                        # MAIN
##-------------------------------------------------------------------------------
message("Starting JunctionSeq analysis...")

args = commandArgs(trailingOnly = TRUE)
juncSeqPipe(args[1], args[2], args[3], args[4])

save.image()

##-------------------------------------------------------------------------------
                                        # USAGE
##-------------------------------------------------------------------------------

## $1: Table with headers:
## col1: sample names(same order as in countFolder)
## col2: conditions
## $2: Flattened (see JunctionSeq manual) GFF file
## $3: Path to folder with counts (see JunctionSeq for counting instructions)
## The path SHOULD NOT end up with "/"
## $4: Prefix for the results

##-------------------------------------------------------------------------------
                                        # EXAMPLE
##-------------------------------------------------------------------------------
