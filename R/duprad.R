## Has to use an old version of picard tools (1.119) because the last
## one is only a single .jar file (i.e no more MarkDuplicate.jar)

## For more check:
#https://www.bioconductor.org/packages/devel/bioc/vignettes/dupRadar/inst/doc/dupRadar.html#fraction-of-multimappers-per-gene

library(dupRadar)

nthreads=5 # not sure gaining a lot by threading
maxmem="40g" # nor by adding a a lot of memory
gtf <- "/home/ekornobis/data/genomes/homo_sapiens/ensembl_h19/Homo_sapiens.GRCh37.82.gtf"
stranded <- 0		# '0' (unstranded), '1' (stranded) and '2' (reversely stranded)
paired <- TRUE

mybam <- "/home/ekornobis/analysis/muchardt/2016_07_12/v3/star/bams/SRR3089084_1_pass2_Aligned.sortedByCoord.out.bam"

#bamDuprm <- markDuplicates(dupremover="picard", bam=mybam, rminput=FALSE, path="/opt/picard-tools-1.119", threads=nthreads, maxmem=maxmem)

# Seems slower than picard
#bamDuprm <- markDuplicates(dupremover="bamutil", bam=mybam, rminput=FALSE, path="/usr/local/bin")

dm <- analyzeDuprates(bamDuprm, gtf, stranded, paired, nthreads)
# Used for comparison with a bad data set
dm.bad <- dupRadar_examples$dm.bad


par(mfrow=c(1,2))
duprateExpDensPlot(DupMat=dm)
duprateExpDensPlot(DupMat=dm.bad)

duprateExpBoxplot(DupMat=dm)
duprateExpBoxplot(DupMat=dm.bad)

expressionHist(dm)
expressionHist(dm.bad)

# calculate the fraction of multimappers per gene
dm$mhRate <- (dm$allCountsMulti - dm$allCounts) / dm$allCountsMulti
# how many genes are exclusively covered by multimappers
sum(dm$mhRate == 1, na.rm=TRUE)

# calculate the fraction of multimappers per gene
dm$mhRate <- (dm$allCountsMulti - dm$allCounts) / dm$allCountsMulti
## Frequency of multimapping
hist(dm$mhRate, 
     breaks=50, 
     col="red", 
     main="Frequency of multimapping rates", 
     xlab="Multimapping rate per gene", 
     ylab="Frequency")
