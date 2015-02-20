#-------------------------------------------------------------------------------
expr.sum <- function(xprs)
    {
                                        # xprs is a express outfile
                                        # return a table without the bundle column
      data = read.csv(xprs, header = T, sep = "\t", row.names=2)
      data = data[,-1]
      
      message(c("Total number of sequences: ", nrow(data)))

      message(c("Min FPKM: ", min(data$fpkm)))
      message(c("Max FPKM: ", max(data$fpkm)))
      message(c("Number of FPKMs > 1000: "), sum(data$fpkm > 10000))
      
                                        # Plots
      layout(matrix(c(1,2),2,2,byrow=T))
      
                                        # Density curves of transcripts length and fpkms
      plot(density(data$length), main="Transcripts length density")
      plot(density(log10(data$fpkm)), xlim=c(0,5), "Log10 FPKMs density")
      
#      return(data)

    }                                                                           

#-------------------------------------------------------------------------------
fast.deseq <- function(eff.counts.table)
  {
#    From a effective counts tables (from express and following the guidelines in
#    MY_BEST_PRACTICES. Get the DR analysis
    library(DESeq)

    eff.counts = read.table(eff.counts.table,header=T, row.names = 1)
    cds = newCountDataSet(round(eff.counts), colnames(eff.counts))
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds, method="blind",sharingMode="fit-only")
    #OR to avoid warning message
    cds = estimateDispersions(cds, method="blind",sharingMode="fit-only", fitType="local")

    #Only 1 comparison !
    res = nbinomTest(cds, colnames(eff.counts)[1], colnames(eff.counts)[2])
    res[res$padj < .1]
    
  }

#-------------------------------------------------------------------------------
counts.plot <-
  function ()
  {
    
#    counts = read.table("../pig_trans_counts.txt", header=T, row.names = 1, sep ="\t")
    counts = read.table("~/Work/Junco/analysis/expression/all_counts", header=T, row.names = 1, sep ="\t")
    annot = read.table("~/Work/Junco/analysis/annotations/transcripts_annot.txt", header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
                                        # norm.f is from edgeR output (the .library.tab)
    norm.f = read.table("../../expression/normalization.txt", header=T, row.names = 1)

                                        # Normalize counts
    for (i in colnames(counts)) {
      realname = substring(i,2) # to remove the "X" added to R colnames
      counts[, i] = counts[, i] * norm.f[realname, 3 ]
    }
       
                                        # setting factors for summing expression
    grps = as.factor(sapply(colnames(counts),
      function(x) paste(substr(x,2,3),substr(x,11,11), sep=""), simplify=T))
    print(grps)
    
                                        # transposing the table to be able to use rowsum by groups
    counts = t(counts)
    c.sum = rowsum(counts, grps)

    par(mfrow=c(2,3))

    # COMMENT / UNCOMMENT for subset
    ####################
    subset=c("ENSTGUT00000003832",
      "ENSTGUT00000003914",
      "ENSTGUT00000004909",
      "ENSTGUT00000010894",
      "ENSTGUT00000013433",
      "ENSTGUT00000019087")
    c.sum = c.sum[,colnames(c.sum)%in%subset]
    ####################

    for (i in 1:ncol(c.sum))
      {
        barplot(c.sum[,i],
                main =  colnames(c.sum)[i],
                col=c("chocolate4","chocolate2","black","white","grey","grey","grey","white"),
                xlab = "Tissues",
                ylab = "Normalized counts")

        # For subtitle
        mtext(paste(annot[colnames(c.sum)[i],1], annot[colnames(c.sum)[i],2], sep= " "))
      }
  }
