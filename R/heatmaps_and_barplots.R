counts.plot <-
  function (meth="sum", by="tissue")
  {
                                        # meth = "sum" or "mean"
                                        # by = "tissue" or "color"
    
                                        # THIS FUNCTION GOT TMP SETTINGS ALWAYS IN, YOU SHOULD CHECK AND EDIT BEFORE USE
    
    library(pheatmap)

                                        #-------------------------------------------------------------------------------
                                        #IO :
    counts = read.table("~/Work/Junco/analysis/expression/all_counts",
      header=T, row.names = 1, sep ="\t")
    annot = read.table("~/Work/Junco/analysis/annotations/transcripts_annot.txt",
      header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
                                        # norm.f is from edgeR output (the .library.tab)
    norm.f = read.table("~/Work/Junco/analysis/expression/normalization.txt",
      header=T, row.names = 1)
                                        # A file with one Transcript name by line for the subset
    subset = scan("~/Work/Junco/analysis/annotations/unchar_prot_DEx4", what="raw")

                                        #-------------------------------------------------------------------------------
                                        #Subsetting : 
    counts = counts[rownames(counts)%in%subset,]
                                        # Removing porblematic ones (too large diff for heatmap)
    remove = c("ENSTGUT00000016277")
    counts = counts[!rownames(counts)%in%remove,]

                                        #-------------------------------------------------------------------------------
                                        #Normalization :
    for (i in colnames(counts)) {
      realname = substring(i,2) # to remove the "X" added to R colnameso
      counts[, i] = counts[, i] * norm.f[realname, 3 ]
    }

                                        #-------------------------------------------------------------------------------
                                        # setting factors for summing expression
    if (by == "tissue" ) {
      message("Grouping by tissues")
      grps = list(sapply(colnames(counts),
        function(x) paste(substr(x,2,3),substr(x,11,11), sep=""), simplify=T))
      mycol=c("chocolate4","chocolate2","black","white","grey","grey","grey","white")
    }
    else if (by == "color") {
      message("Grouping by colors")
      grps = list(c("brown","Lbrown","black","white","brown","Lbrown",
        "black","white","brown","Lbrown","black","brown","Lbrown","black","white",
        "grey","grey","white","grey","grey","grey","grey","grey","grey","white",
        "grey","grey","grey","white"))
      mycol=c("black","chocolate4","grey","chocolate2","white")
    }
    
                                        # transposing the table to be able to aggregate
    counts = t(counts)
                                        # Sum and Mean by groups 
    c.mean = aggregate(counts, grps, mean ) 
    rownames(c.mean) = c.mean[,1]
    c.mean[,1] = NULL # aggregate outputhave the first colum with the groups
    
    c.sum = aggregate(counts, grps, sum)
    rownames(c.sum) = c.sum[,1]
    c.sum[,1] = NULL # aggregate outputhave the first colum with the groups

    print(rownames(c.sum))
    par(mfrow=c(2,3))

                                        #-------------------------------------------------------------------------------
                                        # Barplots
                                        # Plotting Means or sums ?
    if (meth == "sum"){
      message("Plotting sums of counts")
      pdata = c.sum
    }
    else if (meth == "mean") {
      message("Plotting means of counts")
      pdata = c.mean
    }
    
    for (i in 1:ncol(pdata))
                                        # The first colum of an aggregate is the groups
      {
        barplot(pdata[,i],
                main =  colnames(pdata)[i],
                col = mycol, 
                names.arg = rownames(pdata),
                xlab = "Tissues",
                ylab = "Mean counts")

                                        # For subtitle
        mtext(paste(annot[colnames(pdata)[i],1], annot[colnames(pdata)[i],2], sep= " "))
      }

                                        #-------------------------------------------------------------------------------
                                        # Heatmap

                                        # Normalize genes counts by median counts (like cichlid article)
    for (i in 1:ncol(pdata)) {
      pdata[,i] = pdata[,i] / median(pdata[,i])
    }
                                        # New window
    dev.new()
    pdata = t(as.matrix(pdata))
    pheatmap(pdata)

  }
