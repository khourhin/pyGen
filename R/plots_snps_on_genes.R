plot.exons <- function(x, rescale, snps){
    print(c(x[4],x[5]))
    if (x[3] == "gene"){
        print("gene")
        rect(as.numeric(x[4])-rescale,-0.25, as.numeric(x[5])-rescale,0.25, col="lightgrey", border="transparent")
    }
    if (x[3] == "exon"){
        print("exon")
        start=as.numeric(x[4])
        end=as.numeric(x[5])
        exon.snps=snps[snps %in% start:end]
        rect(start-rescale,-1, end-rescale,1, col="black")
        points(exon.snps-rescale, rep(2,length(exon.snps)),pch=25, bg="black")
    }    
}

plotSNPs <- function(geneIDs, gtfFile, snpsFile, annotFile=NULL )
{
    # geneIDs: a vector of geneIDs as in the gtf files
    # gtfFile: the corresponding gtf file from Ensembl
    # snpsPos: the file out of SNPs_annotator (SNP position in $3 and Geneid in $4)
    # AnnotFile: a file with unique annotations for genes: i.e ~/Work/Junco/analysis/Zfinch/genes_ensembl.tab

    message("Importing gtf and snps file")
    gtf <- read.table(gtfFile, sep="\t", header=F)
    snps.f<- read.table(snpsFile, sep="\t", header=T)
    
    for (geneID in geneIDs){
        message(paste("Preparing for:", geneID, sep=" "))
                                        # Select all snps for the considered geneID
        snps <- snps.f[snps.f$GID == geneID, 3 ]
        snps <- snps[!is.na(snps)]
                                        # Get the annotation of the geneID
        subset <- grepl(geneID, gtf[,9] )
        gAnnot <- gtf[subset,]
        
        glen = gAnnot[gAnnot[,3] == "gene", 5 ] - gAnnot[gAnnot[,3] == "gene",4 ] + 1
        # For correct position of the SNPs (SNPs pos - gene start +1)
        rescale = gAnnot[gAnnot[,3] == "gene", 4 ]

        pdf(paste(geneID, ".pdf",sep=""))
        plot(c(-10, glen +10 ), c(-10, 10), type = "n")

                                        # The whole gene:
        apply(gAnnot,1, plot.exons, rescale=rescale, snps=snps)

        points(snps-rescale, rep(2,length(snps)),pch=25)
        # Add title
        if (! is.null(annotFile)){
            annot = read.table(annotFile, header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
            mtext(paste(annot[geneID,1], annot[geneID,2], sep= " "))
        }
        dev.off()
    }
}

