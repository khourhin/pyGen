#!/usr/bin/Rscript

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
    # snpsFile: the file out of SNPs_annotator (SNP position in $3 and Geneid in $4)
    # AnnotFile: a file with unique annotations for genes: i.e ~/Work/Junco/analysis/Zfinch/genes_ensembl.tab

    message("Importing gtf and snps file")
    gtf <- read.table(gtfFile, sep="\t", header=F)
    snps.f<- read.table(snpsFile, sep="\t", header=T)
    
    for (geneID in geneIDs){
        message(paste("Preparing for:", geneID, sep=" "))
                                        # Select all snps for the considered geneID
        snps <- snps.f[snps.f$GID == geneID, 3 ]
        snps <- snps[!is.na(snps)]

        # For plotting only genes with SNPs
        if (length(snps) != 0) {
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
}
    

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
        # USAGE
    # (vector) geneIDs: a vector of geneIDs as in the gtf files
    # $1 gtfFile: the corresponding gtf file from Ensembl
    # $2 snpsFile: the file out of SNPs_annotator (SNP position in $3 and Geneid in $4)
    # $3 AnnotFile: a file with unique annotations for genes: i.e ~/Work/Junco/analysis/Zfinch/genes_ensembl.tab

    # Now only plot genes which got snps !!!

# Example:
# /plots_snps_on_genes.R ../demo_data/Zfinch/Taeniopygia_guttata.taeGut3.2.4.77.gtf ../demo_data/Zfinch/SNPs_annotated.out ../demo_data/Zfinch/genes_ensembl.tab

message("WARNING ! THIS SCRIPT HAS TO BE EDITED BEFORE USAGE !!!")

geneIDs=c("ENSTGUG00000000056","ENSTGUG00000000084","ENSTGUG00000000255","ENSTGUG00000000284","ENSTGUG00000000541","ENSTGUG00000000676","ENSTGUG00000000701","ENSTGUG00000001948","ENSTGUG00000002034","ENSTGUG00000002328","ENSTGUG00000002370","ENSTGUG00000002685","ENSTGUG00000003414","ENSTGUG00000003683","ENSTGUG00000003763","ENSTGUG00000003786","ENSTGUG00000003872","ENSTGUG00000003895","ENSTGUG00000004704","ENSTGUG00000006146","ENSTGUG00000006195","ENSTGUG00000006508","ENSTGUG00000006559","ENSTGUG00000006635","ENSTGUG00000006654","ENSTGUG00000007369","ENSTGUG00000007686","ENSTGUG00000007761","ENSTGUG00000008004","ENSTGUG00000008024","ENSTGUG00000008285","ENSTGUG00000008467","ENSTGUG00000008489","ENSTGUG00000008872","ENSTGUG00000008890","ENSTGUG00000009279","ENSTGUG00000009485","ENSTGUG00000009543","ENSTGUG00000009875","ENSTGUG00000010051","ENSTGUG00000010290","ENSTGUG00000010434","ENSTGUG00000010468","ENSTGUG00000010754","ENSTGUG00000010914","ENSTGUG00000011227","ENSTGUG00000011381","ENSTGUG00000011447","ENSTGUG00000011904","ENSTGUG00000012079","ENSTGUG00000012574","ENSTGUG00000012899","ENSTGUG00000012909","ENSTGUG00000013043","ENSTGUG00000013835","ENSTGUG00000015204","ENSTGUG00000015245","ENSTGUG00000016832","ENSTGUG00000017017")

args = commandArgs(trailingOnly = TRUE)
plotSNPs(geneIDs, args[1], args[2], args[3])
