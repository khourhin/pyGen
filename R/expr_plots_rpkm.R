calc.gene.length <- function(annotFile, genes.with.counts){

    annotF=read.csv(annotFile, sep="\t", header=T, row.names=1)
    gene.length= annotF[4] - annotF[3] + 1
    colnames(gene.length)="length"
    gene.length = gene.length[genes.with.counts,, drop=FALSE ]
    return(gene.length)
}

expr.plots.rpkm <- function(countFile, groupsFile, groupChoice, subsetFile=NULL, annotFile=NULL)
    {
                                        # countFile: the files with counts in, tab delimited
                                        # groupsFile: tab delim file with C1: sample name; C2: factor for grouping1;
                                        # C3: factor for grouping2 etc...
                                        # groupChoice: the number of the group column to use from groupsFile
                                        # (ex: 2 for C2, 3 for C3)
                                        # subsetFile: A file with one Transcript name by line to filter in
                                        # At least MINIMUM !!! two transcript names should be specified
                                        # annotFile (tab separated):
                                        # $1: seqid, $2: gene name $3 description
                                        # $4: transcript start $5: transcript end 
                                        # this file can be created in
                                        # ensembl biomart with only transcript id
        library(edgeR)

        counts <- read.delim(countFile, row.names=1, check.names=F)

        group <- factor(read.delim( groupsFile, header=F )[,groupChoice])
        #annot <- read.table("~/Work/Junco/analysis/annotations/transcripts_annot.txt",
        #                    header=T, sep="\t", quote='"', row.names=1, fill=TRUE)

        if (! is.null(subsetFile)){
            if (file.exists(subsetFile[1])){
                message("Using file for subsetting")
                subset = scan(subsetFile, what="raw")
            }
            else {
                message(paste("Using this vector as subset:", annotFile, sep=" "))
                subset = subsetFile
            }
        }

        message(c("Using theses groups: ", paste(levels(group), collapse=" ")))
        message(c("Seems to have ", length(group), " libraries in the analysis."))
        

        y <- DGEList(counts=counts, group=group )
        y <- calcNormFactors(y)

        # Calculate gene length from annotFile amd RPKMs
        gene.length = calc.gene.length(annotFile, rownames(y$count))
        y$genes = gene.length
        rpkms = rpkm(y)

        # Filter in the subset and transpose for aggregate
        rpkms = rpkms[rownames(rpkms)%in%subset,]
        rpkms = t(rpkms)
        print(subset)
        print(rpkms)
        print(group)

        out.d="hist_heat_out/"
        dir.create(out.d)

        pdf(paste(out.d,"test_out.pdf",sep=""))
        par(mfrow=c(1,1))
        for (i in 1:ncol(rpkms))
            {
                # Uncomment for 1 pdf by gene
#                pdf( paste(c(out.d, colnames(rpkms)[i], ".pdf"), collapse=""))
                boxplot(rpkms[,i]~group,
                        # grouped by colors
#                        col = c("grey15","chocolate4","grey","chocolate","white"),
                        # grouped by tissues
                        col = c("chocolate4","chocolate","grey15","white","grey","grey","grey","white"), 
                        names.arg = rownames(rpkms),
                        ylab = "RPKM",
                        las=3)

                # Subtitle with no "0"
                #mtext(gsub("0","",colnames(rpkms)[i]))
                
                                        # For title
                if (! is.null(annotFile)){
                    annot = read.table(annotFile, header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
                    gName=as.character(annot[colnames(rpkms)[i],1])
                    gDescr=as.character(annot[colnames(rpkms)[i],2])
                    gDescr=strsplit(gDescr, "[", fixed=T)[[1]][1]
                    gLen=as.numeric(annot[colnames(rpkms)[i],4]) - as.numeric(annot[colnames(rpkms)[i],3]) + 1 
                    title(main=paste(gName,gDescr,"\n",gLen,"bp", sep= " "))
                }
                 # Uncomment for 1 pdf by gene
#                dev.off()
            }
        dev.off()
    }

