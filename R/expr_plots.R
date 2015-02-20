expr.plots <- function(countFile, groupsFile, groupChoice, subsetFile=NULL, annotFile=NULL)
    {
                                        # countFile: the files with counts in, tab delimited
                                        # groupsFile: tab delim file with C1: sample name; C2: factor for grouping1;
                                        # C3: factor for grouping2 etc...
                                        # groupChoice: the number of the group column to use from groupsFile
                                        # (ex: 2 for C2, 3 for C3)
                                        # subsetFile: A file with one Transcript name by line to filter in
                                        # At least MINIMUM !!! two transcript names should be specified
                                        # annotFile: with $1: seqid, $2: annotation # this file can be created in
                                        # ensembl biomart with only transcript id
        library(edgeR)
        library(pheatmap)

        counts <- read.delim(countFile, row.names=1, check.names=F)
        group <- factor(read.delim( groupsFile, header=F )[,groupChoice])
        #annot <- read.table("~/Work/Junco/analysis/annotations/transcripts_annot.txt",
        #                    header=T, sep="\t", quote='"', row.names=1, fill=TRUE)

        if (! is.null(subsetFile)){
            if (file.exists(subsetFile[1])){
                message("Using file for subsetting")
                subset = scan(subsetFile, what="raw")
                print(subset)
            }
            else{
                message(paste("Using this vector as subset:", annotFile, sep=" "))
                subset = subsetFile
            }
        }

        message(c("Using theses groups: ", paste(levels(group), collapse=" ")))
        message(c("Seems to have ", length(group), " libraries in the analysis."))
        

        y <- DGEList(counts=counts, group=group )
        y <- calcNormFactors(y)
#        print(y)


        # FIRST: without log2 transformation
        cpms = cpm(y)
        # Filter in the subset and transpose for aggregate
        cpms = cpms[rownames(cpms)%in%subset,]
        cpms = t(cpms)

        # Before for the bar plot by means
        #cpms = aggregate(cpms, by=list(group), mean)
        # Aggreagte add a fisrt columns with the groups names > let's remove it
        #rownames(cpms) = cpms[,1]
        #cpms = cpms[,-1]

        out.d="hist_heat_out/"
        dir.create(out.d)
        
        par(mfrow=c(2,3))
        head(cpms)
        for (i in 1:ncol(cpms))
            {
                pdf( paste(c(out.d, colnames(cpms)[i], ".pdf"), collapse=""))

                boxplot(cpms[,i]~group,
#                        col = c("grey15","chocolate4","grey","chocolate","white"),
                        col = c("chocolate4","chocolate","grey15","white", "grey", "grey", "grey", "white"),
                        names.arg = rownames(cpms),
                        xlab = "Plumage color",
                        ylab = "CPM",
                        las=3)
                                        # Subtitle with no "0"
                mtext(gsub("0","",colnames(cpms)[i]))

                
                # old barplot
                #barplot(cpms[,i],
                 #       main =  colnames(cpms)[i],
                  #      col = c("black","chocolate4","grey","chocolate","white"), 
                   #     names.arg = rownames(cpms),
                    #    xlab = "Tissues",
                     #   ylab = "Mean cpms")
                
                                        # For subtitle
                if (! is.null(annotFile)){
                    annot = read.table(annotFile, header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
                    gName=as.character(annot[colnames(cpms)[i],1])
                    gDescr=as.character(annot[colnames(cpms)[i],2])
                    gDescr=strsplit(gDescr, "[", fixed=T)[[1]][1]
                    title(main=paste(gName,gDescr, sep= " "))
                }
                dev.off()
            }

        # SECOND: with log2 transformation
        # From edgeR manual (to calculate logCPM)
        cpms = cpm(y, prior.count=2, log=T)

        # Filter in the subset and transpose for aggregate
        cpms = cpms[rownames(cpms)%in%subset,]
        cpms = t(cpms)
        
        cpms = aggregate(cpms, by=list(group), mean)
        # Aggreagte add a fisrt columns with the groups names > let's remove it
        rownames(cpms) = cpms[,1]
        cpms = cpms[,-1]

        if (! is.null(annotFile)){
            annot = read.table(annotFile, header=T, sep="\t", quote='"', row.names=1, fill=TRUE)
            gName=as.character(annot[colnames(cpms),1])

            # Get gene names from ensemvbl codes
            colnames(cpms) = gName
        }
        
        for (i in 1:ncol(cpms)) { cpms[,i] = cpms[,i] / median(cpms[,i]) }

                                        # New window
        cpms = t(as.matrix(cpms))
        pdf( paste(out.d, "heatmap.pdf", collapse=""))
        pheatmap(cpms)
        dev.off()
#        return(cpms)
    }

