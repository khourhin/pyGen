#! /usr/bin/Rscript

##-------------------------------------------------------------------------------
exact.test.fun <- function(y, i, j){
                                        # Perform DE exact test    
            y <- calcNormFactors(y)
                                        # plotMDS(y)
                                        # plotBCV(y)
            y <- estimateCommonDisp(y)
            y <- estimateTagwiseDisp(y)

            et <- exactTest(y, pair=c(i,j))
            return(et)
}

##-------------------------------------------------------------------------------
glm.test.fun <- function(y, i, j, group){

    design = model.matrix(~0+group, data=y$samples)
    expr = paste(c("group", i, "-", "group", j), collapse="")
    my.cont = makeContrasts(contrasts = expr, levels=design)
    print(my.cont)

    y <- calcNormFactors(y)
    y <- estimateGLMCommonDisp(y,design)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y, design)

    lrt <- glmLRT(fit, contrast=my.cont)

    return(lrt)
}

##-------------------------------------------------------------------------------
all.groups.test.fun <- function(counts, group, test){
                                        # Launch tests for each group comparison
    gpl = levels(group)
    all.deg = data.frame()
            
    for (i in gpl){
        gpl = gpl[-1]

        for (j in gpl){
            comp <- paste(c(i,"_vs_",j), collapse="")
            message(c("Performing ", test, " test for: ", comp ))

            y <- DGEList(counts=counts, group=group )

# Filter minimum number of cpm (sum of counts > 1 should be superior to min number of replicates for each group
            keep <- rowSums(cpm(y) > 1) >= min(length(group[group==i]), length(group[group==j]))
            y = y[keep,]
            y$samples$lib.size <- colSums(y$counts)

            if (test == "exact"){
                                        # Do the exact test
                et = exact.test.fun(y, i, j)
            } else if (test == "glm"){
                                        # Do glm test
                et = glm.test.fun(y, i, j, group)
            } else {
                message("ERROR with test selection")
            }
            
            summary(de <- decideTestsDGE(et, p=0.05))
            de.table <- et$table[as.logical(de),]

            if (any(as.logical(de))){ #Check if there is any DEG
                de.table <- cbind(comp, row.names(de.table), de.table)
                all.deg <- rbind(all.deg, de.table)

            } else {
                message(c(comp, " did not give any DEG"))
            }
        }
    }
                                        # Add the seqid name to column for reshaping
    colnames(all.deg)[2] = "seq_id"
    return(all.deg)
}

DE.fun <- function(countFile,groupsFile, groupChoice, test="exact", outFile, logFCthres=0, annotFile=NULL){
    
    library(edgeR)
    library(reshape)
    
    counts <- read.delim(countFile, row.names=1, check.names=F)
    group <- factor(read.delim( groupsFile, header=F )[,groupChoice])

    message(c("Using theses groups: ", paste(levels(group), collapse=" ")))
    message(c("Seems to have ", length(group), " samples in the analysis."))

                                        # Perform all tests:
    all.deg <- all.groups.test.fun(counts, group, test)

                                        # To have the table by seqID with comparisons in columns
                                        # Cutof logFC theshold
    sum.table <- all.deg[abs(all.deg[,"logFC"]) >= logFCthres, ]
                                        # Reshape, keeping comp, seqId and FC
    tmp1 <- melt(sum.table[,1:3], c("seq_id","comp"))
    sum.table <- cast(tmp1,  seq_id ~ comp ) 

    if (! is.null(annotFile)){
        annots <-  read.csv(annotFile, row.names=1, header=T, sep="\t")
        # Get annots for the all DEG
        annots1 <-  annots[as.character(all.deg$seq_id), ]
        all.deg <- cbind(all.deg, annots1)
        annots2 <-  annots[as.character(sum.table$seq_id),]
        sum.table <- cbind(sum.table, annots2)
    }
    
     # Output
    write.table(all.deg, paste(c(outFile,"_total.csv"),collapse=""), quote=F, sep="\t", row.names=F)

    write.table(sum.table, paste(c(outFile, "_FC_", logFCthres, "_summary.csv"), collapse=""),
                row.names=F, quote=F, sep="\t", na="")
    return(sum.table)
#    return(all.deg)    
}

##-------------------------------------------------------------------------------
                                        # MAIN
##-------------------------------------------------------------------------------
                                        # USAGE

## Differential expression calculation with edgeR

## $1, countFile: the file with counts in, tab delimited
## $2, groupsFile: tab delim file with C1: sample name; C2...CN: factor for grouping 
###### WARN: groups order should correspond to the order of the ids in count file
## $3, groupChoice: the number of the group column to use from groupsFile
## (ex: 2 for C2, 3 for C3)
## $4, Test= exact or glm
## $5, outFile: the outfile with the DEG
## $6, logFCthres: a logFC value for cutof (for ex: 2)
## This option is only used for the last summary table
## $7, annotFile a file with the annotations (tab delim):
## C1: seq (or clus) id
## C2: annots
                                        # EXAMPLE
## In 
#DE_with_edgeR.R all_counts_junco junco_sample_grps.tab 4 exact myout1 2 transcripts_ensembl.tab
                                        # JOB
args = commandArgs(trailingOnly = TRUE)
DE.fun(args[1], args[2], as.numeric(args[3]), args[4], args[5], as.numeric(args[6]), args[7])
