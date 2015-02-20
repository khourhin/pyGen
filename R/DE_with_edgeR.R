#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
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
    colnames(all.deg)[2]="seq_id"
    return(all.deg)
}

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
DE.fun <- function(countFile,groupsFile, groupChoice, test="exact", outFile, logFCthres=0, annotFile=NULL){
                                        # Differential expression calculation with edgeR
                                        # countFile: the files with counts in, tab delimited
                                        # groupsFile: tab delim file with C1: sample name; C2: factor for grouping1;
                                        # C3: factor for grouping2 etc...
                                        # groupChoice: the number of the group column to use from groupsFile
                                        # (ex: 2 for C2, 3 for C3)
                                        # Test= exact or glm
                                        # outFile: the outfile with the DEG
                                        # logFCthres: optional a logFC value for cutof (for ex: 2)
                                        # This option is only used for the last summary table
                                        # annotFile a file with the annotations (tab delim):
                                        # $1: seq (or clus) id
                                        # $2: annots

    
    library(edgeR)
    
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
    sum.table <- reshape(sum.table[,1:3], idvar="seq_id", timevar="comp", direction="wide")

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
