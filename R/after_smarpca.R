after.smartpca <- function(evecFile, popFile){
    #evecFile: The file .evec after smartpca (eigensoft)
    #popFile: the file create before smartpca: $1=ID, $2=group
    popColUsed = 2 # Set which column to use for the groups colorings/labelling

    
    evec = read.table(evecFile, row.names=1)
    groups = read.table(popFile, row.names=1)

    pdf(paste(evecFile, ".pca.pdf", sep=""))
    # Plot first components and attach the group name
    par(mfrow = c(2,2))

    plot(evec[,1],evec[,2], col = as.factor(groups[rownames(evec),popColUsed]), xlab="comp1", ylab="comp2")
#    text(evec[,1],evec[,2], groups[rownames(evec),popColUsed], pos=3)

    plot(evec[,4],evec[,3], col = as.factor(groups[rownames(evec),popColUsed]), xlab="comp4", ylab="comp3")
#    text(evec[,1],evec[,3], groups[rownames(evec),popColUsed], pos=3)

    plot(evec[,2],evec[,3], col = as.factor(groups[rownames(evec),popColUsed]), xlab="comp2", ylab="comp3")
#    text(evec[,2],evec[,3], groups[rownames(evec),popColUsed], pos=3)

    plot(evec[,2],evec[,4], col = as.factor(groups[rownames(evec),popColUsed]), xlab="comp2", ylab="comp4")
#    text(evec[,2],evec[,4], groups[rownames(evec),popColUsed], pos=3)

    dev.off()
    
    return(evec)
    
    
}

# To launch with Rscript
evecFile = commandArgs(trailingOnly = TRUE)[1]
popFile = commandArgs(trailingOnly = TRUE)[2]
after.smartpca(evecFile, popFile)
