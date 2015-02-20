COUNTS.file = commandArgs(trailingOnly = TRUE)[1]
Col.Group = commandArgs(trailingOnly = TRUE)[2]
INFO.file = commandArgs(trailingOnly = TRUE)[3] ## Output basename
library(edgeR)
library(locfit)

COUNTS_RAW=read.table(COUNTS.file, header=T, check.names=FALSE)
INFO=read.table('stdin', header=T)
row.names(INFO)<-INFO$BARCODE
INFO<-INFO[-1]
COUNTS=(COUNTS_RAW[names(COUNTS_RAW) %in% row.names(INFO)])
INFO_IN_COUNTS=as.data.frame(t(as.data.frame(t(INFO))[row.names(INFO) %in% names(COUNTS)]))
group=factor(INFO_IN_COUNTS[,names(INFO_IN_COUNTS)==Col.Group])
y=DGEList(counts=COUNTS,group=group[names(COUNTS)])
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(group)
levels(y$samples$group)
table(INFO_IN_COUNTS[,names(INFO_IN_COUNTS)==Col.Group])

keep <- rowSums(cpm(y)>1) >= min(c(table(INFO_IN_COUNTS[,names(INFO_IN_COUNTS)==Col.Group])[[1]][1],table(INFO_IN_COUNTS[,names(INFO_IN_COUNTS)==Col.Group])[[2]][1]))
#keep <- rowSums(cpm(y)>1) >= min(c(dim(y$samples[y$samples$group==unique(y$samples$group)[1],])[1], dim(y$samples[y$samples$group==unique(y$samples$group)[2],])[1]))
y <- y[keep,]

dim(y)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, contrast=c(-1,1))
Pv=topTags(lrt, n=Inf)

warnings()

pdf(paste(INFO.file, ".cbv.pdf", sep=""))
plotBCV(y)
dev.off()

pdf(paste(INFO.file, ".mds.pdf", sep=""))
plotMDS(y, col=as.numeric(group), labels=colnames(y))
dev.off()

write.table(summary(de <- decideTestsDGE(lrt, p=0.05)), paste(INFO.file, ".summary.tab", sep=""),quote=F)

de <- decideTestsDGE(lrt, p=0.05)
detags <- rownames(y)[as.logical(de)]
detags=rownames(topTags(lrt, n=Inf))
cpm=cpm(y) [detags,]
#cpm2=cpm(y,normalized.lib.size=T) [detags,]
library=y$samples
write.table(library, file=paste(INFO.file, ".library.tab", sep=""), quote=F)

tcpm=cbind(library, t(cpm))
tcpm_sort=tcpm[order(tcpm$group),]
write.table( Pv , file=paste(INFO.file, ".main.tab", sep=""), quote=F)
write.table(tcpm_sort, file=paste(INFO.file, ".cpm.tab", sep=""), quote=F)
