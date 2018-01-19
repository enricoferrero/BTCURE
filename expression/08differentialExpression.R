pval <- 0.05
top <- 24
library(DESeq2)
library(parallel)
library(BiocParallel)
register(MulticoreParam(workers=detectCores()))
library(gtools)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(org.Hs.eg.db)
library(genefilter)
library(gridExtra)

# load counts
cnts <- read.delim("../../data/analysis/RNAseq/07countReads/counts.txt")
names(cnts) <- sub("\\.", "-", names(cnts))
cnts <- cnts[, order(names(cnts))]

# load samples
smpls <- read.delim("../../data/analysis/RNAseq/00samples/samples.txt")
smpls$Name <- paste(smpls$Sample, smpls$Swelling, sep=".")
smpls <- smpls[order(smpls$Sample), ]

# match samples and counts
smpls <- smpls[smpls$Sample %in% names(cnts), ]
cnts <- cnts[, names(cnts) %in% smpls$Sample]

# create DESeqDataSet
dds <- DESeqDataSetFromMatrix(cnts, smpls, design = ~ Sex + Trial + Medication + Swelling)
colnames(dds) <- smpls$Name

# reorder factor levels
dds$Swelling <- factor(dds$Swelling, levels=c("None", "Low", "High"))
dds$Activity <- factor(dds$Activity, levels=c("Low", "Moderate", "High"))
dds$Stage <- factor(dds$Stage, levels=c("Pre", "Early", "Late"))

## differential expression analysis
# pairwise comparisons: perform Wald test
dds <- DESeq(dds, test="Wald", parallel=TRUE)

# export normalised counts
ncnts <- counts(dds, normalized=TRUE)
write.table(ncnts, "../../data/analysis/RNAseq/07countReads/normCounts.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

## data exploration/visualisation
# pairwise comparisons
# dispersion plot
png("../../data/analysis/RNAseq/08differentialExpression/dispersionPlot.png", height=10*150, width=10*150, res=150)
plotDispEsts(dds)
dev.off()

# rlog transformation
rld <- rlog(dds, blind=TRUE)

# clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
sidecols <- c("#DADAEC", "#807DBA", "#663893")[rld$Swelling]

# samples
sampleDists <- as.matrix(dist(t(assay(rld))))
png("../../data/analysis/RNAseq/08differentialExpression/sampleDists.hm.png", height=10*300, width=10*300, res=300)
heatmap.2(sampleDists, trace="none", col=rev(hmcol), ColSideColors=sidecols, RowSideColors=sidecols, tracecol="red3", margin=c(10,10))
dev.off()

sampleCors <- cor(assay(rld))
png("../../data/analysis/RNAseq/08differentialExpression/sampleCors.hm.png", height=11*300, width=10*300, res=300)
heatmap.2(sampleCors, trace="none", col=hmcol, ColSideColors=sidecols, RowSideColors=sidecols, tracecol="red3", margin=c(10,10))
dev.off()

# genes
# grab indices of top 50 genes with highest expression
topGenes <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:50]
png("../../data/analysis/RNAseq/08differentialExpression/topGenes.hm.png", height=10*300, width=10*300, res=300)
heatmap.2(assay(rld)[topGenes,], dendrogram="column", Rowv=FALSE, scale="none", col=hmcol, ColSideColors=sidecols, trace="none", tracecol="red3", margin=c(10,10))
dev.off()

# grab indices of top 50 genes with highest variance
varGenes <- head(order(-rowVars(assay(rld))), 50)
png("../../data/analysis/RNAseq/08differentialExpression/varGenes.hm.png", height=10*300, width=10*300, res=300)
heatmap.2(assay(rld)[varGenes,], scale="row", col=hmcol, ColSideColors=sidecols, trace="none", tracecol="red3", margin=c(10,10))
dev.off()

# pca Batch
pcaData <- plotPCA(rld, intgroup="Batch", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleBatch.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Batch)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca Lane
pcaData <- plotPCA(rld, intgroup="Lane", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleLane.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Lane)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca Sex
pcaData <- plotPCA(rld, intgroup="Sex", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleSex.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Sex)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca Trial
pcaData <- plotPCA(rld, intgroup="Trial", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleTrial.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Trial)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca Patient
pcaData <- plotPCA(rld, intgroup="Patient", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/samplePatient.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Patient)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca SJC
pcaData <- plotPCA(rld, intgroup="SJC", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleSJC.pca.png", height=10*150, width=10*150, res=150)
print(
    ggplot(pcaData, aes(PC1, PC2, fill=SJC)) +
    geom_point(size=8, shape=21, colour="black") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_fill_gradientn(colours=c("#DADAEC", "#807DBA", "#663893"))
    )
dev.off()

# pca Swelling
pcaData <- plotPCA(rld, intgroup="Swelling", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleSwelling.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Swelling)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893"))
	)
dev.off()

# pca DAS28
pcaData <- plotPCA(rld, intgroup="DAS28", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleDAS28.pca.png", height=10*150, width=10*150, res=150)
print(
    ggplot(pcaData, aes(PC1, PC2, fill=DAS28)) +
    geom_point(size=8, shape=21, colour="black") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_fill_gradientn(colours=c("#DADAEC", "#807DBA", "#663893"))
    )
dev.off()

# pca Activity
pcaData <- plotPCA(rld, intgroup="Activity", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleActivity.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Activity)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893"))
	)
dev.off()

# pca Duration
pcaData <- plotPCA(rld, intgroup="Duration", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleDuration.pca.png", height=10*150, width=10*150, res=150)
print(
    ggplot(pcaData, aes(PC1, PC2, fill=Duration)) +
    geom_point(size=8, shape=21, colour="black") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    scale_fill_gradientn(colours=c("#DADAEC", "#807DBA", "#663893"))
    )
dev.off()

# pca Stage
pcaData <- plotPCA(rld, intgroup="Stage", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleStage.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Stage)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893"))
	)
dev.off()

# pca Medication
pcaData <- plotPCA(rld, intgroup="Medication", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleMedication.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Medication)) +
	geom_point(size=8, shape=21, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance"))
	)
dev.off()

# pca Swelling + Sex
pcaData <- plotPCA(rld, intgroup=c("Swelling", "Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleSwellingSex.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Swelling, shape=Sex)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

# pca Swelling + Medication
pcaData <- plotPCA(rld, intgroup=c("Swelling", "Medication"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleSwellingMedication.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Medication, shape=Swelling)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

# pca Activity + Sex
pcaData <- plotPCA(rld, intgroup=c("Activity", "Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleActivitySex.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Activity, shape=Sex)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

# pca Activity + Medication
pcaData <- plotPCA(rld, intgroup=c("Activity", "Medication"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleActivityMedication.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Medication, shape=Activity)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

# pca Stage + Sex
pcaData <- plotPCA(rld, intgroup=c("Stage", "Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleStageSex.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Stage, shape=Sex)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

# pca Stage + Medication
pcaData <- plotPCA(rld, intgroup=c("Stage", "Medication"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("../../data/analysis/RNAseq/08differentialExpression/sampleStageMedication.pca.png", height=10*150, width=10*150, res=150)
print(
	ggplot(pcaData, aes(PC1, PC2, fill=Medication, shape=Stage)) +
	geom_point(size=8, colour="black") +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	scale_shape_manual(values=21:25) +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   )
	)
dev.off()

## get results
# Low vs None
res1 <- results(dds, contrast=c("Swelling", "Low", "None"), parallel=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/Low.vs.None.MA.png", height=10*150, width=10*150, res=150)
plotMA(res1, alpha=pval, main=paste0(sum(res1$padj < pval, na.rm=TRUE), " DEGs (FDR < ", pval, ")"))
dev.off()
countData <- plotCounts(dds, gene=which.min(res1$pvalue), intgroup=c("Swelling", "Sex"), returnData=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/Low.vs.None.topGene.png", height=10*150, width=10*150, res=150)
print(
	ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
	ggtitle(rownames(res1[which.min(res1$pvalue),]))
	)
dev.off()
plots <- list()
for (i in 1:top) {
	idx <- which(res1$pvalue <= sort(res1$pvalue)[top])[i]
	countData <- plotCounts(dds, idx, intgroup=c("Swelling","Sex"), returnData=TRUE)
	plots[[i]] <- ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
		ggtitle(rownames(res1)[idx])
}
plots <- marrangeGrob(grobs=plots, nrow=2, ncol=2)
ggsave("../../data/analysis/RNAseq/08differentialExpression/Low.vs.None.topGenes.pdf", plots, height=14, width=14)
# annotate
res1$foldChange <- logratio2foldchange(res1$log2FoldChange) 
ann <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(res1), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL"))
res1 <- merge(ann, as.data.frame(res1), by.x="ENSEMBL", by.y="row.names")
write.csv(res1, "../../data/analysis/RNAseq/08differentialExpression/LN.deg.csv", row.names=FALSE, quote=FALSE)

# High vs None
res2 <- results(dds, contrast=c("Swelling", "High", "None"), parallel=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/High.vs.None.MA.png", height=10*150, width=10*150, res=150)
plotMA(res2, alpha=pval, main=paste0(sum(res2$padj < pval, na.rm=TRUE), " DEGs (FDR < ", pval, ")"))
dev.off()
countData <- plotCounts(dds, gene=which.min(res2$pvalue), intgroup=c("Swelling", "Sex"), returnData=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/High.vs.None.topGene.png", height=10*150, width=10*150, res=150)
print(
	ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
	ggtitle(rownames(res2[which.min(res2$pvalue),]))
	)
dev.off()
plots <- list()
for (i in 1:top) {
	idx <- which(res2$pvalue <= sort(res2$pvalue)[top])[i]
	countData <- plotCounts(dds, idx, intgroup=c("Swelling","Sex"), returnData=TRUE)
	plots[[i]] <- ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
		ggtitle(rownames(res2)[idx])
}
plots <- marrangeGrob(grobs=plots, nrow=2, ncol=2)
ggsave("../../data/analysis/RNAseq/08differentialExpression/High.vs.None.topGenes.pdf", plots, height=14, width=14)
# annotate
res2$foldChange <- logratio2foldchange(res2$log2FoldChange) 
ann <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(res2), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL"))
res2 <- merge(ann, as.data.frame(res2), by.x="ENSEMBL", by.y="row.names")
write.csv(res2, "../../data/analysis/RNAseq/08differentialExpression/HN.deg.csv", row.names=FALSE, quote=FALSE)

# High vs Low
res3 <- results(dds, contrast=c("Swelling", "High", "Low"), parallel=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/High.vs.Low.MA.png", height=10*150, width=10*150, res=150)
plotMA(res3, alpha=pval, main=paste0(sum(res3$padj < pval, na.rm=TRUE), " DEGs (FDR < ", pval, ")"))
dev.off()
countData <- plotCounts(dds, gene=which.min(res3$pvalue), intgroup=c("Swelling", "Sex"), returnData=TRUE)
png("../../data/analysis/RNAseq/08differentialExpression/High.vs.Low.topGene.png", height=10*150, width=10*150, res=150)
print(
	ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
	ggtitle(rownames(res3[which.min(res3$pvalue),]))
	)
dev.off()
plots <- list()
for (i in 1:top) {
	idx <- which(res3$pvalue <= sort(res3$pvalue)[top])[i]
	countData <- plotCounts(dds, idx, intgroup=c("Swelling","Sex"), returnData=TRUE)
	plots[[i]] <- ggplot(countData, aes(Swelling, count)) +
	geom_smooth(aes(group=1), colour="red4", size=2, se=TRUE, method="loess") +
	geom_point(aes(fill=Swelling, shape=Sex), size=6, colour="black", position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=c("#DADAEC", "#807DBA", "#663893")) +
	scale_shape_manual(values=21:25) +
	scale_y_log10() +
	ylab("Log(count)") +
	guides(
		   fill = guide_legend(override.aes=list(shape=21)),
		   shape = guide_legend(override.aes=list(colour="black", fill="white"))
		   ) +
		ggtitle(rownames(res3)[idx])
}
plots <- marrangeGrob(grobs=plots, nrow=2, ncol=2)
ggsave("../../data/analysis/RNAseq/08differentialExpression/High.vs.Low.topGenes.pdf", plots, height=14, width=14)
# annotate
res3$foldChange <- logratio2foldchange(res3$log2FoldChange) 
ann <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(res3), keytype="ENSEMBL", columns=c("ENTREZID", "SYMBOL"))
res3 <- merge(ann, as.data.frame(res3), by.x="ENSEMBL", by.y="row.names")
write.csv(res3, "../../data/analysis/RNAseq/08differentialExpression/HL.deg.csv", row.names=FALSE, quote=FALSE)
