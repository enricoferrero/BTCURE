library(parallel)
library(Rsubread)

# samples
smpls <- read.delim("../../data/analysis/RNAseq/00samples/samples.txt")
smpls <- smpls$Sample

# alignments
bams <- list.files("../../data/analysis/RNAseq/06postAlignment", pattern="\\.bam$", full.names=TRUE)
names(bams) <- sub("Sample_(RA-.+).bam", "\\1", basename(bams))

# match samples and alignments
bams <- bams[names(bams) %in% smpls]
bams <- bams[order(names(bams))]

# get counts
cnts <- featureCounts(bams, annot.ext="../../data/analysis/RNAseq/00genome/Homo_sapiens.GRCh38.78.gtf", isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, strandSpecific=2, countChimericFragments=TRUE, nthreads=detectCores())
print(cnts$stat)
cnts <- cnts$counts
colnames(cnts) <- names(bams)
write.table(cnts, "../../data/analysis/RNAseq/07countReads/counts.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

