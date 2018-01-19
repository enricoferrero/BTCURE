library(bsseq)

ncores <- parallel::detectCores()
indir <- file.path("../../data/analysis/WGBS/06methylationCalling")
outdir <- file.path("../../data/analysis/WGBS/07differentialMethylation")
dir.create(outdir, showWarnings=FALSE)

smpls <- read.delim("../../data/analysis/WGBS/00samples/samples.txt")
smpls$Name <- paste(smpls$Sample, smpls$Swelling, sep=".")
smpls <- smpls[order(smpls$Sample),]
smpls2mtch <- paste(smpls$Sample, collapse="|")

covs <- list.files(indir, pattern="\\.bismark.cov$", full.names=TRUE, recursive=TRUE)
covs <- covs[grep(smpls2mtch, covs, ignore.case=TRUE)]
names(covs) <- sub("-trimmed-pair1.fastq.gz_bismark_bt2_pe.bismark.cov$", "", basename(covs)) 
covs <- covs[order(names(covs))]

## create BSseq object
if (file.exists(file.path(outdir, "bsobj.rds"))) {
    bsobj <- readRDS(file.path(outdir, "bsobj.rds"))
} else {
    bsobj <- read.bismark(files=covs, sampleNames=names(covs), fileType="cov")
    pData(bsobj) <- smpls
    saveRDS(bsobj, file.path(outdir, "bsobj.rds"))
}

## smoothing
if (file.exists(file.path(outdir, "fit.rds"))) {
    fit <- readRDS(file.path(outdir, "fit.rds"))
} else {
    fit <- BSmooth(bsobj, mc.cores=ncores, parallelBy="sample", verbose=TRUE)
    pData(fit) <- smpls
    saveRDS(fit, file.path(outdir, "fit.rds"))
}

## remove low coverage CpG
if (file.exists(file.path(outdir, "fit2.rds"))) {
    fit2 <- readRDS(file.path(outdir, "fit2.rds"))
} else {
    covr <- getCoverage(fit)
    toKeep <- which(rowSums(covr[, bsobj$Swelling == "High"] >= 2) >= 2 &
                    rowSums(covr[, bsobj$Swelling == "Low"] >= 2) >= 2 &
                    rowSums(covr[, bsobj$Swelling == "None"] >= 2) >= 2)
    fit2 <- fit[toKeep, ]
    pData(fit2) <- smpls
    saveRDS(fit2, file.path(outdir, "fit2.rds"))
}

## export smoothed methylation levels
if (file.exists(file.path(outdir, "meth.rds"))) {
    meth <- readRDS(file.path(outdir, "meth.rds"))
} else {
    meth <- granges(fit2)
    mcols(meth) <- getMeth(fit2, type="smooth")
    names(mcols(meth)) <- smpls$Name
    saveRDS(meth, file.path(outdir, "meth.rds"))
}

## t stat
# High vs Low
if (file.exists(file.path(outdir, "HL.tstat.rds"))) {
    HL.tstat <- readRDS(file.path(outdir, "HL.tstat.rds"))
} else {
    HL.tstat <- BSmooth.tstat(fit2,
                        group1=which(bsobj$Swelling == "High"),
                        group2=which(bsobj$Swelling == "Low"),
                        estimate.var="same",
                        local.correct=TRUE,
                        mc.cores=ncores,
                        verbose=TRUE)
    saveRDS(HL.tstat, file.path(outdir, "HL.tstat.rds"))
}

# High vs None
if (file.exists(file.path(outdir, "HN.tstat.rds"))) {
    HN.tstat <- readRDS(file.path(outdir, "HN.tstat.rds"))
} else {
    HN.tstat <- BSmooth.tstat(fit2,
                        group1=which(bsobj$Swelling == "High"),
                        group2=which(bsobj$Swelling == "None"),
                        estimate.var="same",
                        local.correct=TRUE,
                        mc.cores=ncores,
                        verbose=TRUE)
    saveRDS(HN.tstat, file.path(outdir, "HN.tstat.rds"))
}

# Low vs None
if (file.exists(file.path(outdir, "LN.tstat.rds"))) {
    LN.tstat <- readRDS(file.path(outdir, "LN.tstat.rds"))
} else {
    LN.tstat <- BSmooth.tstat(fit2,
                        group1=which(bsobj$Swelling == "Low"),
                        group2=which(bsobj$Swelling == "None"),
                        estimate.var="same",
                        local.correct=TRUE,
                        mc.cores=ncores,
                        verbose=TRUE)
    saveRDS(LN.tstat, file.path(outdir, "LN.tstat.rds"))
}


## DMR
# set params
cutoff <- 4.6
nloci <- 3
mdiff <- 0.1

# set plotting params
fit2$Swelling <- factor(fit2$Swelling, levels=c("None", "Low", "High"))
fit2$col <- c("#DADAEC", "#807DBA", "#663893")[fit2$Swelling]
fit2$lwd <- 2

# High vs Low
if (file.exists(file.path(outdir, "HL.dmr.rds"))) {
    HL.dmr <- readRDS(file.path(outdir, "HL.dmr.rds"))
} else {
    HL.tstat@stats <- apply(HL.tstat@stats, 2, as.numeric)
    HL.dmr <- dmrFinder(HL.tstat, cutoff=c(-cutoff, cutoff))
    HL.dmr <- subset(HL.dmr, n >= nloci & abs(meanDiff) >= mdiff)
    HL.dmr <- HL.dmr[order(abs(HL.dmr$areaStat), decreasing=TRUE), ]
    saveRDS(HL.dmr, file.path(outdir, "HL.dmr.rds"))
}
write.csv(HL.dmr, file.path(outdir, "HL.dmr.csv"), quote=FALSE, row.names=FALSE)
top <- ifelse(nrow(HL.dmr) > 100, 100, nrow(HL.dmr))
pdf(file.path(outdir, "HL.dmr.pdf"), width=10, height=5)
plotManyRegions(fit2, HL.dmr[1:top, ], extend=2500, addRegions=HL.dmr)
dev.off()

# High vs None
if (file.exists(file.path(outdir, "HN.dmr.rds"))) {
    HN.dmr <- readRDS(file.path(outdir, "HN.dmr.rds"))
} else {
    HN.tstat@stats <- apply(HN.tstat@stats, 2, as.numeric)
    HN.dmr <- dmrFinder(HN.tstat, cutoff=c(-cutoff, cutoff))
    HN.dmr <- subset(HN.dmr, n >= nloci & abs(meanDiff) >= mdiff)
    HN.dmr <- HN.dmr[order(abs(HN.dmr$areaStat), decreasing=TRUE), ]
    saveRDS(HN.dmr, file.path(outdir, "HN.dmr.rds"))
}
write.csv(HN.dmr, file.path(outdir, "HN.dmr.csv"), quote=FALSE, row.names=FALSE)
top <- ifelse(nrow(HN.dmr) > 100, 100, nrow(HN.dmr))
pdf(file.path(outdir, "HN.dmr.pdf"), width=10, height=5)
plotManyRegions(fit2, HN.dmr[1:top, ], extend=2500, addRegions=HN.dmr)
dev.off()

# Low vs None
if (file.exists(file.path(outdir, "LN.dmr.rds"))) {
    LN.dmr <- readRDS(file.path(outdir, "LN.dmr.rds"))
} else {
    LN.tstat@stats <- apply(LN.tstat@stats, 2, as.numeric)
    LN.dmr <- dmrFinder(LN.tstat, cutoff=c(-cutoff, cutoff))
    LN.dmr <- subset(LN.dmr, n >= nloci & abs(meanDiff) >= mdiff)
    LN.dmr <- LN.dmr[order(abs(LN.dmr$areaStat), decreasing=TRUE), ]
    saveRDS(LN.dmr, file.path(outdir, "LN.dmr.rds"))
}
write.csv(LN.dmr, file.path(outdir, "LN.dmr.csv"), quote=FALSE, row.names=FALSE)
top <- ifelse(nrow(LN.dmr) > 100, 100, nrow(LN.dmr))
pdf(file.path(outdir, "LN.dmr.pdf"), width=10, height=5)
plotManyRegions(fit2, LN.dmr[1:top, ], extend=2500, addRegions=LN.dmr)
dev.off()
