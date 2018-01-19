## Importing the DMRs
DMR_path <- "../0_Original_project/data/analysis/WGBS/07differentialMethylation/"

HvL_dmrs <- read.csv(paste0(DMR_path, "HL.dmr.csv"), stringsAsFactors = F)
HvN_dmrs <- read.csv(paste0(DMR_path, "HN.dmr.csv"), stringsAsFactors = F)
LvN_dmrs <- read.csv(paste0(DMR_path, "LN.dmr.csv"), stringsAsFactors = F)

require(GenomicRanges)
HvL_dmrs_gr <- makeGRangesFromDataFrame(df = HvL_dmrs, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "end")
HvN_dmrs_gr <- makeGRangesFromDataFrame(df = HvN_dmrs, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "end")
LvN_dmrs_gr <- makeGRangesFromDataFrame(df = LvN_dmrs, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "end")

#Genomic annotations
require(rtracklayer)
require(ensembldb)
require(AnnotationDbi)
require(GenomicFeatures)
require(org.Hs.eg.db)
require(ChIPseeker)

hg19tohg38 <- import.chain("../0_Original_project/data/analysis/WGBS/00correlations/hg19ToHg38.over.chain")
anno_folder <- "/home/al580162/data"
anno_GRCh37 <- file.path(anno_folder, "Homo_sapiens.GRCh37.75.txdb.db")
anno_GRCh38 <- file.path(anno_folder, "Homo_sapiens.GRCh38.89.txdb.db")

if(file.exists(anno_GRCh38)){
	txdb_GRCh38 <- loadDb(anno_GRCh38)
} else{
	txdb_GRCh38 <- makeTxDbFromGFF(file.path(anno_folder, "Homo_sapiens.GRCh38.89.gtf"), taxonomyId = 9606)
	txdb_GRCh38 <- renameSeqlevels(txdb_GRCh38, gsub("^([0-9XY]+)$", "chr\\1", seqlevels(txdb_GRCh38)))
	saveDb(x = txdb_GRCh38, file = anno_GRCh38)
}
if(file.exists(anno_GRCh37)){
	txdb_GRCh37 <- loadDb(anno_GRCh37)
} else{
	txdb_GRCh37 <- makeTxDbFromGFF(file.path(anno_folder, "Homo_sapiens.GRCh37.75.gtf"), taxonomyId = 9606)
	txdb_GRCh37 <- renameSeqlevels(txdb_GRCh37, gsub("^([0-9XY]+)$", "chr\\1", seqlevels(txdb_GRCh37)))
	saveDb(x = txdb_GRCh37, file = anno_GRCh37)
}

## Promoters: Find the overlaps between the DMRs and the promoters using ChIPSeeker

# Simply defining promoters as x upstream of TSS and Y downstream of TSS.
#require(ChIPseeker)
#HvL_dmrs_gr_anno <- as.GRanges(annotatePeak(HvL_dmrs_gr, tssRegion = c(-2500, 250), TxDb = txdb, annoDb = "org.Hs.eg.db"))
#HvN_dmrs_gr_anno <- as.GRanges(annotatePeak(HvN_dmrs_gr, tssRegion = c(-2500, 250), TxDb = txdb, annoDb = "org.Hs.eg.db"))
#LvN_dmrs_gr_anno <- as.GRanges(annotatePeak(LvN_dmrs_gr, tssRegion = c(-2500, 250), TxDb = txdb, annoDb = "org.Hs.eg.db"))

# Using Ensembl Regulatory Build
promoter_data <- "EnsemblRegulatoryBuild/v89"
ensreg_promoters <- read.csv(file.path(promoter_data, "GRCh38_promoters_biomart.tsv"), stringsAsFactors = F, header = T, sep = "\t")
ensreg_promoters_df <- data.frame(chr = ensreg_promoters[, 1],
				  start = ensreg_promoters[, 2],
				  end = ensreg_promoters[, 3],
				  feature_type = ensreg_promoters[, 4],
				  feature_status = ensreg_promoters[, 6],
				  celltype = gsub("(^.+)\\(((C|V)B)\\)$", "\\1_\\2", ensreg_promoters[, 5]))
ensreg_promoters_df <- unique(ensreg_promoters_df)
ensreg_actprom_df <- ensreg_promoters_df[which(ensreg_promoters_df$feature_status == "ACTIVE"), ]

ensreg_actprom_df$celltype <- gsub("^CD14\\+.+(_(C|V)B)$", "Mon\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("^CD4\\+.+(_(C|V)B)$", "CD4\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("^CD8\\+.+(_(C|V)B)$", "CD8\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("^M0.+(_(C|V)B)$", "Mac0\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("^M1.+(_(C|V)B)$", "Mac1\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("^M2.+(_(C|V)B)$", "Mac2\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("erythro.+(_(C|V)B)$", "Ery\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df$celltype <- gsub("neutro.+(_(C|V)B)$", "Neu\\1", ensreg_actprom_df$celltype)
ensreg_actprom_df <- unique(ensreg_actprom_df)

ensreg_actprom_df <- aggregate(celltype~., ensreg_actprom_df, FUN = toString)
ensreg_actprom_df$celltype <- gsub(", ", ";", ensreg_actprom_df$celltype)

ensreg_actprom_gr <- keepStandardChromosomes(makeGRangesFromDataFrame(df = ensreg_actprom_df, keep.extra.columns = T), pruning.mode = "coarse")
ensreg_actprom_gr <- ensreg_actprom_gr[ensreg_actprom_gr$feature_type == "Promoter", ]

ensreg_actprom_anno <- as.GRanges(annotatePeak(ensreg_actprom_gr, tssRegion = c(-2500, 250), TxDb = txdb_GRCh38, annoDb = "org.Hs.eg.db"))
write.csv(as.data.frame(ensreg_actprom_anno), file.path(promoter_data, "GRCh38_actprom_blood.csv"))

# Find the overlap between the DMRs and the promoters
dm_promoter_overlap <- function(dmrs, promoters){
	dmrs <- keepStandardChromosomes(dmrs, pruning.mode = "coarse")
	promoters <- keepStandardChromosomes(promoters, pruning.mode = "coarse")

	overlapping_ind <- findOverlaps(promoters, dmrs)

	dmrs_df <- data.frame(dmrs)[subjectHits(overlapping_ind), c("seqnames", "start", "end", "n", "invdensity", "areaStat", "maxStat", "meanDiff", "group1.mean", "group2.mean", "tstat.sd", "direction")]
	colnames(dmrs_df) <- c("dmrChr", "dmrStart", "dmrEnd", "nCpGs", "invdensity", "areaStat", "maxStat", "meanDiff", "group1.mean", "group2.mean", "tstat.sd", "direction")

	promoters_df <- data.frame(promoters)[queryHits(overlapping_ind), c("seqnames", "start", "end", "feature_type", "celltype", "distanceToTSS", "geneChr", "geneStart", "geneEnd", "geneId", "SYMBOL")]
	colnames(promoters_df) <- c("promChr", "promStart", "promEnd", "feature_type", "celltype",  "promDistToTSS", "geneChr", "geneStart", "geneEnd", "geneENS", "geneSYMBOL")

	dm_promoters <- cbind(dmrs_df, promoters_df)

	dm_promoters$dmrChr <- gsub("chr", "", dm_promoters$dmrChr)
	dm_promoters$promChr <- gsub("chr", "", dm_promoters$promChr)
	dm_promoters$geneChr <- gsub("chr", "", dm_promoters$geneChr)

	dm_promoters <- aggregate(celltype~., dm_promoters, FUN = toString)
	dm_promoters$celltype <- gsub(", *", ";", dm_promoters$celltype)

	return(dm_promoters)
}

HvL_dmprom <- dm_promoter_overlap(HvL_dmrs_gr, ensreg_actprom_anno) 
HvN_dmprom <- dm_promoter_overlap(HvN_dmrs_gr, ensreg_actprom_anno)
LvN_dmprom <- dm_promoter_overlap(LvN_dmrs_gr, ensreg_actprom_anno)
write.csv(HvL_dmprom, "HvL_dmpromoters.csv")
write.csv(HvN_dmprom, "HvN_dmpromoters.csv")
write.csv(LvN_dmprom, "LvN_dmpromoters.csv")

HvL_dmprom_cells <- table(unlist(strsplit(as.character(HvL_dmprom$celltype), ";")))
HvN_dmprom_cells <- table(unlist(strsplit(as.character(HvN_dmprom$celltype), ";")))
LvN_dmprom_cells <- table(unlist(strsplit(as.character(LvN_dmprom$celltype), ";")))
BG_prom_cells <- table(unlist(strsplit(as.character(ensreg_actprom_anno$celltype), ";")))

groupwise_fishtest <- function(df){
	stats <- matrix(ncol = 5, nrow = ncol(df))
	colnames(stats) <- c("pval", "padj", "lower_CI95", "upper_CI95", "OR")
	rownames(stats) <- colnames(df)
	for(i in 1:ncol(df)){
		fish_res <- fisher.test(data.frame(case = df[,i], control = rowSums(df[, -i])))
		stats[i, "pval"] <- fish_res$p.value
		stats[i, "lower_CI95"] <- fish_res$conf.int[1]
		stats[i, "upper_CI95"] <- fish_res$conf.int[2]
		stats[i, "OR"] <- fish_res$estimate
	}
	stats[,"padj"] <- p.adjust(stats[,"pval"])	
	return(stats)
}

HvLvBG_prom <- groupwise_fishtest(rbind(HvL_dmprom_cells, BG_prom_cells))
HvNvBG_prom <- groupwise_fishtest(rbind(HvN_dmprom_cells, BG_prom_cells))
write.csv(HvLvBG_prom, "HvLvBG_prom_fisher.csv") 
write.csv(HvNvBG_prom, "HvNvBG_prom_fisher.csv")

# Enhancers: Find the overlaps between the DMRs and the enhancers through ActivePromoterEnhancerLinks.tsv 
#PCHiCvrCHiC <- read.csv("Javierre2016/PCHiC_vs_rCHiC_peak_matrix.tsv", sep = "\t")
#PCHiC_sig <- read.csv("Javierre2016/PCHiC_peak_matrix_cutoff5.tsv", sep = "\t")
ActiveProm <- read.csv("Javierre2016/ActivePromoterEnhancerLinks.tsv", sep = "\t")
ActiveProm$baitChr <- gsub("chr", "", ActiveProm$baitChr)
ActiveProm$oeChr <- gsub("chr", "", ActiveProm$oeChr)

celltypes_assayed <- unique(unlist(strsplit(x = as.character(unique(ActiveProm$cellType.s.)), split = ",")))

# We are currently not sure which cell type is the most prominent cell type in the biopsy, we do suspect an infiltration of immune cells and in particular of T-cells. To find out whether a particular immune cell is enriched in our dataset, we will perform a hypergeometric test to see whether our DMRs are enriched in enhancers of a particular immune cell. 
ActivePromOE_gr <- makeGRangesFromDataFrame(df = ActiveProm, keep.extra.columns = T, seqnames.field = "oeChr", start.field = "oeSt", end.field = "oeEnd")
ActivePromOE_anno <- as.GRanges(annotatePeak(ActivePromOE_gr, tssRegion = c(-2500, 250), TxDb = txdb_GRCh37, annoDb = "org.Hs.eg.db"))
ActivePromOE_anno_df <- as.data.frame(ActivePromOE_anno)
ActivePromOE_anno_df_sub <- ActivePromOE_anno_df[,c("seqnames", "start", "end", "baitChr", "baitSt", "baitEnd", "annotation", "distanceToTSS", "geneId", "SYMBOL", "cellType.s.")]
colnames(ActivePromOE_anno_df_sub) <- c("oeChr_GRCh38", "oeStart_GRCh38", "oeEnd_GRCh38", "baitChr_GRCh37", "baitStart_GRCh37", "baitEnd_GRCh37", "oeAnno", "oeDistTSS", "oeENS", "oeSYMBOL", "celltype")

ActivePromBait_gr <- makeGRangesFromDataFrame(df = ActivePromOE_anno_df_sub, keep.extra.columns = T, seqnames.field = "baitChr_GRCh37", start.field = "baitStart_GRCh37", end.field = "baitEnd_GRCh37")
ActivePromBait_anno <- as.GRanges(annotatePeak(ActivePromBait_gr, tssRegion = c(-2500, 250), TxDb = txdb_GRCh37, annoDb = "org.Hs.eg.db"))
ActivePromBait_anno_df <- as.data.frame(ActivePromBait_anno)
ActivePromBait_anno_df_sub <- ActivePromBait_anno_df[,c("seqnames", "start", "end", "annotation", "distanceToTSS", "geneId", "SYMBOL","oeChr_GRCh38", "oeStart_GRCh38", "oeEnd_GRCh38", "oeAnno", "oeDistTSS", "oeENS", "oeSYMBOL", "celltype")]
colnames(ActivePromBait_anno_df_sub) <- c("baitChr", "baitStart", "baitEnd", "baitAnno", "baitDistTSS", "baitENS", "baitSYMBOL", "oeChr", "oeStart", "oeEnd", "oeAnno", "oeDistTSS", "oeENS", "oeSYMBOL", "celltype")

write.csv(ActivePromBait_anno_df_sub, "ActiveProm_anno_GRCh38.csv")

ActiveProm_GRCh38_OE <- makeGRangesFromDataFrame(df = ActivePromBait_anno_df_sub, keep.extra.columns = T, seqnames.field = "oeChr", start.field = "oeStart", end.field = "oeEnd")

#ActiveProm_GRCh38_OE_red <- reduce(ActiveProm_GRCh38_OE)
#ActiveProm_overlap <- data.frame(findOverlaps(ActiveProm_GRCh38_OE_red, ActiveProm_GRCh38_OE))
#ActiveProm_GRCh38_OE_split <- split(ActiveProm_overlap, ActiveProm_overlap$queryHits)


#Find the overlaps between the DMRs and the enhancers

dm_enhancer_overlap <- function(dmrs, enhancers){
	dmrs <- keepStandardChromosomes(dmrs, pruning.mode = "coarse")
	enhancers <- keepStandardChromosomes(enhancers, pruning.mode = "coarse")

	overlapping_ind <- findOverlaps(enhancers, dmrs)
	dmrs_df <- data.frame(dmrs)[subjectHits(overlapping_ind), c("seqnames", "start", "end", "n", "invdensity", "areaStat", "maxStat", "meanDiff", "group1.mean", "group2.mean", "tstat.sd", "direction")]
	colnames(dmrs_df) <- c("dmrChr", "dmrStart", "dmrEnd", "nCpGs", "invdensity", "areaStat", "maxStat", "meanDiff", "group1.mean", "group2.mean", "tstat.sd", "direction")

	enhancers_df <- data.frame(enhancers)[queryHits(overlapping_ind), c("seqnames", "start", "end", "oeAnno", "oeDistTSS", "oeENS", "oeSYMBOL", "baitChr", "baitStart", "baitEnd", "baitAnno", "baitDistTSS", "baitENS", "baitSYMBOL", "celltype")]
	colnames(enhancers_df) <- c("oeChr", "oeStart", "oeEnd", "oeAnno", "oeDistTSS", "oeENS", "oeSYMBOL", "baitChr", "baitStart", "baitEnd", "baitAnno", "baitDistTSS", "baitENS", "baitSYMBOL", "celltype")

	dm_enhancers <- cbind(dmrs_df, enhancers_df)

	dm_enhancers$dmrChr <- gsub("chr", "", dm_enhancers$dmrChr)
	dm_enhancers$oeChr <- gsub("chr", "", dm_enhancers$oeChr)
	dm_enhancers$baitChr <- gsub("chr", "", dm_enhancers$baitChr)

	dm_enhancers <- aggregate(celltype~., dm_enhancers, FUN = toString)
	dm_enhancers$celltype <- gsub(", *", ";", dm_enhancers$celltype)

	return(dm_enhancers)
}

HvL_dmen <- dm_enhancer_overlap(HvL_dmrs_gr, ActiveProm_GRCh38_OE) 
HvN_dmen <- dm_enhancer_overlap(HvN_dmrs_gr, ActiveProm_GRCh38_OE)
LvN_dmen <- dm_enhancer_overlap(LvN_dmrs_gr, ActiveProm_GRCh38_OE)
write.csv(HvL_dmen, "HvL_dmenhancers.csv")
write.csv(HvN_dmen, "HvN_dmenhancers.csv")
write.csv(LvN_dmen, "LvN_dmenhancers.csv")

#Find cell enhancer the enrichments per contrast
HvL_dmen_cells <- table(unlist(strsplit(as.character(HvL_dmen$celltype), ",")))
HvN_dmen_cells <- table(unlist(strsplit(as.character(HvN_dmen$celltype), ",")))
LvN_dmen_cells <- table(unlist(strsplit(as.character(LvN_dmen$celltype), ",")))
LvN_dmen_cells <- c(0, 0, 0, 2, 2, 0, 0, 0, 0)
BG_dmen_cells <- table(unlist(strsplit(as.character(ActiveProm_GRCh38_OE$celltype), ",")))

HvLvBG_perc <- HvL_dmen_cells/BG_dmen_cells*100
HvNvBG_perc <- HvN_dmen_cells/BG_dmen_cells*100
LvNvBG_perc <- LvN_dmen_cells/BG_dmen_cells*100

write.csv(HvLvBG_perc, "HvLvBG_perc.csv")
write.csv(HvNvBG_perc, "HvNvBG_perc.csv")
write.csv(LvNvBG_perc, "LvNvBG_perc.csv")

#Hypergeometric test
HvLvBG_dmen <- rbind(HvL_dmen_cells, BG_dmen_cells)
HvNvBG_dmen <- rbind(HvN_dmen_cells, BG_dmen_cells)
LvNvBG_dmen <- rbind(LvN_dmen_cells, BG_dmen_cells)

HvLvBG_dmen_fisher <- groupwise_fishtest(HvLvBG_dmen)
HvNvBG_dmen_fisher <- groupwise_fishtest(HvNvBG_dmen)
LvNvBG_dmen_fisher <- groupwise_fishtest(LvNvBG_dmen)

write.csv(HvLvBG_dmen_fisher, "HvLvBG_dmen_fisher.csv")
write.csv(HvNvBG_dmen_fisher, "HvNvBG_dmen_fisher.csv")
write.csv(LvNvBG_dmen_fisher, "LvNvBG_dmen_fisher.csv")
