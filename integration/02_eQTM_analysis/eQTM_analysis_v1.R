# Import the DMR enhancer annotationf
require(GenomicRanges)
enhancer_folder <- "DMenhancers/170704"

HvL_dmrs_enhancers <- read.csv(file.path(enhancer_folder, "HvL_dmenhancers.csv"), stringsAsFactors = F)[, -1] 
HvL_dmrs_enhancers$dmrChr <- gsub("chr", "", as.character(HvL_dmrs_enhancers$dmrChr))
HvL_dmrs_enhancers$contrast <- "HvL"
HvL_dmrs_enhancers_gr <- makeGRangesFromDataFrame(HvL_dmrs_enhancers, keep.extra.columns = T, seqnames.field = "dmrChr", start.field = "dmrStart", end.field = "dmrEnd")

LvN_dmrs_enhancers <- read.csv(file.path(enhancer_folder, "LvN_dmenhancers.csv"), stringsAsFactors = F)[, -1]
LvN_dmrs_enhancers$dmrChr <- gsub("chr", "", as.character(LvN_dmrs_enhancers$dmrChr))
LvN_dmrs_enhancers$contrast <- "LvN"
LvN_dmrs_enhancers_gr <- makeGRangesFromDataFrame(LvN_dmrs_enhancers, keep.extra.columns = T, seqnames.field = "dmrChr", start.field = "dmrStart", end.field = "dmrEnd")

HvN_dmrs_enhancers <- read.csv(file.path(enhancer_folder, "HvN_dmenhancers.csv"), stringsAsFactors = F)[, -1]
HvN_dmrs_enhancers$dmrChr <- gsub("chr", "", as.character(HvN_dmrs_enhancers$dmrChr))
HvN_dmrs_enhancers$contrast <- "HvN"

dm_enhancers <- rbind(HvL_dmrs_enhancers, LvN_dmrs_enhancers, HvN_dmrs_enhancers) 
dm_enhancers$celltype <- gsub(", *", ";", dm_enhancers$celltype)
#dm_enhancers <- aggregate(contrast~dmrChr+dmrStart+dmrEnd+nCpGs+baitChr+baitStart+baitEnd+baitAnno+baitDistTSS+baitENS+baitSYMBOL+oeChr+oeStart+oeEnd+oeAnno+oeDistTSS+oeENS+oeSYMBOL+celltype, dm_enhancers, FUN = toString)
#dm_enhancers$contrast <- gsub(", ", ";", dm_enhancers$contrast)
dm_enhancers_gr <- keepStandardChromosomes(makeGRangesFromDataFrame(dm_enhancers, keep.extra.columns = T, seqnames.field = "dmrChr", start.field = "dmrStart", end.field = "dmrEnd"))


#greduction <- function(gr){
#	gr_red <- reduce(gr)
#	gr_overlap <- data.frame(findOverlaps(gr_red, gr))
#}
#dm_enhancers_gr_red <- reduce(dm_enhancers_gr)

# Import the DMR promoter annotations
promoter_folder <- "DMpromoters/170704"

HvL_dmrs_promoters <- read.csv(file.path(promoter_folder, "HvL_dmpromoters.csv"), stringsAsFactors = F)[, -1]
HvL_dmrs_promoters$dmrChr <- gsub("chr", "", as.character(HvL_dmrs_promoters$dmrChr))
HvL_dmrs_promoters$contrast <- "HvL"

HvN_dmrs_promoters <- read.csv(file.path(promoter_folder, "HvN_dmpromoters.csv"), stringsAsFactors = F)[, -1]
HvN_dmrs_promoters$dmrChr<- gsub("chr", "", as.character(HvN_dmrs_promoters$dmrChr))
HvN_dmrs_promoters$contrast <- "HvN"

dm_promoters <- rbind(HvL_dmrs_promoters, HvN_dmrs_promoters)
dm_promoters$celltype <- gsub(", *", ";", dm_promoters$celltype)
#dm_promoters <- aggregate(contrast~seqnames+start+end+n+distanceToTSS+SYMBOL+feature_type+celltype+geneId, dm_promoters, FUN = toString)
#dm_promoters$contrast <- gsub(", ", ";", dm_promoters$contrast)
dm_promoters_gr <- keepStandardChromosomes(makeGRangesFromDataFrame(dm_promoters, keep.extra.columns = T, seqnames.field = "dmrChr", start.field = "dmrStart", end.field = "dmrEnd"))

# Import the DMR peak counts
require(bsseq)
dmr_folder <- "/home/al580162/projects/BTCURE/0_Original_project/data/analysis/WGBS/07differentialMethylation"

fit2 <- readRDS(file.path(dmr_folder, "fit2.rds"))
fit2 <- keepStandardChromosomes(fit2, pruning.mode = "coarse")

dmen_fit2 <- subsetByOverlaps(fit2, dm_enhancers_gr)
dmen_betas <- getMeth(dmen_fit2, type = "raw")
rownames(dmen_betas) <- paste0(seqnames(dmen_fit2), "_", start(dmen_fit2))
colnames(dmen_betas) <- pData(dmen_fit2)$Sample

dmprom_fit2 <- subsetByOverlaps(fit2, dm_promoters_gr)
dmprom_betas <- getMeth(dmprom_fit2, type = "raw")
rownames(dmprom_betas) <- paste0(seqnames(dmprom_fit2), "_", start(dmprom_fit2))
colnames(dmprom_betas) <- pData(dmprom_fit2)$Sample

#Averaging the DMRs
require(parallel)
no_cores <- detectCores()-1

dmr_beta_mean <- function(dmr_gr, cpg_beta){
	cpg_chr <- gsub("([0-9XY]+)_[0-9]+", "\\1", rownames(cpg_beta))
	cpg_start <- gsub("[0-9XY]+_([0-9]+)", "\\1", rownames(cpg_beta))
	dmr_chr <- gsub("chr", "", as.character(seqnames(dmr_gr)))
	betas <- cpg_beta[which(cpg_chr == dmr_chr & cpg_start >= start(dmr_gr) & cpg_start <= end(dmr_gr)),]
	mean_beta <- colMeans(betas, na.rm = T)
	return(mean_beta)
}

#Enhancers
cl <- makeCluster(no_cores)
clusterExport(cl, c("dmr_beta_mean", "dmen_betas"))
dmen_betas_mean <- parLapply(cl = cl, X = dm_enhancers_gr, fun = function(dmr){dmr_beta_mean(dmr_gr = dmr, cpg_beta = dmen_betas)})
stopCluster(cl)

names(dmen_betas_mean) <- paste0(seqnames(dm_enhancers_gr), ":", start(dm_enhancers_gr), "-", end(dm_enhancers_gr))
dmen_betas_mean_df <- do.call(rbind, dmen_betas_mean)

na_dmen_indices <- which(is.na(dmen_betas_mean_df), arr.ind = T)[,1]
if(length(na_dmen_indices) != 0){
	dmen_betas_mean_df <- dmen_betas_mean_df[-na_dmen_indices, ]
	dm_enhancers_gr <- dm_enhancers_gr[-na_dmen_indices, ]
}
dmen_beta_gene <- data.frame(dmrcoords = as.character(rownames(dmen_betas_mean_df)), 
			     enhcoords = as.character(paste0(dm_enhancers_gr$oeChr, ":", dm_enhancers_gr$oeStart, "-", dm_enhancers_gr$oeEnd)),
			     ensembl = as.character(dm_enhancers_gr$baitENS), 
			     symbol = as.character(dm_enhancers_gr$baitSYMBOL), 
			     contrast = as.character(dm_enhancers_gr$contrast), 
			     celltype = as.character(dm_enhancers_gr$celltype))
dmen_beta_gene <- unique(dmen_beta_gene)
dmen_beta_gene <- aggregate(contrast~., dmen_beta_gene, FUN = toString)
dmen_beta_gene$contrast <- gsub(", ", ";", dmen_beta_gene$contrast)

#Promoters
cl <- makeCluster(no_cores)
clusterExport(cl, c("dmr_beta_mean", "dmprom_betas"))
dmprom_betas_mean <- parLapply(cl = cl, X = dm_promoters_gr, fun = function(dmr){dmr_beta_mean(dmr_gr = dmr, cpg_beta = dmprom_betas)})
stopCluster(cl)

names(dmprom_betas_mean) <- paste0(seqnames(dm_promoters_gr), ":", start(dm_promoters_gr), "-", end(dm_promoters_gr))
dmprom_betas_mean_df <- do.call(rbind, dmprom_betas_mean)

na_prom_indices <- which(is.na(dmprom_betas_mean_df), arr.ind = T)[,1]
if(length(na_prom_indices) != 0){
	dmprom_betas_mean_df <- dmprom_betas_mean_df[-na_prom_indices, ]
	dm_promoters_gr <- dm_promoters_gr[-na_prom_indices, ]
}
dmprom_beta_gene <- data.frame(dmrcoords = as.character(rownames(dmprom_betas_mean_df)), 
			       promcoords = as.character(paste0(dm_promoters_gr$promChr, ":", dm_promoters_gr$promStart, "-", dm_promoters_gr$promEnd)),
			       ensembl = as.character(dm_promoters_gr$geneENS), 
			       symbol = as.character(dm_promoters_gr$geneSYMBOL), 
			       contrast = as.character(dm_promoters_gr$contrast), 
			       celltype = as.character(dm_promoters_gr$celltype))
dmprom_beta_gene <- unique(dmprom_beta_gene)
dmprom_beta_gene <- aggregate(contrast~., dmprom_beta_gene, FUN = toString)
dmprom_beta_gene$contrast <- gsub(", ", ";", dmprom_beta_gene$contrast)
dmprom_beta_gene$ensembl <- as.character(dmprom_beta_gene$ensembl)
# Remove the overlapping DMRs
dmprom_beta_gene_split <- split(dmprom_beta_gene, dmprom_beta_gene$ensembl)

# Import the RNA counts
exprs_folder <- "/home/al580162/projects/BTCURE/0_Original_project/data/analysis/RNAseq/07countReads" 

expr_data <- read.csv(file.path(exprs_folder, "counts.txt"), sep = "\t")
colnames(expr_data) <- gsub("\\.", "-", colnames(expr_data))

# Correlate the DMRs with the enhancers
eqtm_func <- function(dmr_data, expr_data, dmr_gene_annotation){
	overlapping_samples <- intersect(colnames(expr_data), colnames(dmr_data)) 
	overlapping_genes <- intersect(rownames(expr_data), as.character(dmr_gene_annotation$ensembl)) 
	dmr_gene_annotation_culled <- unique(dmr_gene_annotation[which(dmr_gene_annotation$ensembl %in% overlapping_genes), ])
	
	correlation_list <- apply(dmr_gene_annotation_culled, MARGIN = 1, FUN = function(dmreg){
		meth_row <- unlist(dmr_data[unlist(dmreg["dmrcoords"]), overlapping_samples])
		expr_row <- unlist(expr_data[unlist(dmreg["ensembl"]), overlapping_samples])
		cor_stats <- cor.test(x = expr_row, y = meth_row, method = "spearman")

		return(c(pval = cor_stats$p.value, rho = cor_stats$estimate))
	})
	correlation_mat <- t(correlation_list)
	colnames(correlation_mat) <- c("pval", "rho")	
	correlation_mat <- data.frame(dmr_gene_annotation_culled, correlation_mat)
	# Found the following on https://stackoverflow.com/questions/31021868/merge-the-rows-in-r-with-the-same-row-name-concatenating-the-content-in-the-colu. The only issue is that it seprarates entries with a comma, which is unhandy if I want to save the table as a .csv.
	correlation_mat <- aggregate(contrast~., correlation_mat, FUN = toString)
	correlation_mat$contrast <- gsub(", ", ";", correlation_mat$contrast)
	correlation_mat <- correlation_mat[order(correlation_mat$pval),]

	return(correlation_mat)
}

dmen_corr <- eqtm_func(dmen_betas_mean_df, expr_data, dmen_beta_gene)
dmen_corr$feature <- "enhancer"
dmprom_corr <- eqtm_func(dmprom_betas_mean_df, expr_data, dmprom_beta_gene)
dmprom_corr$feature <- "promoter"

write.csv(dmen_corr, "dmenhancers_eqtm.csv")
write.csv(dmprom_corr, "dmpromoters_eqtm.csv")

eqtm_plot <- function(dmr_data, expr_data, coordinates, ensembl_id, factors, plot_name, legend_name){
	#dmr_data: A dataframe containing the methylation data with individual methylation features as rows and samples as columns.
	#expr_data: A dataframe containing the expression data with individual expression features as rows and samples as columns.
	#factors: A vector containing the factors used to color the points with. Names must be the same as the ones used for the methylation and expression data.
	overlapping_samples <- intersect(colnames(expr_data), colnames(dmr_data)) 

	factors <- factors[overlapping_samples]
	dmr_data <- unlist(dmr_data[coordinates, overlapping_samples])
	expr_data <- unlist(expr_data[ensembl_id, overlapping_samples])

	plot_df <- data.frame(methylation = dmr_data, transcription = expr_data, groups = factors) 

	require(ggplot2)
	p <- ggplot(plot_df, aes(x = methylation, y = transcription)) + 
		geom_point(aes(fill = groups), color = "black", shape = 21, size = 5, guide = guide_legend(title = legend_name)) + 
		theme_bw() +
		xlab("Methylation") + 
		ylab("Expression") + 
		xlim(0,1) +
		ggtitle(plot_name) + 
		theme(plot.title = element_text(face = "bold"),
			axis.title = element_text(size = 14, face = "bold"),
			axis.text = element_text(size = 12), 
			legend.title = element_text(size = 14, face = "bold"),
			legend.text = element_text(size = 12),
			legend.position = "bottom")

	return(p)
}

named_factors <- pData(fit2)$Swelling
names(named_factors) <- pData(fit2)$Sample
named_sjc <- pData(fit2)$SJC
names(named_sjc) <- pData(fit2)$Sample

# Cairo(file = "ENSG00000153832.pdf", type = "pdf", dpi = 90, bg = "white", width = 800, height = 800, units = "px")
png("ENSG0000015382.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmen_betas_df, expr_data, "2:231545170-231546134", "ENSG00000153832", named_sjc, "ENSG00000153832"))
dev.off()

eqtm_top_plot <- function(dmr_data, expr_data, eqtm_top, factors, number = 1){
	if(length(number) == 1){
		name_plot <- paste0(eqtm_top[number, "coordinates"], "_", eqtm_top[number, "ensembl"], " (", eqtm_top[number, "symbol"], ")")
		png(file = paste0(eqtm_top[number, "ensembl"], ".png"), width = 800, height = 800, res = 150)
		print(eqtm_plot(dmr_data = dmr_data, expr_data = expr_data, coordinates = as.character(eqtm_top[number, "coordinates"]), factors = factors, ensembl_id = as.character(eqtm_top[number, "ensembl"]), name = name_plot))
		dev.off()
	} else{
		lapply(number, function(iterator){ 
			name_plot <- paste0(eqtm_top[iterator, "coordinates"], "_", eqtm_top[iterator, "ensembl"], " (", eqtm_top[iterator, "symbol"], ")")
			png(file = paste0(iterator, "_", eqtm_top[iterator, "ensembl"], ".png"), width = 800, height = 800, res = 150)
			print(eqtm_plot(dmr_data = dmr_data, expr_data = expr_data, coordinates = as.character(eqtm_top[iterator, "coordinates"]), factors = factors, ensembl_id = as.character(eqtm_top[iterator, "ensembl"]), name = name_plot))
			dev.off()
		})
	}
}

eqtm_top_plot(dmr_data = dmen_betas_mean_df, expr_data = expr_data, eqtm_top = dmen_corr, factors = named_sjc, number = 1:10)
eqtm_top_plot(dmr_data = dmprom_betas_mean_df, expr_data = expr_data, eqtm_top = dmprom_corr, factors = named_sjc, number = 1:10)

rho_sig_plot <- function(corr_df, celltypes){
	require(ggExtra)
	require(ggplot2)
	rho_sig_plotter <- function(cell, corr_df){
			cell_df <- corr_df[grep(paste0(cell), corr_df$celltype), c("pval", "rho")]
#			plot_obj <- ggplot(cell_df, aes(x = rho, y = pval)) +
#				geom_point() + 
#				theme_bw() +
#				ggtitle(cell) + 
#				xlim(-1.2, 1.2) + 
#				ylim(-0.2, 1.2) + 
#				geom_rug()
#
#			dir.create("PvRho")
#			png(file = paste0("PvRho/", cell, ".png"), width = 800, height = 800, res = 150)
#			print(plot_obj)
#			dev.off()

			dir.create("Rho")
			hist_obj <- ggplot(cell_df, aes(x = rho)) +
				xlim(-1, 1) + 
				ylim(0, 2.5) +
				geom_histogram(aes(y=..density..),
						binwidth = 0.05,
						colour = "black", 
						fill = "white") +
				geom_density(alpha = 0.2) +
				ggtitle(cell) +
				theme_bw() 
			png(file = paste0("Rho/", cell, ".png"), width = 800, height = 800, res = 150)
			print(hist_obj)
			dev.off()
	}
	lapply(celltypes, function(cell){rho_sig_plotter(cell = cell, corr_df = corr_df)})
}
enhancer_celltypes <- unique(unlist(strsplit(as.character(dmen_corr$celltype), ",")))
rho_sig_plot(corr_df = dmen_corr, celltypes = enhancer_celltypes)
promoter_celltypes <- unique(gsub("_(C|V)B", "", unlist(strsplit(as.character(dmprom_corr$celltype), ";"))))
rho_sig_plot(corr_df = dmprom_corr, celltypes = promoter_celltypes)

plot_df <- dmen_corr %>% separate_rows(celltype, sep = "(,|;)")
png(file = "enhancer_rho_dist.png", width = 800, height = 800, res = 150)
ggplot(plot_df, aes(x = rho)) +
	geom_density(aes(fill = factor(celltype), alpha = 0.05)) +
	xlim(-1,1) +
	theme(legend.position = "bottom") +
	ggtitle("Enhancer rho distribution")
dev.off()

promoter_celltypes <- unique(unlist(strsplit(as.character(dmprom_corr$celltype), ";")))
rho_sig_plot(corr_df = dmprom_corr, celltypes = promoter_celltypes)

#rho_sig_margplot <- function(corr_df){
#	require(tidyr)
#	source("~/scripts/R_marginal_plot/marginal_plot.R")
#	plot_df <- corr_df %>% separate_rows(celltype, sep = "(,|;)")
#	png(file = "Correlation distributioni.png", width = 800, height = 800, res = 150)
#	marginal_plot(data = plot_df, x = rho, y = pval,group = celltype, lm_formula = NULL, pch = 15, cex = 0.5, bw = "nrd")
#	dev.off()
#}
#
#rho_sig_margplot(dmen_corr)

## Which genes are regulated by several promoters/enhancers?
dmprom_corr2 <- dmprom_corr
dmen_corr2 <- dmen_corr
colnames(dmprom_corr2)[2] <- colnames(dmen_corr2)[2] <- "regcoords" 
dmreg_corr <- rbind(dmprom_corr2, dmen_corr2)
table(table(dmreg_corr$symbol))

#Some genes are present twice because the (slightly) same dmr was found twice
dmreg_corr_split <- split(dmreg_corr, dmreg_corr$symbol)
dmreg_corr_split1 <- dmreg_corr_split[lapply(dmreg_corr_split, nrow) > 1]

dmreg_corr_split_td <- lapply(dmreg_corr_split1, function(entry){
	       reg_ens <- paste0(entry$regcoords, entry$ensembl)
	       deduped <- entry[!duplicated(reg_ens), ]
	       deregduped <- aggregate(regcoords~., deduped, FUN = toString)
	       return(deregduped)
})
dmreg_corr_split_td <- dmreg_corr_split_td[lapply(dmreg_corr_split_td, nrow) > 1]

dmreg_corr_split_df <- do.call(rbind, dmreg_corr_split_td)
dmreg_corr_split_df <- aggregate(regcoords~., dmreg_corr_split_df, FUN = toString)
write.csv(dmreg_corr_split_df, "gene_multireg.csv")

#PAPD5
png("PAPD5_1.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmr_data = dmen_betas_mean_df, expr_data = expr_data, coordinates = "16:50203607-50205157", "ENSG00000121274", factors = named_sjc, plot_name = "PAPD5_enhancer", legend_name = "SJC"))
dev.off()
png("PAPD5_2.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmr_data = dmprom_betas_mean_df, expr_data = expr_data, coordinates = "16:50154665-50154991", "ENSG00000121274", factors = named_sjc, plot_name = "PAPD5_promoter", legend_name = "SJC"))
dev.off()

png("MR1.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmr_data = dmen_betas_mean_df, expr_data = expr_data, coordinates = "1:181131123-181132025", "ENSG00000153029", factors = named_sjc, plot_name = "MR1", legend_name = "SJC"))
dev.off()
png("IER5.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmr_data = dmen_betas_mean_df, expr_data = expr_data, coordinates = "1:181131123-181132025", "ENSG00000162783", factors = named_sjc, plot_name = "IERS", legend_name = "SJC"))
dev.off()
#dmreg_corr$chr <- gsub("^(.+):[0-9]+-[0-9]+", "\\1", dmreg_corr$coordinates)
#dmreg_corr$start <- gsub("^.+:([0-9]+)-[0-9]+", "\\1", dmreg_corr$coordinates)
#dmreg_corr$end <- gsub("^.+:[0-9]+-([0-9]+)", "\\1", dmreg_corr$coordinates)
#
#dmreg_corr_gr <- makeGRangesFromDataFrame(df = dmreg_corr, 
#					  keep.extra.columns = T, 
#					  seqnames.field = "chr", 
#					  start.field = "start",
#					  end.field = "end")
#
#prom_genes <- sort(names(which(table(dmreg_corr$symbol) > 1)))
#prom_genes_list <- lapply(prom_genes, function(gene){dmreg_corr[which(dmreg_corr$symbol == gene),]})
#names(prom_genes_list) <- prom_genes
#prom_genes_list[order(unlist(lapply(prom_genes_list, nrow)))]

## Which regulatory regions contain DMRs and affect multiple genes?
dmreg_corr_multigene_split <- split(dmreg_corr, dmreg_corr$regcoords)
dmreg_corr_multigene_split <- dmreg_corr_multigene_split[which(lapply(dmreg_corr_multigene_split, nrow) > 1)]
dmreg_corr_multigene_split[order(unlist(lapply(dmreg_corr_multigene_split, nrow)))]

dmreg_corr_multigene_df <- do.call(rbind, dmreg_corr_multigene_split)
write.csv(dmreg_corr_multigene_df, "regulator_multigene.csv") 

#SLC7A7
png("SLC7A7.png", width = 800, height = 800, res = 150)
print(eqtm_plot(dmr_data = dmen_betas_mean_df, expr_data = expr_data, coordinates = "14:23366500-23366804", "ENSG00000155465", factors = named_sjc, plot_name = "SLC7A7", legend_name = "SJC"))
dev.off()

