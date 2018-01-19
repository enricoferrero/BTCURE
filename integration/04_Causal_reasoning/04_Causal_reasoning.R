# parameters
ncores <- parallel::detectCores()
data_dir <- file.path("../../data")
output_dir <- file.path("../../output")
eQTM_analyses_dir <- file.path(output_dir, "02_eQTM_analyses")

# load libraries
options(stringsAsFactors=FALSE)
library(CausalR)
#library(dplyr)
options(java.parameters = "-Xmx8g")
#library(xlsx)

### functions
# plot network
plotCausalNetwork <- function(input, res, netw, alpha=1e-5, dlt=1, dset=dataset, path=outpath) {

    require(dplyr)
    require(igraph)
    require(scales)
    
    # process degs
    input <- input %>%
        filter(VALUE != 0) %>%
        rename(name=SYMBOL) %>%
        mutate(shape="circle",
               size=5,
               color=ifelse(VALUE == 1, "brown2",
                            ifelse(VALUE == -1, "cornflowerblue", "grey70"))) %>%
        select(name, shape, size, color)

    # process regulators
    res <- res %>%
        mutate(name=sub("[+-]$", "", rownames(res))) %>%
        filter(`p-value` < alpha,
               `Enrichment p-value` < alpha) %>%
        mutate(shape="square",
               size=round(rescale(Score, c(5, 15))),
               color=ifelse(Regulation == 1, "darkgoldenrod1",
                            ifelse(Regulation == -1, "mediumpurple", "grey70"))) %>%
        select(name, shape, size, color)

    # process network
    netw <- netw %>%
        inner_join(input, by=c("target" = "name")) %>%
        inner_join(res, by=c("source" = "name")) %>%
        select(source, interaction, target)

    if (nrow(input) > 0 & nrow(res) > 0 & nrow(netw) > 0) {

        # create vertices
        vrts <- full_join(res, input) %>%
            filter(name %in% unique(c(netw$source, netw$target))) %>%
            distinct(name)

        # create edges
        edgs <- netw %>%
            filter(source %in% vrts$name & target %in% vrts$name) %>%
            mutate(color=ifelse(interaction == "Activation", "green4",
                                ifelse(interaction == "Inhibition", "red3", "grey40"))) %>%
            select(source, target, color)

        # create graph
        grph <- graph.data.frame(edgs, directed=TRUE, vertices=vrts)
        grph <- simplify(grph, remove.multiple=FALSE, remove.loops=TRUE)
        grph <- delete.vertices(grph, V(grph)[degree(grph) < 1])

        # plot graph
        n <- 35
        png(paste0(path, dset, ".delta", dlt, ".causalNetwork.png"), width=n*150, height=n*150, res=150)
        try(
            plot.igraph(grph,
                        vertex.label.family="sans",
                        vertex.label.font=2,
                        vertex.label.cex=2,
                        vertex.label.color="black",
                        edge.width=7,
                        edge.arrow.size=2,
                        layout=layout.auto),
            silent=FALSE)
        dev.off()

    } else {
        warning("No significant results!")
        return()
    }

}


# save causal network
saveCausalNetwork <- function(input, res, netw, alpha=1e-5, dlt=1, dset=dataset, path=outpath) {

    require(dplyr)

    # process regulators
    res <- res %>%
        mutate(name=sub("[+-]$", "", rownames(res))) %>%
        filter(`p-value` < alpha,
               `Enrichment p-value` < alpha)

    # process degs
    input <- input %>%
        filter(VALUE != 0) %>%
        rename(name=SYMBOL)

    # process network
    netw <- netw %>%
        inner_join(input, by=c("target" = "name")) %>%
        inner_join(res, by=c("source" = "name")) %>%
        select(source, interaction, target)

    if (nrow(input) > 0 & nrow(res) > 0 & nrow(netw) > 0) {

        # create vertices
        vrts <- unique(c(res$name, input$name))

        # create edges
        edgs <- netw %>%
            filter(source %in% vrts & target %in% vrts) %>%
            mutate(interaction=ifelse(interaction == "Activation", "activates",
                                ifelse(interaction == "Inhibition", "inhibits", NA))) %>%
            arrange(source, target)

        # write to file
        write.csv(edgs, paste0(path, dset, ".delta", dlt, ".causalNetwork.csv"), row.names=FALSE)

    } else {
        warning("No significant results!")
        return()
    }

}


# plot all networks
plotAllCausalNetworks <- function(input, res, netw, alpha=1e-5, dlt=1, top=25, dset=dataset, path=outpath) {

    require(dplyr)
    require(igraph)
    require(scales)
    
    # process degs
    input <- input %>%
        filter(VALUE != 0) %>%
        rename(name=SYMBOL) %>%
        mutate(shape="circle",
               size=8,
               color=ifelse(VALUE == 1, "brown2",
                            ifelse(VALUE == -1, "cornflowerblue", "grey70"))) %>%
        select(name, shape, size, color)

    # process regulators
    res <- res %>%
        mutate(name=sub("[+-]$", "", rownames(res))) %>%
        filter(`p-value` < alpha,
               `Enrichment p-value` < alpha) %>%
        mutate(shape="square",
               size=round(rescale(Score, c(8, 20))),
               color=ifelse(Regulation == 1, "darkgoldenrod1",
                            ifelse(Regulation == -1, "mediumpurple", "grey70"))) %>%
        select(name, shape, size, color)

    # process network
    netw <- netw %>%
        inner_join(input, by=c("target" = "name")) %>%
        inner_join(res, by=c("source" = "name")) %>%
        select(source, interaction, target)

    top <- ifelse(nrow(res) > top, top, nrow(res))
    for (i in 1:top) {

        # process network
        netw.i <- netw %>%
            filter(source==res[i, 1])

        # select degs
        input.i <- input %>%
            inner_join(netw.i, by=c("name" = "target")) %>%
            select(name, shape, size, color)

        if (nrow(input) > 0 & nrow(res) > 0 & nrow(netw) > 0) {

            # create vertices
            vrts <- full_join(res[i,], input.i) %>%
                filter(name %in% unique(c(netw.i$source, netw.i$target))) %>%
                distinct(name)

            # create edges
            edgs <- netw.i %>%
                filter(source %in% vrts$name & target %in% vrts$name) %>%
                mutate(color=ifelse(interaction == "Activation", "green4",
                                    ifelse(interaction == "Inhibition", "red3", "grey40"))) %>%
                select(source, target, color)

            # create graph
            grph <- graph.data.frame(edgs, directed=TRUE, vertices=vrts)
            grph <- simplify(grph, remove.multiple=FALSE, remove.loops=TRUE)
            grph <- delete.vertices(grph, V(grph)[degree(grph) < 1])

            # plot graph
            n <- 12
            png(paste0(path, dset, ".delta", dlt, ".", "#", i, ".", res[i, "name"], ".causalNetwork.png"), width=n*150, height=n*150, res=150)
            try(
                plot.igraph(grph,
                            vertex.label.family="sans",
                            vertex.label.font=2,
                            vertex.label.cex=1,
                            vertex.label.color="black",
                            edge.width=4,
                            edge.arrow.size=1,
                            layout=layout.auto),
                silent=FALSE)
            dev.off()

        } else {
            warning("No significant results!")
            return()
        }
    }

}


# export network for Cytoscape
exportCausalNetwork <- function(input, res, netw, alpha=1e-5, dlt=1, dset=dataset, path=outpath) {

    require(dplyr)
    require(igraph)
    require(scales)
    
    # process regulators
    res <- res %>%
        mutate(name=sub("[+-]$", "", rownames(res))) %>%
        filter(`p-value` < alpha,
               `Enrichment p-value` < alpha) %>%
        left_join(input, by=c("name"="SYMBOL")) %>%
        mutate(shape="square",
               size=round(rescale(Score, c(5, 15))),
               color=ifelse(Regulation == 1, "darkgoldenrod",
                            ifelse(Regulation == -1, "mediumpurple", "grey70"))) %>%
        select(name, shape, size, color)
    
    # process degs
    input <- input %>%
        filter(VALUE != 0) %>%
        rename(name=SYMBOL) %>%
        mutate(shape="circle",
               size=5,
               color=ifelse(VALUE == 1, "brown2",
                            ifelse(VALUE == -1, "cornflowerblue", "grey70"))) %>%
        select(name, shape, size, color)

    # process network
    netw <- netw %>%
        inner_join(input, by=c("target" = "name")) %>%
        inner_join(res, by=c("source" = "name")) %>%
        select(source, interaction, target)

    if (nrow(input) > 0 & nrow(res) > 0 & nrow(netw) > 0) {

        # create vertices
        vrts <- full_join(res, input) %>%
            filter(name %in% unique(c(netw$source, netw$target))) %>%
            distinct(name)

        # create edges
        edgs <- netw %>%
            filter(source %in% vrts$name & target %in% vrts$name) %>%
            mutate(color=ifelse(interaction == "Activation", "green4",
                                ifelse(interaction == "Inhibition", "red3", "grey40"))) %>%
            select(source, target, color)

        write.table(edgs, paste0(path, dset, ".delta", dlt, ".cytoscapeNetwork.txt"), sep="\t", quote=FALSE, row.names=FALSE)
        write.table(vrts, paste0(path, dset, ".delta", dlt, ".cytoscapeTable.txt"), sep="\t", quote=FALSE, row.names=FALSE)

        } else {
            warning("No significant results!")
            return()
        }

}

### analysis
# import causal network
netw_dir <- file.path("/GWD/bioinfo/projects/RD-TSci-CommonData/all_datasets/interactions/causal/live")

netw <- read.csv(file.path(netw_dir, "causal.H.2.sif"), sep = "\t", header = F)
names(netw) <- c("source", "interaction", "target")
ccg <- CreateCCG(file.path(netw_dir, "causal.H.2.sif"))

# read experimental data in
degs_merged_ccg <- ReadExperimentalData(file.path(eQTM_analyses_dir, "DEG_summary/DEG_truthtable.txt"), ccg)

# causal reasoning
degs_merged_rank <- RankTheHypotheses(ccg, degs_merged_ccg, delta = 1, doParallel = T, numCores = 7)

# save the output
causal_output_dir <- file.path(output_dir, "04_Causal_reasoning_analyses")
dir.create(causal_output_dir)

write.csv(degs_merged_rank, file.path(causal_output_dir, "DEGs_Causal_inference.csv"))
saveRDS(degs_merged_rank, file.path(causal_output_dir, "degs_merged_rank.rds"))
