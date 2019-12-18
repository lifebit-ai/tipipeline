#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library("scater")
library("DropletUtils")
selected.clusters <- strsplit(args[2],",")[[1]]
list.part1 <- readRDS(args[1])
cluster <- list.part1[[3]]
keep <- cluster %in% selected.clusters
filt.cluster.selected.count <- counts(list.part1[[1]][,keep])
norm.cluster.selected.count <- list.part1[[2]][,keep]
saveRDS(filt.cluster.selected.count,"clusters_filtered.rds")
saveRDS(norm.cluster.selected.count,"clusters_norm.rds")