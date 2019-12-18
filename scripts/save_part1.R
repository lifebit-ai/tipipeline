#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

filt <- readRDS(args[1])
norm <- readRDS(args[2])
clusters.seurat <- readRDS(args[3])

saveRDS(list(filt,norm,clusters.seurat),"part1.rds")