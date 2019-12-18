#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(DropletUtils)
# Main
sce.counts <- DropletUtils::read10xCounts(args[1],col.names = T)
# getting rid of ADT data
is.RNA <- rowData(sce.counts)$Type=="Gene Expression" 
#sce.counts <- sce.counts[is.RNA,] # only Gene Expression data
# getting proper names for genes
rownames(sce.counts) <- scater::uniquifyFeatureNames(rowData(sce.counts)$ID, rowData(sce.counts)$Symbol)
saveRDS(sce.counts,file="preprocess.rds")