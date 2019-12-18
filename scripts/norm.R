#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(scran)

#counts.filtered <- readRDS(args[1])
sce <- readRDS(args[1])
#sce <- SingleCellExperiment(list(counts=counts.filtered))
sce <- computeSumFactors(sce)
sce <- scater::normalize(sce)
counts.normalized <- logcounts(sce)

saveRDS(counts.normalized,file = "norm.rds")