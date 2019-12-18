#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(dynwrap)
counts.normalized <- readRDS(args[1])
counts.filtered <- readRDS(args[2])

#counts.norm.rd <- dyndimred::dimred(counts.normalized, method = "landmark_mds", ndim = 5)
#counts.filt.rd <- dyndimred::dimred(counts.filtered, method = "landmark_mds", ndim = 5)

dynready_data <- wrap_expression( counts = counts.filtered, expression = counts.normalized) 
dynready_data <- add_prior_information(dynready_data,start_id = sample(rownames(counts.filtered), 1))

dynutils::write_h5(dynready_data,"ti_ready_dataset.h5")
