library(cellforestR)

args <- commandArgs(trailingOnly = TRUE)

root_dir <- commandArgs(trailingOnly = TRUE)[1]
spec <- args[2] # TO-DO: Add support for passing spec through Rscript
process <- args[3]

seurat_obj <- cellforest_load(root_dir, spec, process)
