library(cellforestR)
library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- commandArgs(trailingOnly = TRUE)[1]
input_rds_path <- args[2]
output_rds_path <- args[3]
min_genes <- as.numeric(args[4])
max_genes <- as.numeric(args[5])
min_cells <- as.numeric(args[6])
perc_mito_cutoff <- as.numeric(args[7])
#r_functions_filepath <- args[8]

verbose <- as.logical(args[8])
nfeatures <- as.numeric(args[9])

#source(r_functions_filepath)


print("creating Seurat object"); print(date())
srat <- readRDS(input_rds_path)
print("reading metadata"); print(date())
meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1, quote="")
print("metadata filter"); print(date())
srat <- metadata_filter_objs(meta, srat)
print("filtering cells"); print(date())
srat <- add_genes_perc_meta(srat, patterns = list("percent.mito" = "^MT-", "percent.ribo" = "^RP[LS]", "percent.hsp" = "^HSP"))
srat <- filter_cells(srat, min_genes, max_genes, perc_mito_cutoff)
print("normalizing"); print(date())
srat <- NormalizeData(srat, verbose = verbose)
print("finding variable features"); print(date())
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = nfeatures)
print("scaling data"); print(date())
srat <- ScaleData(srat, verbose = verbose)
print("Saving output object"); print(date())
saveRDS(srat, file = output_rds_path)
print("default Seurat normalization DONE"); print(date())
