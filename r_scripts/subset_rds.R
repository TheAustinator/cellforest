# Args are: tenx_directory_path, output_rds_path, seurat_metadata_path, min_cells
args <- commandArgs(trailingOnly = TRUE)

library(Seurat)

input_rds_path <- args[1]
output_rds_path <- args[2]
seurat_metadata_path <- args[3]

input_object <- readRDS(input_rds_path)
metadata <- read.table(seurat_metadata_path, sep = "\t", row.names = 1, header = TRUE)
output_object <- subset(input_object, cells=rownames(metadata))
saveRDS(output_object, output_rds_path)

