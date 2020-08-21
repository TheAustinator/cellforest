library(future)
library(dplyr)
plan("multiprocess", workers = 16)
options(future.globals.maxSize = 16000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)

input_rds_path <- commandArgs(trailingOnly = TRUE)[1]
rdata_output_dir <- args[2]

min_detected_genes <- as.numeric(args[3])
max_detected_genes <- as.numeric(args[4])
perc_mito_cutoff <- as.numeric(args[5])
num_features <- as.numeric(args[6])
num_pcs <- as.numeric(args[7])

r_functions_filepath <- args[8]
source(r_functions_filepath)

output_rds_path <- paste(rdata_output_dir, "output.rds", sep = "/")
filtered_barcodes_path <- paste(rdata_output_dir, "filtered_cell_barcodes.tsv", sep = "/")

print("Reading RDS")
seurat_object <- readRDS(input_rds_path)
print("Filtering cells")
seurat_object <- filter_cells(seurat_object, min_detected_genes, max_detected_genes, perc_mito_cutoff)
print("Running standard workflow")
seurat_object <- run_standard_workflow(seurat_object, nfeatures = num_features, npcs = num_pcs, verbose = TRUE)
print("Saving RDS")
#saveRDS(seurat_object, output_rds_path)
write.table(colnames(seurat_object), filtered_barcodes_path, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
print("Writing embeddings")
save_embeddings(seurat_object, rdata_output_dir)
