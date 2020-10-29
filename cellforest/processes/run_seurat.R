library(cellforestR)

args <- commandArgs(trailingOnly <- TRUE)
print(args)
tenx_directory_path <- commandArgs(trailingOnly = TRUE)[1]
output_rds_path <- args[2]
seurat_metadata_path <- args[3]
min_cells <- args[4]
min_detected_genes <- args[5]
max_detected_genes <- args[6]
percent_mito_cutoff <- args[7]
corrected_umi_output_path <- args[8]
pearson_residual_file_path <- args[9]
output_embeddings_path <- args[10]
output_loadings_path <- args[11]
num_pcs <- args[12]
resolution <- as.numeric(args[13])

#r_functions_filepath <- args[14]
#source(r_functions_filepath)

seurat_object <- srat_from_tenx(tenx_directory_path, seurat_metadata_path, min_cells)

seurat_object <- filter_cells(seurat_object, min_detected_genes, max_detected_genes, percent_mito_cutoff)
seurat_object <- run_sctransform(seurat_object, corrected_umi_output_path, pearson_residual_file_path)
seurat_object <- run_pca(seurat_object, output_embeddings_path, output_loadings_path, output_stdev_path, npcs = num_pcs)
seurat_object <- find_clusters(seurat_object, num_pcs, resolution)

print("Saving output object")
print(date())
saveRDS(seurat_object, file = output_rds_path)