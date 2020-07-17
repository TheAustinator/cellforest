args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_rds_path <- args[2]
output_rds_path <- args[3]
output_embeddings_path <- args[4]
output_loadings_path <- args[5]
npcs <- as.numeric(args[6])

r_functions_filepath <- args[7]
source(r_functions_filepath)

filter_outputs <- metadata_filter(input_metadata_path, input_rds_path)
seurat_object <- filter_outputs$seurat_object
seurat_object <- run_pca(seurat_object, output_embeddings_path, output_loadings_path)
print("Saving seurat object object"); print(date())
saveRDS(seurat_object, file = output_rds_path)
print("PCA DONE"); print(date())
