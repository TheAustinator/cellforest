# Args are: tenx_directory_path, output_rds_path, seurat_metadata_path, min_cells
args <- commandArgs(trailingOnly = TRUE)

tenx_directory_path <- commandArgs(trailingOnly = TRUE)[1]
output_rds_path <- args[2]
seurat_metadata_path <- args[3]
min_cells <- args[4]

r_functions_filepath <- args[5]
source(r_functions_filepath)

seurat_object <- create_seurat_object(tenx_directory_path, seurat_metadata_path, min_cells)
print(date())
print("Saving Seurat object")
saveRDS(seurat_object, output_rds_path)