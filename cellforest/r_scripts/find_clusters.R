args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_rds_path <- args[2]
output_rds_path <- args[3]
output_clusters_path <- args[4]
num_pcs <- as.numeric(args[5])
resolution <- as.numeric(args[6])
eps <- as.numeric(args[7])
# This is janky, but R sucks and I can't find a better way
r_functions_filepath <- args[8]
print("eps"); print(eps)
print("functions"); print(r_functions_filepath)
source(r_functions_filepath)

print("metadata filter")
filter_outputs <- metadata_filter(input_metadata_path, input_rds_path)
seurat_object <- filter_outputs$seurat_object
seurat_object <- find_clusters(seurat_object, output_clusters_path, num_pcs, resolution, eps)
print("Saving RDS"); print(date())
saveRDS(seurat_object, file = output_rds_path)
