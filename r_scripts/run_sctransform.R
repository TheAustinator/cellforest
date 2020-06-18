args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_tenx_directory_path <- args[2]
output_rds_path <- args[3]
min_genes <- as.numeric(args[4])
max_genes <- as.numeric(args[5])
min_cells <- as.numeric(args[6])
perc_mito_cutoff <- as.numeric(args[7])
r_functions_filepath <- args[8]

output_corrected_umi_path <- args[9]
output_pearson_residual_path <- args[10]

source(r_functions_filepath)

seurat_object <- create_seurat_object(input_tenx_directory_path, input_metadata_path, min_cells)
print("reading metadata"); print(date())
metadata <- read_tsv(input_metadata_path)
print("metadata filter"); print(date())
filter_outputs <- metadata_filter_objs(metadata, seurat_object)
seurat_object <- filter_outputs$seurat_object
print("filtering cells"); print(date())
seurat_object <- filter_cells(seurat_object, min_genes, max_genes, perc_mito_cutoff)
print("running sctransform"); print(date())
seurat_object <- run_sctransform(seurat_object, output_corrected_umi_path, output_pearson_residual_path)
print("Saving output object"); print(date())
saveRDS(seurat_object, file = output_rds_path)
print("sctransform DONE"); print(date())