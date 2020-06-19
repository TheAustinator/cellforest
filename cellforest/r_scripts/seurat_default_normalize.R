args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_tenx_directory_path <- args[2]
output_rds_path <- args[3]
min_genes <- as.numeric(args[4])
max_genes <- as.numeric(args[5])
min_cells <- as.numeric(args[6])
perc_mito_cutoff <- as.numeric(args[7])
r_functions_filepath <- args[8]

verbose <- as.logical(args[9])
nfeatures <- as.numeric(args[10])

source(r_functions_filepath)

print(input_metadata_path)

print("creating Seurat object")
seurat_object <- create_seurat_object(input_tenx_directory_path, input_metadata_path, min_cells)
print("reading metadata"); print(date())
metadata <- read_tsv(input_metadata_path)
print("metadata filter"); print(date())
filter_outputs <- metadata_filter_objs(metadata, seurat_object)
seurat_object <- filter_outputs$seurat_object
print("filtering cells"); print(date())
seurat_object <- filter_cells(seurat_object, min_genes, max_genes, perc_mito_cutoff)
print("normalizing"); print(date())
seurat_object <- NormalizeData(seurat_object, verbose=verbose)
print("finding variable features"); print(date())
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures=nfeatures)
print("scaling data"); print(date())
seurat_object <- ScaleData(seurat_object, verbose=verbose)
print("Saving output object"); print(date())
saveRDS(seurat_object, file = output_rds_path)
print("default Seurat normalization DONE"); print(date())

