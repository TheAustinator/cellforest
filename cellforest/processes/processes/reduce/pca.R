library(cellforestR)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- commandArgs(trailingOnly = TRUE)[1]
input_rds_path <- args[2]
output_embeddings_path <- args[3]
output_loadings_path <- args[4]
output_stdev_path <- args[5]
npcs <- as.numeric(args[6])

#r_functions_filepath <- args[7]
#source(r_functions_filepath)

print("creating Seurat object"); print(date())
srat <- readRDS(input_rds_path)
print("reading metadata"); print(date())
meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1, quote="")
print("metadata filter"); print(date())
srat <- metadata_filter_objs(meta, srat)
print("running pca"); print(date())
srat <- run_pca(srat, output_embeddings_path, output_loadings_path, output_stdev_path, npcs = npcs)
#print("Saving seurat object object"); print(date())
#saveRDS(srat, file = output_rds_path)
print("PCA DONE"); print(date())
