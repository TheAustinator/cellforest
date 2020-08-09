args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
output_clusters_path <- args[2]
root_dir <- args[3]
spec_str <- args[4]
num_pcs <- as.numeric(args[5])
resolution <- as.numeric(args[6])
eps <- as.numeric(args[7])

r_functions_filepath <- args[8]
print("eps"); print(eps)
print("functions"); print(r_functions_filepath)
source(r_functions_filepath)
library(cellforestR)

print("loading metadata"); print(date())
meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1)
print("cellforestR loading seurat object"); print(date())
srat <- cellforest_load(root_dir, spec_str, "cluster")
print("metadata filter"); print(date())
srat <- metadata_filter_objs(meta, srat)
print("clustering"); print(date())
srat <- find_clusters(srat, output_clusters_path, num_pcs, resolution, eps)
print("clustering DONE"); print(date())
