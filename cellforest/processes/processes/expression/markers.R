library(future)
library(dplyr)
library(parallel)
library(Seurat)
library(cellforestR)

options(future.globals.maxSize = 8000 * 1024^2)
plan("multiprocess", workers = detectCores() - 1)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- commandArgs(trailingOnly = TRUE)[1]
output_markers_dir <- args[2]
root_dir <- args[3]
spec_str <- args[4]
logfc_thresh <- as.numeric(args[5])
test <- args[6]


print("loading metadata"); print(date())
meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1, quote="")
print("cellforestR loading seurat object"); print(date())
srat <- cellforest_load(root_dir, spec_str, "cluster")
print("metadata filter"); print(date())
srat <- metadata_filter_objs(meta, srat)
print("Finding cluster markers"); print(date())
cluster_ids <- unique(srat@meta.data$cluster_id)
print(paste0("cluster_ids: ", cluster_ids)); print(date())
#markers_df_list <- list()
markers <- FindAllMarkers(srat, test.use = test, logfc.threshold = logfc_thresh)
#for (value in cluster_ids) {
#  out <- tryCatch(
#    {
#      markers <- FindMarkers(srat, ident.1 = value, test.use = test, logfc.threshold = logfc_thresh)
#      print(paste0("saving markers for cluster ", value)); print(date())
#      filepath <- paste0(output_markers_dir, "/markers_", value, ".tsv")
#      write.table(markers, sep = "\t", file = filepath, quote = FALSE)
#    },
#    error = function(e) {
#      print(paste0("error in cluster ", value))
#      print(e)
#    }
#  )
#}
print("Writing markers"); print(date())
filepath <- paste0(output_markers_dir, "/markers.tsv")
write.table(markers, sep = "\t", file = filepath, quote = FALSE)
print("DONE"); print(date())
