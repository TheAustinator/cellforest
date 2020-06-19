library(future)
library(dplyr)

plan("multiprocess", workers = 16)
options(future.globals.maxSize = 8000 * 1024^2)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_rds_path <- args[2]
output_markers_dir <- args[3]
test <- args[4]
logfc_thresh <- as.numeric(args[5])

# This is janky, but R sucks and I can't find a better way
r_functions_filepath <- args[6]
source(r_functions_filepath)

filter_outputs <- metadata_filter(input_metadata_path, input_rds_path)
metadata <- filter_outputs$metadata
seurat_object <- filter_outputs$seurat_object

print("Finding cluster markers"); print(date())
cluster_ids <- unique(seurat_object@meta.data$seurat_clusters)
print(paste0("cluster_ids: ", print(cluster_ids))); print(date())
markers_df_list <- list()
for (value in cluster_ids) {
  out <- tryCatch(
    markers <- FindMarkers(seurat_object, ident.1 = value, test.use = test, logfc.threshold = logfc_thresh),
    error = function(e) {
      print("error in cluster " value)
    }
  )

  filepath <- paste0(output_markers_dir, "/markers_", value, ".tsv")
  markers_df_list
  write.table(markers, sep = "\t", file = filepath, quote = FALSE)
}
print("DONE")
print(date())
