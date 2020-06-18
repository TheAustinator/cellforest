library(future)
library(dplyr)
library(Seurat)
plan("multiprocess", workers = 6)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_rds_path <- args[2]
output_diffexp_path <- args[3]
test <- args[4]
ident1 <- args[5]
ident2 <- args[6]
groupby <- args[7]
logfc_thresh <- as.numeric(args[8])
# This is janky, but R sucks and I can't find a better way
r_functions_filepath <- args[9]
source(r_functions_filepath)

filter_outputs <- metadata_filter(input_metadata_path, input_rds_path)
metadata <- filter_outputs$metadata
seurat_object <- filter_outputs$seurat_object

# Make sure the rows in metadata[groupby] are pulled out in the same order as the cell labels
# in the seurat_object. Remember, there's no guarantee they're in the same order and metadata[groupby]
# has no concept of cell labels
print(paste0("grouping by ", groupby)); print(date())
rownames(metadata) <- metadata$cell_id
seurat_object[[groupby, , ]] <- metadata[colnames(seurat_object), groupby, drop = F]

print("axding cluster_id"); print(date())
print(colnames(metadata))
seurat_object$cluster_id <- metadata$cluster_id
seurat_object[["cluster_id", , ]] <- metadata[colnames(seurat_object), "cluster_id", drop = F]

Idents(seurat_object) <- metadata$cluster_id
print(groupby)
cluster_ids <- unique(Idents(seurat_object))
print(paste0(length(cluster_ids),  " clusters exist: ", cluster_ids)); print(date())

print(paste0("Identifying markers with logfc_thresh: ", logfc_thresh)); print(date())

datalist <- list()
i <- 1
for (value in cluster_ids) {
  print("cluster ", value)

  out <- tryCatch(
    {
    markers <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2, group.by = groupby, subset.ident = value, test.use = test, logfc.threshold = logfc_thresh)
    markers$cluster <- value
    markers$gene_symbol <- rownames(markers)
    rownames(markers) <- NULL
    datalist[[i]] <- markers
    i <- i + 1
  },
    error = function(e) e
    # {
    # print(paste0("ERROR on cluster ", cond)); print(date())
    # markers <- data.frame()
    # print("error")
    # return(NA)
    # }
  )
}

all_markers <- do.call(rbind, datalist)
write.table(all_markers, sep="\t", file=output_diffexp_path, quote=FALSE, row.names=FALSE)
