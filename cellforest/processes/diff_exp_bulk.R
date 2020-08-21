library(future)
library(dplyr)
library(Seurat)
plan("multiprocess", workers = 6)

args = commandArgs(trailingOnly = TRUE)

input_metadata_path <- commandArgs(trailingOnly = TRUE)[1]
input_rds_path <- args[2]
output_diffexp_path = args[3]
test <- args[4]
ident1 <- args[5]
ident2 <- args[6]
groupby <- args[7]
logfc_thresh <- as.numeric(args[8])

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

Idents(seurat_object) <- groupby
print(groupby)
unique_idents <- unique(Idents(seurat_object))
print(paste0(length(unique_idents), " groups exist: ", unique_idents)); print(date())

print(paste0("Identifying markers with logfc_thresh: ", logfc_thresh)); print(date())
markers <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2, test.use = test, logfc.threshold = logfc_thresh) # group.by=groupby, subset.ident=value,
print("writing markers"); print(date())
write.table(markers, sep = "\t", file = output_diffexp_path, quote = FALSE, row.names = FALSE)
print("diffexp_bulk DONE"); print(date())
