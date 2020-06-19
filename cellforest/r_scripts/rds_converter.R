library(Matrix)

args <- commandArgs(trailingOnly = TRUE)

rds_path <- args[1]
output_matrix_path <- args[2]
output_cell_ids_path <- args[3]
output_genes_path <- args[4]
output_metadata_path <- args[5]



seurat_object <- readRDS(rds_path)
writeMM(seurat_object@assays$RNA@data, output_matrix_path)

write.table(seurat_object@assays$RNA@data@Dimnames[[1]],file=output_genes_path,col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(seurat_object@assays$RNA@data@Dimnames[[2]], file=output_cell_ids_path, col.names=FALSE, row.names=FALSE, quote=FALSE)

if ("meta.data" %in% slotNames(seurat_object)) {
  df <- seurat_object@meta.data
  cell_ids <- as.data.frame(rownames(seurat_object@meta.data), col.names=c("cell_id"))
  colnames(cell_ids) <- "cell_id"
  df <- cbind(cell_ids, df)
  write.table(df, file=output_metadata_path, col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}