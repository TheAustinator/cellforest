# The __main__ construct from python doesn't seem to exist in R, so it's
# easier to stuff the functions in a sourceable file and keep the executable
# scripts stand alone

library(Seurat)
library(dplyr)
library(future)
library(sctransform)
library(stringr)
library(Matrix)
library(readr)

#' @export
metadata_filter_paths <- function(input_metadata_path, input_rds_path) {
  print("reading metadata"); print(date())
  meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1)
  print("reading rds"); print(date())
  seurat_object <- readRDS(input_rds_path)
  return(metadata_filter_objs(metadata, seurat_object))
}

#' @importFrom Seurat AddMetaData
#' @export
metadata_filter_objs <- function(meta, srat) {
  print(paste0("filtering cells by metadata and keeping ", nrow(meta), " / ", nrow(srat@meta.data))); print(date())
  srat <- subset(srat, cells = rownames(meta))
  srat <- AddMetaData(srat, meta)
  if (nrow(srat@meta.data) != nrow(meta)) {
    stop("Seurat object must contain all cell_ids present in metadata for filtering")
  }
  return(srat)
}

# TODO: only used in one process -- maybe move
#' @importFrom Seurat FindNeighbors FindClusters Idents
#' @export
find_clusters <- function(srat, output_clusters_path, num_pcs, resolution, nn_eps) {
  num_pcs <- 1:num_pcs
  print("Finding Neighbors"); print(date())
  srat <- FindNeighbors(object = srat, dims = num_pcs, verbose = TRUE, nn.eps = nn_eps, assay = "pca", graph.name = "pca_snn")
  print("Finding Clusters"); print(date())
  srat <- FindClusters(object = srat, resolution = resolution, verbose = TRUE, n.start = 10, graph.name = "pca_snn")
  print("Writing Clusters"); print(date())
  write.table(Idents(srat), sep = "\t", quote = FALSE, file = output_clusters_path, col.names = "cluster_id")
  return(srat)
}

#' @importFrom Seurat CreateSeuratObject AddMetaData Read10X
#' @export
srat_from_tenx <- function(input_tenx_directory_path, input_metadata_path, min_cells = 3) {
  print("Load 10X data"); print(date())
  metadata_df = read_tsv(input_metadata_path)
  metadata_df <- as.data.frame(metadata_df)
  rownames(metadata_df) <- metadata_df[, 1]
  tenx_data = Read10X(data.dir = input_tenx_directory_path, gene.column = 2)
  print("Create Seurat object"); print(date())
  tenx_subset = tenx_data[, colnames(tenx_data) %in% rownames(metadata_df)]
  seurat_object = CreateSeuratObject(counts = tenx_subset, min.cells = min_cells, meta.data = metadata_df)
  print("calculating percent mito"); print(date())
  mito.genes = grep(pattern = "^MT-", x = rownames(x = seurat_object$RNA@data), value = TRUE)
  percent.mito = Matrix::colSums(seurat_object$RNA@counts[mito.genes,]) / Matrix::colSums(seurat_object$RNA@counts)
  seurat_object = AddMetaData(object = seurat_object, metadata = percent.mito, col.name = "percent.mito")
  return(seurat_object)
}

#' @export
save_embeddings <- function(seurat_object, output_dir) {
  write.table(seurat_object$pca@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "pca_embeddings.tsv"), row.names = TRUE)
  write.table(seurat_object$pca@feature.loadings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "pca_loadings.tsv"), row.names = TRUE)
  write.table(seurat_object$umap@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "umap_embeddings.tsv"), row.names = TRUE)
}

#' @importFrom Seurat RunPCA
#' @export
run_pca <- function(seurat_object, output_embeddings_path, output_loadings_path, output_stdev_path, npcs = 30) {
  print("Running PCA"); print(date())
  seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(seurat_object), verbose = TRUE, npcs = npcs)

  # Should add this in later and this script will be a decision point
  # JackStrawPlot(object = seurat_object, PCs = 1:12)
  print("Writing PCA results"); print(date())
  write.table(seurat_object$pca@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = output_embeddings_path, row.names = TRUE)
  write.table(seurat_object$pca@feature.loadings, sep = "\t", col.names = TRUE, quote = FALSE, file = output_loadings_path, row.names = TRUE)
  write.table(seurat_object$pca@stdev, sep = "\t", col.names = TRUE, quote = FALSE, file = output_stdev_path, row.names = TRUE)
  return(seurat_object)
}

#' @importFrom Seurat PercentageFeatureSet
#' @export
add_genes_perc_meta <- function(seurat_object, patterns = list("percent.mito" = "^MT-", "percent.ribo" = "^RP[LS]", "percent.hsp" = "^HSP")) {
  # add columns of percentage feature set to metadata from named list of patterns
  for (key in names(patterns)) {
    seurat_object[[key]] <- PercentageFeatureSet(seurat_object, pattern = patterns[[key]])
  }

  return(seurat_object)
}

#' @importFrom Seurat FetchData
#' @export
filter_cells <- function(seurat_object, min_detected_genes, max_detected_genes, percent_mito_cutoff) {

  print("Filter cells")
  print(date())
  num_cells_before = length(seurat_object@active.ident)
  print(num_cells_before)
  filters = FetchData(seurat_object, vars = c("percent.mito", "nFeature_RNA"))

  seurat_object = seurat_object[, which(x = filters$nFeature_RNA > min_detected_genes &
    filters$nFeature_RNA < max_detected_genes &
    filters$percent.mito < percent_mito_cutoff)]
  num_cells_after = length(seurat_object@active.ident)

  print(num_cells_after)
  print(paste0("Filtered cells:", num_cells_before, "->", num_cells_after))

  return(seurat_object)
}

