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


metadata_filter_paths <- function(input_metadata_path, input_rds_path) {
  print("reading metadata"); print(date())
  meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1)
  print("reading rds"); print(date())
  seurat_object <- readRDS(input_rds_path)
  return(metadata_filter_objs(metadata, seurat_object))
}

metadata_filter_objs <- function(meta, srat) {
  print(paste0("filtering cells by metadata and keeping ", nrow(meta), " / ", nrow(srat@meta.data))); print(date())
  srat <- subset(srat, cells = rownames(meta))
  srat <- AddMetaData(srat, meta)
  if (nrow(srat@meta.data) != nrow(meta)) {
    stop("Seurat object must contain all cell_ids present in metadata for filtering")
  }
  return(srat)
}


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


find_cluster_markers <- function(seurat_object, output_diffexp_path) {
  cluster_ids = unique(seurat_object@meta.data$cluster_id)
  for (value in cluster_ids) {
    markers = FindMarkers(seurat_object, ident.1 = value)
    write.table(markers, sep = "\t", file = output_markers_path, quote = FALSE)
  }
  return(seurat_object)
}


diff_exp <- function(seurat_object, output_diffexp_path, test, ident1, ident2, groupby, logfc_thresh) {
  print("getting cluster identities")
  groups <- unique(Idents(seurat_object))
  print("groups")
  print(groups)
  datalist <- list()
  i <- 1
  for (value in cluster_ids) {
    out <- tryCatch(
      {
      markers <- FindMarkers(seurat_object, ident.1 = ident1, ident.2 = ident2, group.by = groupby, subset.ident = value, test.use = test, logfc.threshold = logfc_threshold)
      markers$cluster = value
      markers$gene_symbol = rownames(markers)
      rownames(markers) <- NULL
      datalist[[i]] <- markers
      i <- i + 1
    },
      error = function(cond) {
        markers <- data.frame()
        return(NA)
      })
  }
  all_markers <- do.call(rbind, datalist)
  write.table(all_markers, sep = "\t", file = output_path, quote = FALSE, row.names = FALSE)
}


create_seurat_object <- function(input_tenx_directory_path, input_metadata_path, min_cells = 3) {
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

run_standard_workflow <- function(seurat_object, nfeatures = 2000, npcs = 30, verbose = TRUE) {
  print(date())
  seurat_object <- NormalizeData(seurat_object, verbose = verbose)
  print("Done normalizing")
  print(date())
  seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = nfeatures)
  print("Done finding variable")
  print(date())
  seurat_object <- ScaleData(seurat_object, verbose = verbose)
  print("Done scaling")
  print(date())
  seurat_object <- RunPCA(object = seurat_object, verbose = verbose, npcs = npcs)
  print("Done PCA")
  print(date())
  seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:npcs, verbose = verbose)
  print("Done UMAP")
  print(date())
  return(seurat_object)
}

save_embeddings <- function(seurat_object, output_dir) {
  write.table(seurat_object$pca@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "pca_embeddings.tsv"), row.names = TRUE)
  write.table(seurat_object$pca@feature.loadings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "pca_loadings.tsv"), row.names = TRUE)
  write.table(seurat_object$umap@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = file.path(output_dir, "umap_embeddings.tsv"), row.names = TRUE)
}

run_integration <- function(seurat_object, normalize = "sct") {
  runlane_list <- SplitObject(seurat_object, split.by = "run_lane")

  if (normalize == "std") {
    for (i in 1:length(runlane_list)) {
      runlane_list[[i]] <- NormalizeData(runlane_list[[i]], verbose = TRUE)
      runlane_list[[i]] <- FindVariableFeatures(runlane_list[[i]], selection.method = "vst",
                                                nfeatures = 2000, verbose = TRUE)
    }
    anchors <- FindIntegrationAnchors(object.list = runlane_list, dims = 1:30)
    runlane_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

    DefaultAssay(runlane_integrated) <- "integrated"

    # Run the standard workflow for visualization and clustering
    runlane_integrated <- ScaleData(runlane_integrated, verbose = TRUE)
  }
  else if (normalize == "sct") {
    for (i in 1:length(runlane_list)) {
      runlane_list[[i]] <- SCTransform(runlane_list[[i]], verbose = FALSE)
    }
  }

  runlane_integrated <- RunPCA(runlane_integrated, npcs = 30, verbose = TRUE)
  runlane_integrated <- RunUMAP(runlane_integrated, reduction = "pca", dims = 1:30)

  return(runlane_integrated)

}

run_pca <- function(seurat_object, output_embeddings_path, output_loadings_path, npcs = 30) {
  print("Running PCA"); print(date())
  seurat_object <- RunPCA(object = seurat_object, verbose = TRUE, npcs = npcs)

  # Should add this in later and this script will be a decision point
  # JackStrawPlot(object = seurat_object, PCs = 1:12)
  print("Writing PCA results"); print(date())
  write.table(seurat_object$pca@cell.embeddings, sep = "\t", col.names = TRUE, quote = FALSE, file = output_embeddings_path, row.names = TRUE)
  write.table(seurat_object$pca@feature.loadings, sep = "\t", col.names = TRUE, quote = FALSE, file = output_loadings_path, row.names = TRUE)
  return(seurat_object)
}

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

run_sctransform <- function(seurat_object, corrected_umi_output_path, pearson_residual_output_path) {
  print("sctransform")
  print(date())
  seurat_object = SCTransform(object = seurat_object, verbose = TRUE)

  writeMM(seurat_object$SCT@counts, file = corrected_umi_output_path)
  write.table(seurat_object$SCT@scale.data, file = pearson_residual_output_path, sep = "\t", quote = FALSE)

  return(seurat_object)
}
