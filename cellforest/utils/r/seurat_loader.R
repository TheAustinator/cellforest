library(Seurat)
library(reticulate)

PCA_EMBED_KEY = "pca_cf"
UMAP_EMBED_KEY = "umap_cf"

#' Load RDS matrix and embeddings at process run
#'
#' @param cf_branch Cellforest branch at selected process
#' 
#' @examples
#' library(reticulate)
#' cellforest <- import("cellforest")
#' 
#' source_python("playground/spec.py")
#' source("cellforest/utils/r/seurat_loader.R")
#' 
#' root_dir <- "tests/data/example_usage/root"
#' cf_branch <- cellforest$load(root_dir, spec=spec)
#' cf_branch$goto_process("reduce")
#' seurat_obj <- get_seurat_object(cf_branch)  # loads RDS and adds embeddings
#' 
#' DimPlot(seurat_obj, reduction = "pca_cf")
get_seurat_object <- function(cf_branch) {
  current_process <- cf_branch$current_process
  current_path_map <- cf_branch[current_process]$path_map
  rds_path <- toString(current_path_map$rna_r)
  seurat_object <- readRDS(file = rds_path)
  sprintf('Creating Seurat object at process "%s"', current_process); print(date())

  # check if rds path prefix matches process run path prefix
  meta_path <- toString(current_path_map$meta)
  rds_path_prefix <- substr(rds_path, start = 1, stop = tail(which(strsplit(rds_path, "")[[1]] == "/"), n = 1))
  meta_path_prefix <- substr(meta_path, start = 1, stop = tail(which(strsplit(meta_path, "")[[1]] == "/"), n = 1))

  if (rds_path_prefix != meta_path_prefix) {
    precursors <- cf_branch$spec$get_precursors_lookup(incl_current = TRUE)[[current_process]]
    for (process in precursors) {
      if (process == "reduce") {  # how to check for actual process name rather than alias?
        seurat_object <- add_dim_reduc_embed(
          seurat_object,
          current_path_map
        )
      }
    }
  }
  sprintf('Seurat object at process "%s" created.', current_path_map); print(date())

  return(seurat_object)
}

add_dim_reduc_embed <- function(seurat_object, path_map, dim_reduc_funcs = c("pca", "umap")) {
  if ("pca" %in% dim_reduc_funcs) {
    print("Loading PCA embeddings and loadings"); print(date())
    seurat_object <- add_pca_embed(seurat_object, path_map)
    sprintf("PCA embeddings and loadings loaded. Access them at seurat_object$%s", PCA_EMBED_KEY); print(date())
  }
  
  if ("umap" %in% dim_reduc_funcs) {
    print("Loading UMAP embeddings"); print(date())
    seurat_object <- add_umap_embed(seurat_object, path_map)
    sprintf("UMAP embeddings loaded. Access them at seurat_object$%s", UMAP_EMBED_KEY); print(date())
  }

  return(seurat_object)
}

add_pca_embed <- function(seurat_object, path_map) {
  input_embeddings_path <- toString(path_map$pca_embeddings)
  input_loadings_path <- toString(path_map$pca_loadings)
  input_embeddings <- data.matrix(read.table(input_embeddings_path, sep = "\t", header = TRUE, row.names = 1))
  input_loadings <- data.matrix(read.table(input_loadings_path, sep = "\t", header = TRUE, row.names = 1))

  seurat_object[[PCA_EMBED_KEY]] <- CreateDimReducObject(
    embeddings = input_embeddings,
    loadings = input_loadings,
    key = "PC_",
    assay = DefaultAssay(seurat_object)
  )

  return(seurat_object)
}

add_umap_embed <- function(seurat_object, path_map) {
  input_embeddings_path <- toString(path_map$meta)  # umap embeddings path should point to meta
  input_embeddings <- data.matrix(read.table(input_embeddings_path, sep = "\t", header = TRUE, row.names = 1))

  seurat_object[[UMAP_EMBED_KEY]] <- CreateDimReducObject(
    embeddings = input_embeddings,
    key = "UMAP_",
    assay = DefaultAssay(seurat_object)
  )

  return(seurat_object)
}