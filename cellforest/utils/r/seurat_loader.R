library(Seurat)
library(reticulate)
library(glue)

PCA_EMBED_KEY = "pca"
UMAP_EMBED_KEY = "umap"

#' Load RDS matrix and embeddings at process run
#'
#' @param cf_branch Cellforest branch at selected process
#' 
#' @return Seurat object with cached dimensionality reduction embeddings
#' 
#' @examples
#' library(reticulate)
#' cellforest <- import("cellforest")
#' 
#' source("cellforest/utils/r/seurat_loader.R")
#' 
#' root_dir <- "tests/data/example_usage/root"
#' cf_branch <- cellforest$load(root_dir, spec = example_spec)
#' cf_branch$goto_process("reduce")
#' seurat_obj <- get_seurat_object(cf_branch)  # loads RDS and adds embeddings
#' 
#' DimPlot(seurat_obj, reduction = "pca")
get_seurat_object <- function(cf_branch) {
  current_process <- cf_branch$current_process
  current_path_map <- cf_branch[current_process]$path_map
  rds_path <- toString(current_path_map$rna_r)
  seurat_object <- readRDS(file = rds_path)
  print(toString(glue("Creating Seurat object at process '{current_process}'"))); print(date())

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
  print(toString(glue("Seurat object at process '{current_process}' created"))); print(date())

  return(seurat_object)
}

add_dim_reduc_embed <- function(seurat_object, path_map, dim_reduc_funcs = c("pca", "umap")) {
  if ("pca" %in% dim_reduc_funcs) {
    print("Loading PCA embeddings and loadings"); print(date())
    seurat_object <- add_pca_embed(seurat_object, path_map)
    print(toString(glue("PCA embeddings and loadings loaded. Access them at seurat_object${PCA_EMBED_KEY}"))); print(date())
  }
  
  if ("umap" %in% dim_reduc_funcs) {
    print("Loading UMAP embeddings"); print(date())
    seurat_object <- add_umap_embed(seurat_object, path_map)
    print(toString(glue("UMAP embeddings loaded. Access them at seurat_object${UMAP_EMBED_KEY}"))); print(date())
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

# list of names lists will convert to list of dictionaries in Python
example_spec <- type.convert(list(
  list(
    "process" = "normalize",
    "params" = list(
      "min_genes" = 4,
      "max_genes" = 5002,
      "min_cells" = 5,
      "nfeatures" = 30,
      "perc_mito_cutoff" = 20,
      "method" = "seurat_default"
    ),
    "subset" = list(
      "sample" = "sample_1"
    )
  ),
  
  list(
    "process" = "reduce",
    "params" = list(
      "pca_npcs" = 3,
      "umap_n_neighbors" = 3,
      "umap_min_dist" = 0.1,
      "umap_n_components" = 2,
      "umap_metric" = "euclidean"
    )
  )
))