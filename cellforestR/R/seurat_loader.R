#' @include example_spec.R
#' @include helper_functions.R

PCA_EMBED_KEY = "pca"
UMAP_EMBED_KEY = "umap"

#' @title Interface for load RDS matrix and embeddings at process run
#'
#' @description Wrapper around get_seurat_object(), main interface for creating
#' Seurat object from a CellForest process run
#'
#' @importFrom reticulate import
#'
#' @section Installation:
#' 1. `setwd(PATH_TO_CELLFOREST)`
#' 2. `install(cellforestR)`
#'
#' @param root_dir CellForest root directory
#' @param spec Specification for CellBranch in Python or R (see `example_spec.R`)
#' @param process Process name (alias, if exists) at which Seurat object should be loaded
#'
#' @return Seurat object with cached dimensionality reduction embeddings
#'
#' @export
#'
#' @examples
#' root_dir <- "tests/data/example_usage/root"
#' seurat_obj <- cellforest_load(root_dir, example_spec_r, "reduce")
#'
#' DimPlot(seurat_obj, reduction = "umap")
cellforest_load <- function(root_dir, spec, process) {
  cellforest <- import("cellforest")
  cf_branch <- cellforest$load(root_dir, spec = spec)
  cf_branch$goto_process(process)
  seurat_obj <- get_seurat_object(cf_branch)

  return(seurat_obj)
}

#' @title Load RDS matrix and embeddings at process run
#'
#' @description Loads RDS from most recent layer (including current layer).
#' If RDS comes from parental folder, append embeddings from the current layer.
#'
#' @importFrom Seurat CreateDimReducObject
#' @importFrom glue glue
#'
#' @param cf_branch CellForest branch at selected process
#'
#' @return Seurat object with cached dimensionality reduction embeddings
#'
#' @export
#'
#' @examples
#' library(reticulate)
#' library(Seurat)
#'
#' library(cellforestR)
#' cellforest <- import("cellforest")
#'
#' root_dir <- "tests/data/example_usage/root"
#' cf_branch <- cellforest$load(root_dir, spec = example_spec_py)
#' cf_branch$goto_process("reduce")  # put process or alias (if exists)
#' seurat_obj <- get_seurat_object(cf_branch)
get_seurat_object <- function(cf_branch) {
  current_process <- cf_branch$current_process
  current_path_map <- cf_branch[current_process]$path_map
  rds_path <- toString(current_path_map$rna_r)
  seurat_object <- readRDS(file = rds_path)
  print(toString(glue("Creating Seurat object at process '{current_process}'")))

  # check if rds is located in the same folder as metadata (if not -> needs update)
  meta_path <- toString(current_path_map$meta)
  rds_path_prefix <- get_prefix_from_path(rds_path)
  meta_path_prefix <- get_prefix_from_path(meta_path)

  if (rds_path_prefix != meta_path_prefix) {
    all_precursors <- cf_branch$spec$get_precursors_lookup(incl_current = TRUE)
    precursors <- all_precursors[[current_process]]

    for (process_name in precursors) {
      process <- cf_branch[process_name]$process
      if (process == "reduce") {
        seurat_object <- add_dim_reduc_embed(
          seurat_object,
          current_path_map
        )
      }

      # TO-DO: Add loading for new processes
    }
  }
  print(toString(glue("Seurat object at process '{current_process}' created")))

  return(seurat_object)
}

add_dim_reduc_embed <- function(seurat_object, path_map, dim_reduc_funcs = c("pca", "umap")) {
  if ("pca" %in% dim_reduc_funcs) {
    print("Loading PCA embeddings and loadings")
    seurat_object <- add_pca_embed(seurat_object, path_map)
    print(toString(glue("PCA embeddings and loadings loaded. Access them at seurat_object${PCA_EMBED_KEY}")))
  }

  if ("umap" %in% dim_reduc_funcs) {
    print("Loading UMAP embeddings")
    seurat_object <- add_umap_embed(seurat_object, path_map)
    print(toString(glue("UMAP embeddings loaded. Access them at seurat_object${UMAP_EMBED_KEY}")))
  }

  return(seurat_object)
}

#' @title Append PCA embeddings and loadings to RDS matrix
#'
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#'
#' @param seurat_object Seurat object
#' @param path_map ProcessRun's path map (paths of relevant files)
#'
#' @return Seurat object with PCA embeddings and loadings
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

#' @title Append UMAP embeddings to RDS matrix
#'
#' @importFrom Seurat CreateDimReducObject DefaultAssay
#'
#' @param seurat_object Seurat object
#' @param path_map ProcessRun's path map (paths of relevant files)
#'
#' @return Seurat object with UMAP embeddings
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