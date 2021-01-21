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
#' 2. `install("cellforestR")`
#'
#' @param root_dir CellForest root directory
#' @param branch_spec Specification for CellBranch in Python or R (see `example_spec.R`)
#' @param process Process name (alias, if exists) at which Seurat object should be loaded
#'
#' @return Seurat object with cached dimensionality reduction embeddings
#'
#' @export
#'
#' @examples
#' library(cellforestR)
#'
#' root_dir <- "tests/data/example_usage/root"
#' seurat_obj <- cellforest_load(root_dir, example_spec_r, "reduce")  # given that processes have been run
#'
#' DimPlot(seurat_obj, reduction = "umap")
cellforest_load <- function(root_dir, branch_spec, process, subset = NULL, filter_ = NULL) {
  cellforest <- import("cellforest")
  cf_branch <- cellforest$load(root_dir, branch_spec = branch_spec)
  cf_branch$goto_process(process)
  seurat_obj <- get_seurat_object(cf_branch)
  # TODO: use Kristin's idea of adding a new column e.g. MYFLAG and using that
  #if (!is.null(subset) && !is.na(subset)) {
  #  temp_env <- new.env()
  #
  #  eval(parse(text = subset), envir = temp_env)
  #  AddMetaData()
  #  seurat_obj <- subset(seurat_obj, subset = subset[1] == subset[2])
  #}
  #if (!is.null(filter_) && !is.na(filter_)) {
  #  seurat_obj <- subset(seurat_obj, subset = filter_[1] != filter_[2])
  #}
  return(seurat_obj)
}

#' @title Load RDS matrix and embeddings at process run
#'
#' @description Loads RDS from most recent layer (including current layer).
#' If RDS comes from parental folder, append embeddings from the current layer.
#'
#' @importFrom Seurat CreateDimReducObject AddMetaData DefaultAssay
#' @importFrom glue glue
#'
#' @param cf_branch CellForest branch at selected process
#' @param dim_reduc_funs List of dimensionality reduction embeddings to append
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
#' cf_branch <- cellforest$load(root_dir, branch_spec = example_spec_py)
#' cf_branch$goto_process("reduce")  # put process or alias (if exists)
#' seurat_obj <- get_seurat_object(cf_branch)
get_seurat_object <- function(cf_branch) {
  current_process <- cf_branch$current_process
  current_path_map <- cf_branch[current_process]$path_map
  rds_path <- toString(current_path_map$rna_r)
  seurat_obj <- readRDS(file = rds_path)
  DefaultAssay(seurat_obj) <- "RNA"

  print(toString(glue("Creating Seurat object at process '{current_process}'")))
  if (current_process == "root") {
    meta_path <- paste0(dirname(rds_path), "/meta.tsv")
    meta <- read.table(meta_path, sep = "\t", header = TRUE, row.names = 1, quote="")
    seurat_obj <- AddMetaData(seurat_obj, meta)
    return(seurat_obj)
  }

  spec <- cf_branch$spec
  rds_process <- basename(dirname(dirname(rds_path)))
  precursors <- spec$get_precursors_lookup(incl_root = TRUE, incl_current = TRUE)[[current_process]]
  processes_to_load <- precursors[-(1:match(rds_process, precursors))]

  # check if rds is located in the same folder as metadata (if not -> needs update)
  if (!is.null(processes_to_load)) {
    for (process_name in processes_to_load) {
      process_path_map <- cf_branch[process_name]$path_map
      meta <- cf_branch$meta
      cols <- setdiff(names(meta), names(seurat_obj[[]]))
      if (length(cols) > 0) {
        seurat_obj <- AddMetaData(seurat_obj, meta[cols])
      }

      process <- cf_branch[process_name]$process
      if (process == "reduce") {
        seurat_obj <- add_dim_reduc_embed(
          seurat_obj,
          process_path_map,
        )
      }
      if (process == "cluster") {
        Idents(seurat_obj) <- seurat_obj[["cluster_id"]]
      }
      # TO-DO: Add loading for new processes
    }
  }
  print(toString(glue("Seurat object created at process '{current_process}'")))

  return(seurat_obj)
}

add_dim_reduc_embed <- function(seurat_object, path_map, dim_reduc_funcs = c("pca", "umap")) {
  if ("pca" %in% dim_reduc_funcs) {
    print("Loading PCA embeddings and loadings")
    seurat_object <- add_pca_embed(seurat_object, path_map)
    print(toString(glue("PCA embeddings and loadings loaded. Access them at <SEURAT_OBJECT>${PCA_EMBED_KEY}")))
  }

  if ("umap" %in% dim_reduc_funcs) {
    print("Loading UMAP embeddings")
    seurat_object <- add_umap_embed(seurat_object, path_map)
    print(toString(glue("UMAP embeddings loaded. Access them at <SEURAT_OBJECT>${UMAP_EMBED_KEY}")))
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
  input_embeddings_path <- toString(path_map$pca_embeddings)  # TO-DO: In the future, fetch from meta
  input_loadings_path <- toString(path_map$pca_loadings)
  input_stdev_path <- toString(path_map$pca_stdev)

  input_embeddings <- data.matrix(read.table(input_embeddings_path, sep = "\t", header = TRUE, row.names = 1))
  input_loadings <- data.matrix(read.table(input_loadings_path, sep = "\t", header = TRUE, row.names = 1))
  input_stdev <- as.list(read.table(input_stdev_path, sep = "\t", header = TRUE, row.names = 1))$x  # TODO-QC: is there a better way to get a list?

  seurat_object[[PCA_EMBED_KEY]] <- CreateDimReducObject(
    embeddings = input_embeddings,
    loadings = input_loadings,
    stdev = input_stdev,
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
  input_embeddings_path <- toString(path_map$meta)
  input_embeddings <- data.matrix(read.table(input_embeddings_path, sep = "\t", header = TRUE, row.names = 1))

  seurat_object[[UMAP_EMBED_KEY]] <- CreateDimReducObject(
    embeddings = input_embeddings,
    key = "UMAP_",
    assay = DefaultAssay(seurat_object)
  )

  return(seurat_object)
}

