#' @title R spec definition
#'
#' @description List of names lists will convert to list of dictionaries in Python
#'
#' @export
example_spec_r <- type.convert(list(
  list(
    "_PROCESS_" = "normalize",
    "_PARAMS_" = list(
      "min_genes" = 5,
      "max_genes" = 5000,
      "min_cells" = 5,
      "nfeatures" = 30,
      "perc_mito_cutoff" = 20,
      "method" = "seurat_default"
    ),
    "_SUBSET_" = list(
      "sample" = "sample_1"
    )
  ),

  list(
    "_PROCESS_" = "reduce",
    "_PARAMS_" = list(
      "pca_npcs" = 3,
      "umap_n_neighbors" = 3,
      "umap_min_dist" = 0.1,
      "umap_n_components" = 2,
      "umap_metric" = "euclidean"
    )
  )
))

#' @importFrom reticulate py_run_string py
py_run_string('example_spec_py = [
  {
      "_PROCESS_": "normalize",
      "_PARAMS_": {
          "min_genes": 5,
          "max_genes": 5000,
          "min_cells": 5,
          "nfeatures": 30,
          "perc_mito_cutoff": 20,
          "method": "seurat_default",
      },
      "_SUBSET_": {"sample": "sample_1"},
  },
  {
      "_PROCESS_": "reduce",
      "_ALIAS_": "pca+umap",
      "_PARAMS_": {
          "pca_npcs": 3,
          "umap_n_neighbors": 3,
          "umap_min_dist": 0.1,
          "umap_n_components": 2,
          "umap_metric": "euclidean",
      },
  },
]')

#' @title Python spec definition
#'
#' @description Converts of Python list of dicts into a reticulate object
#'
#' @export
example_spec_py <- py$example_spec_py
