#' @title R spec definition
#'
#' @description List of names lists will convert to list of dictionaries in Python
#'
#' @export
example_spec_r <- type.convert(list(
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

#' @importFrom reticulate py_run_string py
py_run_string('example_spec_py = [
  {
      "process": "normalize",
      "params": {
          "min_genes": 4,
          "max_genes": 5002,
          "min_cells": 5,
          "nfeatures": 30,
          "perc_mito_cutoff": 20,
          "method": "seurat_default",
      },
      "subset": {"sample": "sample_1"},
  },
  {
      "process": "reduce",
      "params": {
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
