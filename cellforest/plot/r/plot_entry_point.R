# source this file to the plotting code to pass parameters
# and write plots

library(Seurat)
library(ggplot2)
library(cellforestR)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)

r_plot_scripts_path <- args[1]
root_dir <- args[2]
spec_str <- args[3]
subset_key <- args[4]
subset_val <- args[5]
current_process <- args[6]
plot_filepath <- args[7]
plot_width_px <- as.integer(args[8])
plot_height_px <- as.integer(args[9])
py_run_string(args[10])  # kwargs = {...}
kwargs <- py$kwargs

#source(r_functions_filepath)
if (!is.null(subset_key) && !all(is.na(subset_key))) {
  subset <- c(subset_key, subset_val)
} else {
  subset <- NULL
}
seurat_obj <- cellforest_load(root_dir, spec_str, current_process, subset)


# Following keyword arguments are avalable for usage in plotting functions
group.by <- if (is.null(kwargs$group.by)) {NULL} else {kwargs$group.by}  # stratify
alpha <- ifelse(is.null(kwargs$alpha), 0.2, as.numeric(kwargs$alpha))  # point transparency
size <- ifelse(is.null(kwargs$size), 0.2, as.numeric(kwargs$size))  # point size
npcs <- ifelse(is.null(kwargs$npcs), 5, as.integer(as.numeric(kwargs$npcs)))  # number of principal components (usable for PCA and UMAP)
