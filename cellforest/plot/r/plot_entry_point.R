# source this file to the plotting code to pass parameters
# and write plots

args <- commandArgs(trailingOnly = TRUE)

r_plot_scripts_path <- args[1]
root_dir <- args[2]
path_to_temp_spec <- args[3]
current_process <- args[4]
plot_filepath <- args[5]
plot_width_px <- as.integer(args[6])
plot_height_px <- as.integer(args[7])
r_functions_filepath <- args[8]

library(Seurat)
library(ggplot2)
library(reticulate)
library(cellforestR)

py_run_string(args[9])  # kwargs = {...}
kwargs <- py$kwargs

source(r_functions_filepath)
spec <- py_load_object(path_to_temp_spec)
seurat_obj <- cellforest_load(root_dir, spec, current_process)

# Following keyword arguments are avalable for usage in plotting functions
group.by <- if (is.null(kwargs$group.by)) {NULL} else {kwargs$group.by}  # stratify
alpha <- ifelse(is.null(kwargs$alpha), 0.2, as.numeric(kwargs$alpha))  # point transparency
size <- ifelse(is.null(kwargs$size), 0.2, as.numeric(kwargs$size))  # point size
npcs <- ifelse(is.null(kwargs$npcs), 5, as.integer(as.numeric(kwargs$npcs)))  # number of principal components (usable for PCA and UMAP)
