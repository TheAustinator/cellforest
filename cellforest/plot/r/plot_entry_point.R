# source this file to the plotting code to pass parameters
# and write plots

args <- commandArgs(trailingOnly = TRUE)

root_dir <- args[1]
path_to_temp_spec <- args[2]
current_process <- args[3]
plot_filepath <- args[4]
r_plot_scripts_path <- args[5]
plot_width_px <- as.integer(args[6])
plot_height_px <- as.integer(args[7])
r_functions_filepath <- args[8]
kwargs <- args[9]  # stratify, alpha

library(Seurat)
library(reticulate)
library(cellforestR)

source(r_functions_filepath)
spec <- py_load_object(path_to_temp_spec)
seurat_obj <- cellforest_load(root_dir, spec, current_process)

if (plot_width_px > 750 & plot_height_px > 750) {  # big plot
    png(filename = plot_filepath, width = plot_width_px, height = plot_height_px, res = 150)
} else {  # default plot
    png(filename = plot_filepath, width = plot_width_px, height = plot_height_px)
}
