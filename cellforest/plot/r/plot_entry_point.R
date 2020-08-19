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

library(Seurat)
library(ggplot2)
library(reticulate)
library(cellforestR)

py_run_string(args[9])  # kwargs = {...}
kwargs <- py$kwargs

source(r_functions_filepath)
spec <- py_load_object(path_to_temp_spec)
seurat_obj <- cellforest_load(root_dir, spec, current_process)

if (plot_width_px > 750 & plot_height_px > 750) {  # big plot
    png(filename = plot_filepath, width = plot_width_px, height = plot_height_px, res = 150)
} else {  # default plot
    png(filename = plot_filepath, width = plot_width_px, height = plot_height_px)
}

# TODO-QC: nicer way to batch-handle kwargs
if (!is.null(kwargs$npcs)) {
    kwargs$npcs <- as.integer(as.numeric(kwargs$npcs))
}
if (!is.null(kwargs$size)) {
    kwargs$size <- as.numeric(kwargs$size)
}
if (!is.null(kwargs$alpha)) {
    kwargs$alpha <- as.numeric(kwargs$alpha)
}