args <- commandArgs(trailingOnly = TRUE)

root_dir <- args[1]
path_to_temp_spec <- args[2]
current_process <- args[3]
plot_filepath <- args[4]
r_plot_scripts_path <- args[5]
r_functions_filepath <- args[6]

library(Seurat)
library(reticulate)
library(cellforestR)

source(r_functions_filepath)
spec <- py_load_object(path_to_temp_spec)
seurat_obj <- cellforest_load(root_dir, spec, current_process)

png(filename=plot_filepath)
# after this comes the actual plot code