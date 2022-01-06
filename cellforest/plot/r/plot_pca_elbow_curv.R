r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

# TODO-QC: Add Kristin's elbow plot threshold picker
ElbowPlot(seurat_obj)

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))