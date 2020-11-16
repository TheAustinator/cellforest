r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

# TODO-QC: can we add an explicit xlab here (for cluster)?
VlnPlot(seurat_obj, features = "percent.ribo", group.by = group.by) +
    NoLegend()

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))