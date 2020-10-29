r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

library(ggplot2)

seurat_obj <- CalculateBarcodeInflections(seurat_obj, threshold.low = 500)
BarcodeInflectionsPlot(seurat_obj) + ylab("log10(nCount_RNA)") + NoLegend()

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))