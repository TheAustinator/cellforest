r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

library(scater)

pbmc.sce <- as.SingleCellExperiment(seurat_obj)
plotHighestExprs(pbmc.sce, exprs_values = "counts", colour_cells_by = group.by)

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))