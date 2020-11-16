r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

default_reduction <- "umap"
reduction = ifelse(is.null(kwargs$reduction), default_reduction, kwargs$reduction)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(seurat_obj, reduction = reduction) +
    theme(legend.position = "bottom")

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))