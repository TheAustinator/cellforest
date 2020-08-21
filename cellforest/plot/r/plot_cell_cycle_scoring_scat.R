source('cellforest/plot/r/plot_entry_point.R')

default_reduction <- "umap"
reduction = ifelse(is.null(kwargs$reduction), default_reduction, kwargs$reduction)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(seurat_obj, reduction = reduction) +
    theme(legend.position = "bottom")

source('cellforest/plot/r/plot_exit_point.R')