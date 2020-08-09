source('cellforest/plot/r/plot_entry_point.R')

library(scater)

pbmc.sce <- as.SingleCellExperiment(seurat_obj)
plotHighestExprs(pbmc.sce, exprs_values = "counts", colour_cells_by = "nFeature_RNA")

dev.off()