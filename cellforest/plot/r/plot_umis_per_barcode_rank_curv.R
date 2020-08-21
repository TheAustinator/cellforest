source('cellforest/plot/r/plot_entry_point.R')

library(ggplot2)

seurat_obj <- CalculateBarcodeInflections(seurat_obj, threshold.low = 100)
BarcodeInflectionsPlot(seurat_obj) + ylab("log10(nCount_RNA)") + NoLegend()

source('cellforest/plot/r/plot_exit_point.R')