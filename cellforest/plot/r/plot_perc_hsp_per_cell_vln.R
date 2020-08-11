source('cellforest/plot/r/plot_entry_point.R')

VlnPlot(seurat_obj, features = "percent.hsp") + NoLegend()

dev.off()