source('cellforest/plot/r/plot_entry_point.R')

VlnPlot(seurat_obj, features = "percent.ribo", group.by = kwargs$group.by) + NoLegend()

dev.off()