source('cellforest/plot/r/plot_entry_point.R')

FeatureScatter(seurat_obj, "percent.ribo", "nCount_RNA", group.by=kwargs$group.by) + NoLegend()

dev.off()