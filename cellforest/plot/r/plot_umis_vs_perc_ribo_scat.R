source('cellforest/plot/r/plot_entry_point.R')

FeatureScatter(seurat_obj, "percent.ribo", "nCount_RNA", group.by = group.by) + 
    theme(legend.position = "bottom")

source('cellforest/plot/r/plot_exit_point.R')