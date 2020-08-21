source('cellforest/plot/r/plot_entry_point.R')

# TODO-QC: can we add an explicit xlab here (for cluster)?
VlnPlot(seurat_obj, features = "percent.hsp", group.by = group.by) +
    NoLegend()

source('cellforest/plot/r/plot_exit_point.R')