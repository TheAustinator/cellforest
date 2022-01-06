r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

FeatureScatter(seurat_obj, "percent.mito", "nCount_RNA", group.by = group.by) + 
    theme(legend.position = "bottom")
    
source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))