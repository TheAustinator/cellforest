r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

library(ggforce)

pca_npcs <- min(ncol(seurat_obj$pca), npcs)
loadings <- data.frame(
    seurat_obj$pca@feature.loadings[, c(1:pca_npcs)]
)

ggplot(loadings, aes(x = .panel_x, y = .panel_y)) + 
    geom_point(
        alpha = ifelse(is.null(kwargs$alpha), 0.2, kwargs$alpha),
        size = ifelse(is.null(kwargs$size), 0.2, kwargs$size)
    ) + 
    geom_autodensity() +
    geom_density2d() +
    facet_matrix(vars(everything()), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE) +
    theme(legend.position = "bottom")

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))