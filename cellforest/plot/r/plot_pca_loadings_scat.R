source('cellforest/plot/r/plot_entry_point.R')

library(ggforce)

loadings <- data.frame(
    seurat_obj$pca@feature.loadings[, c(1:ifelse(is.null(kwargs$npcs), 5, kwargs$npcs))]
)

ggplot(loadings, aes(x = .panel_x, y = .panel_y)) + 
    geom_point(
        alpha = ifelse(is.null(kwargs$alpha), 0.2, kwargs$alpha),
        size = ifelse(is.null(kwargs$size), 0.2, kwargs$size)
    ) + 
    geom_autodensity() +
    geom_density2d() +
    facet_matrix(vars(everything()), layer.diag = 2, layer.upper = 3, grid.y.diag = FALSE)

dev.off()