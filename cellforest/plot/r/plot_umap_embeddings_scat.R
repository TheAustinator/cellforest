source('cellforest/plot/r/plot_entry_point.R')

library(ggforce)

default_npcs <- 2

embeddings <- data.frame(seurat_obj$umap@cell.embeddings[, c(1:ifelse(is.null(kwargs$npcs), default_npcs, kwargs$npcs))])
if (!is.null(kwargs$group.by)) {
    embeddings[[kwargs$group.by]] <- seurat_obj@meta.data[match(row.names(embeddings), row.names(seurat_obj@meta.data)), kwargs$group.by]
}

ggplot(embeddings, aes_string(x = ".panel_x", y = ".panel_y", fill = kwargs$group.by, colour = kwargs$group.by)) + 
    geom_point(
        alpha = ifelse(is.null(kwargs$alpha), 0.2, kwargs$alpha),
        size = ifelse(is.null(kwargs$size), 0.2, kwargs$size),
        position = "auto"
    ) + 
    geom_autodensity(alpha = 0.3, colour = NA, position = "identity") +
    facet_matrix(
        vars(names(embeddings[c(1:ifelse(is.null(kwargs$npcs), default_npcs, kwargs$npcs))])),
        layer.diag = 2
    )

dev.off()