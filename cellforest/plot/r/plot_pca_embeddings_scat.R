source('cellforest/plot/r/plot_entry_point.R')

library(ggforce)

embeddings <- data.frame(seurat_obj$pca@cell.embeddings[, c(1:ifelse(is.null(kwargs$npcs), 5, kwargs$npcs))])
if (!is.null(kwargs$group.by)) {
    embeddings[[kwargs$group.by]] <- seurat_obj@meta.data[match(row.names(embeddings), row.names(seurat_obj@meta.data)), kwargs$group.by]
}

ggplot(embeddings, aes_string(x = ".panel_x", y = ".panel_y", fill = kwargs$group.by, colour = kwargs$group.by)) + 
    geom_point(shape = 16, size = 0.5, position = "auto") + 
    geom_autodensity(alpha = 0.3, colour = NA, position = "identity") +
    facet_matrix(
        vars(names(embeddings[c(1:ifelse(is.null(kwargs$npcs), 5, kwargs$npcs))])),
        layer.diag = 2
    )

dev.off()