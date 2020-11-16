r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

library(ggforce)

pca_npcs <- min(ncol(seurat_obj$pca), npcs)
embeddings <- data.frame(seurat_obj$pca@cell.embeddings[, c(1:pca_npcs)])
if (!is.null(kwargs$group.by)) {
    embeddings[[kwargs$group.by]] <- seurat_obj@meta.data[match(row.names(embeddings), row.names(seurat_obj@meta.data)), kwargs$group.by]
}

ggplot(embeddings, aes_string(x = ".panel_x", y = ".panel_y", fill = kwargs$group.by, colour = kwargs$group.by)) + 
    geom_point(
        alpha = alpha,
        size = ifelse(is.null(kwargs$size), 0.2, kwargs$size),
        position = "auto"
    ) + 
    geom_autodensity(alpha = 0.3, colour = NA, position = "identity") +
    facet_matrix(
        vars(names(embeddings[c(1:pca_npcs)])),
        layer.diag = 2
    ) +
    theme(legend.position = "bottom")

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))