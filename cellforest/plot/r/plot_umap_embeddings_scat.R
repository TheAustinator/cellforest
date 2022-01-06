r_plot_scripts_path <- commandArgs(trailingOnly = TRUE)[1]
source(paste0(r_plot_scripts_path, "/plot_entry_point.R"))

library(ggforce)

umap_npcs <- min(ncol(seurat_obj$umap), npcs)
embeddings <- data.frame(seurat_obj$umap@cell.embeddings[, c(1:umap_npcs)])
if (!is.null(group.by)) {
    embeddings[[group.by]] <- seurat_obj@meta.data[match(row.names(embeddings), row.names(seurat_obj@meta.data)), group.by]
}

ggplot(embeddings, aes_string(x = ".panel_x", y = ".panel_y", fill = group.by, colour = group.by)) + 
    geom_point(
        alpha = ifelse(is.null(alpha), 0.2, alpha),
        size = ifelse(is.null(size), 0.2, size),
        position = "auto"
    ) + 
    geom_autodensity(alpha = 0.3, colour = NA, position = "identity") +
    facet_matrix(
        vars(names(embeddings[c(1:umap_npcs)])),
        layer.diag = 2
    ) +
    theme(legend.position = "bottom")

source(paste0(r_plot_scripts_path, "/plot_exit_point.R"))