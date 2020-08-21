source('cellforest/plot/r/plot_entry_point.R')

library(data.table)

output_markers_dir <- args[10]

n_clusters <- length(unique(Idents(pbmc)))
list_of_marker_tables = list()
for (i in 1:n_clusters) {
    marker_genes <- read.table(file = paste0(output_markers_dir, "/marker_", i, ".tsv"), sep = "\t", quote = FALSE)
    list_of_marker_tables <- c(list_of_marker_tables, marker_genes)
}

all_marker_genes <- rbindlist(list_of_marker_tables)
plot_data <- data.frame(all_marker_genes %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% tally())

ggplot(data = plot_data, aes(x = cluster, y = n, fill = cluster)) +
  geom_bar(stat="identity") +
  labs(title = paste0('Number of marker genes per cluster (max_pval < 0.05)')) +
  xlab('cluster') +
  ylab('number of marker genes') +
  NoLegend() +
  theme_minimal()

source('cellforest/plot/r/plot_exit_point.R')