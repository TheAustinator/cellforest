library(reticulate)
# TODO: move this to config
use_python("/Users/austinmckay/code/cellforest/venv/bin/python")
library(Seurat)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
pickle_path <- args[1]
meta_path <- args[2]
output_rds_path <- args[3]

store <- py_load_object(pickle_path)
counts <- t(store$matrix)
colnames(counts) <- as.vector(store$cell_ids)
rownames(counts) <- as.vector(t(store$features["genes"]))
srat <- CreateSeuratObject(counts)
meta <- read.table(meta_path, sep = "\t", header = TRUE, row.names = 1)
# TODO: move mito to original python
mito.genes <- grep(pattern = "^MT-", x = rownames(x = srat$RNA@data), value = TRUE)
percent.mito <- Matrix::colSums(srat$RNA@counts[mito.genes,]) / Matrix::colSums(srat$RNA@counts)
seurat_object <- AddMetaData(object = srat, metadata = percent.mito, col.name = "percent.mito")
saveRDS(seurat_object, output_rds_path)
