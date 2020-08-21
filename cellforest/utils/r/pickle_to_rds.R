args <- commandArgs(trailingOnly = TRUE)
pickle_path <- commandArgs(trailingOnly = TRUE)[1]
meta_path <- args[2]
output_rds_path <- args[3]
python_executable <- args[4]

library(reticulate)
# TODO: move this to config
use_python(python_executable)
library(Seurat)
library(Matrix)


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
