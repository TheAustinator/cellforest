library("reticulate")
use_python("/Users/austinmckay/code/cellforest/venv/bin/python")
library("Matrix")

args <- commandArgs(trailingOnly = TRUE)

input_rds_path <- args[1]
output_meta_path <- args[2]
output_pickle_path <- args[3]
counts_store_module <- args[4]

source_python(counts_store_module)


srat <- readRDS(input_rds_path)
matrix <- t(srat@assays$RNA@data)
cell_ids <- as.data.frame(rownames(matrix))
colnames(cell_ids) <- as.numeric(0)    # default python column name when header not stored
# TODO: check whether seurat has ensgs, b/c otherwise, stuck with just gene names
features <- as.data.frame(cbind(ensgs = "None", genes = colnames(matrix), mode = "Gene Expression"))
build_counts_store(matrix, cell_ids, features, output_pickle_path)
write.table(srat[[]], output_meta_path)