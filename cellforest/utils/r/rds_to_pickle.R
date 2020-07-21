args <- commandArgs(trailingOnly = TRUE)

input_rds_path <- args[1]
output_meta_path <- args[2]
output_pickle_path <- args[3]
counts_store_module <- args[4]
python_executable <- args[5]

library("reticulate")
use_python(python_executable)
library("Matrix")

source_python(counts_store_module)


srat <- readRDS(input_rds_path)
matrix <- t(srat@assays$RNA@data)
cell_ids <- as.data.frame(rownames(matrix))
colnames(cell_ids) <- as.numeric(0)    # default python column name when header not stored
# TODO: check whether seurat has ensgs, b/c otherwise, stuck with just gene names
features <- as.data.frame(cbind(ensgs = "None", genes = colnames(matrix), mode = "Gene Expression"))
# sourced from python
build_counts_store(matrix, cell_ids, features, output_pickle_path)
meta <- srat[[]]
# TODO: update colnames here

write.table(meta, output_meta_path)