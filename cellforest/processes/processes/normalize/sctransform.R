library(future)
library(parallel)
library(cellforestR)
library(Seurat)

options(future.globals.maxSize = 8000 * 1024^2)
plan("multiprocess", workers = detectCores() - 1)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- commandArgs(trailingOnly = TRUE)[1]
input_tenx_directory_path <- args[2]
output_rds_path <- args[3]
min_genes <- as.numeric(args[4])
max_genes <- as.numeric(args[5])
min_cells <- as.numeric(args[6])
perc_mito_cutoff <- as.numeric(args[7])
#r_functions_filepath <- args[8]

output_corrected_umi_path <- args[8]
output_pearson_residual_path <- args[9]

#source(r_functions_filepath)

run_sctransform <- function(seurat_object, corrected_umi_output_path, pearson_residual_output_path) {
  print("sctransform")
  print(date())
  seurat_object = SCTransform(object = seurat_object, verbose = TRUE)

  writeMM(seurat_object$SCT@counts, file = corrected_umi_output_path)
  write.table(seurat_object$SCT@scale.data, file = pearson_residual_output_path, sep = "\t", quote = FALSE)

  return(seurat_object)
}

seurat_object <- srat_from_tenx(input_tenx_directory_path, input_metadata_path, min_cells)
print("reading metadata"); print(date())
metadata <- read_tsv(input_metadata_path)
print("metadata filter"); print(date())
filter_outputs <- metadata_filter_objs(metadata, seurat_object)
seurat_object <- filter_outputs$seurat_object
print("filtering cells"); print(date())
seurat_object <- add_genes_perc_meta(seurat_object, patterns = list("percent.mito" = "^MT-", "percent.ribo" = "^RP[LS]", "percent.hsp" = "^HSP"))
seurat_object <- filter_cells(seurat_object, min_genes, max_genes, perc_mito_cutoff)
print("running sctransform"); print(date())
seurat_object <- run_sctransform(seurat_object, output_corrected_umi_path, output_pearson_residual_path)
print("Saving output object"); print(date())
saveRDS(seurat_object, file = output_rds_path)
print("sctransform DONE"); print(date())