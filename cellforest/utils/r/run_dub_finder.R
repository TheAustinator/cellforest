library(DoubletFinder)
library(Seurat)
library(R.utils)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
srat_path <- args[1]
output_path <- args[2]
dub_rate_per_1k <- as.numeric(args[3])    # 0.008
n_pcs <- as.numeric(args[4])    # 15
nfeatures <- as.numeric(args[5])    # 2000

print(paste0("srat_path: ", srat_path))
print(paste0("output_path: ", output_path))
print(paste0("dub_rate_per_1k: ", dub_rate_per_1k))
print(paste0("n_pcs: ", n_pcs))
print(paste0("nfeatures: ", nfeatures))


# TEMP
# srat <- readRDS(srat_path)
# if (is.null(srat@assays$originalexp)){srat@assays$originalexp=srat@assays$RNA}
# if (is.null(srat@assays$RNA)){srat@assays$RNA=srat@assays$originalexp}
# srat@assays$RNA = NULL
# srat$nFeature_RNA = srat$nFeature_originalexp

srat <- CreateSeuratObject(Read10X_h5(srat_path)$`Gene Expression`, project="placeholder", min.cells = 3,min.genes=200)


srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = nfeatures)
srat <- ScaleData(srat)
srat <- RunPCA(srat)
srat <- RunUMAP(srat, dims = 1:n_pcs)
# pK identification
print("sweep")
sweep_res_list <- paramSweep_v3(srat, PCs = 1:n_pcs, sct = FALSE)
print("summarize")
sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
print("pK")
bcmvn <- find.pK(sweep_stats)
pK <- as.numeric(bcmvn[which.max(bcmvn$BCmetric),][["pK"]])
nExp_poi <- round(dub_rate_per_1k / 1000 * nrow(srat@meta.data) * nrow(srat@meta.data))
# doubletFinder
print("running doublet finder")
srat <- doubletFinder_v3(srat, PCs = 1:n_pcs, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
dub_class_col <- colnames(srat@meta.data)[[which(startsWith(colnames(srat@meta.data), "DF"))]]
# write
print("writing labels")
write.csv(srat@meta.data[dub_class_col], file = output_path, col.names = F)
