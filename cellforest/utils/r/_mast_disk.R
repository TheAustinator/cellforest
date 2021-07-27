library(MAST)
library(Matrix)
library(scater)
library(parallel)

options(mc.cores = detectCores() - 1)
args <- commandArgs(trailingOnly = TRUE)

FORMULA <- args[1]

sce <- readRDS("/tmp/ad_sce.rds")
gc()
sca <- SceToSingleCellAssay(sce)
gc()
zlm_res <- zlm( as.formula(FORMULA), sca, parallel=TRUE)
gc()
df <- summary(zlm_res)$datatable
write.csv(df, "/tmp/df_zlm.csv")