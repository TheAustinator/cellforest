library(MAST)
library(Matrix)
library(scater)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

FORMULA <- args[1]
N_CORES <- as.integer(args[2])

if (N_CORES == -1) {
  options(mc.cores = detectCores() - 1)
} else {
  options(mc.cores = N_CORES)
}


sce <- readRDS("/tmp/ad_sce.rds")
gc()
sca <- SceToSingleCellAssay(sce)
gc()
zlm_res <- zlm( as.formula(FORMULA), sca, parallel=TRUE)
gc()
df <- summary(zlm_res)$datatable
write.csv(df, "/tmp/df_zlm.csv")