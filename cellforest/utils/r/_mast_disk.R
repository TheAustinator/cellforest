library(MAST)
library(Matrix)
library(scater)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

RUN_ID <- args[1]
FORMULA <- args[2]
CORES <- as.integer(args[3])

print(paste0("MAST: run_id=",  RUN_ID, "; formula: ", FORMULA, "; cores=", CORES))

if (CORES == 0) {
  CORES <- detectCores() - 1
}
options(mc.cores = CORES)

sce <- readRDS(paste0("/tmp/cf_ad_sce_", RUN_ID, ".rds"))
gc()
sca <- SceToSingleCellAssay(sce)
gc()
zlm_res <- zlm( as.formula(FORMULA), sca, parallel=TRUE)
gc()
df <- summary(zlm_res,doLRT=TRUE)$datatable
write.csv(df, paste0("/tmp/cf_df_zlm_", RUN_ID, ".csv"))
