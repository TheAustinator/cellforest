args <- commandArgs(trailingOnly = TRUE)

input_dir_10x <- args[1]
input_path_clusters <- args[2]
output_path_dx_est <- args[3]
output_path_dx_h5 <- args[4]
output_path_sx_est <- args[5]
output_path_sx_prof <- args[6]
output_path_sx_h5 <- args[7]


run_decontx <- function(input_path_outs, input_path_clusters = "", output_path_dx_est, output_path_dx_h5) {
  input_path_h5 <- file.path(input_path_outs, "filtered_feature_bc_matrix.h5")
  sce <- DropletUtils::read10xCounts(input_path_h5)
  if (length(input_path_clusters) > 4) {
    clusters <- read.csv(input_path_clusters)
  } else {
    clusters <- NULL
  }
  sce <- celda::decontX(sce, z = clusters)
  write.csv(SummarizedExperiment::colData(sce), output_path_dx_est)
  dxc <- sce@assays@data$decontXcounts
  # TODO: currently outputs h5 with no gene names but ensembl ids
  DropletUtils::write10xCounts(output_path_dx_h5, dxc, barcodes = sce@colData$Barcode, overwrite = T)
  return(sce)
}

run_soupx <- function(input_dir_10x, clusters, output_path_sx_est, output_path_sx_prof, output_path_sx_h5) {
  sx <- SoupX::load10X(input_dir_10x)
  sx <- SoupX::setClusters(sx, clusters)
  sx <- SoupX::autoEstCont(sx)
  write.csv(sx$metaData, output_path_sx_est)
  write.csv(sx$soupProfile, output_path_sx_prof)
  sx_counts <- SoupX::adjustCounts(sx)
  DropletUtils::write10xCounts(output_path_sx_h5, sx_counts, overwrite = T)
  # return(c(sx, sx_counts))
}

run_all <- function(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_dx_h5, output_path_sx_est, output_path_sx_prof, output_path_sx_h5) {
  sce_dx <- run_decontx(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_dx_h5)
  run_soupx(input_dir_10x, sce_dx$decontX_clusters, output_path_sx_est, output_path_sx_prof, output_path_sx_h5)   # list[sx, sx_counts] <-
  # return(c(sce_dx, sx, sx_counts))
}

run_all(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_dx_h5, output_path_sx_est, output_path_sx_prof, output_path_sx_h5)
