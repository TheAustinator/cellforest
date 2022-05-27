args <- commandArgs(trailingOnly = TRUE)

library("DropletUtils")

input_dir_10x <- args[1]
input_path_clusters <- args[2]
output_path_dx_est <- args[3]
output_path_sx_est <- args[4]
output_path_sx_prof <- args[5]

print(paste0("input_dir_10x: ", input_dir_10x))
print(paste0("input_path_clusters: ", input_path_clusters))
print(paste0("output_path_dx_est: ", output_path_dx_est))
print(paste0("output_path_sx_est: ", output_path_sx_est))
print(paste0("output_path_sx_prof: ", output_path_sx_prof))


run_decontx <- function(input_path_outs, input_path_clusters = "", output_path_dx_est) {
  input_path_h5 <- file.path(input_path_outs, "filtered_feature_bc_matrix.h5")
  print(paste0("decontx: reading counts from ", input_path_h5))
  sce <- DropletUtils::read10xCounts(input_path_h5)
  print("input_path_clusters: ", )
  if (length(input_path_clusters) > 4) {
    clusters <- read.csv(input_path_clusters)
  } else {
    clusters <- NULL
  }
  sce <- celda::decontX(sce, z = clusters)
  print("decontx: writing estimates")
  write.csv(SummarizedExperiment::colData(sce), output_path_dx_est)
  dxc <- sce@assays@data$decontXcounts
  # TODO: currently outputs h5 with no gene names but ensembl ids
  print("decontx: writing matrix")
  DropletUtils::write10xCounts("/data/decontx_counts.h5", dxc, barcodes = sce@colData$Barcode, overwrite = T)
  return(sce)
}

run_soupx <- function(input_dir_10x, clusters, output_path_sx_est, output_path_sx_prof) {
  sx <- SoupX::load10X(input_dir_10x)
  sx <- SoupX::setClusters(sx, clusters)
  sx <- SoupX::autoEstCont(sx)
  print(paste0("soupx: saving estimates to ", output_path_sx_est))
  write.csv(sx$metaData, output_path_sx_est)
  print(paste0("soupx: saving profile to ", output_path_sx_prof))
  write.csv(sx$soupProfile, output_path_sx_prof)
  print("soupx: correcting")
  sx_counts <- SoupX::adjustCounts(sx)
  print("soupx: writing matrix")
  DropletUtils::write10xCounts("/data/soupx_counts.h5", sx_counts, overwrite = T)
  # return(c(sx, sx_counts))
}

run_decontx <- function(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_sx_est, output_path_sx_prof) {
  sce_dx <- run_decontx(input_dir_10x, input_path_clusters, output_path_dx_est)
}

run_soupx <- function(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_sx_est, output_path_sx_prof) {
  run_soupx(input_dir_10x, input_path_clusters, output_path_sx_est, output_path_sx_prof)   # list[sx, sx_counts] <-
}


tryCatch( {
            run_decontx(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_sx_est, output_path_sx_prof)
          },
          error = function(e) {
            print("Decontx Failed")
          }
        )

tryCatch( {
            run_soupx(input_dir_10x, input_path_clusters, output_path_dx_est, output_path_sx_est, output_path_sx_prof)
          },
          error = function(e) {
            print("Soupx Failed")
          }
        )

