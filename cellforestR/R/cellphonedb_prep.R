#' Scale down metadata by a downsample factor
#'
#' @param seurat_obj Seurat object
#' @param downsample_factor By how many times to reduce number of cells per group
#' @param rounding Rounding after devision by `downsample_factor` ("up", "down" or NULL for standard rounding)
#' @param group_by Down-sample per defined group (usually, column name that indicates cell type)
#'
#' @return Row-wise subset of the given metadata dataframe
#'
#' @export
downsample_metadata <- function(seurat_obj, downsample_factor = 100, rounding = "up", group_by = "Subclass_Cell_Identity") {
  downsampled_cell_type_counts <- summary(seurat_obj@meta.data[[group_by]]) / downsample_factor
  if (rounding == "up") {
    num_samples_per_cell_type <- ceiling(downsampled_cell_type_counts)
  } else if (rouding == "down") {
    num_samples_per_cell_type <- floor(downsampled_cell_type_counts)
  } else {
    num_samples_per_cell_type <- round(downsampled_cell_type_counts)
  }

  list_of_downsampled_grouped_meta <- list()
  i <- 1
  for (group_name in names(num_samples_per_cell_type)) {
    grouped_meta <- seurat_obj@meta.data %>% filter(seurat_obj@meta.data[[group_by]] == group_name)

    n_select <- num_samples_per_cell_type[group_name]
    downsampled_grouped_meta <- grouped_meta[sample(nrow(grouped_meta), n_select), ]

    list_of_downsampled_grouped_meta[[i]] <- downsampled_grouped_meta
    i <- i + 1
  }

  downsampled_meta <- do.call("rbind", list_of_downsampled_grouped_meta)
  downsampled_meta <- downsampled_meta[group_by]

  return(downsampled_meta)
}

#' Downsample counts based on downsampled meta
#'
#' @param seurat_obj Seurat object
#' @param downsampled_meta Metadata dataframe, output of downsample_metadata()
#' @param features Dataframe of mapping from Ensembl ID to gene name
#'
#' @return Column-wise (barcode-wise) subset of counts matrix
#'
#' @export
downsample_counts <- function(seurat_obj, downsampled_meta, features = NULL) {
  barcodes <- rownames(downsampled_meta)
  downsampled_seur_obj <- subset(seurat_obj, cells = barcodes)

  downsampled_counts <- as.matrix(GetAssayData(object = downsampled_seur_obj, slot = "counts"))  # make dense matrix
  rownames(downsampled_counts) <- mapvalues(rownames(downsampled_counts), from = features$V2, to = features$V1)  # get ENSEMBL IDs mapping

  return(downsampled_counts)
}

#' Save prepared counts matrix
#'
#' @param downsampled_counts Counts matrix
#' @param save_dir Folder into which to save counts.tsv
#'
#' @export
save_counts <- function(downsampled_counts, save_dir) {
  write.table(downsampled_counts, paste0(save_dir, "/counts.tsv"),
              sep = "\t", row.names = TRUE, col.names=NA, quote = F)
}

#' Save prepared metadata
#'
#' @param downsampled_metadata metadata dataframe
#' @param save_dir Folder into which to save meta.tsv

#' @export
save_meta <- function(downsampled_metadata, save_dir) {
  write.table(downsampled_metadata, paste0(save_dir, "/meta.tsv"),
              sep = "\t", row.names = TRUE, col.names=NA, quote = F)
}

