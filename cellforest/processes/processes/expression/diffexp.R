library(future)
library(dplyr)
library(Seurat)
plan("multiprocess", workers = 6)

args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
output_diffexp_path <- args[2]
root_dir <- args[3]
spec_str <- args[4]
test <- args[5]
logfc_thresh <- as.numeric(args[6])
ident1 <- args[7]
ident2 <- args[8]
groupby <- args[9]

r_functions_filepath <- args[10]
source(r_functions_filepath)
library("cellforestR")


print("loading metadata"); print(date())
meta <- read.table(input_metadata_path, sep = "\t", header = TRUE, row.names = 1)
print("cellforestR loading seurat object"); print(date())
srat <- cellforest_load(root_dir, spec_str, "cluster")
print("metadata filter"); print(date())
srat <- metadata_filter_objs(meta, srat)

Idents(srat) <- srat[[groupby]]
group_names <- unique(Idents(srat))
print(paste0(length(group_names), " clusters exist: ", group_names)); print(date())

print(paste0("Identifying markers with logfc_thresh: ", logfc_thresh)); print(date())

datalist <- list()
i <- 1
for (value in group_names) {
  #print("cluster ", value)

  out <- tryCatch(
    {
      markers <- FindMarkers(
        srat,
        ident.1 = ident1,
        ident.2 = ident2,
        group.by = groupby,
        test.use = test,
        logfc.threshold = logfc_thresh
        # subset.ident = value
      )
      markers$group <- value
      markers$gene_symbol <- rownames(markers)
      rownames(markers) <- NULL
      datalist[[i]] <- markers
      i <- i + 1
    },
    error = function(e) e
    # {
    # print(paste0("ERROR on cluster ", cond)); print(date())
    # markers <- data.frame()
    # print("error")
    # return(NA)
    # }
  )
}

all_markers <- do.call(rbind, datalist)
write.table(all_markers, sep = "\t", file = output_diffexp_path, quote = FALSE, row.names = FALSE)
