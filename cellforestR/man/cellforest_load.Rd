% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat_loader.R
\name{cellforest_load}
\alias{cellforest_load}
\title{Interface for load RDS matrix and embeddings at process run}
\usage{
cellforest_load(root_dir, branch_spec, process)
}
\arguments{
\item{root_dir}{CellForest root directory}

\item{branch_spec}{Specification for CellBranch in Python or R (see \code{example_spec.R})}

\item{process}{Process name (alias, if exists) at which Seurat object should be loaded}
}
\value{
Seurat object with cached dimensionality reduction embeddings
}
\description{
Wrapper around get_seurat_object(), main interface for creating
Seurat object from a CellForest process run
}
\section{Installation}{

\enumerate{
\item \code{setwd(PATH_TO_CELLFOREST)}
\item \code{install("cellforestR")}
}
}

\examples{
library(cellforestR)

root_dir <- "tests/data/example_usage/root"
seurat_obj <- cellforest_load(root_dir, example_spec_r, "reduce")  # given that processes have been run

DimPlot(seurat_obj, reduction = "umap")
}
