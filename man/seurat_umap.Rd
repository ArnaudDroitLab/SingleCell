% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat.R
\name{seurat_umap}
\alias{seurat_umap}
\title{Compute the UMAP using RunUMAP from Seurat}
\usage{
seurat_umap(seurat, n.neighbors = 30, dims = 1:20)
}
\arguments{
\item{seurat}{The Seurat object.}

\item{n.neighbors}{Number of neighbors to use when computing the UMAP. Default 30}

\item{dims}{Number of dimension used. Default 1:20}
}
\value{
Seurat object with UMAP.
}
\description{
Compute the UMAP using RunUMAP from Seurat
}
