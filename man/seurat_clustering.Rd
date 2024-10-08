% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seurat.R
\name{seurat_clustering}
\alias{seurat_clustering}
\title{Find the clusters using FindClusters from Seurat}
\usage{
seurat_clustering(seurat, resolution = 1, prefix = "RNA_snn")
}
\arguments{
\item{seurat}{The Seurat object.}

\item{resolution}{Which final resolution to use. Typically resolutions range between 0.1 and 2. Default 1}

\item{prefix}{The prefix to use. Default to "RNA_snn"}
}
\value{
Seurat object with clusters
}
\description{
Find the clusters using FindClusters from Seurat
}
