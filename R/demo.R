#' Get demo seurat object 
#'
#' @return Seurat object
#'
#' @examples
#' pbmc <- get_demo_seurat_object()
#'
#' @export
get_demo_seurat_object <- function() {
  pbmc <- readRDS(system.file("extdata", "pbmc.rds", package = "SingleCell"))
  pbmc[["percent_mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-", assay = "RNA")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 700 & percent_mt < 5)
  return(pbmc)
}

get_demo_seurat_object_pca <- function() {
  pbmc <- get_demo_seurat_object()
  
  pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::RunPCA(pbmc)
  
  return(pbmc)
}

get_demo_seurat_object_clust <- function() {
  pbmc <- get_demo_seurat_object()
  
  pbmc <- Seurat::NormalizeData(pbmc)
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::RunPCA(pbmc)
  pbmc <- Seurat::FindNeighbors(pbmc, k.param = 20, graph.name = "RNA_snn")
  pbmc <- Seurat::FindClusters(seurat, resolution = 0.5)
  
  return(pbmc)
}

get_demo_seurat_object_umap <- function() {
  pbmc <- get_demo_seurat_object()
  
  pbmc <- Seurat::NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::RunPCA(pbmc)
  pbmc <- Seurat::FindNeighbors(pbmc, k.param = 20, graph.name = "RNA_snn")
  pbmc <- Seurat::FindClusters(pbmc, resolution = 0.5, prefix = "RNA_snn")
  pbmc <- Seurat::RunUMAP(pbmc, n.neighbors = 30, n.components = 2, dims = 1:20)
  
  return(pbmc)
  
}