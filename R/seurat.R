#' Check that assay is in seurat
#'
#' Check if object is of class seurat, check that assay is a string, and check that assay is present in seurat
#'
#' @param seurat an object to check
#' @param assay an assay to chack
#'
#' @return nothing
check_assay <- function(seurat, assay="RNA") {
  checkmate::assert_class(seurat, "Seurat")
  checkmate::assert_string(assay)
  if (!assay %in% names(seurat@assays)) {
    stop(paste0("Assay ", assay, " not in seurat object."))
  }
}

check_reduction <- function(seurat, reduction = "pca") {
  checkmate::assert_class(seurat, "Seurat")
  checkmate::assert_string(reduction)
  if (!reduction %in% names(seurat@reductions)) {
    stop(paste0("Reduction ", reduction, " not in seurat object."))
  }
}

#' load Cellranger data for Seurat
#'
#' @param path_to_matrix The path to a cellranger matrix directory
#' @param project The project name, put here the name of the sample currently being analyzed
#'
#' @return A Seurat object containing the data from the matrix
#' @export
seurat_load_10x <- function(path_to_matrix, project = "seurat") {
  checkmate::assert_directory_exists(path_to_matrix)
  checkmate::assert_string(project)

  mtx <- Seurat::Read10X(path_to_matrix)
  seurat <- Seurat::CreateSeuratObject(counts = mtx, project = project)

  return(seurat)
}

#' Get the mitochondrial percentage per cell from a list of genes in a Seurat object.
#'
#' If both organism and mitochondrial_genes are empty, no mitochondrial QC will be performed.
#'
#' @param seurat A Seurat object
#' @param assay Which assay to use. Default RNA
#' @param organism An organism supported by the library. Currently supported : human and mouse.
#' "" to skip automatic mitochondrial detection. Default ""
#' @param mitochondrial_genes A vector of mitochondrial genes. If a supported organism is requested,
#' it will be skipped. Default NULL
#'
#' @return A Seurat object
#' @export
seurat_compute_mt <- function(seurat, assay = "RNA", organism="", mitochondrial_genes=NULL) {
  check_assay(seurat, assay)
  checkmate::assert_string(organism, na.ok = TRUE)
  checkmate::assert_vector(mitochondrial_genes, null.ok = TRUE)
  if (!organism %in% c("human", "mouse") & (is.null(mitochondrial_genes))) {
    print("Cannot find mitochondrial genes, organism is not human and mouse, and no list has been given")
    return(seurat)
  } else if (organism == "human") {
    mitochondrial_genes = grep("^MT-", row.names(seurat), value = TRUE)
  } else if (organism == "mouse") {
    mitochondrial_genes = grep("^mt-", row.names(seurat), value = TRUE)
  }

  seurat <- Seurat::PercentageFeatureSet(seurat, features = mitochondrial_genes,
                                         col.name = "percent_mt", assay = assay)
  return(seurat)
}

#' Apply filters on Cells and genes in a Seurat object
#'
#' If no percent_mt is found in meta data, min_mt and max_mt will be skipped.
#'
#' @param seurat The seurat object
#' @param assay Which assay int he seurat object to filter on. Default RNA
#' @param min_genes All cells having a lower number of genes expressed will be filtered out. Default 0
#' @param max_genes All cells having a higher number of genes expressed will be filtered out. Default Inf
#' @param min_cells All genes having a lower number of cells expressing it will be filtered out. Default 0
#' @param max_cells All genes having a higher number of cells expressing it will be filtered out. Default Inf
#' @param min_counts All cells having a lower number of counts will be filtered out. Default 0
#' @param max_counts All cells having a higher number of counts will be filtered out. Default Inf
#' @param min_mt All cells having a lower percentage of mitochondria will be filtered out. Default 0
#' @param max_mt All cells having a higher percentage of mitochondria will be filtered out. Default Inf
#'
#' @return A Seurat object with its count matrix filtered.
#' @export
seurat_filter <- function(seurat, assay = "RNA", min_genes = 100, min_counts = 100,
                          min_cells = 1, min_mt = 0, max_genes = Inf, max_counts=Inf,
                          max_cells = Inf, max_mt = Inf) {
  check_assay(seurat, assay)
  checkmate::assert_numeric(min_genes)
  checkmate::assert_numeric(min_counts)
  checkmate::assert_numeric(min_cells)
  checkmate::assert_numeric(min_mt)
  checkmate::assert_numeric(max_genes)
  checkmate::assert_numeric(max_counts)
  checkmate::assert_numeric(max_cells)
  checkmate::assert_numeric(max_mt)

  n_cells = length(Seurat::Cells(seurat))
  n_genes = nrow(seurat@assays[[assay]]@counts)

  if (min_cells > 0 | max_cells < length(SeuratObject::Cells(seurat))) {
    counts <- as.matrix(SeuratObject::GetAssayData(seurat, slot = "counts", assay = assay))
    genes_count <- rowSums(counts > 0)
    genes_filter <- names(genes_count[genes_count >= min_cells & genes_count <= max_cells])
    if (Seurat::DefaultAssay(seurat) != assay) {
      temp_assay <- Seurat::DefaultAssay(seurat)
      Seurat::DefaultAssay(seurat) <- assay
      seurat = subset(seurat, features = genes_filter)
      Seurat::DefaultAssay(seurat) <- temp_assay
    } else {
      seurat = subset(seurat, features = genes_filter)
    }
  }

  seurat = subset(seurat, subset = (nFeature_RNA >= min_genes & nFeature_RNA <= max_genes &
                                      nCount_RNA >= min_counts & nCount_RNA <= max_counts))

  if ("percent_mt" %in% colnames(seurat@meta.data)) {
    seurat = subset(seurat, subset = (percent_mt >= min_mt & percent_mt <= max_mt))
  }

  n_cells2 = length(Seurat::Cells(seurat))
  n_genes2 = nrow(seurat@assays$RNA@counts)
  print(data.frame(Removed = c(n_cells-n_cells2, n_genes-n_genes2), Kept = c(n_cells2, n_genes2),
                   row.names = c("Cells", "Genes")))
  return(seurat)
}

#' Seurat NormalizeData
#'
#' @param seurat Seurat object to use.
#' @param assay Assay to normalize. Default RNA
#'
#' @return A Seurat object with normalized data.
#' @export
seurat_normalize <- function(seurat, assay = "RNA") {
  check_assay(seurat, assay)
  seurat <- Seurat::NormalizeData(seurat, assay = assay)
  return(seurat)
}

#' Run Seurat::FindVariableFeatures with nfeatures and selected method
#'
#' @param seurat Seurat object.
#' @param assay Which assay to use.
#' @param nfeatures Number of variable genes to select. Default 2000
#' @param method What to use for features selection. Default vst
#'
#' @return A Seurat object with variable features selected.
#' @export
seurat_features <- function(seurat, assay = "RNA", nfeatures = 2000, method = "vst") {
  check_assay(seurat, assay)
  checkmate::assert_int(nfeatures, lower = 0)
  checkmate::assert_string(method)
  seurat <- Seurat::FindVariableFeatures(seurat, selection.method = method, nfeatures = nfeatures, assay = assay)
  return(seurat)
}

#' Run Seurat::ScaleData
#'
#' @param seurat The Seurat object.
#' @param assay Assay to scale. Default RNA
#' @param features Features to scale. Default all.genes
#'
#' @return A seurat object with scaled features.
#' @export
seurat_scale <- function(seurat, assay = "RNA", features = NULL) {
  check_assay(seurat, assay)
  checkmate::assert_character(features, null.ok = TRUE)
  seurat <- Seurat::ScaleData(seurat, features = features, assay = assay)
  return(seurat)
}

#' Compute the PCA using RunPCA from Seurat
#'
#' @param seurat The Seurat object.
#' @param assay Assay to use. Default RNA
#' @param npcs Number of components to compute. Default 50
#'
#' @return A seurat object with PCA.
#' @export
seurat_pca <- function(seurat, assay = "RNA", npcs = 50) {
  check_assay(seurat, assay)
  checkmate::assert_int(npcs, lower = 2)
  seurat <- Seurat::RunPCA(seurat, assay = assay, npcs = npcs)
}


#' Integrate a list of Seurat objects together.
#'
#' @param analysis_list A list of several objects to integrate together. These objects must all be of the same type.
#' @param nfeatures The number of features to use in integration selection. Default 5000
#' @param assay Which assay to use. Default "RNA"
#'
#' @return A Seurat object with all the analysis from the list integrated together.
#' @export
seurat_integrate <- function(seurat_list, nfeatures = 5000, assay = "RNA") {
  checkmate::assert_int(nfeatures)
  for (i in seurat_list) {
    check_assay(i, assay)
  }
  seurat_list <- lapply(seurat_list, function(x) {
    Seurat::RenameCells(x, new.names = paste0(x@meta.data[["orig.ident"]], "_", Seurat::Cells(x)))
  })
  features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nfeatures)
  anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, assay = rep_len(assay, length(seurat_list)))

  seurat_integrated <- Seurat::IntegrateData(anchorset = anchors)
  return(seurat_integrated)
}


























