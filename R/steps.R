#' Load from a tool to a R package.
#'
#' Currently only supports from Cellranger to Seurat.
#'
#' @param path_to_sample Path to the sample directory from the counting tool. Path to the matrix inside the sample directory will be detected automatically for Cellranger.
#' @param sample Name of the sample. Will be added in the metadata of the returned object.
#' @param from Which tool does the count matrix come from. Currently only Cellranger supported.
#' @param to To which package should the matrix be converted to. Currently only Seurat supported.
#'
#' @return An object from the package in the argument to. By default will return a Seurat object.
#' @export
load_data <- function(path_to_sample, sample, from = "Cellranger", to = "Seurat") {
  checkmate::assert_directory_exists(path_to_sample)
  checkmate::assert_string(sample)
  checkmate::assert_string(from)
  checkmate::assert_string(to)
  if (from == "Cellranger" & to == "Seurat") {
    path_to_matrix = file.path(path_to_sample, "outs", "filtered_feature_bc_matrix")
    analysis = seurat_load_10x(path_to_matrix, project=sample)
  } else {
    stop(paste0("Unsupported conversion from ", from, " to ", to, "."))
  }

  return(analysis)
}


#' Filter out genes and cells outside the min and max limits for several parameters.
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to save the images.
#' @param method The analysis method. Default Seurat
#' @param organism An organism supported by the library for automatic mitochondrial detection. Currently supported : human and mouse.
#' Default "" to skip automatic mitochondrial detection.
#' @param mitochondrial_genes A vector of mitochondrial genes. If a supported organism is requested, it will be skipped. Default empty
#' @param min_genes All cells having a lower number of genes expressed will be filtered out. Default 0
#' @param max_genes All cells having a higher number of genes expressed will be filtered out. Default Inf
#' @param min_cells All genes having a lower number of cells expressing it will be filtered out. Default 0
#' @param max_cells All genes having a higher number of cells expressing it will be filtered out. Default Inf
#' @param min_counts All cells having a lower number of counts will be filtered out. Default 0
#' @param max_counts All cells having a higher number of counts will be filtered out. Default Inf
#' @param min_mt All cells having a lower percentage of mitochondria will be filtered out. Default 0
#' @param max_mt All cells having a higher percentage of mitochondria will be filtered out. Default Inf
#' @param plots_dir The path to save the plots. Keep empty to skip plot saving. Default ""
#'
#' @return An object of the analysis type filtered.
#' @export
filter_data <- function(analysis, sample, assay = "RNA", method = "Seurat", organism = "",
                        mitochondrial_genes = c(), min_genes = 100, min_counts = 0,
                        min_cells = 1, min_mt = 0, max_genes = Inf, max_counts=Inf,
                        max_cells = Inf, max_mt = 20, plots_dir = "") {

  checkmate::assert_string(sample)
  checkmate::assert_string(method)
  checkmate::assert_string(organism, na.ok = TRUE)
  checkmate::assert_vector(mitochondrial_genes, null.ok = TRUE)
  checkmate::assert_numeric(min_genes)
  checkmate::assert_numeric(min_counts)
  checkmate::assert_numeric(min_cells)
  checkmate::assert_numeric(min_mt)
  checkmate::assert_numeric(max_genes)
  checkmate::assert_numeric(max_counts)
  checkmate::assert_numeric(max_cells)
  checkmate::assert_numeric(max_mt)
  checkmate::assert_string(plots_dir, na.ok = TRUE)

  if (method == "Seurat") {
    check_assay(analysis, assay)
    analysis <- seurat_compute_mt(analysis, assay, organism, mitochondrial_genes)
    df_plot <- data.frame(samples = analysis@meta.data$orig.ident,
                          nCount_RNA = analysis@meta.data$nCount_RNA,
                          nFeature_RNA = analysis@meta.data$nFeature_RNA)
    list_plot <- list()
    list_plot[["plot_count"]] <- plot_filter(df_plot, x_name = "samples", y_name = "nCount_RNA",
                                             low = min_counts, high = max_counts)
    list_plot[["plot_feature"]] <- plot_filter(df_plot, x_name = "samples", y_name = "nFeature_RNA",
                                               low = min_genes, high = max_genes)
    if ("percent_mt" %in% colnames(analysis@meta.data)) {
      df_plot$percent_mitochondria <- analysis@meta.data$percent_mt
      list_plot[["plot_mitochondria"]] <- plot_filter(df_plot, x_name = "samples", y_name = "percent_mitochondria",
                                                      low = min_mt, high = max_mt)
    }
    if (checkmate::check_directory_exists(plots_dir) == TRUE) {
      filename = file.path(plots_dir, paste0(sample, "_filter_plots.pdf"))
      pdf(filename, onefile = TRUE)
      for (i in names(list_plot)) {
        print(list_plot[[i]])
      }
      dev.off()
    } else {
      print(paste0("Directory ", plots_dir, " does not exist, not saving images and showing them to screen."))
      for (i in names(list_plot)) {
        print(list_plot[[i]])
      }
    }
    analysis <- seurat_filter(analysis, assay, min_genes = min_genes, min_counts = min_counts,
                              min_cells = min_cells, min_mt = min_mt, max_genes = max_genes,
                              max_counts = max_counts, max_cells = max_cells, max_mt = max_mt)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Perform the data normalization using normalize, variable features and scaling.
#'
#' @param analysis
#' @param method
#' @param assay
#' @param nfeatures
#' @param selection_method
#' @param features
#'
#' @return
#' @export
normalize_data <- function(analysis, method = "Seurat", assay = "RNA", nfeatures = 2000,
                           selection_method = "vst", features = NULL) {
  checkmate::assert_string(method)
  checkmate::assert_int(nfeatures, lower = 0)
  checkmate::assert_string(method)
  checkmate::assert_vector(features, null.ok = TRUE)
  if (method == "Seurat") {
    check_assay(analysis, assay)
    analysis <- seurat_normalize(analysis, assay)
    analysis <- seurat_features(analysis, assay, nfeatures, selection_method)
    analysis <- seurat_scale(analysis, assay, features)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}


