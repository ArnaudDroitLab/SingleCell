#' Load from a tool to a R package.
#'
#' Currently only supports from Cellranger to Seurat.
#'
#' @param path_to_count The path to where the directories for each sample are.
#' Path to the matrix inside the sample directory will be detected automatically for Cellranger.
#' @param sample Name of the sample. Will be added in the metadata of the returned object.
#' @param from Which tool does the count matrix come from. Currently only Cellranger supported.
#' @param to To which package should the matrix be converted to. Currently only Seurat supported.
#'
#' @return An object from the package in the argument to. By default will return a Seurat object.
#' @export
load_data <- function(path_to_count, sample, from = "Cellranger", to = "Seurat") {
  checkmate::assert_directory_exists(path_to_count)
  checkmate::assert_string(sample)
  checkmate::assert_string(from)
  checkmate::assert_string(to)
  if (from == "Cellranger" & to == "Seurat") {
    path_to_matrix = file.path(path_to_count, sample, "outs", "filtered_feature_bc_matrix")
    analysis = seurat_load_10x(path_to_matrix, project=sample)
  } else {
    stop(paste0("Unsupported conversion from ", from, " to ", to, "."))
  }
  
  checkmate::assert_class(analysis, to)
  return(analysis)
}

#' Filter out genes and cells outside the min and max limits for several parameters.
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to name the images files.
#' @param method The analysis method. Default "Seurat"
#' @param assay Which assay to make the filters on. Default "RNA"
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
#' @param results_dir The path to save the dataframe. Keep empty to skip plot saving. Default ""
#'
#' @return An object of the analysis type filtered.
#' @importFrom ggplot2 ggsave
#' @importFrom ggpubr ggarrange
#' @export
filter_data <- function(analysis, sample = "", method = "Seurat", assay = "RNA", organism = "",
                        mitochondrial_genes = c(), min_genes = 100, min_counts = 100,
                        min_cells = 1, min_mt = 0, max_genes = Inf, max_counts=Inf,
                        max_cells = Inf, max_mt = 20, plots_dir = "", results_dir = "") {
  
  checkmate::assert_string(sample)
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
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
  checkmate::assert_directory(results_dir)
  
  if (method == "Seurat") {
    check_assay(analysis, assay)
    analysis <- seurat_compute_mt(analysis, assay, organism, mitochondrial_genes)
    df_plot <- data.frame(samples = analysis@meta.data$orig.ident,
                          nCount_RNA = analysis@meta.data$nCount_RNA,
                          nFeature_RNA = analysis@meta.data$nFeature_RNA)
    
    list_plot <- list()
    list_plot[["count"]] <- plot_filter(df_plot, x_name = "samples", y_name = "nCount_RNA",
                                        low = min_counts, high = max_counts)
    list_plot[["feature"]] <- plot_filter(df_plot, x_name = "samples", y_name = "nFeature_RNA",
                                          low = min_genes, high = max_genes)
    if ("percent_mt" %in% colnames(analysis@meta.data)) {
      df_plot$percent_mitochondria <- analysis@meta.data$percent_mt
      list_plot[["mitochondria"]] <- plot_filter(df_plot, x_name = "samples", y_name = "percent_mitochondria",
                                                 low = min_mt, high = max_mt)
    }
    
    if (checkmate::check_directory_exists(plots_dir) == TRUE) {
      filename = file.path(plots_dir, paste0(sample, "_filter_plots.pdf"))
      # pdf(filename, onefile = TRUE)
      for (i in names(list_plot)) {
        ggplot2::ggsave(paste(sample, i,  "filter_plot.png", sep = "_"), plot = list_plot[[i]],
                        device = "png", path = plots_dir, dpi = 200, width = 500 + 300*length(unique(df_plot$samples)), # 800 px width is ok for 1 sample
                        height = 1200, units = "px", limitsize = FALSE)
        # print(list_plot[[i]])
      }
      # dev.off()
    } else {
      print(paste0("Directory ", plots_dir, " does not exist, not saving images and showing them to screen."))
      for (i in names(list_plot)) {
        print(list_plot[[i]])
      }
    }
    
    if (length(list_plot) == 2) {
      plot_filter_complete <- ggpubr::ggarrange(list_plot[["count"]], list_plot[["feature"]], common.legend = TRUE, ncol = 2, legend = "right")
      ggplot2::ggsave(paste(sample, "complete_filter_plot.png", sep = "_"), plot_filter_complete,
                      device = "png", path = plots_dir, dpi = 200, width = 500 + 300*2, 
                      height = 1200, units = "px", limitsize = FALSE)
    } else if (length(list_plot) == 3) {
      plot_filter_complete <- ggpubr::ggarrange(list_plot[["count"]], list_plot[["feature"]], list_plot[["mitochondria"]], common.legend = TRUE, ncol = 3, legend = "right")
      ggplot2::ggsave(paste(sample, "complete_filter_plot.png", sep = "_"), plot_filter_complete,
                      device = "png", path = plots_dir, dpi = 200, width = 500 + 300*3, 
                      height = 1200, units = "px", limitsize = FALSE)
    } else {
      warning(paste0("There are no filter plots?"))
    }
    
    # print(results_dir)
    analysis <- seurat_filter(sample = sample, seurat = analysis, assay = assay, min_genes = min_genes, min_counts = min_counts,
                              min_cells = min_cells, min_mt = min_mt, max_genes = max_genes,
                              max_counts = max_counts, max_cells = max_cells, max_mt = max_mt, results_dir = results_dir,
                              plots_dir = plots_dir)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Perform the data normalization using normalize, variable features and scaling.
#'
#' @param analysis The analysis object.
#' @param method The analysis method. Default "Seurat"
#' @param assay Which assay to make the filters on. Default "RNA"
#' @param nfeatures Number of variable genes to select. Default 2000
#' @param selection_method What method to use for features selection. Default vst
#' @param features Vector of features to scale. Keep NULL to scale only variable features. Default NULL
#' @param perform_normalization If normalization should be perform. If integration was performed prior, this should be FALSE. If this feature is runned after integration, this would automatically by FALSE. 
#'
#' @return An analysis object of type method normalized.
#' @export
normalize_data <- function(analysis, method = "Seurat", assay = "RNA", nfeatures = 2000,
                           selection_method = "vst", features = NULL, perform_normalization = FALSE, step = "loading_data") {
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
  checkmate::assert_int(nfeatures, lower = 0)
  checkmate::assert_vector(features, null.ok = TRUE)
  if (method == "Seurat") {
    check_assay(analysis, assay)
    if (step == "normalizing") {
      
      if (perform_normalization) {
        analysis <- seurat_normalize(analysis, assay)
        analysis <- seurat_features(analysis, assay, nfeatures, selection_method)
        analysis <- seurat_scale(analysis, assay, features)
      } else {
        analysis <- seurat_features(analysis, assay, nfeatures, selection_method)
        analysis <- seurat_scale(analysis, assay, features)
      }
    } else {
      analysis <- seurat_normalize(analysis, assay)
      analysis <- seurat_features(analysis, assay, nfeatures, selection_method)
      analysis <- seurat_scale(analysis, assay, features)
    }
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Compute the PCA
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to name the images files.
#' @param method The analysis method. Default "Seurat"
#' @param assay Which assay to compute the PCA on. Default "RNA"
#' @param npcs Number of components to compute. Default 50
#' @param plots_dir The path to save the plots. Keep empty to skip plot saving. Default ""
#'
#' @return An analysis object of type method with PCA.
#' @importFrom ggplot2 ggsave
#' @export
pca <- function(analysis, sample = "", method = "Seurat", assay = "RNA", npcs = 50, plots_dir = "") {
  checkmate::assert_string(method)
  checkmate::assert_string(sample)
  checkmate::assert_class(analysis, method)
  checkmate::assert_int(npcs, lower = 2)
  if (method == "Seurat") {
    check_assay(analysis, assay)
    analysis <- seurat_pca(analysis, assay = assay, npcs = npcs)
    list_plot <- list()
    list_plot[["Elbow"]] <- plot_seurat_elbow(analysis, reduction = "pca", npc = npcs)
    list_plot[["PCA"]] <- plot_seurat_dim(analysis, reduction = "pca", colour_by = "orig.ident")
    if (checkmate::check_directory_exists(plots_dir) == TRUE) {
      for (i in names(list_plot)) {
        ggplot2::ggsave(paste(sample, i,  "pca_plot.png", sep = "_"), plot = list_plot[[i]],
                        device = "png", path = plots_dir, dpi = 200, width = 1500,
                        height = 1200, units = "px", limitsize = FALSE)
      }
    }
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Integrate several analysis together.
#'
#' @param analysis_list A list of several objects to integrate together. These objects must all be of the same type.
#' @param method The analysis method. Default "Seurat"
#' @param nfeatures The number of features to use in integration selection. Default 5000
#' @param assay Which assay to use. Default "RNA"
#'
#' @return An analysis object with all the analysis from the list integrated together.
#' @export
integrate_data <- function(analysis_list, method = "Seurat", nfeatures = 5000, assay = "RNA",
                           k.weight = k.weight, k.filter = k.filter) {
  checkmate::assert_string(method)
  checkmate::assert_int(nfeatures)
  if (method == "Seurat") {
    for (i in analysis_list) {
      check_assay(i, assay)
    }
    analysis <- seurat_integrate(analysis_list, nfeatures = nfeatures, assay = assay, 
                                 k.weight = k.weight, k.filter = k.filter)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Compute neighbors
#'
#' @param analysis The analysis object.
#' @param method The analysis method. Default "Seurat"
#' @param k.param Defines k for the k-nearest neighbor algorithm. Default 20
#'
#' @return An analysis object of type method with neighbors.
#' @export
neighbors <- function(analysis, method = "Seurat", k.param = 20) {
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
  checkmate::assert_int(k.param)
  if (method == "Seurat") {
    analysis <- seurat_neighbors(analysis, k.param = k.param)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Perform the clustering step, and optionally a clustree graph
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to name the images files.
#' @param method The analysis method. Default "Seurat"
#' @param res_clustree Vector of resolutions to use in the clustree.
#' Typically resolutions range between 0.1 and 2. Keep empty to skip. Default c()
#' @param resolution Which final resolution to keep (can be a resolution not in `res_clustree`).
#' Typically resolutions range between 0.1 and 2. Default 1
#' @param plots_dir The path to save the plots. Keep empty to skip plot saving. Default ""
#'
#' @return An analysis object of type method with clusters.
#' @importFrom ggplot2 ggsave
#' @import ggraph 
#' @export
clustering <- function(analysis, sample = "", method = "Seurat", res_clustree = c(), resolution = 1, plots_dir = "") {
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
  checkmate::assert_double(res_clustree, lower = 0, null.ok = TRUE)
  checkmate::assert_double(resolution, len = 1, lower = 0)
  
  if (method == "Seurat") {
    if (length(res_clustree) > 0) {
      for (res in res_clustree) {
        analysis <- seurat_clustering(analysis, resolution = res, prefix = "RNA_snn")
      }
      if (checkmate::check_directory_exists(plots_dir) == TRUE) {
        clustree_plot <- plot_seurat_clustree(analysis, prefix = "RNA_snn_res.")
        ggplot2::ggsave(paste(sample, "clustree.png", sep = "_"), plot = clustree_plot,
                        device = "png", path = plots_dir, dpi = 200, width = 1500,
                        height = 2000, units = "px", limitsize = FALSE)
      }
    }
    analysis <- seurat_clustering(analysis, resolution = resolution)
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}


#' Compute and plot the umap.
#'
#' Will plot a first umap showing the orig.ident,
#' a second one using the column in plot_clustering.
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to name the images files.
#' @param method The analysis method. Default "Seurat"
#' @param n_neighbors number of neighbors to use when computing the UMAP. Default 30
#' @param plot_clustering Which clustering to colour the UMAP with. Default "RNA_snn_res.1"
#' @param plots_dir The path to save the plots. Keep empty to skip plot saving. Default ""
#'
#' @return an analysis object.
#' @importFrom ggplot2 ggsave
#' @importFrom Seurat NoLegend
#' @export
#'
#' @examples
umap <- function(analysis, sample = "", method = "Seurat", n_neighbors = 30, plot_clustering = "RNA_snn_res.1", plots_dir = "") {
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
  checkmate::assert_string(sample)
  checkmate::assert_string(plots_dir)
  checkmate::assert_double(n_neighbors, len = 1, lower = 0)
  if (plots_dir != "") {checkmate::assert_directory(plots_dir)}
  
  if (method == "Seurat") {
    analysis <- seurat_umap(analysis, n.neighbors = n_neighbors)
    if (checkmate::check_directory_exists(plots_dir) == TRUE) {
      list_plot <- list()
      list_plot[["sample"]] <- plot_seurat_dim(analysis, reduction = "umap", colour_by = "orig.ident")
      list_plot[["clusters"]] <- plot_seurat_dim(analysis, reduction = "umap", colour_by = plot_clustering)
      list_plot[["clusters_numbers"]] <- plot_label_umap(analysis, colour_by = plot_clustering, pt.size = 1.5)
      list_plot[["nCount"]] <- plot_seurat_dim(analysis, reduction = "umap", colour_by = "nCount_RNA")
      
      if ("percent_mt" %in% colnames(analysis@meta.data)) {
        list_plot[["mitochondria"]] <- plot_seurat_dim(analysis, reduction = "umap", colour_by = "percent_mt")
      } else {
        print("Mitochondria column in meta.data does not exist")
      }
      
      for (i in names(list_plot)) {
        if (i != "clusters") {
          ggplot2::ggsave(paste(sample, i,  "umap_plot.png", sep = "_"), plot = list_plot[[i]],
                          device = "png", path = plots_dir, dpi = 200, width = 1500,
                          height = 1200, units = "px", limitsize = FALSE)
        }
      }
    }
    
  } else {
    stop(paste0(method, " is an unsupported method."))
  }
  return(analysis)
}

#' Find DE between all clusters
#'
#' @param analysis The analysis object.
#' @param sample Sample name, will be used to name the table files.
#' @param method The analysis method. Default "Seurat"
#' @param test Which statistical tests to use for DE. To see available options, see the documentation for Seurat::FindAllMarkers. Default "wilcox"
#' @param logfc_threshold Filter genes based on a minimum logged fold change. Default 0.25
#' @param pvalue_threshold Filter genes based on a maximum pvalue. Default 0.05
#' @param results_dir The path to save the tables. Keep empty to skip plot saving. Default ""
#' @param variable what to isolate from the FetchData. Default "seurat_clusters"
#'
#' @return A dataframe that contains the different gene expression for each cluster against all other clusters.
#' @export
#'
#' @examples
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @import readr
#' @importFrom dplyr group_by
#' @importFrom dplyr slice_min
#' @importFrom dplyr slice_max
#' @importFrom readr write_csv
#' @importFrom dplyr arrange
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom dplyr rename
#' @importFrom Seurat FetchData

find_all_DE <- function(analysis, sample = "", method = "Seurat",
                        test = "wilcox", logfc_threshold = 0.25, pvalue_threshold = 0.05,
                        results_dir = "", variable = "seurat_clusters", force_DE = FALSE) {
  
  checkmate::assert_string(method)
  checkmate::assert_class(analysis, method)
  checkmate::assert_string(variable)
  checkmate::assert_string(sample)
  checkmate::assert_string(results_dir)
  checkmate::assert_string(test)
  checkmate::assert_number(logfc_threshold, lower = 0)
  checkmate::assert_number(pvalue_threshold, lower = 0, upper = 1)
  checkmate::assert_logical(force_DE)
  
  if (results_dir != "") {checkmate::assert_directory(results_dir)}
  
  path_DE_saved <- file.path(results_dir, paste0(sample, "_DE.csv"))
  if (!force_DE) {
    print("force_DE = FALSE, skipping finding Markers. Congratulations.")
    return()
  }
  
  if (method == "Seurat") {
    df_de <- seurat_all_DE(analysis, sample = sample, assay = "RNA", slot = "data", method = method, logfc_threshold = logfc_threshold,
                           pvalue_threshold = pvalue_threshold, min_pct = 0.1, only.pos = FALSE)
  }
  
  if (results_dir != "") {
    top_de <- df_de %>% dplyr::group_by(cluster) %>% dplyr::slice_min(p_val, n = 10) %>% dplyr::slice_max(avg_log2FC, n = 10)
    
    readr::write_csv(df_de, file.path(results_dir, paste0(sample, "_DE.csv")))
    readr::write_csv(top_de, file.path(results_dir, paste0(sample, "_top10_DE.csv")))
    
    df_stats <- Seurat::FetchData(analysis, vars = variable) %>%
      dplyr::group_by(!!sym(variable)) %>%
      dplyr::summarise(CellCount = n(), .groups = "drop") %>%
      dplyr::arrange(!!sym(variable)) %>%
      dplyr::mutate(Differentially_Expressed = df_de %>%
                      dplyr::group_by(cluster) %>%
                      dplyr::summarise(n = n(), .groups = "drop") %>%
                      dplyr::arrange(cluster) %>%
                      dplyr::pull(n),
                    UpRegulated = df_de %>%
                      dplyr::filter(avg_log2FC > 0) %>%
                      dplyr::group_by(cluster) %>%
                      dplyr::summarise(n = n(), .groups = "drop") %>%
                      dplyr::arrange(cluster) %>%
                      dplyr::pull(n),
                    DownRegulated = df_de %>%
                      dplyr::filter(avg_log2FC < 0) %>%
                      dplyr::group_by(cluster) %>%
                      dplyr::summarise(n = n(), .groups = "drop") %>%
                      dplyr::arrange(cluster) %>%
                      dplyr::pull(n)) %>%
      dplyr::rename(Cluster = variable)
    readr::write_csv(df_stats, file.path(results_dir, paste0(sample,"_summary_per_clusters.csv")))
    
  }
  return(df_de)
}




annotate_clusters <- function(analysis, sample = "", method = "manual", ..., group.by = "orig.ident", results_dir = "", plots_dir = "") {
  checkmate::assert_string(method)
  checkmate::assert_string(sample)
  checkmate::assert_string(group.by, null.ok = FALSE)
  group_by = features_in_seurat(seurat, group.by)
  if (is.null(group_by)) {
    stop("group.by ", group.by, " not found in seurat object.")
  }
  checkmate::assert_string(results_dir)
  if (results_dir != "") {checkmate::assert_directory(results_dir)}
  checkmate::assert_string(plots_dir)
  if (plots_dir != "") {checkmate::assert_directory(plots_dir)}
  
  if (method == "manual") {
    clust_genes <- list(...)
    checkmate::assert_list(clust_genes, any.missing = FALSE, min.len = 1, null.ok = TRUE, types = "character")
    if (is.null(clust_genes)) {
      warning("No signature provided for manual cluster annotation.")
      return()
    }

    for (cell_type in clust_genes) {
      markers <- clust_genes[[cell_type]]
      markers <- features_in_seurat(seurat, markers)
      if (is.null(markers)) {
        warning("None of the genes for cell type ", cell_type, " were found in the Seurat object, skipping violin plots.")
        next
      }
      p <- plot_seurat_violin(seurat, features = markers, group.by = group.by)
      ggplot2::ggsave(paste(sample, cell_type,  "anotation_violin.png", sep = "_"), plot = p,
                       device = "png", path = plots_dir, dpi = 200, width = 1500,
                       height = 75 + 300 * length(markers), units = "px")
    }
  }
}





