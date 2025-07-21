#' Integrate data from several samples.
#'
#' Takes a list of samples and either a path to a cellranger directory with all the samples within, or a list of Seurat object.
#'
#' @param samples A list of the sample names, used for matrix detection and plots saving. This parameter is mandatory
#' @param path_to_count Function need either that if step is 1 or `analysis_list` if step is 2 to 5.
#' The path to where the directories for each sample are. The function will automatically
#' detect the structure in each directory to find the matrix, using the `from` parameter. Default ""
#' @param analysis_list Function need either that if step is 2 to 5 or `analysis_list` if step is 1.
#' A list of all the objects used for the integration. Default list()
#' @param step Steps ranging from `loading_data` to `integrating`. The following steps are : `loading_data`, `filtering`, `normalizing`, `PCA`, and `integrating`. 
#' Choosing `loading_data` will do the entire process, while `integrating` will only do the integration and assume that every previous step 
#' needed has already been done by the user. Default 1
#' @param from From which tools the count matrix comes from. Currently supports : Cellranger. Default "Cellranger"
#' @param method The analysis method used through the analysis, as well as the objects type in `analysis_list`.
#' Currently supports : Seurat. Default "Seurat"
#' @param assay Which to use for all functions asking for it. Default "RNA"
#' @param save_path In which directory the plots, results and report are saved. Empty to not save anything. Default "."
#' Keep empty to just return the analysis and not save it. Default "."
#' @param file_name What name to give to the final file. Will be skipped if `save_path` is empty. Default "analysis.rds"
#' @param make_report Whether to make a report with the plot generated or not.
#' If `save_path` is empty, this parameter will be skipped. Default TRUE
#' @param organism The organism is used for mitochondrial genes detection. Supports mouse and human.
#' Keep empty to skip mitochondrial detection. Default ""
#' @param mitochondrial_genes List of mitochondrial genes. This argument will be completely ignored
#' if a supported organism is provided. Leave `organism` and `mitochondrial_genes` empty to skip mitochondrial filtering.
#' @param min_genes All cells having a lower number of genes expressed will be filtered out. Default 0
#' @param max_genes All cells having a higher number of genes expressed will be filtered out. Default Inf
#' @param min_cells All genes having a lower number of cells expressing it will be filtered out. Default 0
#' @param max_cells All genes having a higher number of cells expressing it will be filtered out. Default Inf
#' @param min_counts All cells having a lower number of counts will be filtered out. Default 0
#' @param max_counts All cells having a higher number of counts will be filtered out. Default Inf
#' @param min_mt All cells having a lower percentage of mitochondria will be filtered out. Default 0
#' @param max_mt All cells having a higher percentage of mitochondria will be filtered out. Default Inf
#' @param nfeatures_normalize Number of variable genes to select for scaling. Default 2000
#' @param selection_method_normalize What method to use for features selection during scaling. Default "vst"
#' @param npcs_pca Number of principal components to compute in the PCA. Default 50
#' @param nfeatures_integration Number of features to use in the integration process. Default 5000
#' @param force_report Whether to overwrite the report if it already exists. If this is FALSE, and a file `integration.Rmd` already exists,
#' report creation will be skipped. Default FALSE
#' @param k.filter How many neighbors (k) to use when filtering anchors
#' @param k.weight Number of neighbors to consider when weighting anchors
#'
#' @return an analysis object containing the integrated data.
#' @export
#'
#' @importFrom magrittr %>%
integrate <- function(samples,
                      path_to_count = "",
                      analysis_list = list(),
                      step = "loading_data",
                      from = "Cellranger",
                      method = "Seurat",
                      assay = "RNA",
                      save_path = ".",
                      file_name = "integrated.rds",
                      make_report = TRUE,
                      organism = "",
                      mitochondrial_genes = c(),
                      min_genes = 0,
                      min_counts = 0,
                      min_cells = 1,
                      min_mt = 0,
                      max_genes = Inf,
                      max_counts = Inf,
                      max_cells = Inf,
                      max_mt = Inf,
                      nfeatures_normalize = 2000,
                      selection_method_normalize = "vst",
                      npcs_pca = 50,
                      nfeatures_integration = 5000,
                      force_report = FALSE, 
                      k.weight = 100, 
                      k.filter = 100) {
  
  checkmate::assert_character(samples, min.len = 2)
  checkmate::assert_string(step)
  if (!step %in% c("loading_data", "filtering_list", "normalizing_list", "PCA_list", "integrated")) {stop("The step chosen is not in the given list of steps for integration.")}
  checkmate::assert_string(path_to_count)
  if (path_to_count == "" & step == "loading_data") {stop("Path_to_count is mandatory when using step loading_data.")}
  if (path_to_count != "") {checkmate::assert_directory(path_to_count)}
  if (!from %in% c("Cellranger")) {stop(paste0(from, " is an unsupported matrix detection method."))}
  checkmate::assert_character(method)
  if (!method %in% c("Seurat")) {stop(paste0(method, " is an unsupported analysis method."))}
  if (step >= "filtering") {checkmate::assert_list(analysis_list, types = method, null.ok = FALSE, len = length(samples))}
  checkmate::assert_string(assay)
  checkmate::assert_string(save_path)
  
  checkmate::assert_string(organism)
  checkmate::assert_character(mitochondrial_genes, null.ok = TRUE)
  checkmate::assert_int(min_genes)
  checkmate::assert_int(min_counts)
  checkmate::assert_int(min_cells)
  checkmate::assert_int(min_mt)
  checkmate::assert_number(max_genes)
  checkmate::assert_number(max_counts)
  checkmate::assert_number(max_cells)
  checkmate::assert_number(max_mt)
  
  checkmate::assert_int(nfeatures_normalize)
  checkmate::assert_string(selection_method_normalize)
  
  checkmate::assert_int(npcs_pca, lower = 2)
  
  checkmate::assert_int(nfeatures_integration)
  
  if (save_path != "") {
    checkmate::assert_directory(save_path)
    plots_dir <- file.path(save_path, "plots")
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }
    results_dir <- file.path(save_path, "results")
    if (!dir.exists(results_dir)) {
      dir.create(results_dir)
    }
  } else {
    plots_dir <- ""
    results_dir <- ""
  }
  
  
  
  if (step == "loading_data") {
    analysis_list <- lapply(samples, function(x) {load_data(path_to_count,
                                                            x,
                                                            from = from,
                                                            to = method)})
    
    print("finished loading")
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))
  
  if (step == "fitlering_list") {
    analysis_list <- lapply(1:length(samples),
                            function(x) {filter_data(analysis_list[[x]],
                                                     samples[x],
                                                     method = method,
                                                     assay = assay,
                                                     organism = organism,
                                                     mitochondrial_genes = mitochondrial_genes,
                                                     min_genes = min_genes,
                                                     min_counts = min_counts,
                                                     min_cells = min_cells,
                                                     min_mt = min_mt,
                                                     max_genes = max_genes,
                                                     max_counts=max_counts,
                                                     max_cells = max_cells,
                                                     max_mt = max_mt,
                                                     plots_dir = plots_dir,
                                                     results_dir = results_dir)})
    step <- "normalizing_list"
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))
  
  if (step == "normalizing_list") {
    analysis_list <- lapply(1:length(samples),
                            function(x) {normalize_data(analysis_list[[x]],
                                                        method = method,
                                                        assay = assay,
                                                        nfeatures = nfeatures_normalize,
                                                        selection_method = selection_method_normalize)})
    perform_normalization <- FALSE
    step <- "PCA_list"
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))
  
  if (step == "PCA_list") {
    analysis_list <- lapply(1:length(samples),
                            function(x) {pca(analysis_list[[x]],
                                             samples[x],
                                             method = method,
                                             assay = assay,
                                             npcs = npcs_pca,
                                             plots_dir = plots_dir)})
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))
  
  if (step == "integrating") {
    analysis <- integrate_data(analysis_list,
                               method = method,
                               nfeatures = nfeatures_integration,
                               assay = assay, 
                               k.weight = k.weight, k.filter = k.filter)
    step <- "filtering"
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))
  
  if (results_dir != "" & file_name != "") {
    file_path = file.path(results_dir, file_name)
    saveRDS(analysis, file = file_path)
  }
  
  make_integration_report(samples = samples, report_path = save_path, report_name = "integration.Rmd", plots_relative_path = "plots", data_relative_path = "results", force = force_report)
  return(analysis)
}


#' Analyze and report the Seurat Single Cell pipeline.
#'
#' @param analysis A Seurat object to use consisting of multiple assays.
#' @param sample Sample name. Default "integrated".
#' @param step Steps ranging from `filtering` to `UMAP`. The following steps are : `filtering`, `normalizing`, `PCA`, `finding_neighbors`, `finding_clusters`, and `UMAP`. 
#' Choosing `filtering` will do the entire process, while `integrating` will only do the integration and assume that every previous step 
#' needed has already been done by the user. Default "filtering"
#' @param method The analysis method used through the analysis, as well as the objects type in `analysis_list`.
#' Currently supports : Seurat. Default "Seurat"
#' @param assay Which to use for all functions asking for it. Default "RNA"
#' @param save_path In which directory the final analysis is saved, in rds format.
#' Keep empty to just return the analysis and not save it. Default "."
#' @param file_name What name to give to the final file. Will be skipped if `save_path` is empty. Default "analysis.rds"
#' @param organism The organism is used for mitochondrial genes detection. Supports mouse and human.
#' Keep empty to skip mitochondrial detection. Default ""
#' @param mitochondrial_genes List of mitochondrial genes. This argument will be completely ignored
#' @param min_genes All cells having a lower number of genes expressed will be filtered out. Default 0
#' @param max_genes All cells having a higher number of genes expressed will be filtered out. Default Inf
#' @param min_cells All genes having a lower number of cells expressing it will be filtered out. Default 0
#' @param max_cells All genes having a higher number of cells expressing it will be filtered out. Default Inf
#' @param min_counts All cells having a lower number of counts will be filtered out. Default 0
#' @param max_counts All cells having a higher number of counts will be filtered out. Default Inf
#' @param min_mt All cells having a lower percentage of mitochondria will be filtered out. If no organism or
#' mitochondrial genes list is offered, this argument will be ignored. Default 0
#' @param max_mt All cells having a higher percentage of mitochondria will be filtered out. If no organism or
#' mitochondrial genes list is offered, this argument will be ignored. Default Inf
#' @param nfeatures_normalize Number of variable genes to select for scaling. Default 2000
#' @param perform_normalization If normalisation should be performed at this stage of the analysis. Default TRUE. 
#' @param selection_method_normalize What method to use for features selection during scaling. Default "vst"
#' @param npcs_pca Number of principal components to compute in the PCA. Default 50
#' @param resolutions_clustree Which resolutions to compute and run the clustree on.
#' If this parameter is NULL, clustree will be skipped. Clustree is an optional part of step finding_clusters. Default c(1:10/10, 5:8/4)
#' @param resolution_clustering Which final resolution to keep (can be a resolution not in `resolutions_clustree`).
#' Typically resolutions range between 0.1 and 2. Default 1
#' @param k.param_neighbors Number of neighbors to use when computing UMAP. Default 20.
#' @param n_neighbors Number of neighbors to use when computing the UMAP. Default 30
#' @param de_test Which statistical tests to use for DE. To see available options, see the documentation for Seurat::FindAllMarkers. Default "wilcox"
#' @param de_logfc Filter genes based on a minimum logged fold change. Default 0.25
#' @param de_pvalue Filter genes based on a maximum pvalue. Default 0.05
#' @param force_report Whether to overwrite the report if it already exists. If this is FALSE, and a file `analysis.Rmd` already exists,
#' report creation will be skipped. Default FALSE
#' @param variable The column name that is going to be used for extracting the cluster number per barcode. Default "seurat_clusters"
#' @param finding_DEG Whether to force the DE to be recomputed if a DE file already exists. Default FALSE
#' @param skip Which steps to skip, keep empty to do all the steps. Default c()
#'
#' @return An analysis object.
#' @export
#'
#' @examples
#'
#' @importFrom magrittr %>%
analyze_integrated <- function(analysis,
                               sample = "integrated",
                               step = "filtering",
                               method = "Seurat",
                               assay = "RNA",
                               save_path = ".",
                               file_name = "analysis.rds",
                               organism = "",
                               mitochondrial_genes = c(),
                               min_genes = 100,
                               min_counts = 100,
                               min_cells = 1,
                               min_mt = 0,
                               max_genes = Inf,
                               max_counts = Inf,
                               max_cells = Inf,
                               max_mt = 20,
                               nfeatures_normalize = 2000,
                               perform_normalization = TRUE,
                               selection_method_normalize = "vst",
                               npcs_pca = 50,
                               k.param_neighbors = 20,
                               resolutions_clustree = c(1:10/10, 5:8/4),
                               resolution_clustering = 1,
                               n_neighbors = 30,
                               de_test = "wilcox",
                               de_logfc = 0.25,
                               de_pvalue = 0.05,
                               force_report = FALSE,
                               variable = "seurat_clusters",
                               finding_DEG = FALSE,
                               skip = NULL) {
  
  checkmate::assert_string(step)
  if (!step %in% c("filtering", "normalizing", "PCA", "finding_neighbors", "finding_clusters", "UMAP")) {stop("The step chosen is not in the given list of steps for clusterisation.")}
  
  checkmate::assert_string(method)
  if (!method %in% c("Seurat")) {stop(paste0(method, " is an unsupported analysis method."))}
  checkmate::assert_class(analysis, method)
  checkmate::assert_string(assay)
  checkmate::assert_string(save_path)
  if (save_path != "") {checkmate::assert_directory(save_path)}
  
  checkmate::assert_string(organism)
  checkmate::assert_character(mitochondrial_genes, null.ok = TRUE)
  checkmate::assert_int(min_genes)
  checkmate::assert_int(min_counts)
  checkmate::assert_int(min_cells)
  checkmate::assert_int(min_mt)
  checkmate::assert_number(max_genes)
  checkmate::assert_number(max_counts)
  checkmate::assert_number(max_cells)
  checkmate::assert_number(max_mt)
  
  checkmate::assert_int(nfeatures_normalize)
  checkmate::assert_string(selection_method_normalize)
  
  checkmate::assert_int(npcs_pca, lower = 2)
  
  checkmate::assert_int(k.param_neighbors, lower = 2)
  
  checkmate::assert_vector(resolutions_clustree, min.len = 1)
  checkmate::assert_number(resolution_clustering, lower = 0)
  checkmate::assert_logical(finding_DEG)
  checkmate::assert_string(skip, null.ok = TRUE)
  
  if (save_path != "") {
    checkmate::assert_directory(save_path)
    plots_dir <- file.path(save_path, "plots")
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }
    results_dir <- file.path(save_path, "results")
    if (!dir.exists(results_dir)) {
      dir.create(results_dir)
    }
  } else {
    plots_dir <- ""
    results_dir <- ""
  }
  
  if (!is.null(step) && is.character(step) && step == "filtering") {
    analysis <- filter_data(analysis,
                            sample,
                            method = method,
                            assay = assay,
                            organism = organism,
                            mitochondrial_genes = mitochondrial_genes,
                            min_genes = min_genes,
                            min_counts = min_counts,
                            min_cells = min_cells,
                            min_mt = min_mt,
                            max_genes = max_genes,
                            max_counts=max_counts,
                            max_cells = max_cells,
                            max_mt = max_mt,
                            plots_dir = plots_dir,
                            results_dir = results_dir)
    checkmate::assert_class(analysis, method) 
    
    step <- "normalizing"
  }
  
  if (!is.null(step) && is.character(step) && step == "normalizing") {
    analysis <- normalize_data(analysis,
                               method = method,
                               assay = assay,
                               nfeatures = nfeatures_normalize,
                               selection_method = selection_method_normalize, 
                               perform_normalization = perform_normalization)
    checkmate::assert_class(analysis, method)
    
    step <- "PCA"
  }
  
  if (!is.null(step) && is.character(step) && step == "PCA") {
    analysis <- pca(analysis,
                    sample,
                    method = method,
                    assay = assay,
                    npcs = npcs_pca,
                    plots_dir = plots_dir)
    checkmate::assert_class(analysis, method)
    
    step <- "finding_neighbors"
  }
  
  if (!is.null(step) && is.character(step) && step == "finding_neighbors") {
    analysis <- neighbors(analysis,
                          method = method,
                          k.param = k.param_neighbors)
    checkmate::assert_class(analysis, method)
    
    step <- "finding_clusters"
  }
  
  if (!is.null(step) && is.character(step) && step == "finding_clusters") {
    analysis <- clustering(analysis,
                           sample = sample,
                           method = method,
                           res_clustree = resolutions_clustree,
                           resolution = resolution_clustering,
                           plots_dir = plots_dir)
    checkmate::assert_class(analysis, method)
    
    step <- "UMAP"
  }
  
  if (!is.null(step) && is.character(step) && step == "UMAP") {
    analysis <- umap(analysis,
                     sample = sample,
                     method = method,
                     n_neighbors = n_neighbors,
                     plots_dir = plots_dir,
                     plot_clustering = paste0("RNA_snn_res.", resolution_clustering))
  }
  
  if (finding_DEG) {
    de_genes <- find_all_DE(analysis,
                            sample = sample,
                            method = method,
                            test = de_test,
                            logfc_threshold = de_logfc,
                            pvalue_threshold = de_pvalue,
                            results_dir = results_dir,
                            variable = variable)
  }
  
  if (save_path != "") {
    file_path = file.path(results_dir, file_name)
    saveRDS(analysis, file = file_path)
  }
  make_analysis_report(sample = sample, report_path = save_path, report_name = "analysis.Rmd", plots_relative_path = "plots", data_relative_path = "results", force = force_report)
  
  return(analysis)
}
