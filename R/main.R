#' Integrate data from a list and samples and their matrix path.
#'
#' @param samples A list of the sample names, used for matrix detection and plots saving. This parameter is mandatory
#' @param path_to_count Function need either that if step is 1 or `analysis_list` if step is 2 to 5.
#' The path to where the directories for each sample are. The function will automatically
#' detect the structure in each directory to find the matrix, using the `from` parameter. Default ""
#' @param analysis_list Function need either that if step is 2 to 5 or `analysis_list` if step is 1.
#' A list of all the objects used for the integration. Default list()
#' @param step Steps range from 1 to 5, with each number making the function start at one of the following steps :
#' 1 : `load data`, 2 : `filter data`, 3 : `normalize data`, 4 : `compute pca`, 5 : `integrate`. Choosing 1 will do the entire process,
#' while 5 will only do the integration and assume that every previous step needed has already been done by the user. Default 1
#' @param from From which tools the count matrix comes from. Currently supports : Cellranger. Default "Cellranger"
#' @param method The analysis method used through the analysis, as well as the objects type in `analysis_list`.
#' Currently supports : Seurat. Default "Seurat"
#' @param assay Which to use for all functions asking for it. Default "RNA"
#' @param plots_dir In which directory the plots are saved. Empty to skip plots saving. Default "."
#' @param save_path In which directory the final analysis is saved, in rds format.
#' Keep empty to just return the analysis and not save it. Default ""
#' @param file_name What name to give to the final file. Will be skipped if `save_path` is empty. Default "analysis.rds"
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
#'
#' @return an analysis object containing the integrated data.
#' @export
integrate <- function(samples,
                           path_to_count = "",
                           analysis_list = list(),
                           step = 1,
                           from = "Cellranger",
                           method = "Seurat",
                           assay = "RNA",
                           plots_dir = ".",
                           save_path = "",
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
                           selection_method_normalize = "vst",
                           npcs_pca = 50,
                           nfeatures_integration = 5000) {

  checkmate::assert_character(samples, min.len = 2)
  checkmate::assert_int(step, lower = 1, upper = 5)
  checkmate::assert_string(path_to_count)
  if (path_to_count == "" & step == 1) {stop("Path_to_count is mandatory when using step 1.")}
  if (path_to_count != "") {checkmate::assert_directory(path_to_count)}
  checkmate::assert_string(from)
  if (!from %in% c("Cellranger")) {stop(paste0(from, " is an unsupported matrix detection method."))}
  checkmate::assert_character(method)
  if (!method %in% c("Seurat")) {stop(paste0(method, " is an unsupported analysis method."))}
  checkmate::assert_list(analysis_list, types = method, null.ok = TRUE)
  if (is.null(analysis_list) & step >= 2) {stop("Analysis_list is mandatory when using steps 2 to 5.")}
  checkmate::assert_string(assay)
  checkmate::assert_string(plots_dir)
  if (plots_dir != "") {checkmate::assert_directory(plots_dir)}
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

  checkmate::assert_int(nfeatures_integration)

  if (step<=1) {
    analysis_list <- lapply(samples, function(x) {load_data(path_to_count,
                                                            x,
                                                            from = from,
                                                            to = method)})
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))

  if (step<=2) {
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
                                                     plots_dir = plots_dir)})
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))

  if (step<=3) {
    analysis_list <- lapply(1:length(samples),
                            function(x) {normalize_data(analysis_list[[x]],
                                                        method = method,
                                                        assay = assay,
                                                        nfeatures = nfeatures_normalize,
                                                        selection_method = selection_method_normalize)})
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))

  if (step<=4) {
    analysis_list <- lapply(1:length(samples),
                            function(x) {pca(analysis_list[[x]],
                                             samples[x],
                                             method = method,
                                             assay = assay,
                                             npcs = npcs_pca,
                                             plots_dir = plots_dir)})
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))

  if (step<=5) {
    analysis <- integrate_data(analysis_list,
                          method = method,
                          nfeatures = nfeatures_integration,
                          assay = assay)
  }
  checkmate::assert_list(analysis_list, types = method, len = length(samples))

  if (save_path != "") {
    file_path = file.path(save_path, file_name)
    saveRDS(analysis, file = file_path)
  }
  return(analysis)
}









