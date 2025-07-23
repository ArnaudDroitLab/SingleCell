#' Title
#'
#' @param name 
#' @param report_part 
#' @param force 
#'
#' @return
#' @export
#'
#' @examples
initialize_report_qmd <- function(name, report_path = ".", force = FALSE) {
  
  checkmate::assert_string(name, pattern = "^[a-zA-Z0-9\\-_.]+$")
  checkmate::assert_directory_exists(report_path)
  checkmate::assert_logical(force, len = 1)
  
  template_path <- system.file("qmd/header.qmd", package = "SingleCell")
  internal_files <- system.file("qmd/.report", package = "SingleCell")
  
  internal_path <- file.path(report_path, "_report")
  
  if (!dir.exists(internal_path)) {
    dir.create(internal_path)
  }
  
  report <- file.path(report_path, name)
  
  if (file.exists(report) & force) {
    warning(report, " already exists, overwriting it.")
  } else if (file.exists(report) & !force) {
    stop(report, " already exists. Set `force` to TRUE if you want to overwrite it.")
  }
  
  file.copy(template_path, to = report, overwrite = force)
  file.copy(internal_files, to = report_path, overwrite = FALSE, recursive = TRUE)
}

#' Title
#'
#' @param analysis 
#' @param samples 
#' @param report_steps 
#' @param plots_relative_path 
#' @param report_path 
#' @param file_name 
#' @param data_relative_path 
#' @param force 
#'
#' @return
#' @export
#'
#' @examples
make_integrate_report_qmd <- function(analysis, samples, report_steps = report_steps, plots_relative_path = "plots", report_path = ".", file_name = "integrated", 
                                      data_relative_path = "results", force = FALSE) {
  
  report_name <- paste0(file_name, ".qmd")
  report <- file.path(report_path, report_name)
  
  initialize_report_qmd(report_name, report_path, force = force)
  
  lines <- readLines(report)
  lines[2] <- paste0("title : \"", file_name, "\"")
  lines[3] <- paste0("subtitle: \"Integration report\"")
  lines[60] <- paste0("The goal of this report is to portrait the statistical insights of the single-cell RNA-seq analysis to better diagnose the SingleCell::integrate function that was use on the given Seurat objects. ")
  writeLines(lines, report)
  
  fig_num <- 1
  table_num <- 1
  cat("## Results\n\n", file = report, sep = "", append = TRUE)
  
  ## Filtering part
  
  if (report_steps[[1]]) {
    cat("### Filtering\n\n", file = report, sep = "", append = TRUE)
    cat("::: panel-tabset\n\n", file = report, sep = "", append = TRUE)
    
    for (sample in samples) {
      cat("## ", sample, "\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: Filtering statistics for ", sample, ".](", plots_relative_path, "/", sample, "_complete_filter_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_complete_filter_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
      cat("```{r, warning=FALSE}\n\ntable <- read.csv('", data_relative_path, "/", sample, "_filtering_stats.csv')\n\nDT::datatable(\n\ttable,\n\trownames = FALSE,\n\textensions = 'Buttons',\n\toptions = list(\n\t\tdom = 'Bfrtip',\n\t\tbuttons = list(\n\t\t\tlist(\n\t\t\t\textend = 'colvis'\n\t\t\t)\n\t\t)\n\t),\n\tcaption = htmltools::tags$caption(\n\t\tstyle = 'caption-side: top; text-align: center;', \n\t\thtmltools::tags$b('Table ", table_num, ":'), 'Cell filtrated and proportions after filtration for ", sample,".'\n\t)\n)\n\n```\n\n", file = report, sep = "", append = TRUE)
      fig_num <- fig_num + 1
      table_num <- table_num + 1
    }
    
    cat(":::\n\n", file = report, sep = "", append = TRUE)
  }
  
  if (report_steps[[2]]) {
    
    cat("### Elbow plots\n\n", file = report, sep = "", append = TRUE)
    cat("::: panel-tabset\n\n", file = report, sep = "", append = TRUE)
    
    for (sample in samples) {
      cat("## ", sample, "\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: Elbow plot for choosing dimensions for UMAPs and FindNeighbors in ", sample, ".](", plots_relative_path, "/", sample, "_Elbow_pca_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_Elbow_pca_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
      fig_num <- fig_num + 1
    }
    
    cat(":::\n\n", file = report, sep = "", append = TRUE)
    
    cat("### PCA plots\n\n", file = report, sep = "", append = TRUE)
    
    cat("::: panel-tabset\n\n", file = report, sep = "", append = TRUE)
    
    for (sample in samples) {
      cat("## ", sample, "\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: PCA plot of ", sample, ".](", plots_relative_path, "/", sample, "_PCA_pca_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_PCA_pca_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
      fig_num <- fig_num + 1
    }
    
    cat(":::\n\n", file = report, sep = "", append = TRUE)
  }
  
  ## Conclusion part
  cat("## Conclusion\n\n", file = report, sep = "", append = TRUE)
  cat("We have successfully integrate the different samples for ", file_name, ". We can now assign cell types to identified clusters based on a list of markers.", file = report, sep = "", append = TRUE)
  
  quarto::quarto_render(report)
}

#' Title
#'
#' @param analysis 
#' @param sample 
#' @param report_path 
#' @param file_name 
#' @param report_steps 
#' @param plots_relative_path 
#' @param data_relative_path 
#' @param force 
#'
#' @return
#' @export
#'
#' @examples
make_analysis_report_qmd <- function(analysis, sample, report_path, file_name, report_steps, 
                                     plots_relative_path = "plots", data_relative_path = "results", force = FALSE) {
  
  report_name <- paste0(file_name, ".qmd")
  report <- file.path(report_path, report_name)
  
  initialize_report_qmd(report_name, report_path, force = force)
  
  lines <- readLines(report)
  lines[2] <- paste0("title : \"", file_name, "\"")
  lines[3] <- paste0("subtitle: \"Analyze integrated report for clusterisation\"")
  lines[60] <- paste0("The goal of this report is to portrait the statistical insights of the single-cell RNA-seq analysis to better diagnose the SingleCell::analyze_integrated function that was use on the given Seurat object. ")
  writeLines(lines, report)
  
  fig_num <- 1
  table_num <- 1
  cat("## Results\n\n", file = report, sep = "", append = TRUE)
  
  ## Filtering part
  if (report_steps[[1]]) {
    cat("### Filtering\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: Filtering statistics for ", sample, ".](", plots_relative_path, "/", sample, "_complete_filter_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_complete_filter_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    cat("```{r, warning=FALSE}\n\ntable <- read.csv('", data_relative_path, "/", sample, "_filtering_stats.csv')\n\nDT::datatable(\n\ttable,\n\trownames = FALSE,\n\textensions = 'Buttons',\n\toptions = list(\n\t\tdom = 'Bfrtip',\n\t\tbuttons = list(\n\t\t\tlist(\n\t\t\t\textend = 'colvis'\n\t\t\t)\n\t\t)\n\t),\n\tcaption = htmltools::tags$caption(\n\t\tstyle = 'caption-side: top; text-align: center;', \n\t\thtmltools::tags$b('Table ", table_num, ":'), 'Cell filtrated and proportions after filtration for ", sample,".'\n\t)\n)\n\n```\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
    table_num <- table_num + 1
  }
  
  ## PCA part
  if (report_steps[[2]]) {
    cat("### PCA diagnosis\n\n", file = report, sep = "", append = TRUE)
    cat("::: panel-tabset\n\n", file = report, sep = "", append = TRUE)
    cat("## Elbow plot\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: Elbow plot for choosing dimensions for UMAPs and FindNeighbors in ", sample, ".](", plots_relative_path, "/", sample, "_Elbow_pca_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_Elbow_pca_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
    cat("## PCA plot\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: PCA plot of ", sample, ".](", plots_relative_path, "/", sample, "_PCA_pca_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_PCA_pca_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    cat(":::\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
  }
  
  ## Finding clusters part
  if (report_steps[[3]]) {
    cat("### Finding clusters\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: All clusters for ", sample, " based on different resolutions.](", plots_relative_path, "/", sample, "_clustree.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_clustree.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
  }
  
  ## UMAPs part
  if (report_steps[[4]]) {
    cat("### UMAPs\n\n", file = report, sep = "", append = TRUE)
    cat("::: panel-tabset\n\n", file = report, sep = "", append = TRUE)
    
    if (length(unique(analysis@meta.data$orig.ident)) == 1) {
      cat("## Sample identification\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: UMAP based on the optimized number of dimensions of showing the sample for ", sample, ".](", plots_relative_path, "/", sample, "_sample_umap_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_sample_umap_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    } else {
      cat("## Samples identification\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: UMAP based on the optimized number of dimensions of showing the samples for ", sample, ".](", plots_relative_path, "/", sample, "_sample_umap_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_sample_umap_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    }
    fig_num <- fig_num + 1
    
    cat("## Base clusters\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: UMAP categorizing the different clusters for ", sample, ". These base clusters are from `seurat_clusters` meta.data column.](", plots_relative_path, "/", sample, "_clusters_numbers_umap_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_clusters_numbers_umap_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
    
    cat("## Number of features\n\n", file = report, sep = "", append = TRUE)
    cat("![**Figure ", fig_num, "**: UMAP portaying the number of read counts per cells in ", sample, " scaled in log2. These results are based on the `nCount_RNA` column in the meta.data.](", plots_relative_path, "/", sample, "_nCount_umap_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_nCount_umap_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
    fig_num <- fig_num + 1
    
    if ("percent_mt" %in% colnames(analysis@meta.data)) {
      cat("## Percentage of mitochondria\n\n", file = report, sep = "", append = TRUE)
      cat("![**Figure ", fig_num, "**: UMAP portaying the mitochondria content per cells in ", sample, " in percentage. These results are based on the `percent_mt` column in the meta.data.](", plots_relative_path, "/", sample, "_mitochondria_umap_plot.png)\n\n:::{.callout-note collaspe='true'}\nYou can find the figure here : `", plots_relative_path, "/", sample, "_mitochondria_umap_plot.png`\n:::\n\n", file = report, sep = "", append = TRUE)
      fig_num <- fig_num + 1
    }
    
    cat(":::\n\n", file = report, sep = "", append = TRUE)
  }
  
  ## Conclusion part
  cat("## Conclusion\n\n", file = report, sep = "", append = TRUE)
  if (length(unique(analysis@meta.data$orig.ident)) == 1) {
    cat("We have successfully analyze the sample for ", sample, ". We can now assign cell types to identified clusters based on a list of markers.", file = report, sep = "", append = TRUE)
  } else {
    cat("We have successfully analyze the different samples for ", sample, ". We can now assign cell types to identified clusters based on a list of markers.", file = report, sep = "", append = TRUE)
  }
  
  quarto::quarto_render(report)
}

