#' Initialize a Rmd report from the predefined template.
#'
#' @param name Name of the file
#' @param report_path Path to where the report should be located. Must be the upper directory of the plots directory (if /path/to/plot is the plot directory, then report should be in /path/to). Default to "."
#' @param force Set to true if you want to overwrite the the report that could already exists`report_path`. Default FALSE
#'
#' @return Nothing
#' @export
#'
#' @examples
initialize_report <- function(name, report_path = ".", force = FALSE) {
  checkmate::assert_string(name, pattern = "^[a-zA-Z0-9\\-_.]+$")
  checkmate::assert_directory_exists(report_path)
  checkmate::assert_logical(force, len = 1)

  template_path <- system.file("rmd/template.Rmd", package = "SingleCell")
  internal_files <- system.file("rmd/.report", package = "SingleCell")
  internal_path <- file.path(report_path, ".report")
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

#' Add a markdown title to a report
#'
#' @param report Path to the report.
#' @param title Content of the title. Default to ""
#' @param level A number between 1 and 5, 1 being a title and going up means a lower level of subtitle. Default to 1
#'
#' @return Nothing
#' @export
#'
#' @examples
add_title_report <- function(report, title = "", level = 1) {
  checkmate::assert_file_exists(report)
  checkmate::assert_string(title)
  checkmate::assert_number(level, lower = 1, upper = 5)

  cat(rep("#", times = level), " ", title, "\n\n", file = report, sep = "", append = TRUE)
}

#' Add some text to a report
#'
#' @param report Path to the report.
#' @param txt Text to add. Default to ""
#'
#' @return Nothing
#' @export
#'
#' @examples
add_text_report <- function(report, txt = "") {
  checkmate::assert_file_exists(report)
  checkmate::assert_string(txt)

  cat(txt, "\n\n", file = report, sep = "", append = TRUE)
}

#' Add an image to a report
#'
#' @param report Path to the report.
#' @param image_relative_path Path to the image, relative to the report.
#' @param title Title to give the image
#'
#' @return Nothing
#' @importFrom tools file_path_sans_ext
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_to_title
#' @export
#'
#' @examples
add_image_markdown_report <- function(report, image_relative_path, title = "") {
  checkmate::assert_file_exists(report)
  checkmate::assert_file(file.path(dirname(report), image_relative_path))
  checkmate::assert_string(title)

  if (title == "") {
    title <- basename(image_relative_path) %>%
      tools::file_path_sans_ext() %>%
      stringr::str_replace_all(pattern = "_", replacement = " ") %>%
      stringr::str_to_title()
  }

  cat("![", title, "](", image_relative_path, " '", title, "')\n\n", file = report, sep = "", append = TRUE)
}

add_images_knit_report <- function(report, images_relative_path, ncol = 2) {
  checkmate::assert_file_exists(report)
  checkmate::assert_number(ncol, lower = 1)

  for (image_relative_path in images_relative_path) {
    checkmate::assert_file(file.path(dirname(report), images_relative_path))
  }
  
  if (length(images_relative_path) == 0) {
    cat("")
  } else {
    cat("```{r, fig.show='hold', out.width='", 100%/%ncol, "%'}\n",
        "knitr::include_graphics(c('", paste0(images_relative_path, collapse = "','"), "'))\n",
        "```\n\n", file = report, sep = "", append = TRUE)
  }

}

#' Add a dataframe to a report.
#'
#' Will read a csv file and display its content on the report.
#'
#' @param report Path to the report.
#' @param csv_relative_path Path to the csv, relative to the report.
#'
#' @return Nothing
#' @export
add_df_report <- function(report, csv_relative_path) {
  checkmate::assert_file_exists(report)
  checkmate::assert_file(file.path(dirname(report), csv_relative_path))

  cat("```{r}\n",
      "df <- read.csv('", csv_relative_path, "', header = TRUE, sep = ',', row.names = 1)\n",
      "df\n",
      "```\n\n", file = report, sep = "", append = TRUE)
}


#' Add a dataframe to a report.
#'
#' Will read a csv file and display its content on the report.
#'
#' @param report Path to the report.
#' @param csv_relative_path Path to the csv, relative to the report.
#'
#' @return Nothing
#' @export
add_df_report_no_rownames <- function(report, csv_relative_path) {
  checkmate::assert_file_exists(report)
  checkmate::assert_file(file.path(dirname(report), csv_relative_path))

  cat("```{r}\n",
      "df <- read.csv('", csv_relative_path, "', header = TRUE, sep = ',')\n",
      "df\n",
      "```\n\n", file = report, sep = "", append = TRUE)
}

add_check_box <- function(report, checklist) {
  checkmate::assert_file_exists(report)
  
  cat("- [ ] &nbsp;&nbsp;", checklist)
}


#' Add a chunk of code to read content of file in a report.
#'
#' @param report Path to the report.
#' @param relative_textfile Path to the text file, relative to the report.
#'
#' @return Nothing
#' @export
add_textfile_report <- function(report, relative_textfile) {
  checkmate::assert_file_exists(report)
  checkmate::assert_file(file.path(dirname(report), relative_textfile))

  cat("```{r comment=''}\n",
      "cat(readLines('", relative_textfile, "'), sep = '\n')\n",
      "```\n\n", file = report, sep = "", append = TRUE)
}


#' Render the markdown document.
#'
#' @param report Path to the report.
#' @param title Title to give to the report
#'
#' @return Nothing
#' @importFrom rmarkdown render
#' @export
#'
#' @examples
build_report <- function(report, title) {
  checkmate::assert_file_exists(report)
  checkmate::assert_string(title)

  rmarkdown::render(
    input = report,
    params = list(title = title),
    output_format = "all"
  )
}

#' Make the integration report, with the files to include hardcoded inside.
#'
#' @param samples List of samples to include in the report.
#' @param report_path Path to where the report should be located.
#' @param report_name Name to give to the report file.
#' @param plots_relative_path Path to the plots directory, relative to the `report_path`.
#' Must be a subdirectory in `report_path`. Default "plots"
#' @param data_relative_path Path to the csv dataframes directory, relative to the `report_path`.
#' Must be a subdirectory in `report_path`. Default "results"
#' @param title The title of the report. Default to "Integration report"
#' @param force Set to true if you want to overwrite the report that could already exists. Default FALSE
#'
#' @return Path to the report as a string.
#' @export
#'
#' @examples
make_integration_report <- function(samples, report_path, report_name = "integration.Rmd", plots_relative_path = "plots", data_relative_path = "results", title = "Integration report", force = FALSE) {

  # Every file will be sample_name.ext, you should report here name.ext, sample_ will be automatically added. Currently only accepts .csv and .png
  steps_files <- list(Filtering = list(plots = c("count_filter_plot.png", "feature_filter_plot.png", "mitochondria_filter_plot.png"), df = c("filtering_stats.csv")),
                      PCA = list(plots = c("Elbow_pca_plot.png", "PCA_pca_plot.png"))
                      )
  report <- file.path(report_path, report_name)

  if (file.exists(report)) {
    if (force) {
      warning("Report ", report, " already exists, overwriting it.")
    } else {
      warning("Report ", report, " already exists, skipping report creation.")
      return(FALSE)
    }
  }

  initialize_report(report_name, report_path, force = force)

  for (step in names(steps_files)) {
    add_title_report(report, toString(step), 1)
    textfile <- file.path(".report", paste0("integration_", step, ".txt"))
    if (file.exists(file.path(report_path, textfile))) {
      add_textfile_report(report, textfile)
    }
    for (sample in samples) {
      add_title_report(report, sample, 2)
      # add_title_report(report, sample, 2)
      for (type in names(steps_files[[step]])) {
        filepaths <- steps_files[[step]][[type]]
        if (length(filepaths) == 0) {
          next
        }

        if (type == "plots") {
          filepaths <- lapply(filepaths, function(x) file.path(plots_relative_path, paste0(sample, "_", x)))
          filepaths <- filepaths[file.exists(file.path(report_path, filepaths))]
          # print(length(filepaths))
          add_images_knit_report(report, filepaths, ncol = if (length(filepaths)%%3 == 0) 3 else 2)

        } else if (type == "df") {
          filepaths <- lapply(filepaths, function(x) file.path(data_relative_path, paste0(sample, "_", x)))
          filepaths <- filepaths[file.exists(file.path(report_path, filepaths))]
          for (df in filepaths) add_df_report(report, df)
        }
      }
    }
  }
  build_report(report, title)
  return(report)


}

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
      cat("```{r, warning=FALSE}\n\ntable <- read.csv('", data_relative_path, "/", sample, "_filtering_stats.csv')\n\nDT::datatable(\n\ttable,\n\trownames = FALSE,\n\textensions = 'Buttons',\n\toptions = list(\n\t\tdom = 'Bfrtip',\n\t\tbuttons = list(\n\t\t\tlist(\n\t\t\t\textend = 'colvis'\n\t\t\t)\n\t\t)\n\t),\n\tcaption = htmltools::tags$caption(\n\t\tstyle = 'caption-side: top; text-align: center;', \n\t\thtmltools::tags$b('Table ", table_num, ":'), 'Cell counts and proportions.'\n\t)\n)\n\n```\n\n", file = report, sep = "", append = TRUE)
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
    cat("```{r, warning=FALSE}\n\ntable <- read.csv('", data_relative_path, "/", sample, "_filtering_stats.csv')\n\nDT::datatable(\n\ttable,\n\trownames = FALSE,\n\textensions = 'Buttons',\n\toptions = list(\n\t\tdom = 'Bfrtip',\n\t\tbuttons = list(\n\t\t\tlist(\n\t\t\t\textend = 'colvis'\n\t\t\t)\n\t\t)\n\t),\n\tcaption = htmltools::tags$caption(\n\t\tstyle = 'caption-side: top; text-align: center;', \n\t\thtmltools::tags$b('Table ", table_num, ":'), 'Cell counts and proportions.'\n\t)\n)\n\n```\n\n", file = report, sep = "", append = TRUE)
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

