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

  cat("```{r, fig.show='hold', out.width='", 100%/%ncol, "%'}\n",
      "knitr::include_graphics(c('", paste0(images_relative_path, collapse = "','"), "'))\n",
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
        if (length(filepaths) < 0) {
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


#' Make the integration report, with the files to include hardcoded inside.
#'
#' @param sample List of samples to include in the report.
#' @param report_path Path to where the report should be located.
#' @param report_name Name to give to the report file.
#' @param plots_relative_path Path to the plots directory, relative to the `report_path`.
#' Must be a subdirectory in `report_path`. Default "plots"
#' @param data_relative_path Path to the csv dataframes directory, relative to the `report_path`.
#' Must be a subdirectory in `report_path`. Default "results"
#' @param title The title of the report. Default to "Analysis report"
#' @param force Set to true if you want to overwrite the report that could already exists. Default FALSE
#'
#' @return Path to the report as a string.

#' @export
#'
#' @examples
make_analysis_report <- function(sample, report_path, report_name, plots_relative_path = "plots", data_relative_path = "results", title = "Analysis report", force = FALSE) {

  # Every file will be sample_name.ext, you should report here name.ext, sample_ will be automatically added. Currently only accepts .csv and .png
  steps_files <- list(Filtering = list(plots = c("count_filter_plot.png", "feature_filter_plot.png", "mitochondria_filter_plot.png"), 
                                       df = c("filtering_stats.csv")),
                      PCA = list(plots = c("Elbow_pca_plot.png", "PCA_pca_plot.png")),
                      "Clustering tree" = list(plots = "clustree.png"),
                      UMAP = list(plots = c("sample_umap_plot.png", "clusters_numbers_umap_plot.png", 
                                            "mitochondria_umap_plot.png", "nCount_umap_plot.png")),
                      DE = list(df = c("top10_DE.csv", "summary_per_clusters.csv")))
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
  
  checklist_steps <- list(Filtering = "Filtering", 
                          Elbow = "Elbow", 
                          Tree = "Clustering Tree", 
                          Table = "Cluster Table", 
                          UMAP = "UMAP")

  for (step in names(steps_files)) {
    add_title_report(report, step, 1)
    textfile <- file.path(".report", paste0("integration_", step, ".txt"))
    if (file.exists(file.path(report_path, textfile))) {
      add_textfile_report(report, textfile)
    }
    for (type in names(steps_files[[step]])) {
      filepaths <- steps_files[[step]][[type]]
      if (length(filepaths) == 0) {
        next
      }

      if (type == "plots") {
        filepaths <- lapply(filepaths, function(x) file.path(plots_relative_path, paste0(sample, "_", x)))
        filepaths <- filepaths[file.exists(file.path(report_path, filepaths))]
        if (step %in% c("Clustering tree", "UMAP", "Filtering_stats")) {
          add_images_knit_report(report, filepaths, ncol = 1)
        } 
        else {
          add_images_knit_report(report, filepaths, ncol = if (length(filepaths)%%3 == 0) 3 else 2)
        }

      } else if (type == "df") {
        filepaths <- lapply(filepaths, function(x) file.path(data_relative_path, paste0(sample, "_", x)))
        filepaths <- filepaths[file.exists(file.path(report_path, filepaths))]
        if (step == "DE") {
          for (df in filepaths) {add_df_report_no_rownames(report, df)}
        } else {
          for (df in filepaths) {add_df_report(report, df)}
        }

      }
    }
  }
  build_report(report, title)
  return(report)

}

