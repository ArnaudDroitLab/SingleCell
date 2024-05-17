theme_white_bg <- function(){
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line())
}

#' Plot qc as violin plot, adding cutoffs lines.
#'
#' @param df A data frame containing x_name and y_name.
#' @param x_name The column of df to use as groups in the violin plot, one violin per group.
#' @param y_name The column of df to use as y axis.
#' @param low The lower cutoff to put on the plot. Defaults to 0, which will not put a lower cutoff.
#' @param high The higher cutoff to put on the plot. Defaults to Inf, which will not put a higher cutoff.
#'
#' @return A ggplot object.
#' @export
plot_filter <- function(df, x_name = "x", y_name = "y", low = 0, high = Inf) {
  checkmate::assert_data_frame(df)
  checkmate::assert_string(x_name)
  if (!x_name %in% names(df)) {stop(paste0(x_name, " not a column in df."))}
  checkmate::assert_string(y_name)
  if (!y_name %in% names(df)) {stop(paste0(y_name, " not a column in df."))}
  checkmate::assert_numeric(low)
  checkmate::assert_numeric(high)
  if (low > high) {stop("Low cannot be superior to high.")}

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], fill = .data[[x_name]])) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(size = 0.2, color = ggplot2::alpha("black", 0.4), fill = ggplot2::alpha("black", 0.4)) +
    theme_white_bg() +
    ggplot2::theme(axis.title.x=ggplot2::element_blank())

  low = max(low, min(0, min(df[[y_name]])))
  high = min(high, max(df[[y_name]]) + 1)

  p <- p + ggplot2::geom_hline(yintercept=low, linetype="dashed", color = "red")
  p <- p + ggplot2::geom_hline(yintercept=high, linetype="dashed", color = "red")
  return(p)
}

#' Make the Elbow plot of a dimension reduction
#'
#' The elbow plot is the standard deviation in y axis and the component number in x axis.
#' Most of the time it has the shape of elbow, hence why it is called elbow plot.
#' The goal of this plot is to help choose the correct number of components for further analysis.
#'
#' @param seurat A Seurat object
#' @param reduction Which reduction to us. Must be present in the Seurat object. Default "pca"
#' @param npc Number of components to plot. Default 50
#'
#' @return The elbow plot as a ggplot object
#' @export
plot_seurat_elbow <- function(seurat, reduction = "pca", npc = 50) {
  checkmate::assert_class(seurat, "Seurat")
  check_reduction(seurat, reduction)
  checkmate::assert_int(npc)
  if (length(seurat@reductions[[reduction]]@stdev)<npc) {stop(paste0(reduction, " does not have ", npc, " components."))}

  df <- data.frame("Standard_Deviation" = seurat@reductions[[reduction]]@stdev[1:npc], PC = 1:npc)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["PC"]], y = .data[["Standard_Deviation"]])) +
    ggplot2::geom_point() +
    theme_white_bg()
  return(p)
}

#' Plot the cells on a reduction, with colours from the object.
#'
#' X axis is the first component, Y axis second component. Colour_by is fetched
#' from the seurat object, and can be a gene or a column from the metadata.
#'
#' @param seurat The Seurat object
#' @param reduction Which reduction to use. Default "pca"
#' @param colour_by Which data to use to colour cells. Default "orig.ident"
#' @param assay Which assay to use (if applicable). Default "RNA"
#' @param slot Which slot to use within the assay (if applicable). Default "data"
#'
#' @return A ggplot object
#' @export
plot_seurat_dim <- function(seurat, reduction = "pca", colour_by = "orig.ident", assay = "RNA", slot = "data") {
  check_assay(seurat, assay)
  check_reduction(seurat, reduction)
  checkmate::assert_string(colour_by)
  checkmate::assert_string(slot)

  x = paste0(reduction,"1")
  y = paste0(reduction,"2")
  df <- data.frame(x = seurat@reductions[[reduction]]@cell.embeddings[,1],
                   y = seurat@reductions[[reduction]]@cell.embeddings[,2],
                   colour_by = Seurat::FetchData(seurat, colour_by, assay = assay, slot = slot))
  colnames(df) <- c(x, y, colour_by)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]],
                                        fill = .data[[colour_by]], colour = .data[[colour_by]])) +
    ggplot2::geom_point() +
    theme_white_bg()
  return(p)
}
