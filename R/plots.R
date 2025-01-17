theme_white_bg <- function(){
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line())
}

colorblind_palette_7 <- function() {
  c("#F0E442", "#56B4E9", "#E69F00", "#009E73", "#D55E00", "#0072B2", "#CC79A7")
}

get_colorblind_x <- function(x) {
  colorblind_palette_7()[0:(x-1)%%7+1]
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
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_jitter
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 alpha
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 geom_hline
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
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x=ggplot2::element_blank())
  
  low = max(low, min(0, min(df[[y_name]])))
  high = min(high, max(df[[y_name]]) + 1)
  
  p <- p + ggplot2::geom_hline(yintercept=low, linetype="dashed", color = "red")
  p <- p + ggplot2::geom_hline(yintercept=high, linetype="dashed", color = "red")
  return(p)
}

#' Barplots that represent the filtering statistics for the seurat filtration
#'
#' @param df Data frame as the output by seurat_filter that has the filtering statistics
#' @param x_name The first column that needs to be work with for the first barplot. Default "Genes"
#' @param y_name The second column that needs to be work with for the second barplot. Default "Cells"
#'
#' @return Barplots as a ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 theme_bw
#' @export
#'
#' @examples

plot_filtering_stats <- function(df, x_name = "Genes", y_name = "Cells") {
  
  checkmate::assert_data_frame(df)
  
  
  if (!x_name %in% colnames(df)) {
    stop(paste0("First column is missing from the data frame."))
  }
  
  if (!y_name %in% colnames(df)) {
    stop(paste0("Second column is missing from the data frame."))
  }
  
  
  if (df[[x_name]][1] - df[[x_name]][2] != df[[x_name]][3]) {stop("Filtration cannot add features.")}
  if (df[[y_name]][1] - df[[y_name]][2] != df[[y_name]][3]) {stop("Filtration cannot add features.")}
  
  if (isFALSE(identical(rownames(df), c("Before", "After", "Filtered_out", "Percentage")))) {
    stop(paste0("Dataframe should be ordered like this = 'Before', 'After', 'Filtered_out'."))
  }
  
  p <- ggplot2::ggplot(data = df, aes(x = row.names(df), .data[[x_name]])) +
    ggplot2::geom_bar(stat = "identity", color = "#D44C7E", fill = "#F39BBC", width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data[[x_name]]), vjust = 1.5, size = 3) +
    ggplot2::scale_x_discrete(name = "", limits = c("Before", "After", "Filtered_out"),
                              labels = c("Before", "After", "Filtered")) +
    ggplot2::theme_bw()
  q <- ggplot2::ggplot(data = df, aes(x = row.names(df), .data[[y_name]])) +
    ggplot2::geom_bar(stat = "identity", color = "#FFC107", fill = "#FFE493", width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data[[y_name]]), vjust = 1.5, size = 3) +
    ggplot2::scale_x_discrete(name = "", limits = c("Before", "After", "Filtered_out"),
                              labels = c("Before", "After", "Filtered")) +
    ggplot2::theme_bw()
  return(p + q)
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
#' @param k.param.neighbors Number of k.param for neighbors joining. Default 20
#'
#' @return The elbow plot as a ggplot object
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_y_continuous
#' @export
plot_seurat_elbow <- function(seurat, reduction = "pca", npc = 50, k.param.neighbors = 20) {
  checkmate::assert_class(seurat, "Seurat")
  check_reduction(seurat, reduction)
  checkmate::assert_int(npc)
  checkmate::assert_int(k.param.neighbors)
  if (length(seurat@reductions[[reduction]]@stdev)<npc) {stop(paste0(reduction, " does not have ", npc, " components."))}
  
  df <- data.frame("Standard_Deviation" = seurat@reductions[[reduction]]@stdev[1:npc], PC = 1:npc)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["PC"]], y = .data[["Standard_Deviation"]])) +
    ggplot2::geom_point() + ggplot2::scale_y_continuous("Standard Deviation") +
    ggplot2::theme_bw() + ggplot2::geom_vline(aes(xintercept = k.param.neighbors), colour = "steelblue3")
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
#' @return A ggplot graph
#' @importFrom Seurat FetchData
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_colour_gradient
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 element_rect
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
  
  if (colour_by == "nCount_RNA") {
    df <- df %>% mutate(!!sym(colour_by) := log2(!!sym(colour_by))) %>% dplyr::arrange(!!sym(colour_by))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_point(aes(colour = .data[[colour_by]])) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_gradient(high = "#429DEF", low = "#ECECEC") +
      ggplot2::labs(color = paste0("log2(RNA count)")) +
      ggplot2::theme(legend.position = "top")
  } else if (colour_by == "percent_mt") {
    df <- df %>% dplyr::arrange(!!sym(colour_by))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_point(aes(colour = .data[[colour_by]])) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_gradient(high = "#F39243", low = "#ECECEC") +
      ggplot2::labs(color = "Mitochondrial percentage") +
      ggplot2::theme(legend.position = "top")
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]],
                                          colour = .data[[colour_by]])) +
      ggplot2::geom_point() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "top", legend.key = element_rect(fill = "white", colour = "black")) +
      ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) +
      ggplot2::theme(legend.position = "top")
  }
  return(p)
}

#' Plot a clustree
#'
#' Clustree is a graph that shows the moving of cells between the clusters in each resolution. The clustering numbering between each resolution
#'
#' @param seurat The Seurat object
#' @param prefix The prefix to use to gather all the resolutions in the graph.
#' This should be the whole prefix before the resolution number. Default "RNA_snn_res."
#'
#' @return a ggplot graph
#' @importFrom clustree clustree
#' @export
plot_seurat_clustree <- function(seurat, prefix = "RNA_snn_res.") {
  checkmate::assert_class(seurat, "Seurat")
  checkmate::assert_character(prefix, max.len = 1)
  greping <- paste0("^", prefix)
  seurat_prefix_column <- grep(greping, colnames(seurat@meta.data), value = TRUE)
  for (col in seurat_prefix_column) {
    checkmate::assert_factor(seurat@meta.data[[col]])
  }
  p <- clustree::clustree(seurat@meta.data, prefix=prefix)
  return(p)
}


plot_seurat_violin <- function(seurat, features, group.by = "orig.ident", assay = "RNA", slot = "data", show_points = TRUE) {
  
  df_large <- Seurat::FetchData(seurat, vars = c(group.by, features), assay = "RNA", slot = "data")
  # Use melt to change data.frame format
  df_long <- reshape2::melt(df_large, id.vars = group.by, measure.vars = features,
                         variable.name = "feature", value.name = "expression")
  colors <- get_color_x(length(unique(df_large[[group.by]])), color_list = colorblind_palette_7())
  
  p <- ggplot(df_long, aes(x = factor(.data[[group.by]]), y = .data[["expression"]], fill = .data[[group.by]])) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(color = "grey3", size = 0.2, alpha = 0.2) +
    ggplot2::scale_y_continuous("Expression level", position="right") +
    facet_grid(rows = vars(feature), scales = "free", switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.spacing = unit(0, "lines"),
                   legend.position = "none",
                   strip.text.y.left = element_text(angle = 0),
                   strip.text = element_text(face = "bold")) +
    ggplot2::scale_fill_manual(values = colors) +
    xlab("Identity")
  return(p)

}














