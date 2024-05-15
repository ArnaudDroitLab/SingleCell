#' Plot qc as violin plot, adding cutoffs lines if needed.
#'
#' @param df A data frame containing x_name and y_name.
#' @param x_name The column of df to use as groups in the violin plot, one violin per group.
#' @param y_name The column of df to use as y axis.
#' @param low The lower cutoff to put on the plot. Defaults to 0, which will not put a lower cutoff.
#' @param high The higher cutoff to put on the plot. Defaults to Inf, which will not put a higher cutoff.
#' @param force_cutoff Whether to force cutoffs to be added, even if they should not (mainly used to add a lower cutoff at 0). Defaults to FALSE.
#'
#' @return A ggplot object.
#' @export
plot_filter <- function(df, x_name = "x", y_name = "y", low = 0, high = Inf, force_cutoff = FALSE) {
  checkmate::assert_data_frame(df)
  checkmate::assert_numeric(low)
  checkmate::assert_numeric(high)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]], fill = .data[[x_name]])) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(size = 0.2, color = ggplot2::alpha("black", 0.4), fill = ggplot2::alpha("black", 0.4)) +
    theme_white_bg() +
    ggplot2::theme(axis.title.x=ggplot2::element_blank())
  if (low != 0 | force_cutoff) {
    p <- p + ggplot2::geom_hline(yintercept=low, linetype="dashed", color = "red")
  }
  if (high < Inf | force_cutoff) {
    p <- p + ggplot2::geom_hline(yintercept=high, linetype="dashed", color = "red")
  }
  return(p)
}


theme_white_bg <- function(){
  ggplot2::theme(panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line())
}

# qc_vlnplot <- function(seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 0, x_col = "orig.ident", show_x = F, show_points = T) {
#   ncol <- ifelse(ncol == 0, length(features), ncol)
#   data <- seurat@meta.data[,c(features, x_col)]
#   plots <- list()
#   for (e in features){
#     p <- ggplot(data, aes(x = as.character(.data[[x_col]]), y = .data[[e]], fill = .data[[x_col]])) + geom_violin() + theme(legend.position = "none") + theme_white_bg()
#     if (show_points) {p <- p + geom_jitter(size = 0.2, color = alpha("black", 0.4), fill = alpha("black", 0.4))}
#     if (!show_x) {p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x = NULL)}
#     plots[[e]] <- p
#   }
#
#   grid.arrange(grobs = plots, ncol = ncol)
#   # return(plots)
# }
