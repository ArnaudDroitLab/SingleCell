theme_white_bg <- function(){
  ggplot2::theme(panel.background = ggplot2::element_blank(),
                 axis.line = ggplot2::element_line())
}

colorblind_palette_7 <- function() {
  c("#F0E442", "#56B4E9", "#E69F00", "#009E73", "#D55E00", "#0072B2", "#CC79A7")
}

get_color_x <- function(x, color_list = colorblind_palette_7()) {
  color_list[0:(x-1)%%7+1]
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
    ggplot2::theme(axis.title.x=ggplot2::element_blank()) + 
    ggplot2::scale_fill_manual(values = "#237BFF")
  
  low = max(low, min(0, min(df[[y_name]])))
  high = min(high, max(df[[y_name]]) + 1)
  
  p <- p + ggplot2::geom_hline(yintercept=low, linetype="dashed", color = "#002f76")
  p <- p + ggplot2::geom_hline(yintercept=high, linetype="dashed", color = "#002f76")
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
    ggplot2::theme_bw() + ggplot2::geom_vline(aes(xintercept = k.param.neighbors), colour = "#002f76")
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
      ggplot2::scale_colour_gradient(high = "#009784", low = "#d6fffa") +
      ggplot2::labs(color = paste0("log2(RNA count)")) +
      ggplot2::theme(legend.position = "top")
  } else if (colour_by == "percent_mt") {
    df <- df %>% dplyr::arrange(!!sym(colour_by))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = .data[[y]])) +
      ggplot2::geom_point(aes(colour = .data[[colour_by]])) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_gradient(high = "#990099", low = "#ffe9ff") +
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
    
    if (length(unique(df$colour_by)) == 1) {
      p <- p + ggplot2::scale_colour_manual(values = "#6B5FAC")
    } else {
      color_palette_UMAP <- c("#FA73E2", "#6B5FAC", "#000023", "#332DA3", "#237BFF", "#92CDF9", "#22aa99")
      color_values <- generate_gradient_palette(seurat, color_palette_UMAP, n_colors_out = 100, n_clusters = colour_by)
      p <- p + ggplot2::scale_colour_manual(values = color_values)
    }
    
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
#' @import ggraph 
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


#' Make a violin plot, using the features (genes or elements of metadata) from the Seurat object.
#'
#' @param seurat The Seurat object from which to
#' @param features
#' @param group.by
#' @param assay
#' @param slot
#' @param show_points
#'
#' @return
#' @export
#'
#' @examples
plot_seurat_violin <- function(seurat, features, group.by = "orig.ident", assay = "RNA", slot = "data", show_points = TRUE, taxonomy = "ensembl", threshold = NULL, alpha = 0.2, box_fill = "lightgrey") {
  
  df_large <- Seurat::FetchData(seurat, vars = c(group.by, features), assay = "RNA", slot = "data")
  # Use melt to change data.frame format
  
  if (taxonomy == "ensembl") {
    gene_list_translated <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = colnames(df_large)[2:length(df_large)],
      column = "SYMBOL", 
      keytype = "ENSEMBL"
    )
    gene_list_names <- names(gene_list_translated)
    gene_list_translated <- unname(gene_list_translated)
  } else {
    gene_list_names <- gene_list
    gene_list_translated <- gene_list
  }
  
  colnames(df_large) <- c(group.by, gene_list_translated)
  
  df_long <- reshape2::melt(df_large, id.vars = group.by, measure.vars = gene_list_translated,
                         variable.name = "feature", value.name = "expression")
  
  if (length(threshold) > 0) {
    df_long <- df_long %>%
      dplyr::filter(expression > threshold)
  }
  
  colors <- get_color_x(length(unique(df_large[[group.by]])), color_list = colorblind_palette_7())
  
  p <- ggplot(df_long, aes(x = factor(.data[[group.by]]), y = .data[["expression"]], fill = .data[[group.by]])) +
    ggplot2::geom_violin() +
    ggplot2::geom_jitter(color = "grey3", size = 0.2, alpha = alpha) +
    ggplot2::scale_y_continuous("Expression level", position="right") +
    facet_grid(rows = vars(feature), scales = "free", switch = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.spacing = unit(0, "lines"),
                   legend.position = "none",
                   strip.text.y.left = element_text(angle = 0),
                   strip.text = element_text(face = "bold"), 
                   strip.background = element_rect(fill = box_fill, color = "black")) +
    ggplot2::scale_fill_manual(values = colors) +
    xlab("Identity")
  return(p)

}

#' Put labels on UMAPs using centroids and ggrepel.
#'
#' @param seurat 
#' @param assays 
#' @param slot 
#' @param colour_by 
#' @param label.size 
#' @param pt.size 
#' @param width 
#' @param height 
#' @param path 
#' @param project_name 
#' @param color_palette_ordered 
#'
#' @return
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @export
#'
#' @examples
plot_label_umap <- function(seurat, assay = "RNA", slot = "data", colour_by = "orig.ident", label.size = 4, pt.size = 0.5, color_palette_ordered = NULL) {
  
  check_assay(seurat, assay)
  check_reduction(seurat, "umap")
  checkmate::assert_string(colour_by)
  checkmate::assert_string(slot)
  
  df <- data.frame(umap1 = seurat@reductions[["umap"]]@cell.embeddings[,1],
                   umap2 = seurat@reductions[["umap"]]@cell.embeddings[,2],
                   color = seurat@meta.data[[colour_by]])
  
  centers <- aggregate(cbind(UMAP_1 = df[,1], UMAP_2 = df[,2]),
                       by = list(group = df$color),
                       FUN = mean)
  
  if (!is.null(color_palette_ordered)) {
    color_values <- generate_gradient_palette(seurat, color_palette_ordered, n_colors_out = 100, n_clusters = colour_by)
    centers$fill_color <- color_values
    centers$text_color <- get_text_contrast_color(centers$fill_color)
  } else {
    centers$text_color <- "black"
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = umap1, 
                               y = umap2, 
                               color = color)) +
    ggplot2::geom_point(size = pt.size, alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::labs(color = colour_by) +
    ggplot2::theme(legend.position = "none") + 
    ggrepel::geom_label_repel(data = centers, 
                              ggplot2::aes(x = UMAP_1, y = UMAP_2, fill = group, label = group),
                              size = label.size, 
                              color = centers$text_color)
  
  if (!is.null(color_palette_ordered)) {
    p <- p + ggplot2::scale_color_manual(values = color_values) + 
      ggplot2::scale_fill_manual(values = color_values)
  }

  return(p)
}


get_text_contrast_color <- function(hex_color) {
  rgb <- grDevices::col2rgb(hex_color) 
  yiq <- (299 * rgb[1, ] + 587 * rgb[2, ] + 114 * rgb[3, ]) / 1000
  ifelse(yiq >= 128, "black", "white")
}

generate_gradient_palette <- function(seurat, color_palette_ordered, n_colors_out = 100, n_clusters = "orig.ident") {
  
  gradient_func <- colorRampPalette(color_palette_ordered)
  full_gradient <- gradient_func(n_colors_out)
  
  n_extract <- length(unique(seurat@meta.data[[n_clusters]]))
  
  idx <- round(seq(1, n_colors_out, length.out = n_extract))
  selected_colors <- full_gradient[idx]
  
  return(selected_colors)
}











