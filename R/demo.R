#' Get demo seurat object 
#'
#' @return Seurat object
#'
#' @examples
#' pbmc <- get_demo_seurat_object()
#'
#' @export
get_demo_seurat_object <- function() {
  readRDS(system.file("extdata/pbmc.rds"), package = "SingleCell")
}