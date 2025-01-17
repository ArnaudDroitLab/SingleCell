library(testthat)
library(Seurat)
library(clustree)

# Set a working seurat object

seurat_obj <- get_demo_seurat_object_clust()

test_that("plot_seurat_clustree handles invalid seurat input", {
  
  # Valid input (Seurat object and correct prefix)
  expect_s3_class(plot_seurat_clustree(seurat_obj), "gg")
  
  # Invalid input: Not a Seurat object
  expect_error(plot_seurat_clustree(list()), "Assertion on 'seurat' failed: Must inherit from class 'Seurat', but has class 'list'.")
})

test_that("plot_seurat_clustree handles 'prefix' argument correctly", {

  # Test for valid prefix
  expect_s3_class(plot_seurat_clustree(seurat_obj, prefix = "RNA_snn_res."), "gg")
  
  # Test for invalid prefix (empty string or multiple strings)
  expect_error(plot_seurat_clustree(seurat_obj, prefix = ""), 
               "Must be of type 'factor', not 'double'.")
  expect_error(plot_seurat_clustree(seurat_obj, prefix = c("RNA_snn_res.", "extra_prefix")), 
               "Assertion on 'prefix' failed: Must have length <= 1, but has length 2.")
})

test_that("plot_seurat_clustree should not handle Seurat objects without clustering information", {
  # Create Seurat object without clustering
  seurat_test <- get_demo_seurat_object_pca()
  
  # Check if error is raised when no clustering is available
  expect_error(plot_seurat_clustree(seurat_test), "Less than two column names matched the prefix: RNA_snn_res.")
})

test_that("plot_seurat_clustree returns a ggplot object", {

  p <- plot_seurat_clustree(seurat_obj)
  
  # Check if the returned object is a ggplot object
  expect_s3_class(p, "gg")
})

test_that("plot_seurat_clustree handles missing clusters in meta.data", {

  # Create a seurat object with false clustering information
  seurat_test <- seurat_obj
  seurat_test[["RNA_snn_res.0.1"]] <- rep("strawberry", length(seurat_test[["nCount_RNA"]]))
  
  # Check if error is raised when no cluster information is available in 'meta.data'
  expect_error(plot_seurat_clustree(seurat_test), "Must be of type 'factor', not 'character'.")
})

test_that("plot_seurat_clustree correctly plots available clusters", {
  
  # Check if clustree function runs without errors
  p <- plot_seurat_clustree(seurat_obj)
  plot <- ggplot2::ggplot_build(p)
  layer_names <- sapply(plot$plot$layers, function(layer) class(layer$geom)[1])
  
  # Ensure the returned plot has a clustree structure (e.g., specific aesthetic mappings)
  expect_true("GeomEdgePath" %in% layer_names)
})