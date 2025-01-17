library(testthat)
library(Seurat)
library(ggplot2)

# Seurat data for testing
seurat <- get_demo_seurat_object_pca()


test_that("plot_seurat_elbow throws an error with invalid Seurat object", {
  
  expect_error(plot_seurat_elbow(1, reduction = "pca", npc = 20, k.param.neighbors = 10), 
               "Assertion on 'seurat' failed: Must inherit from class 'Seurat', but has class 'numeric'.")
})

test_that("plot_seurat_elbow throws an error with invalid reduction", {

  expect_error(plot_seurat_elbow(seurat, reduction = "invalid_reduction", npc = 20, k.param.neighbors = 10),
               "Reduction invalid_reduction not in seurat object.")
})

test_that("plot_seurat_elbow checks for numeric npc and param", {
  
  # Attempt to use a non-numeric value for npc
  expect_error(plot_seurat_elbow(seurat, reduction = "pca", npc = "help", k.param.neighbors = 10),
               "Assertion on 'npc' failed: Must be of type 'single integerish value', not 'character'.")
  
  # Attempt to use a non-numeric value for param
  expect_error(plot_seurat_elbow(seurat, reduction = "pca", npc = 50, k.param.neighbors = "me"),
               "Assertion on 'k.param.neighbors' failed: Must be of type 'single integerish value', not 'character'.")
  
})

test_that("plot_seurat_elbow throws an error with invalid npc", {
  
  # Attempt to use more components than available
  expect_error(plot_seurat_elbow(seurat, reduction = "pca", npc = 80, k.param.neighbors = 10),
               "pca does not have 80 components.")
})

test_that("plot_seurat_elbow outputs a proper ggplot with a vline", {
  
  # Call the function with valid inputs
  p <- plot_seurat_elbow(seurat, reduction = "pca", npc = 20, k.param.neighbors = 10)
  
  # Check if a ggplot object is returned
  expect_s3_class(p, "gg")
  
  # Check if plot contains the expected vertical line (k.param.neighbors)
  expect_true(any(grepl("xintercept", as.character(ggplot_build(p)))))
})

test_that("plot_seurat_elbow generates the plot without errors", {
  
  p <- plot_seurat_elbow(seurat, reduction = "pca", npc = 20, k.param.neighbors = 10)
  plot <- ggplot2::ggplot_build(p)
  
  # Check that the plot was generated and contains expected components
  expect_s3_class(p, "gg")
  expect_true("Standard_Deviation" %in% ggplot2::ggplot_build(p)$plot$labels$y)
})

