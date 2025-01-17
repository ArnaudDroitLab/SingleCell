library(testthat)
library(Seurat)
library(ggplot2)

# Seurat data for testing
seurat_obj <- get_demo_seurat_object_umap()


test_that("plot_seurat_dim checks for valid seurat input", {
  
  expect_s3_class(plot_seurat_dim(seurat_obj), "gg")
  
  expect_error(plot_seurat_dim(list()), "Must inherit from class 'Seurat', but has class 'list'.")
})

test_that("plot_seurat_dim handles 'colour_by' argument correctly", {
  
  # Test for 'orig.ident'
  plot_orig_ident <- plot_seurat_dim(seurat_obj, colour_by = "orig.ident")
  expect_true(inherits(plot_orig_ident, "gg"))
  
  # Test for 'nCount_RNA' (log-transformed)
  plot_nCount_RNA <- plot_seurat_dim(seurat_obj, colour_by = "nCount_RNA")
  expect_true(inherits(plot_nCount_RNA, "gg"))
  
  # Test for 'percent_mt'
  plot_percent_mt <- plot_seurat_dim(seurat_obj, colour_by = "percent_mt")
  expect_true(inherits(plot_percent_mt, "gg"))
})

test_that("plot_seurat_dim returns a ggplot object regardless of mentionning colour_by", {
  p <- plot_seurat_dim(seurat_obj)
  expect_s3_class(p, "gg")
})

test_that("plot_seurat_dim with 'nCount_RNA' transforms data correctly", {

  # Test for log2 transformation
  plot_nCount_RNA <- plot_seurat_dim(seurat_obj, colour_by = "nCount_RNA")
  # Check if the legend title includes the expected text
  plot <- ggplot2::ggplot_build(plot_nCount_RNA)
  expect_true("log2(RNA count)" %in% plot$plot$labels$colour)
})

test_that("plot_seurat_dim handles invalid reduction", {
  expect_error(plot_seurat_dim(seurat_obj, reduction = "tsne"), 
               "Reduction tsne not in seurat object.")
})

test_that("plot_seurat_dim handles missing 'colour_by' correctly", {
  expect_error(plot_seurat_dim(seurat_obj, colour_by = "nonexistent_column"),
               "None of the requested variables were found: nonexistent_column")
})

test_that("plot_seurat_dim handles custom 'assay' and 'slot' parameters", {
  
  # Test for custom assay and slot (e.g., 'RNA' assay and 'data' slot)
  plot_custom_assay <- plot_seurat_dim(seurat_obj, assay = "RNA", slot = "data")
  expect_true(inherits(plot_custom_assay, "gg"))
  
  # Test for invalid slot
  expect_error(plot_seurat_dim(seurat_obj, slot = "invalid_slot"), 
               "'arg' should be one of 'data', 'scale.data', 'counts'")
  
  # Test for invalid assay
  expect_error(plot_seurat_dim(seurat_obj, assay = "lol"),
               "Invalid assay 'lol' in Seurat object")
})
