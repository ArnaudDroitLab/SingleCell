library(testthat)
library(Seurat)

# Unit tests for check_reduction function
test_that("check_reduction works correctly", {
  
  # Create a dummy Seurat object
  seurat_obj <- get_demo_seurat_object_umap()
  
  # Add a reduction to the Seurat object (other than the default PCA)
  seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = matrix(1:100, ncol = 10), key = "UMAP_")
  
  # Test that check_reduction works with valid Seurat object and reduction name
  expect_silent(check_reduction(seurat_obj, reduction = "pca"))
  expect_silent(check_reduction(seurat_obj, reduction = "umap"))
  
  # Test that check_reduction throws an error if the reduction is not present
  expect_error(check_reduction(seurat_obj, reduction = "tsne"), 
               "Reduction tsne not in seurat object.")
  
  # Test that check_reduction throws an error if the Seurat object is not of class 'Seurat'
  expect_error(check_reduction(list(), reduction = "pca"), 
               "Assertion on 'seurat' failed: Must inherit from class 'Seurat', but has class 'list'")
  
  # Test that check_reduction throws an error if the reduction is not a string
  expect_error(check_reduction(seurat_obj, reduction = 123), 
               "Assertion on 'reduction' failed: Must be of type 'string', not 'double'.")
})