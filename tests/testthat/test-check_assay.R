library(testthat)
library(Seurat)

# Unit tests for check_assay function
test_that("check_assay works correctly", {
  
  # Create a dummy Seurat object
  seurat_obj <- get_demo_seurat_object()
  
  # Add an assay to the Seurat object (other than the default RNA assay)
  seurat_obj[["assay1"]] <- seurat[["RNA"]]
  
  # Test that check_assay works with valid Seurat object and assay name
  expect_silent(check_assay(seurat_obj, assay = "RNA"))
  expect_silent(check_assay(seurat_obj, assay = "assay1"))
  
  # Test that check_assay throws an error if the assay is not present
  expect_error(check_assay(seurat_obj, assay = "nonexistent_assay"), 
               "Assay nonexistent_assay not in seurat object.")
  
  # Test that check_assay throws an error if the Seurat object is not of class 'Seurat'
  expect_error(check_assay(list(), assay = "RNA"), 
               "Assertion on 'seurat' failed: Must inherit from class 'Seurat', but has class 'list'")
  
  # Test that check_assay throws an error if the assay is not a string
  expect_error(check_assay(seurat_obj, assay = 123), 
               "Assertion on 'assay' failed: Must be of type 'string', not 'double'.")
})