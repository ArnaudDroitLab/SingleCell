library(testthat)
library(Seurat)

# Seurat object

seurat <- get_demo_seurat_object()

# Example setup: Create a mock Seurat object
create_mock_seurat <- function() {
  set.seed(123)
  data <- matrix(rnorm(2000), nrow = 200, ncol = 10)
  rownames(data) <- paste0("Gene", 1:200)
  colnames(data) <- paste0("Cell", 1:10)
  
  seurat_obj <- CreateSeuratObject(counts = data)
  
  # Remove zero variance features
  seurat_obj <- seurat_obj[Matrix::rowSums(seurat_obj@assays$RNA@counts) > 0, ]
  
  # Identify variable features
  seurat_obj <- FindVariableFeatures(seurat_obj)
  
  # Check how many variable features were found
  variable_features <- VariableFeatures(seurat_obj)
  if (length(variable_features) == 0) {
    stop("No variable features found; ensure your data has sufficient variability.")
  }
  
  # Scale the data
  seurat_obj <- ScaleData(seurat_obj)
  
  # Set npc safely
  npc <- min(30, length(variable_features), ncol(seurat_obj) - 1)
  
  # Check for enough cells
  if (ncol(seurat_obj) < 2) {
    stop("Not enough cells to perform PCA.")
  }
  
  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = npc)
  
  return(seurat_obj)
}



# Create a mock Seurat object for testing
mock_seurat <- create_mock_seurat()

test_that("plot_seurat_elbow handles correct inputs", {
  expect_s3_class(plot_seurat_elbow(mock_seurat, npc = min(30, length(VariableFeatures(mock_seurat)), ncol(mock_seurat) - 1)), "gg")
})

test_that("plot_seurat_elbow checks Seurat object class", {
  expect_error(plot_seurat_elbow("not_a_seurat"), 
               "Assertion on 'seurat' failed: Must inherit from class 'Seurat', but has class 'character'.", fixed = TRUE)
})

test_that("plot_seurat_elbow checks for valid reduction", {
  expect_error(plot_seurat_elbow(mock_seurat, reduction = "invalid_reduction"), 
               "not in seurat object.")
})

test_that("plot_seurat_elbow checks npc against available components", {
  expect_error(plot_seurat_elbow(mock_seurat, npc = 20), 
               "pca does not have 20 components.", fixed = TRUE)
})

test_that("plot_seurat_elbow works with default parameters", {
  plot <- plot_seurat_elbow(mock_seurat)
  expect_s3_class(plot, "gg")
})

test_that("plot_seurat_elbow handles custom npc and k.param.neighbors", {
  plot <- plot_seurat_elbow(mock_seurat, npc = 5, k.param.neighbors = 3)
  expect_s3_class(plot, "gg")
  expect_true("Standard Deviation" %in% colnames(plot$data))
})
