library(testthat)
library(ggplot2)

# Example valid data frame for testing
test_df <- data.frame(
  Genes = c(200, 150, 50),
  Cells = c(100, 75, 25),
  row.names = c("Before", "After", "Filtered_out")
)

# Tests for plot_filtering_stats function
test_that("plot_filtering_stats works with valid input", {
  plot_output <- plot_filtering_stats(test_df, x_name = "Genes", y_name = "Cells")
  
  # Check if the output is a ggplot object
  expect_s3_class(plot_output, "gg")
  
  # Check if the plot contains two layers (the two bar plots)
  expect_equal(length(plot_output[[1]]$layers) + length(plot_output[[2]]$layers), 4)
})

test_that("plot_filtering_stats checks for data frame input", {
  expect_error(plot_filtering_stats("not_a_dataframe", x_name = "Genes", y_name = "Cells"), 
               "Must be of type 'data.frame', not 'character'.", fixed = TRUE)
  
  expect_error(plot_filtering_stats(NULL, x_name = "Genes", y_name = "Cells"), 
               "Assertion on 'df' failed: Must be of type 'data.frame', not 'NULL'.", fixed = TRUE)
  
  expect_error(plot_filtering_stats(42, x_name = "Genes", y_name = "Cells"), 
               "Assertion on 'df' failed: Must be of type 'data.frame', not 'double'.", fixed = TRUE)
})

test_that("plot_filtering_stats handles missing columns", {
  df_missing_genes <- data.frame(
    Cells = c(100, 75, 25),
    row.names = c("Before", "After", "Filtered_out")
  )
  
  expect_error(plot_filtering_stats(df_missing_genes, x_name = "Genes", y_name = "Cells"), 
               paste0("First column is missing from the data frame."), fixed = TRUE)  # Update this message based on your checks
  
  df_missing_cells <- data.frame(
    Genes = c(200, 100, 50),
    row.names = c("Before", "After", "Filtered_out")
  )
  
  expect_error(plot_filtering_stats(df_missing_cells, x_name = "Genes", y_name = "Cells"), 
               paste0("Second column is missing from the data frame."), fixed = TRUE)  # Update this message based on your checks
})

test_that("plot_filtering stats handles values that should not colide with the concept of filtration", {
  truncated_values <- data.frame(
    Genes = c(100, 150, 50), # Here, filtration strangely added features to the filtration.
    Cells = c(75, 50, 25),
    row.names = c("Before", "After", "Filtered_out")
  )
  expect_error(plot_filtering_stats(truncated_values, x_name = "Genes", y_name = "Cells"), 
               "Filtration cannot add features.", fixed = TRUE)
  
  truncated_values <- data.frame(
    Genes = c(150, 50, 50), # Here, filtration was strangely computed. 
    Cells = c(75, 50, 25),
    row.names = c("Before", "After", "Filtered_out")
  )
  expect_error(plot_filtering_stats(truncated_values, x_name = "Genes", y_name = "Cells"), 
               "Filtration cannot add features.", fixed = TRUE)
  
  truncated_values <- data.frame(
    Genes = c(150, 100, 50), 
    Cells = c(50, 75, 25), # Here, filtration strangely added features to the filtration.
    row.names = c("Before", "After", "Filtered_out")
  )
  expect_error(plot_filtering_stats(truncated_values, x_name = "Genes", y_name = "Cells"), 
               "Filtration cannot add features.", fixed = TRUE)
  
  truncated_values <- data.frame(
    Genes = c(150, 100, 50), 
    Cells = c(75, 50, 50), # Here, filtration was strangely computed. 
    row.names = c("Before", "After", "Filtered_out")
  )
  expect_error(plot_filtering_stats(truncated_values, x_name = "Genes", y_name = "Cells"), 
               "Filtration cannot add features.", fixed = TRUE)
})

test_that("plot_filtering_stat should indicate if the row names are not proper.", {
  test_names <- data.frame(
    Genes = c(200, 150, 50),
    Cells = c(50, 50, 0),
    row.names = c("Filtered_out", "Before", "After"))
    
  # Let's be severe on our codes as if this code did not go beyond the SingleCell package. 
  
    expect_error(plot_filtering_stats(test_names, x_name = "Genes", y_name = "Cells"), 
                 paste0("Dataframe should be ordered like this = 'Before', 'After', 'Filtered_out'."))
})