library(testthat)
library(ggplot2)  # Import ggplot2 for plotting
library(SingleCell) 

# Sample data frame for testing
test_df <- data.frame(
  group = rep(c("A", "B", "C"), each = 10),
  value = c(rnorm(10, mean = 5), rnorm(10, mean = 10), rnorm(10, mean = 15))
)

test_that("plot_filter produces a ggplot object", {
  result <- plot_filter(test_df, x_name = "group", y_name = "value")
  expect_s3_class(result, "gg")
})

test_that("plot_filter handles cutoffs correctly", {
  p <- plot_filter(test_df, x_name = "group", y_name = "value", low = 7, high = 12)
  
  # Check if the plot contains the specified cutoff lines
  low_line <- ggplot_build(p)$layout$panel_params[[1]]$y.range[1] <= 7
  high_line <- ggplot_build(p)$layout$panel_params[[1]]$y.range[2] >= 12
  
  expect_true(low_line)
  expect_true(high_line)
})

test_that("plot_filter checks for data frame input", {
  expect_error(plot_filter("not_a_dataframe", x_name = "group", y_name = "value"), 
               "Must be of type 'data.frame', not 'character", fixed = TRUE)
})

test_that("plot_filter checks for valid column names", {
  expect_error(plot_filter(test_df, x_name = "invalid_column", y_name = "value"), 
               "invalid_column not a column in df.")
  
  expect_error(plot_filter(test_df, x_name = "group", y_name = "invalid_column"), 
               "invalid_column not a column in df.")
})

test_that("plot_filter checks for numeric cutoffs", {
  expect_error(plot_filter(test_df, x_name = "group", y_name = "value", low = "low"), 
               "Must be of type 'numeric'", fixed = TRUE)
  
  expect_error(plot_filter(test_df, x_name = "group", y_name = "value", high = "high"), 
               "Must be of type 'numeric'", fixed = TRUE)
})

test_that("plot_filter checks low and high cutoff values", {
  expect_error(plot_filter(test_df, x_name = "group", y_name = "value", low = 10, high = 5), 
               "Low cannot be superior to high.", fixed = TRUE)
})

test_that("plot_filter works with default cutoff values", {
  result <- plot_filter(test_df, x_name = "group", y_name = "value")
  expect_s3_class(result, "gg")
})
