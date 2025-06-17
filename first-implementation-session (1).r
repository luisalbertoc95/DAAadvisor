# First Implementation Session: Getting Started with DAAadvisor
# This is a complete example of how to start implementing the package

# =============================================================================
# STEP 1: Create Package Structure
# =============================================================================

# Run this in your terminal/RStudio
if (!require("devtools")) install.packages("devtools")
if (!require("usethis")) install.packages("usethis")
if (!require("testthat")) install.packages("testthat")

# Create package
usethis::create_package("~/DAAadvisor")

# Set up testing infrastructure
usethis::use_testthat()

# Set up git
usethis::use_git()

# Create main R files
usethis::use_r("data_profiler")
usethis::use_r("method_selector")  
usethis::use_r("method_wrappers")
usethis::use_r("utils")

# =============================================================================
# STEP 2: Create DESCRIPTION File
# =============================================================================

# In DESCRIPTION file:
usethis::use_description(
  list(
    Title = "Intelligent Differential Abundance Analysis with Automatic Method Selection",
    Description = "Provides intelligent selection of differential abundance methods 
    for microbiome data based on comprehensive data profiling. Unlike existing 
    tools that require manual interpretation of method comparisons, DAAadvisor 
    automatically recommends the most suitable methods based on data characteristics.",
    `Authors@R` = 'person("Your", "Name", email = "your.email@example.com", 
                          role = c("aut", "cre"))',
    License = "MIT + file LICENSE",
    Depends = "R (>= 4.0.0)"
  )
)

# Add dependencies
usethis::use_package("ggplot2", "Imports")
usethis::use_package("dplyr", "Imports")
usethis::use_package("stats", "Imports")

# =============================================================================
# STEP 3: Implement First Core Function with Tests
# =============================================================================

# File: R/utils.R
# Let's start with basic utility functions

#' Calculate Gini Coefficient
#'
#' @param x Numeric vector of values
#' @return Gini coefficient between 0 and 1
#' @export
#' @examples
#' calculateGini(c(1, 1, 1, 1))  # Perfect equality = 0
#' calculateGini(c(0, 0, 0, 100))  # Perfect inequality = 1
calculateGini <- function(x) {
  # Input validation
  if (!is.numeric(x)) {
    stop("Input must be numeric")
  }
  
  # Remove negative values
  x <- x[x >= 0]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  # Handle all zeros
  if (sum(x) == 0) {
    return(0)
  }
  
  # Sort values
  x <- sort(x)
  n <- length(x)
  
  # Calculate Gini coefficient
  index <- seq_len(n)
  gini <- 2 * sum(index * x) / (n * sum(x)) - (n + 1) / n
  
  return(gini)
}

#' Calculate Simpson's Diversity Index
#'
#' @param count_matrix Matrix with features as rows, samples as columns
#' @return Vector of Simpson's diversity indices for each sample
#' @export
calculateSimpsonDiversity <- function(count_matrix) {
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("Input must be a matrix or data.frame")
  }
  
  count_matrix <- as.matrix(count_matrix)
  
  # Calculate relative abundances
  rel_abund <- sweep(count_matrix, 2, colSums(count_matrix), "/")
  
  # Simpson's diversity
  simpson <- 1 - colSums(rel_abund^2, na.rm = TRUE)
  
  return(simpson)
}

# =============================================================================
# STEP 4: Write Tests First (TDD Approach)
# =============================================================================

# File: tests/testthat/test-utils.R
library(testthat)
library(DAAadvisor)

test_that("calculateGini works correctly", {
  # Perfect equality
  expect_equal(calculateGini(c(1, 1, 1, 1)), 0)
  
  # Perfect inequality
  expect_equal(calculateGini(c(0, 0, 0, 100)), 0.75)
  
  # Handle zeros
  expect_equal(calculateGini(c(0, 0, 0, 0)), 0)
  
  # Handle negative values (should be removed)
  expect_equal(calculateGini(c(-1, 1, 2, 3)), 
               calculateGini(c(1, 2, 3)))
  
  # Error on non-numeric
  expect_error(calculateGini("not numeric"))
})

test_that("calculateSimpsonDiversity works correctly", {
  # Create test matrix
  test_mat <- matrix(c(10, 0, 0,
                      0, 10, 0,
                      0, 0, 10), 
                    nrow = 3, ncol = 3)
  
  # Each sample has one dominant feature
  expect_equal(calculateSimpsonDiversity(test_mat), c(0, 0, 0))
  
  # Even distribution
  even_mat <- matrix(rep(1, 9), nrow = 3)
  simpson_even <- calculateSimpsonDiversity(even_mat)
  expect_true(all(simpson_even > 0.6))
  
  # Error on wrong input
  expect_error(calculateSimpsonDiversity("not a matrix"))
})

# =============================================================================
# STEP 5: Start Core Data Profiler Implementation
# =============================================================================

# File: R/data_profiler.R

#' S3 Class for Data Profile
#'
#' @description
#' Container for comprehensive microbiome data characteristics
#' 
#' @details
#' A DataProfile object contains:
#' - Basic metrics (samples, features, data type)
#' - Sparsity analysis at multiple scales
#' - Compositional bias metrics
#' - Zero-inflation assessment
#' - Metadata characteristics
#' 
createDataProfile <- function() {
  structure(
    list(
      data_type = NA_character_,
      n_samples = NA_integer_,
      n_features = NA_integer_,
      sparsity = list(),
      zero_inflation = list(),
      compositional_bias = list(),
      metadata_factors = list(),
      timestamp = Sys.time(),
      version = packageVersion("DAAadvisor")
    ),
    class = "DataProfile"
  )
}

#' Print method for DataProfile
#' @export
print.DataProfile <- function(x, ...) {
  cat("DAAadvisor Data Profile\n")
  cat("=======================\n")
  cat(sprintf("Data Type: %s\n", x$data_type))
  cat(sprintf("Samples: %d | Features: %d\n", x$n_samples, x$n_features))
  
  if (!is.null(x$sparsity$overall)) {
    cat(sprintf("Overall Sparsity: %.1f%%\n", x$sparsity$overall * 100))
  }
  
  if (!is.null(x$compositional_bias$score)) {
    cat(sprintf("Compositional Bias: %.2f\n", x$compositional_bias$score))
  }
  
  cat(sprintf("\nProfiled at: %s\n", format(x$timestamp, "%Y-%m-%d %H:%M:%S")))
  invisible(x)
}

#' Profile microbiome data characteristics
#'
#' @param count_table Matrix or data.frame with features as rows, samples as columns
#' @param metadata Data.frame with sample metadata  
#' @param data_type Character: "auto", "16S", "ASV", "gene", "viral"
#' @return DataProfile object
#' @export
#' @examples
#' # Create example data
#' set.seed(123)
#' counts <- matrix(rnbinom(1000, size = 0.1, mu = 10), nrow = 100)
#' metadata <- data.frame(
#'   sample_id = paste0("S", 1:10),
#'   condition = rep(c("A", "B"), each = 5)
#' )
#' 
#' profile <- profileData(counts, metadata)
#' print(profile)
profileData <- function(count_table, metadata, data_type = "auto") {
  
  # Input validation
  if (!is.matrix(count_table) && !is.data.frame(count_table)) {
    stop("count_table must be a matrix or data.frame")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data.frame")
  }
  
  # Convert to matrix
  count_mat <- as.matrix(count_table)
  
  # Check dimensions
  if (ncol(count_mat) != nrow(metadata)) {
    stop("Number of samples in count_table must match rows in metadata")
  }
  
  # Initialize profile
  profile <- createDataProfile()
  
  # Basic information
  profile$n_samples <- ncol(count_mat)
  profile$n_features <- nrow(count_mat)
  
  # Detect data type if auto
  if (data_type == "auto") {
    profile$data_type <- .detectDataType(count_mat)
  } else {
    profile$data_type <- data_type
  }
  
  # Calculate sparsity metrics
  profile$sparsity <- .calculateSparsityMetrics(count_mat)
  
  # Assess zero inflation
  profile$zero_inflation <- .assessZeroInflation(count_mat)
  
  # Detect compositional bias
  profile$compositional_bias <- .detectCompositionalBias(count_mat)
  
  # Analyze metadata
  profile$metadata_factors <- .analyzeMetadata(metadata)
  
  return(profile)
}

# Helper function implementations
.detectDataType <- function(count_mat) {
  # Simple heuristic based on sparsity and feature names
  sparsity <- mean(count_mat == 0)
  
  if (sparsity > 0.9) {
    return("viral")  # Extremely sparse
  } else if (sparsity > 0.7) {
    return("16S")    # Highly sparse
  } else {
    return("gene")   # Moderate sparsity
  }
}

.calculateSparsityMetrics <- function(count_mat) {
  list(
    overall = mean(count_mat == 0),
    by_sample = colMeans(count_mat == 0),
    by_feature = rowMeans(count_mat == 0),
    gini = calculateGini(as.vector(count_mat)),
    structural_zeros = .detectStructuralZeros(count_mat)
  )
}

.detectStructuralZeros <- function(count_mat) {
  # Check if any features are completely absent in some samples
  n_zeros_per_feature <- rowSums(count_mat == 0)
  n_samples <- ncol(count_mat)
  
  # Features missing in >90% of samples but present in some
  potentially_structural <- sum(n_zeros_per_feature > 0.9 * n_samples & 
                               n_zeros_per_feature < n_samples)
  
  return(potentially_structural > 0)
}

.assessZeroInflation <- function(count_mat) {
  # Placeholder - implement full zero-inflation test
  list(
    score = mean(count_mat == 0),
    excessive_zeros = mean(count_mat == 0) > 0.5
  )
}

.detectCompositionalBias <- function(count_mat) {
  # Placeholder - implement full compositional bias detection
  list(
    score = 0.5,
    high_bias = FALSE
  )
}

.analyzeMetadata <- function(metadata) {
  # Get basic metadata characteristics
  list(
    n_factors = ncol(metadata),
    factor_names = names(metadata),
    group_sizes = if(ncol(metadata) > 0) table(metadata[[1]]) else NULL
  )
}

# =============================================================================
# STEP 6: Create Session Context for Next Claude Session
# =============================================================================

# Save this as: session_context.R
session_context <- list(
  package_name = "DAAadvisor",
  current_status = list(
    utils = "basic functions implemented",
    data_profiler = "skeleton created, needs helper functions",
    method_selector = "not started",
    method_wrappers = "not started"
  ),
  next_tasks = c(
    "Implement .detectCompositionalBias() fully",
    "Add more sophisticated zero-inflation testing",
    "Create visualization methods for DataProfile",
    "Start method selector module"
  ),
  design_decisions = list(
    s3_classes = TRUE,
    tidyverse_style = TRUE,
    extensive_validation = TRUE
  ),
  key_functions_completed = c(
    "calculateGini",
    "calculateSimpsonDiversity", 
    "profileData (skeleton)"
  )
)

# Save progress
save(session_context, file = "session_context.RData")

# =============================================================================
# NEXT CLAUDE SESSION PROMPT TEMPLATE
# =============================================================================

next_session_prompt <- "
I'm continuing development of DAAadvisor, an R package for intelligent 
differential abundance analysis that improves on DAtest.

COMPLETED SO FAR:
- Package structure created with usethis
- Basic utility functions (Gini coefficient, Simpson diversity)
- DataProfile S3 class skeleton
- profileData() main function structure
- Unit tests for utility functions

CURRENT STATE:
- data_profiler.R has placeholder functions that need implementation
- Need to implement: .detectCompositionalBias(), .assessZeroInflation()
- Need visualization methods for DataProfile objects

NEXT GOAL:
Implement .detectCompositionalBias() with these requirements:
1. Calculate CLR transformation of count data
2. Find correlation between high-abundance and other features
3. Return proportion of negative correlations
4. Handle edge cases (all zeros, single feature)

Here's the current code structure:
[paste relevant code]

Please help me implement this function with proper error handling and tests.
"

print("Package skeleton created! Run tests with devtools::test()")
print("Next: Use the prompt template above for your next Claude session")