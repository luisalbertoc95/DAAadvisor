#!/usr/bin/env Rscript

# Test different approaches to get actual data objects
library(curatedMetagenomicData)

cat("Testing different API approaches...\n")

# Try different parameters
cat("\nMethod 1: dryrun = FALSE\n")
tryCatch({
    data1 <- curatedMetagenomicData("2021-03-31.HMP_2019_ibdmdb.relative_abundance", 
                                   counts = FALSE, 
                                   dryrun = FALSE,
                                   rownames = "short")
    cat("Method 1 result:", class(data1), length(data1), "\n")
}, error = function(e) cat("Method 1 error:", e$message, "\n"))

# Try the ExperimentHub approach
cat("\nMethod 2: Direct ExperimentHub\n")
tryCatch({
    library(ExperimentHub)
    eh <- ExperimentHub()
    # Query for the specific dataset
    query_result <- query(eh, "HMP_2019_ibdmdb")
    cat("ExperimentHub query results:", length(query_result), "\n")
    if (length(query_result) > 0) {
        # Try to get the first result
        data2 <- eh[[names(query_result)[1]]]
        cat("Method 2 result:", class(data2), "\n")
    }
}, error = function(e) cat("Method 2 error:", e$message, "\n"))

# Try without dryrun parameter
cat("\nMethod 3: No dryrun parameter\n")
tryCatch({
    data3 <- curatedMetagenomicData("2021-03-31.HMP_2019_ibdmdb.relative_abundance", 
                                   counts = FALSE, 
                                   rownames = "short")
    cat("Method 3 result:", class(data3), length(data3), "\n")
    
    # If it's still a string, maybe we need to process it differently
    if (is.character(data3)) {
        cat("Still getting strings. Checking if it's actually stored somewhere...\n")
        # Maybe the data is cached locally?
        cat("Checking global environment...\n")
        ls_result <- ls(.GlobalEnv)
        cat("Global objects:", paste(ls_result, collapse = ", "), "\n")
    }
}, error = function(e) cat("Method 3 error:", e$message, "\n"))

# Check the package version and help
cat("\nPackage version:\n")
packageVersion("curatedMetagenomicData")

cat("\nChecking if dryrun=FALSE actually downloads...\n")
try({
    # Try with explicit dryrun=FALSE 
    test_data <- curatedMetagenomicData("HMP_2019_ibdmdb", 
                                       counts = FALSE, 
                                       dryrun = FALSE, 
                                       rownames = "short")
    cat("With dryrun=FALSE:", class(test_data), "\n")
    
    if (length(test_data) > 0 && !is.character(test_data)) {
        cat("Success! Got actual data object\n")
        if (length(test_data) > 0) {
            first_obj <- test_data[[1]]
            cat("First object class:", class(first_obj), "\n")
        }
    }
})