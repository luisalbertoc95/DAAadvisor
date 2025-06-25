#!/usr/bin/env Rscript

# Test the updated curatedMetagenomicData API
library(curatedMetagenomicData)

cat("Testing curatedMetagenomicData API...\n")

# Check available studies
cat("Getting sample metadata...\n")
meta <- sampleMetadata

# Look for IBD studies
ibd_studies <- meta[grepl("IBD|inflammatory", meta$disease, ignore.case = TRUE), ]
cat("Found", nrow(ibd_studies), "IBD-related samples\n")

if (nrow(ibd_studies) > 0) {
    cat("Available studies:", paste(unique(ibd_studies$study_name), collapse = ", "), "\n")
}

# Test the new API
cat("\nTesting API calls...\n")

tryCatch({
    # Try pattern-based search
    cat("Trying pattern-based search for IBD...\n")
    
    # Method 1: Search by pattern
    ibd_data <- curatedMetagenomicData("IBD", counts = FALSE, rownames = "short")
    cat("Pattern search result: ", length(ibd_data), " datasets\n")
    
    if (length(ibd_data) > 0) {
        cat("First dataset shape:", dim(ibd_data[[1]]), "\n")
    }
    
}, error = function(e) {
    cat("Error with pattern search:", e$message, "\n")
})

tryCatch({
    # Method 2: Search by study name
    cat("Trying specific study names...\n")
    
    # Get unique study names from metadata
    study_names <- unique(meta$study_name)
    ibd_study_names <- study_names[grepl("IBD|inflammatory|bowel", study_names, ignore.case = TRUE)]
    
    cat("IBD study names found:", paste(ibd_study_names, collapse = ", "), "\n")
    
    if (length(ibd_study_names) > 0) {
        # Try first IBD study
        study_data <- curatedMetagenomicData(ibd_study_names[1], counts = FALSE, rownames = "short")
        cat("Study-specific search result:", length(study_data), "datasets\n")
    }
    
}, error = function(e) {
    cat("Error with study search:", e$message, "\n")
})

cat("\nAPI test completed!\n")