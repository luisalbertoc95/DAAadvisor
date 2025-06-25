#!/usr/bin/env Rscript

library(curatedMetagenomicData)

# Test manual download with correct accessing
cat("Testing access methods...\n")

# Download HMP IBD data
cmd_data <- curatedMetagenomicData("HMP_2019_ibdmdb", counts = FALSE, rownames = "short")

cat("Downloaded", length(cmd_data), "datasets\n")

# Check if it's a list or ExperimentList
cat("Class:", class(cmd_data), "\n")

# Try different access methods
if (length(cmd_data) > 0) {
    # Method 1: Access by index and check
    for (i in 1:min(3, length(cmd_data))) {  # Check first 3
        cat("Dataset", i, ":\n")
        dataset <- cmd_data[[i]]
        cat("  Class:", class(dataset), "\n")
        cat("  Dimensions:", dim(dataset), "\n")
        
        # Check if it's relative abundance
        if ("assays" %in% slotNames(dataset)) {
            assay_names <- names(assays(dataset))
            cat("  Assay names:", paste(assay_names, collapse = ", "), "\n")
        }
        
        # Check metadata for this dataset
        if ("colData" %in% slotNames(dataset)) {
            metadata <- colData(dataset) 
            cat("  Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")
            cat("  Sample count:", nrow(metadata), "\n")
        }
        cat("\n")
    }
}