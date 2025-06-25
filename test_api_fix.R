#!/usr/bin/env Rscript

# Test the correct 2-step API approach
library(curatedMetagenomicData)
library(dplyr)
library(readr)

cat("Testing 2-step curatedMetagenomicData API...\n")

# Step 1: Get available dataset names
dataset_names <- curatedMetagenomicData("HMP_2019_ibdmdb", counts = FALSE, rownames = "short")

cat("Available datasets:\n")
print(dataset_names)

# Step 2: Filter for relative abundance
rel_abund_names <- dataset_names[grepl("relative_abundance", dataset_names)]
cat("\nRelative abundance datasets:\n")
print(rel_abund_names)

if (length(rel_abund_names) > 0) {
    # Step 3: Download specific dataset by exact name
    dataset_name <- rel_abund_names[1]
    cat("\nDownloading:", dataset_name, "\n")
    
    # This should return the actual data object
    data_obj <- curatedMetagenomicData(dataset_name, counts = FALSE, rownames = "short")
    
    cat("Downloaded object class:", class(data_obj), "\n")
    cat("Data dimensions:", dim(data_obj), "\n")
    
    # Test data access
    abundance_matrix <- assay(data_obj)
    metadata_df <- as.data.frame(colData(data_obj))
    
    cat("Abundance matrix shape:", dim(abundance_matrix), "\n")
    cat("Metadata shape:", dim(metadata_df), "\n")
    cat("Metadata columns:", paste(colnames(metadata_df), collapse = ", "), "\n")
    
    if ("disease" %in% colnames(metadata_df)) {
        disease_counts <- table(metadata_df$disease)
        cat("Disease distribution:\n")
        print(disease_counts)
    }
    
    cat("\n✅ API fix validated successfully!\n")
} else {
    cat("❌ No relative abundance datasets found\n")
}