#!/usr/bin/env Rscript

library(curatedMetagenomicData)
library(dplyr)
library(readr)

# Test manual download
cat("Testing manual download...\n")

# Download HMP IBD data
cmd_data <- curatedMetagenomicData("HMP_2019_ibdmdb", counts = FALSE, rownames = "short")

cat("Downloaded", length(cmd_data), "datasets\n")
cat("Dataset names:\n")
print(names(cmd_data))

# Look for relative abundance
rel_abund_names <- names(cmd_data)[grepl("relative_abundance", names(cmd_data))]
cat("Relative abundance datasets:\n")
print(rel_abund_names)

if (length(rel_abund_names) > 0) {
    # Use the most recent one
    latest_rel_abund <- rel_abund_names[length(rel_abund_names)]
    cat("Using dataset:", latest_rel_abund, "\n")
    
    data <- cmd_data[[latest_rel_abund]]
    
    cat("Data dimensions:", dim(data), "\n")
    cat("Assay dimensions:", dim(assay(data)), "\n")
    
    # Check metadata
    metadata <- colData(data)
    cat("Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")
    
    if ("disease" %in% colnames(metadata)) {
        unique_diseases <- unique(metadata$disease)
        cat("Unique diseases:", paste(unique_diseases, collapse = ", "), "\n")
        
        disease_counts <- table(metadata$disease)
        print(disease_counts)
    }
}