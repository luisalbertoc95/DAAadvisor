#!/bin/bash
# Script to organize all DAAadvisor artifacts after manual export
# Save this as: organize_artifacts.sh

# Create directory structure
echo "Creating DAAadvisor directory structure..."

mkdir -p DAAadvisor/{R,inst/{templates,examples},vignettes,docs,python_prototype}

# ==============================================================================
# ARTIFACT INDEX - What to save and where
# ==============================================================================

cat > DAAadvisor/ARTIFACTS_INDEX.md << 'EOF'
# DAAadvisor Artifacts Index

## Artifacts to Export (in order of creation):

### 1. Planning & Design Documents
- **Artifact**: daa-tool-implementation-plan
- **Type**: Markdown
- **Save as**: `docs/implementation_plan.md`
- **Description**: Comprehensive implementation plan and architecture

- **Artifact**: implementation-strategy
- **Type**: Markdown  
- **Save as**: `docs/implementation_strategy.md`
- **Description**: Strategy for building with Claude

### 2. R Package Core Files
- **Artifact**: daa-advisor-r-package
- **Type**: Markdown
- **Save as**: `docs/package_design.md`
- **Description**: Main R package structure and core functions

- **Artifact**: daa-method-selection-algorithm
- **Type**: R code
- **Save as**: `R/method_selector_core.R`
- **Description**: Core algorithm for intelligent method selection

- **Artifact**: first-implementation-session
- **Type**: R code
- **Save as**: `R/initial_implementation.R`
- **Description**: Starting implementation with tests

- **Artifact**: daa-advisor-example
- **Type**: R code
- **Save as**: `inst/examples/comparison_with_DAtest.R`
- **Description**: Practical comparison example

### 3. Comparison & Analysis Framework
- **Artifact**: daa-comparison-framework
- **Type**: Markdown
- **Save as**: `docs/comparison_framework.md`
- **Description**: Advanced comparison and assessment methods

### 4. Python Prototype (for reference)
- **Artifact**: daa-tool-starter-code
- **Type**: Python
- **Save as**: `python_prototype/daa_tool.py`
- **Description**: Initial Python prototype

- **Artifact**: daa-visualization-module
- **Type**: Python
- **Save as**: `python_prototype/visualization.py`
- **Description**: Visualization module in Python

- **Artifact**: daa-method-wrapper-template
- **Type**: Python
- **Save as**: `python_prototype/method_wrappers.py`
- **Description**: Method wrapper templates

### 4. Documentation
- **Artifact**: daa-requirements-readme
- **Type**: Markdown
- **Save as**: `README.md` (in root)
- **Description**: Package README and requirements

## Export Instructions:
1. Click on each artifact in Claude
2. Click the download/save button
3. Save with the filename specified above
4. Run this organization script after saving all files
EOF

# ==============================================================================
# File organization script
# ==============================================================================

cat > DAAadvisor/organize_files.R << 'EOF'
# R script to split combined files into proper package structure
# Run after exporting all artifacts

# Read the main R package design
if (file.exists("docs/package_design.md")) {
  content <- readLines("docs/package_design.md")
  
  # Extract R code sections
  in_code <- FALSE
  current_file <- NULL
  current_content <- character()
  
  for (line in content) {
    if (grepl("^```r$", line) && grepl("# (.*\\.R)$", content[which(content == line) - 1])) {
      in_code <- TRUE
      current_file <- gsub(".*# (.*\\.R).*", "\\1", content[which(content == line) - 1])
      current_content <- character()
    } else if (grepl("^```$", line) && in_code) {
      in_code <- FALSE
      if (!is.null(current_file)) {
        writeLines(current_content, current_file)
        cat("Created:", current_file, "\n")
      }
    } else if (in_code) {
      current_content <- c(current_content, line)
    }
  }
}

# Extract code from method selection algorithm
if (file.exists("R/method_selector_core.R")) {
  # This file is already in the right place
  cat("Method selector core already in place\n")
}

# Create proper package structure files
cat("Creating package infrastructure...\n")

# DESCRIPTION file
description <- 'Package: DAAadvisor
Type: Package
Title: Intelligent Differential Abundance Analysis with Automatic Method Selection
Version: 0.1.0
Authors@R: 
    person("Your", "Name", , "your.email@example.com", role = c("aut", "cre"))
Description: Provides intelligent selection of differential abundance methods 
    for microbiome data based on comprehensive data profiling. Unlike existing 
    tools that require manual interpretation of method comparisons, DAAadvisor 
    automatically recommends the most suitable methods based on data characteristics.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.3
Depends: R (>= 4.0.0)
Imports:
    stats,
    methods,
    ggplot2,
    dplyr
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
VignetteBuilder: knitr
'
writeLines(description, "DESCRIPTION")

# Create .Rbuildignore
buildignore <- 'docs/
python_prototype/
organize_files.R
ARTIFACTS_INDEX.md
.git
'
writeLines(buildignore, ".Rbuildignore")

cat("\nPackage structure created successfully!\n")
cat("Next steps:\n")
cat("1. Run devtools::document() to generate documentation\n")
cat("2. Run devtools::test() to run tests\n")
cat("3. Run devtools::check() to check package\n")
EOF

# Make scripts executable
chmod +x DAAadvisor/organize_files.R

# ==============================================================================
# Quick reference for manual export
# ==============================================================================

cat > DAAadvisor/EXPORT_CHECKLIST.md << 'EOF'
# Export Checklist

## Step 1: Export each artifact from Claude
- [ ] daa-tool-implementation-plan → docs/implementation_plan.md
- [ ] daa-advisor-r-package → docs/package_design.md
- [ ] daa-comparison-framework → docs/comparison_framework.md
- [ ] daa-method-selection-algorithm → R/method_selector_core.R
- [ ] first-implementation-session → R/initial_implementation.R
- [ ] daa-advisor-example → inst/examples/comparison_with_DAtest.R
- [ ] implementation-strategy → docs/implementation_strategy.md
- [ ] daa-requirements-readme → README.md
- [ ] daa-tool-starter-code → python_prototype/daa_tool.py
- [ ] daa-visualization-module → python_prototype/visualization.py
- [ ] daa-method-wrapper-template → python_prototype/method_wrappers.py

## Step 2: Run organization
```bash
cd DAAadvisor
Rscript organize_files.R
```

## Step 3: Initialize R package
```r
library(devtools)
setwd("DAAadvisor")
document()
test()
check()
```

## File Structure After Organization:
```
DAAadvisor/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── R/
│   ├── method_selector_core.R
│   ├── initial_implementation.R
│   ├── data_profiler.R
│   ├── method_wrappers.R
│   └── utils.R
├── inst/
│   ├── examples/
│   │   └── comparison_with_DAtest.R
│   └── templates/
├── docs/
│   ├── implementation_plan.md
│   ├── package_design.md
│   ├── comparison_framework.md
│   └── implementation_strategy.md
├── python_prototype/
│   ├── daa_tool.py
│   ├── visualization.py
│   └── method_wrappers.py
└── tests/
    └── testthat/
```
EOF

echo "Organization script created!"
echo ""
echo "INSTRUCTIONS:"
echo "1. Manually export each artifact from Claude using the names in ARTIFACTS_INDEX.md"
echo "2. Save them in a folder called 'DAAadvisor'"
echo "3. Run this script to organize them properly"
echo "4. Follow the checklist in EXPORT_CHECKLIST.md"
echo ""
echo "The artifacts contain:"
echo "- Complete R package implementation plan"
echo "- Core algorithm implementations"
echo "- Comparison with DAtest"
echo "- Python prototype for reference"
echo "- Full documentation and examples"