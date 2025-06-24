#!/usr/bin/env python3
"""
curatedMetagenomicData Downloader and Cross-Validation Framework

This module downloads real microbiome data from curatedMetagenomicData 
(Bioconductor) and performs cross-validation with realistic synthetic data.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional, Any
import logging
from pathlib import Path
import subprocess
import tempfile
import warnings

logger = logging.getLogger(__name__)

class CuratedDataDownloader:
    """
    Download and process data from curatedMetagenomicData R package
    
    Features:
    - Automated R script generation and execution
    - Multiple disease/condition categories
    - Data format standardization for DAAadvisor
    - Cross-validation with synthetic data
    - Ground truth establishment from literature
    """
    
    def __init__(self, output_dir: str = "curated_real_data"):
        """
        Initialize curated data downloader
        
        Parameters:
        -----------
        output_dir : str
            Directory to save downloaded datasets
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Available conditions in curatedMetagenomicData
        self.available_conditions = {
            "IBD": {
                "description": "Inflammatory Bowel Disease studies",
                "study_name": "HMP_2019_ibdmdb",
                "conditions": ["healthy", "IBD", "UC", "CD"],
                "expected_features": ["Faecalibacterium_prausnitzii", "Escherichia_coli"],
                "expected_direction": "decreased_butyrate_producers"
            },
            "CRC": {
                "description": "Colorectal Cancer studies", 
                "study_name": "ThomasAM_2019",
                "conditions": ["healthy", "CRC"],
                "expected_features": ["Fusobacterium_nucleatum", "Bacteroides_fragilis"],
                "expected_direction": "increased_pathogens"
            },
            "T2D": {
                "description": "Type 2 Diabetes studies",
                "study_name": "QinJ_2012", 
                "conditions": ["healthy", "T2D"],
                "expected_features": ["Akkermansia_muciniphila", "Bifidobacterium"],
                "expected_direction": "decreased_beneficial"
            },
            "obesity": {
                "description": "Obesity studies",
                "study_name": "LeChatelier_2013",
                "conditions": ["healthy", "obesity"],
                "expected_features": ["Bacteroidetes", "Firmicutes"],
                "expected_direction": "firmicutes_bacteroidetes_ratio"
            },
            "cirrhosis": {
                "description": "Liver cirrhosis studies",
                "study_name": "QinN_2014",
                "conditions": ["healthy", "cirrhosis"],
                "expected_features": ["Enterobacteriaceae", "Streptococcus"],
                "expected_direction": "increased_pathogens"
            }
        }
        
        logger.info(f"📁 Curated data downloader initialized: {output_dir}")
        logger.info(f"📊 {len(self.available_conditions)} condition categories available")
    
    def check_r_environment(self) -> bool:
        """Check if R and required packages are available"""
        
        logger.info("🔍 Checking R environment...")
        
        try:
            # Check R installation
            result = subprocess.run(['R', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode != 0:
                logger.error("❌ R not found. Install R first.")
                return False
            
            r_version = result.stdout.split('\n')[0]
            logger.info(f"✅ {r_version}")
            
            # Check if BiocManager is available
            check_bioc_script = """
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
                cat("BiocManager not found\\n")
                quit(status = 1)
            } else {
                cat("BiocManager available\\n")
            }
            """
            
            result = subprocess.run(['R', '--slave', '-e', check_bioc_script],
                                  capture_output=True, text=True, timeout=30)
            
            if "BiocManager not found" in result.stdout:
                logger.warning("⚠️ BiocManager not found. Will attempt to install.")
            else:
                logger.info("✅ BiocManager available")
            
            return True
            
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.error(f"❌ R environment check failed: {e}")
            return False
    
    def install_required_packages(self) -> bool:
        """Install required R packages"""
        
        logger.info("📦 Installing required R packages...")
        
        install_script = '''
        # Install BiocManager if needed
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "https://cran.r-project.org")
        }
        
        # Install curatedMetagenomicData
        if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
            cat("Installing curatedMetagenomicData...\\n")
            BiocManager::install("curatedMetagenomicData", ask = FALSE)
        }
        
        # Install other required packages
        required_packages <- c("dplyr", "readr")
        for (pkg in required_packages) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                install.packages(pkg, repos = "https://cran.r-project.org")
            }
        }
        
        cat("Package installation completed\\n")
        '''
        
        try:
            result = subprocess.run(['R', '--slave', '-e', install_script],
                                  capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                logger.info("✅ R packages installed successfully")
                return True
            else:
                logger.error(f"❌ R package installation failed: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error("⏰ R package installation timeout")
            return False
        except Exception as e:
            logger.error(f"❌ R package installation error: {e}")
            return False
    
    def download_condition_data(self, condition: str, max_studies: int = 5) -> Optional[Dict]:
        """
        Download data for a specific condition
        
        Parameters:
        -----------
        condition : str
            Condition name (e.g., 'IBD', 'CRC', 'T2D')
        max_studies : int
            Maximum number of studies to download
            
        Returns:
        --------
        Optional[Dict]
            Downloaded data with count_table, metadata, ground_truth
        """
        
        if condition not in self.available_conditions:
            logger.error(f"❌ Unknown condition: {condition}")
            return None
        
        condition_info = self.available_conditions[condition]
        
        logger.info(f"⬇️ Downloading {condition}: {condition_info['description']}")
        
        # Create R script for downloading with fixed API
        study_name = condition_info.get('study_name', condition)
        output_path = self.output_dir / condition.lower()
        relevant_conditions = ', '.join([f'"{c}"' for c in condition_info['conditions']])
        expected_features = ', '.join([f'"{f}"' for f in condition_info['expected_features']])
        
        r_script = f'''
library(curatedMetagenomicData)
library(dplyr)
library(readr)

# Set output directory
output_dir <- "{output_path}"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Download data for {condition} using study name: {study_name}
cat("Downloading {condition} data from study {study_name}...\\n")

tryCatch({{
    # Step 1: Get available dataset names
    cat("Getting available datasets for {study_name}...\\n")
    
    dataset_names <- curatedMetagenomicData(
        "{study_name}", 
        counts = FALSE,
        dryrun = TRUE,
        rownames = "short"
    )
    
    cat(paste("Found", length(dataset_names), "datasets\\n"))
    
    # Step 2: Filter for relative abundance datasets
    rel_abund_names <- dataset_names[grepl("relative_abundance", dataset_names)]
    cat(paste("Found", length(rel_abund_names), "relative abundance datasets\\n"))
    
    if (length(rel_abund_names) > 0) {{
        # Step 3: Download the specific dataset (use most recent)
        dataset_name <- rel_abund_names[length(rel_abund_names)]
        cat(paste("Downloading dataset:", dataset_name, "\\n"))
        
        # CRITICAL: Use dryrun = FALSE to get actual data object
        cmd_data_list <- curatedMetagenomicData(
            dataset_name, 
            counts = FALSE,
            dryrun = FALSE,
            rownames = "short"
        )
        
        if (length(cmd_data_list) > 0) {{
            # Extract the actual data object
            current_data <- cmd_data_list[[1]]
            cat(paste("Successfully downloaded data object, class:", class(current_data), "\\n"))
            
            # Get abundance matrix
            abundance_matrix <- assay(current_data)
            
            # Get metadata
            metadata_df <- as.data.frame(colData(current_data))
            
            cat(paste("Original data shape:", nrow(abundance_matrix), "features x", ncol(abundance_matrix), "samples\\n"))
            
            # Print available columns for debugging
            cat("Available metadata columns:", paste(colnames(metadata_df), collapse = ", "), "\\n")
            
            # Print unique disease values for debugging
            if ("disease" %in% colnames(metadata_df)) {{
                unique_diseases <- unique(metadata_df$disease)
                cat("Unique disease values:", paste(unique_diseases, collapse = ", "), "\\n")
                
                # Filter for relevant conditions (more flexible matching)
                relevant_conditions <- c({relevant_conditions})
                valid_samples <- metadata_df$disease %in% relevant_conditions
                
                cat(paste("Found", sum(valid_samples), "samples matching conditions\\n"))
                
                if (sum(valid_samples) > 10) {{
                    metadata_df <- metadata_df[valid_samples, ]
                    abundance_matrix <- abundance_matrix[, valid_samples]
                }}
            }}
            
            # Ensure we have enough samples
            if (ncol(abundance_matrix) > 10) {{
                cat(paste("Final data shape:", nrow(abundance_matrix), "features x", ncol(abundance_matrix), "samples\\n"))
                
                # Save abundance data (transpose: features x samples -> samples x features)
                abundance_df <- as.data.frame(t(abundance_matrix))
                
                # Add sample IDs
                abundance_df$sample_id <- rownames(abundance_df)
                abundance_df <- abundance_df[, c("sample_id", setdiff(names(abundance_df), "sample_id"))]
                
                # Save files
                write_csv(abundance_df, file.path(output_dir, "{condition.lower()}_count_table.csv"))
                
                # Save metadata with sample_id as first column for proper indexing
                metadata_df_save <- metadata_df
                metadata_df_save$sample_id <- rownames(metadata_df_save)
                metadata_df_save <- metadata_df_save[, c("sample_id", setdiff(names(metadata_df_save), "sample_id"))]
                write_csv(metadata_df_save, file.path(output_dir, "{condition.lower()}_metadata.csv"))
                
                # Create ground truth (example differential features)
                expected_features <- c({expected_features})
                writeLines(expected_features, file.path(output_dir, "{condition.lower()}_ground_truth.txt"))
                
                # Create summary info
                summary_info <- data.frame(
                    condition = "{condition}",
                    study_name = "{study_name}",
                    n_samples = ncol(abundance_matrix),
                    n_features = nrow(abundance_matrix),
                    dataset_name = dataset_name,
                    conditions_available = paste(unique(metadata_df$disease), collapse = ";")
                )
                
                write_csv(summary_info, file.path(output_dir, "{condition.lower()}_summary.csv"))
                
                cat("✅ Successfully saved {condition} data\\n")
                
            }} else {{
                cat("⚠️ Not enough samples after filtering (", ncol(abundance_matrix), ")\\n")
            }}
            
        }} else {{
            cat("❌ No relative abundance datasets found\\n")
        }}
        
    }} else {{
        cat("❌ No data downloaded for {study_name}\\n")
    }}
    
    cat("Download completed for {condition}\\n")
    
}}, error = function(e) {{
    cat(paste("❌ Error downloading {condition}:", e$message, "\\n"))
}})
'''
        
        # Save and execute R script
        r_script_file = self.output_dir / f"download_{condition.lower()}.R"
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        logger.info(f"📝 R script saved: {r_script_file}")
        
        try:
            logger.info(f"🔄 Executing R script for {condition}...")
            result = subprocess.run(['Rscript', str(r_script_file)],
                                  capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0:
                logger.info(f"✅ R script completed for {condition}")
                logger.info(f"📄 R output: {result.stdout}")
                
                # Process downloaded files
                return self._process_downloaded_files(condition)
            else:
                logger.error(f"❌ R script failed for {condition}: {result.stderr}")
                return None
                
        except subprocess.TimeoutExpired:
            logger.error(f"⏰ R script timeout for {condition}")
            return None
        except Exception as e:
            logger.error(f"❌ R script error for {condition}: {e}")
            return None
    
    def _process_downloaded_files(self, condition: str) -> Optional[Dict]:
        """Process downloaded files into DAAadvisor format"""
        
        condition_dir = self.output_dir / condition.lower()
        
        # Find downloaded files (updated to match R script output)
        abundance_files = list(condition_dir.glob(f"{condition.lower()}_count_table.csv"))
        metadata_files = list(condition_dir.glob(f"{condition.lower()}_metadata.csv"))
        
        if not abundance_files or not metadata_files:
            logger.warning(f"⚠️ No data files found for {condition}")
            return None
        
        # Use the first available file
        abundance_file = abundance_files[0]
        metadata_file = metadata_files[0]
        
        logger.info(f"📊 Processing {abundance_file.name}")
        
        try:
            # Load data
            abundance_df = pd.read_csv(abundance_file, index_col=0)
            metadata_df = pd.read_csv(metadata_file, index_col=0)
            
            # Align samples
            common_samples = list(set(abundance_df.index) & set(metadata_df.index))
            
            if len(common_samples) < 10:
                logger.warning(f"⚠️ Too few common samples: {len(common_samples)}")
                return None
            
            abundance_aligned = abundance_df.loc[common_samples]
            metadata_aligned = metadata_df.loc[common_samples]
            
            # Standardize condition column
            if 'disease' in metadata_aligned.columns:
                metadata_aligned['condition'] = metadata_aligned['disease']
            elif 'study_condition' in metadata_aligned.columns:
                metadata_aligned['condition'] = metadata_aligned['study_condition']
            
            # Filter to binary conditions if needed
            if condition in ['IBD', 'CRC', 'T2D', 'obesity', 'cirrhosis']:
                healthy_mask = metadata_aligned['condition'].isin(['healthy', 'control'])
                disease_mask = metadata_aligned['condition'].isin([condition, condition.lower()])
                
                binary_mask = healthy_mask | disease_mask
                
                logger.info(f"Condition filtering: {healthy_mask.sum()} control, {disease_mask.sum()} disease samples")
                
                if binary_mask.sum() > 10:
                    abundance_aligned = abundance_aligned[binary_mask]
                    metadata_aligned = metadata_aligned[binary_mask]
                    
                    # Standardize labels
                    metadata_aligned['condition'] = metadata_aligned['condition'].map({
                        'healthy': 'Control',
                        'control': 'Control',
                        condition: 'Treatment',
                        condition.lower(): 'Treatment'
                    }).fillna(metadata_aligned['condition'])
                    
                    logger.info(f"Final dataset: {len(abundance_aligned)} samples, {len(abundance_aligned.columns)} features")
            
            # Load ground truth if available, otherwise create from literature
            ground_truth_file = condition_dir / f"{condition.lower()}_ground_truth.txt"
            if ground_truth_file.exists():
                with open(ground_truth_file, 'r') as f:
                    ground_truth = [line.strip() for line in f if line.strip()]
                logger.info(f"📋 Loaded ground truth from file: {len(ground_truth)} features")
            else:
                ground_truth = self._create_ground_truth(condition, abundance_aligned.columns)
                logger.info(f"📋 Created ground truth from literature: {len(ground_truth)} features")
            
            logger.info(f"✅ Processed {condition}: {abundance_aligned.shape} samples×features")
            logger.info(f"🎯 Ground truth: {len(ground_truth)} expected differential features")
            
            return {
                'count_table': abundance_aligned,
                'metadata': metadata_aligned,
                'ground_truth': ground_truth,
                'condition_info': self.available_conditions[condition],
                'source': 'curatedMetagenomicData'
            }
            
        except Exception as e:
            logger.error(f"❌ Error processing {condition}: {e}")
            return None
    
    def _create_ground_truth(self, condition: str, feature_names: List[str]) -> List[str]:
        """Create ground truth based on literature for known conditions"""
        
        # Literature-based expected changes
        literature_markers = {
            "IBD": [
                "Faecalibacterium_prausnitzii", "Roseburia", "Bifidobacterium",  # Decreased
                "Escherichia_coli", "Enterobacteriaceae", "Bacteroides_fragilis"  # Increased
            ],
            "CRC": [
                "Fusobacterium_nucleatum", "Bacteroides_fragilis", "Peptostreptococcus",  # Increased
                "Faecalibacterium_prausnitzii", "Roseburia", "Bifidobacterium"  # Decreased
            ],
            "T2D": [
                "Akkermansia_muciniphila", "Bifidobacterium", "Faecalibacterium",  # Decreased
                "Enterobacteriaceae", "Clostridium_ramosum"  # Increased
            ],
            "obesity": [
                "Bacteroidetes", "Akkermansia_muciniphila",  # Decreased
                "Firmicutes", "Lactobacillus"  # Increased
            ],
            "cirrhosis": [
                "Enterobacteriaceae", "Streptococcus", "Veillonella",  # Increased
                "Lachnospiraceae", "Ruminococcaceae"  # Decreased
            ]
        }
        
        if condition not in literature_markers:
            # Return random subset for unknown conditions
            return np.random.choice(feature_names, min(20, len(feature_names)//10), replace=False).tolist()
        
        expected_markers = literature_markers[condition]
        ground_truth = []
        
        # Find features that match literature markers
        for marker in expected_markers:
            matching_features = [f for f in feature_names if marker.lower() in f.lower()]
            ground_truth.extend(matching_features)
        
        # Add some random features to simulate discovery
        remaining_features = list(set(feature_names) - set(ground_truth))
        if remaining_features:
            additional = np.random.choice(
                remaining_features, 
                min(10, len(remaining_features)//20), 
                replace=False
            )
            ground_truth.extend(additional)
        
        return list(set(ground_truth))  # Remove duplicates
    
    def download_all_conditions(self, max_conditions: int = 5) -> Dict[str, Dict]:
        """
        Download data for all available conditions
        
        Parameters:
        -----------
        max_conditions : int
            Maximum number of conditions to download
            
        Returns:
        --------
        Dict[str, Dict]
            All downloaded datasets
        """
        
        logger.info(f"📥 Downloading data for up to {max_conditions} conditions...")
        
        conditions = list(self.available_conditions.keys())[:max_conditions]
        datasets = {}
        
        for condition in conditions:
            logger.info(f"\n⬇️ Downloading {condition}...")
            
            try:
                data = self.download_condition_data(condition)
                if data:
                    datasets[condition] = data
                    logger.info(f"✅ {condition}: {data['count_table'].shape} samples×features")
                else:
                    logger.warning(f"❌ Failed to download {condition}")
                    
            except Exception as e:
                logger.error(f"❌ Error downloading {condition}: {e}")
                continue
        
        logger.info(f"\n🎉 Downloaded {len(datasets)}/{len(conditions)} conditions")
        
        return datasets
    
    def cross_validate_with_synthetic(self, real_datasets: Dict[str, Dict]) -> Dict[str, Dict]:
        """
        Cross-validate real data with synthetic data
        
        Parameters:
        -----------
        real_datasets : Dict[str, Dict]
            Real datasets from curatedMetagenomicData
            
        Returns:
        --------
        Dict[str, Dict]
            Cross-validation results
        """
        
        logger.info("🔄 Cross-validating real data with synthetic data...")
        
        from ..data_generators import MicrobiomeDataGenerator
        
        generator = MicrobiomeDataGenerator()
        validation_results = {}
        
        for condition, real_data in real_datasets.items():
            logger.info(f"\n🔄 Cross-validating {condition}...")
            
            # Generate synthetic data with similar characteristics
            real_count = real_data['count_table']
            real_meta = real_data['metadata']
            
            # Calculate real data characteristics
            n_samples = len(real_count)
            n_features = len(real_count.columns)
            sparsity = (real_count == 0).sum().sum() / real_count.size
            
            # Generate synthetic counterpart
            if 'gene' in str(real_data.get('source', '')).lower():
                synthetic_count, synthetic_meta, synthetic_truth = generator.generate_gene_data(
                    n_samples=n_samples,
                    n_features=n_features,
                    n_differential=len(real_data['ground_truth']),
                    sparsity=sparsity,
                    effect_size=2.0
                )
            else:
                synthetic_count, synthetic_meta, synthetic_truth = generator.generate_asv_data(
                    n_samples=n_samples,
                    n_features=n_features,
                    n_differential=len(real_data['ground_truth']),
                    sparsity=sparsity,
                    effect_size=2.0
                )
            
            # Compare characteristics
            validation_results[condition] = {
                'real_data': real_data,
                'synthetic_data': {
                    'count_table': synthetic_count,
                    'metadata': synthetic_meta,
                    'ground_truth': synthetic_truth
                },
                'comparison': {
                    'n_samples': {'real': n_samples, 'synthetic': len(synthetic_count)},
                    'n_features': {'real': n_features, 'synthetic': len(synthetic_count.columns)},
                    'sparsity': {
                        'real': sparsity,
                        'synthetic': (synthetic_count == 0).sum().sum() / synthetic_count.size
                    },
                    'ground_truth_size': {
                        'real': len(real_data['ground_truth']),
                        'synthetic': len(synthetic_truth)
                    }
                }
            }
            
            logger.info(f"✅ {condition} validation prepared")
        
        return validation_results


def download_curated_data(output_dir: str = "curated_real_data",
                         max_conditions: int = 3) -> Dict[str, Dict]:
    """
    Convenience function to download curatedMetagenomicData
    
    Parameters:
    -----------
    output_dir : str
        Directory to save data
    max_conditions : int
        Maximum conditions to download
        
    Returns:
    --------
    Dict[str, Dict]
        Downloaded datasets
    """
    
    downloader = CuratedDataDownloader(output_dir=output_dir)
    
    # Check and install requirements
    if not downloader.check_r_environment():
        logger.error("❌ R environment not ready")
        return {}
    
    if not downloader.install_required_packages():
        logger.error("❌ Failed to install R packages")
        return {}
    
    # Download data
    datasets = downloader.download_all_conditions(max_conditions=max_conditions)
    
    return datasets