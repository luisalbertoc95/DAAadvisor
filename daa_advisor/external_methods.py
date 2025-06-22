#!/usr/bin/env python3
"""
External Method Comparison Framework

This module implements comparison with external differential abundance tools
for comprehensive benchmarking against the state-of-the-art.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import logging
import subprocess
import tempfile
from pathlib import Path
import warnings
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)

class ExternalMethod(ABC):
    """Abstract base class for external methods"""
    
    def __init__(self, name: str):
        self.name = name
    
    @abstractmethod
    def is_available(self) -> bool:
        """Check if method is available/installed"""
        pass
    
    @abstractmethod
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """Run the external method"""
        pass

class LEfSeMethod(ExternalMethod):
    """LEfSe (Linear discriminant analysis Effect Size) method"""
    
    def __init__(self):
        super().__init__("LEfSe")
    
    def is_available(self) -> bool:
        """Check if LEfSe is available"""
        try:
            result = subprocess.run(['lefse-format_input.py', '--help'], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except:
            return False
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """Run LEfSe analysis"""
        
        if not self.is_available():
            raise RuntimeError("LEfSe not available")
        
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Create LEfSe input format
            input_file = tmpdir / "input.txt"
            formatted_file = tmpdir / "formatted.in"
            results_file = tmpdir / "results.res"
            
            # Write input file
            with open(input_file, 'w') as f:
                # Header with class information
                f.write("class\t" + "\t".join(metadata.index) + "\n")
                f.write("condition\t" + "\t".join(metadata['condition']) + "\n")
                
                # Feature data
                for feature in count_table.columns:
                    values = count_table[feature].values
                    f.write(f"{feature}\t" + "\t".join(map(str, values)) + "\n")
            
            try:
                # Format input
                subprocess.run([
                    'lefse-format_input.py', str(input_file), str(formatted_file),
                    '-c', '1', '-s', '2', '-u', '3', '-o', '1000000'
                ], check=True, capture_output=True)
                
                # Run LEfSe
                subprocess.run([
                    'lefse-run_lefse.py', str(formatted_file), str(results_file)
                ], check=True, capture_output=True)
                
                # Parse results
                results = []
                if results_file.exists():
                    with open(results_file, 'r') as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 4:
                                results.append({
                                    'feature': parts[0],
                                    'log_max_mean': float(parts[1]) if parts[1] != '-' else 0,
                                    'lda_score': float(parts[2]) if parts[2] != '-' else 0,
                                    'pvalue': float(parts[3]) if parts[3] != '-' else 1.0
                                })
                
                results_df = pd.DataFrame(results)
                if not results_df.empty:
                    results_df['padj'] = results_df['pvalue']  # LEfSe doesn't do FDR correction
                    results_df['log2fc'] = results_df['log_max_mean']
                
                return results_df
                
            except subprocess.CalledProcessError as e:
                logger.error(f"LEfSe failed: {e}")
                return pd.DataFrame()

class MaAsLin2Method(ExternalMethod):
    """MaAsLin2 method using R interface"""
    
    def __init__(self):
        super().__init__("MaAsLin2")
    
    def is_available(self) -> bool:
        """Check if MaAsLin2 is available"""
        try:
            import rpy2.robjects as ro
            from rpy2.robjects.packages import importr
            
            # Check if MaAsLin2 is installed
            utils = importr('utils')
            installed_packages = ro.r('installed.packages()')
            return 'Maaslin2' in list(installed_packages.rx(True, 0))
        except:
            return False
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """Run MaAsLin2 analysis"""
        
        if not self.is_available():
            raise RuntimeError("MaAsLin2 not available")
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            
            pandas2ri.activate()
            
            # Import MaAsLin2
            maaslin2 = importr('Maaslin2')
            
            with tempfile.TemporaryDirectory() as tmpdir:
                output_dir = Path(tmpdir) / "maaslin2_output"
                
                # Prepare data
                r_count_table = pandas2ri.py2rpy(count_table.T)  # Features x Samples
                r_metadata = pandas2ri.py2rpy(metadata)
                
                # Run MaAsLin2
                results = maaslin2.Maaslin2(
                    input_data=r_count_table,
                    input_metadata=r_metadata,
                    output=str(output_dir),
                    fixed_effects='condition',
                    normalization='TSS',  # Total Sum Scaling
                    transform='LOG',
                    analysis_method='LM',  # Linear Model
                    plot_heatmap=False,
                    plot_scatter=False
                )
                
                # Read results
                results_file = output_dir / "significant_results.tsv"
                if results_file.exists():
                    results_df = pd.read_csv(results_file, sep='\t')
                    
                    # Standardize column names
                    results_df = results_df.rename(columns={
                        'feature': 'feature',
                        'coef': 'log2fc',
                        'pval': 'pvalue',
                        'qval': 'padj'
                    })
                    
                    return results_df[['feature', 'log2fc', 'pvalue', 'padj']]
                else:
                    return pd.DataFrame()
                    
        except Exception as e:
            logger.error(f"MaAsLin2 failed: {e}")
            return pd.DataFrame()

class CornCobMethod(ExternalMethod):
    """corncob method using R interface"""
    
    def __init__(self):
        super().__init__("corncob")
    
    def is_available(self) -> bool:
        """Check if corncob is available"""
        try:
            import rpy2.robjects as ro
            from rpy2.robjects.packages import importr
            
            utils = importr('utils')
            installed_packages = ro.r('installed.packages()')
            return 'corncob' in list(installed_packages.rx(True, 0))
        except:
            return False
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, **kwargs) -> pd.DataFrame:
        """Run corncob analysis"""
        
        if not self.is_available():
            raise RuntimeError("corncob not available")
        
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            
            pandas2ri.activate()
            
            # Import required packages
            corncob = importr('corncob')
            phyloseq = importr('phyloseq')
            
            # Create phyloseq object
            r_otu_table = pandas2ri.py2rpy(count_table.T)  # Features x Samples
            r_sample_data = pandas2ri.py2rpy(metadata)
            
            # Run R code to create phyloseq and run corncob
            ro.r('''
            library(corncob)
            library(phyloseq)
            
            run_corncob <- function(otu_table, sample_data) {
                # Create phyloseq object
                ps <- phyloseq(otu_table(otu_table, taxa_are_rows = TRUE),
                              sample_data(sample_data))
                
                # Run corncob
                results <- list()
                taxa_names <- taxa_names(ps)
                
                for (i in 1:min(length(taxa_names), 50)) {  # Limit to first 50 taxa for speed
                    tryCatch({
                        taxon <- taxa_names[i]
                        result <- bbdml(formula = reformulate("condition"),
                                       phi.formula = ~ 1,
                                       data = ps,
                                       taxa = taxon)
                        
                        # Extract results
                        coef_summary <- summary(result)$coefficients
                        if (nrow(coef_summary) > 1) {
                            results[[taxon]] <- data.frame(
                                feature = taxon,
                                log2fc = coef_summary[2, 1],
                                pvalue = coef_summary[2, 4],
                                stringsAsFactors = FALSE
                            )
                        }
                    }, error = function(e) {
                        # Skip problematic taxa
                        NULL
                    })
                }
                
                # Combine results
                if (length(results) > 0) {
                    final_results <- do.call(rbind, results)
                    final_results$padj <- p.adjust(final_results$pvalue, method = "fdr")
                    return(final_results)
                } else {
                    return(data.frame())
                }
            }
            ''')
            
            # Run analysis
            r_results = ro.r['run_corncob'](r_otu_table, r_sample_data)
            
            # Convert back to pandas
            if r_results.nrow[0] > 0:
                results_df = pandas2ri.rpy2py(r_results)
                return results_df[['feature', 'log2fc', 'pvalue', 'padj']]
            else:
                return pd.DataFrame()
                
        except Exception as e:
            logger.error(f"corncob failed: {e}")
            return pd.DataFrame()

class ExternalMethodBenchmark:
    """
    Comprehensive benchmarking against external methods
    """
    
    def __init__(self):
        self.external_methods = {
            'LEfSe': LEfSeMethod(),
            'MaAsLin2': MaAsLin2Method(),
            'corncob': CornCobMethod()
        }
        
        # Check availability
        self.available_methods = {}
        for name, method in self.external_methods.items():
            if method.is_available():
                self.available_methods[name] = method
                logger.info(f"✅ {name} available")
            else:
                logger.warning(f"❌ {name} not available")
    
    def run_external_comparison(self, 
                              count_table: pd.DataFrame,
                              metadata: pd.DataFrame,
                              ground_truth: List[str]) -> Dict[str, Dict]:
        """
        Run comparison with all available external methods
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Count data
        metadata : pd.DataFrame
            Sample metadata
        ground_truth : List[str]
            True differential features
            
        Returns:
        --------
        Dict[str, Dict]
            Results for each external method
        """
        
        logger.info(f"🔄 Running external method comparison ({len(self.available_methods)} methods)")
        
        results = {}
        
        for method_name, method in self.available_methods.items():
            logger.info(f"  Running {method_name}...")
            
            try:
                # Run method
                method_results = method.run(count_table, metadata)
                
                if not method_results.empty:
                    # Calculate metrics
                    metrics = self._calculate_method_metrics(
                        method_results, ground_truth, list(count_table.columns)
                    )
                    
                    results[method_name] = {
                        'results': method_results,
                        'metrics': metrics
                    }
                    
                    logger.info(f"    ✅ {method_name}: {metrics.get('f1_score', 0):.3f} F1 score")
                else:
                    logger.warning(f"    ❌ {method_name}: No results returned")
                    
            except Exception as e:
                logger.error(f"    ❌ {method_name} failed: {e}")
                continue
        
        return results
    
    def _calculate_method_metrics(self, 
                                results: pd.DataFrame,
                                ground_truth: List[str],
                                all_features: List[str],
                                alpha: float = 0.05) -> Dict:
        """Calculate performance metrics for external method"""
        
        try:
            # Get significant features
            padj_col = None
            for col in ['padj', 'qvalue', 'adj_pvalue']:
                if col in results.columns:
                    padj_col = col
                    break
            
            if padj_col is None:
                # Use raw p-values if no adjusted available
                padj_col = 'pvalue'
                alpha = 0.01  # More stringent threshold
            
            significant_features = set(results[results[padj_col] < alpha]['feature'].values)
            
            # Calculate metrics
            tp = len(significant_features.intersection(set(ground_truth)))
            fp = len(significant_features - set(ground_truth))
            fn = len(set(ground_truth) - significant_features)
            tn = len(set(all_features) - set(ground_truth) - significant_features)
            
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
            fdr = fp / (tp + fp) if (tp + fp) > 0 else 0
            
            return {
                'sensitivity': sensitivity,
                'specificity': specificity,
                'precision': precision,
                'f1_score': f1_score,
                'fdr': fdr,
                'n_significant': len(significant_features),
                'n_true_positives': tp,
                'n_false_positives': fp
            }
            
        except Exception as e:
            logger.error(f"Error calculating metrics: {e}")
            return {}

class ComputationalEfficiencyBenchmark:
    """
    Benchmark computational efficiency and resource usage
    """
    
    def __init__(self):
        self.results = {}
    
    def benchmark_scalability(self,
                            sample_sizes: List[int] = [50, 100, 200, 500],
                            feature_sizes: List[int] = [100, 500, 1000, 2000]) -> Dict:
        """
        Benchmark scalability across different data sizes
        
        Parameters:
        -----------
        sample_sizes : List[int]
            Sample sizes to test
        feature_sizes : List[int]
            Feature counts to test
            
        Returns:
        --------
        Dict
            Scalability results
        """
        
        logger.info("⏱️ Running computational efficiency benchmark...")
        
        from .data_generators import MicrobiomeDataGenerator
        from .core import DifferentialAbundanceTool
        import time
        import psutil
        import os
        
        generator = MicrobiomeDataGenerator()
        tool = DifferentialAbundanceTool()
        
        results = []
        
        for n_samples in sample_sizes:
            for n_features in feature_sizes:
                logger.info(f"  Testing {n_samples} samples × {n_features} features...")
                
                try:
                    # Generate test data
                    count_table, metadata, _ = generator.generate_asv_data(
                        n_samples=n_samples,
                        n_features=n_features,
                        n_differential=min(50, n_features // 10)
                    )
                    
                    # Monitor resource usage
                    process = psutil.Process(os.getpid())
                    start_memory = process.memory_info().rss / 1024 / 1024  # MB
                    
                    # Run analysis
                    start_time = time.time()
                    
                    analysis_results = tool.analyze(
                        count_table=count_table,
                        metadata=metadata,
                        use_consensus=True
                    )
                    
                    end_time = time.time()
                    end_memory = process.memory_info().rss / 1024 / 1024  # MB
                    
                    # Calculate metrics
                    runtime = end_time - start_time
                    memory_usage = end_memory - start_memory
                    n_methods_success = len(analysis_results['analyses'])
                    
                    results.append({
                        'n_samples': n_samples,
                        'n_features': n_features,
                        'runtime_seconds': runtime,
                        'memory_usage_mb': memory_usage,
                        'n_methods_success': n_methods_success,
                        'samples_per_second': n_samples / runtime,
                        'features_per_second': n_features / runtime
                    })
                    
                    logger.info(f"    Runtime: {runtime:.2f}s, Memory: {memory_usage:.1f}MB, Methods: {n_methods_success}")
                    
                except Exception as e:
                    logger.error(f"    Failed: {e}")
                    continue
        
        self.results['scalability'] = pd.DataFrame(results)
        return self.results


def run_comprehensive_external_benchmark(count_table: pd.DataFrame,
                                        metadata: pd.DataFrame,
                                        ground_truth: List[str]) -> Dict:
    """
    Run comprehensive benchmark including external methods
    
    Parameters:
    -----------
    count_table : pd.DataFrame
        Count data
    metadata : pd.DataFrame
        Sample metadata  
    ground_truth : List[str]
        True differential features
        
    Returns:
    --------
    Dict
        Complete benchmark results
    """
    
    logger.info("🏆 Running comprehensive external method benchmark...")
    
    # Initialize benchmarks
    external_benchmark = ExternalMethodBenchmark()
    efficiency_benchmark = ComputationalEfficiencyBenchmark()
    
    # Run external method comparison
    external_results = external_benchmark.run_external_comparison(
        count_table, metadata, ground_truth
    )
    
    # Run efficiency benchmark
    efficiency_results = efficiency_benchmark.benchmark_scalability()
    
    # Run DAAadvisor for comparison
    from .core import DifferentialAbundanceTool
    
    tool = DifferentialAbundanceTool()
    daa_results = tool.analyze(count_table, metadata, use_consensus=True)
    
    # Calculate DAAadvisor metrics
    from .publication_benchmark import PublicationBenchmark
    
    benchmark = PublicationBenchmark()
    daa_metrics = benchmark.calculate_comprehensive_metrics(
        daa_results['analyses'],
        ground_truth,
        list(count_table.columns)
    )
    
    # Combine results
    combined_results = {
        'daaadvisor': {method: {'metrics': metrics} for method, metrics in daa_metrics.items()},
        'external_methods': external_results,
        'efficiency': efficiency_results
    }
    
    logger.info("✅ Comprehensive external benchmark completed!")
    
    return combined_results