#!/usr/bin/env python3
"""
R-based method implementations for DAAadvisor
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Optional
import logging
import warnings

from .base import DAAMethod

logger = logging.getLogger(__name__)

# Suppress rpy2 warnings for cleaner output
warnings.filterwarnings("ignore", category=UserWarning, module="rpy2")

try:
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr
    from rpy2.rinterface_lib.embedded import RRuntimeError
    
    R_AVAILABLE = True
    
    # Initialize R interface
    r = robjects.r
    
except ImportError:
    R_AVAILABLE = False
    logger.warning("rpy2 not available, R methods will not be functional")


class RMethodBase(DAAMethod):
    """Base class for R-based methods"""
    
    def __init__(self):
        if not R_AVAILABLE:
            raise ImportError("rpy2 is required for R-based methods")
        self._check_r_packages()
    
    def _check_r_packages(self):
        """Check if required R packages are installed"""
        required_packages = self.get_required_r_packages()
        missing_packages = []
        
        for package in required_packages:
            try:
                importr(package)
                logger.debug(f"R package {package} is available")
            except RRuntimeError:
                missing_packages.append(package)
                logger.warning(f"R package {package} is not installed")
        
        if missing_packages:
            raise ImportError(f"Missing R packages: {missing_packages}")
    
    def get_required_r_packages(self) -> list:
        """Return list of required R packages"""
        return []
    
    def _prepare_data_for_r(self, count_table: pd.DataFrame, metadata: pd.DataFrame):
        """Prepare data for R analysis"""
        # Ensure data types are correct
        count_table = count_table.astype(int)
        
        # Remove samples with all zeros
        sample_sums = count_table.sum(axis=1)
        valid_samples = sample_sums > 0
        
        if not valid_samples.all():
            logger.warning(f"Removing {(~valid_samples).sum()} samples with zero counts")
            count_table = count_table[valid_samples]
            metadata = metadata.loc[count_table.index]
        
        # Remove features with all zeros
        feature_sums = count_table.sum(axis=0)
        valid_features = feature_sums > 0
        
        if not valid_features.all():
            logger.warning(f"Removing {(~valid_features).sum()} features with zero counts")
            count_table = count_table.loc[:, valid_features]
        
        return count_table, metadata
    
    def _convert_to_r(self, data):
        """Convert Python data to R objects"""
        with localconverter(robjects.default_converter + pandas2ri.converter):
            return robjects.conversion.py2rpy(data)
    
    def _convert_from_r(self, r_obj):
        """Convert R objects to Python data"""
        with localconverter(robjects.default_converter + pandas2ri.converter):
            return robjects.conversion.rpy2py(r_obj)


class ALDEx2Method(RMethodBase):
    """ALDEx2 implementation for differential abundance analysis"""
    
    def name(self) -> str:
        return "aldex2"
    
    def get_required_r_packages(self) -> list:
        return ["ALDEx2"]
    
    def run(self, 
            count_table: pd.DataFrame, 
            metadata: pd.DataFrame, 
            formula: Optional[str] = None,
            group_column: Optional[str] = None,
            mc_samples: int = 128,
            denom: str = "all",
            test: str = "t",
            **kwargs) -> pd.DataFrame:
        """
        Run ALDEx2 analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        group_column : str, optional
            Column name for grouping variable
        mc_samples : int
            Number of Monte Carlo samples
        denom : str
            Features to use as denominator ('all', 'iqlr', 'zero', etc.)
        test : str
            Statistical test ('t', 'kw', 'glm')
        """
        
        self.validate_input(count_table, metadata)
        count_table, metadata = self._prepare_data_for_r(count_table, metadata)
        
        # Import ALDEx2
        aldex2 = importr("ALDEx2")
        
        # Determine grouping column
        if group_column is None:
            group_column = metadata.columns[0]
            logger.info(f"Using '{group_column}' as grouping variable")
        
        if group_column not in metadata.columns:
            raise ValueError(f"Group column '{group_column}' not found in metadata")
        
        # Prepare conditions vector
        conditions = metadata[group_column].astype(str)
        unique_conditions = conditions.unique()
        
        if len(unique_conditions) != 2:
            raise ValueError(f"ALDEx2 requires exactly 2 groups, found {len(unique_conditions)}: {unique_conditions}")
        
        logger.info(f"Running ALDEx2 with {len(count_table)} samples and {len(count_table.columns)} features")
        logger.info(f"Groups: {unique_conditions[0]} ({(conditions == unique_conditions[0]).sum()}) vs {unique_conditions[1]} ({(conditions == unique_conditions[1]).sum()})")
        
        try:
            # Convert data to R
            r_counts = self._convert_to_r(count_table.T)  # ALDEx2 expects features x samples
            r_conditions = self._convert_to_r(conditions.values)
            
            # Run ALDEx2 CLR transformation
            logger.info(f"Performing CLR transformation with {mc_samples} Monte Carlo samples...")
            r_clr = aldex2.aldex_clr(r_counts, r_conditions, mc_samples=mc_samples, denom=denom)
            
            # Perform statistical testing
            logger.info(f"Performing {test} test...")
            if test == "t":
                r_results = aldex2.aldex_ttest(r_clr, paired_test=False, hist_plot=False)
            elif test == "kw":
                r_results = aldex2.aldex_kw(r_clr)
            elif test == "glm":
                r_results = aldex2.aldex_glm(r_clr, r_conditions)
            else:
                raise ValueError(f"Unknown test type: {test}")
            
            # Convert results back to Python
            results_df = self._convert_from_r(r_results)
            
            # Standardize column names
            if 'we.ep' in results_df.columns:
                results_df['pvalue'] = results_df['we.ep']
            elif 'kw.ep' in results_df.columns:
                results_df['pvalue'] = results_df['kw.ep']
            else:
                results_df['pvalue'] = results_df.iloc[:, -1]  # Last column is usually p-value
            
            if 'we.eBH' in results_df.columns:
                results_df['padj'] = results_df['we.eBH']
            elif 'kw.eBH' in results_df.columns:
                results_df['padj'] = results_df['kw.eBH']
            else:
                # Calculate FDR if not available
                from statsmodels.stats.multitest import multipletests
                _, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
                results_df['padj'] = padj
            
            if 'diff.btw' in results_df.columns:
                results_df['log2fc'] = results_df['diff.btw']
            elif 'diff.win' in results_df.columns:
                results_df['log2fc'] = results_df['diff.win']
            else:
                results_df['log2fc'] = 0.0
            
            if 'we.estat' in results_df.columns:
                results_df['statistic'] = results_df['we.estat']
            elif 'kw.estat' in results_df.columns:
                results_df['statistic'] = results_df['kw.estat']
            else:
                results_df['statistic'] = 0.0
            
            # Add feature names
            results_df['feature'] = count_table.columns
            results_df = results_df.reset_index(drop=True)
            
            # Sort by p-value
            results_df = results_df.sort_values('pvalue')
            
            logger.info(f"ALDEx2 analysis complete. Found {(results_df['padj'] < 0.05).sum()} significant features (FDR < 0.05)")
            
            return self.standardize_output(results_df)
            
        except Exception as e:
            logger.error(f"ALDEx2 analysis failed: {e}")
            raise RuntimeError(f"ALDEx2 analysis failed: {e}")
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get method parameters"""
        return {
            'group_column': {
                'type': str,
                'default': None,
                'description': 'Column name for grouping variable'
            },
            'mc_samples': {
                'type': int,
                'default': 128,
                'description': 'Number of Monte Carlo samples for CLR transformation'
            },
            'denom': {
                'type': str,
                'default': 'all',
                'choices': ['all', 'iqlr', 'zero'],
                'description': 'Features to use as denominator'
            },
            'test': {
                'type': str,
                'default': 't',
                'choices': ['t', 'kw', 'glm'],
                'description': 'Statistical test to perform'
            }
        }
    
    def cite(self) -> str:
        """Return citation information"""
        return ("Fernandes, A.D., et al. (2014). Unifying the analysis of high-throughput "
                "sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing "
                "and selective growth experiments by compositional data analysis. "
                "Microbiome, 2(1), 15.")


class ANCOMBCMethod(RMethodBase):
    """ANCOM-BC implementation for differential abundance analysis"""
    
    def name(self) -> str:
        return "ancom-bc"
    
    def get_required_r_packages(self) -> list:
        return ["ANCOMBC", "phyloseq", "TreeSummarizedExperiment"]
    
    def run(self, 
            count_table: pd.DataFrame, 
            metadata: pd.DataFrame, 
            formula: Optional[str] = None,
            group_column: Optional[str] = None,
            alpha: float = 0.05,
            **kwargs) -> pd.DataFrame:
        """
        Run ANCOM-BC analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        formula : str, optional
            Formula for analysis (e.g., "~ group + batch")
        group_column : str, optional
            Column name for grouping variable
        alpha : float
            Significance level
        """
        
        self.validate_input(count_table, metadata)
        count_table, metadata = self._prepare_data_for_r(count_table, metadata)
        
        # Import required packages
        ancombc = importr("ANCOMBC")
        tse = importr("TreeSummarizedExperiment")
        
        # Determine grouping column and formula
        if group_column is None:
            group_column = metadata.columns[0]
            logger.info(f"Using '{group_column}' as grouping variable")
        
        if formula is None:
            formula = f"~ {group_column}"
            logger.info(f"Using formula: {formula}")
        
        logger.info(f"Running ANCOM-BC with {len(count_table)} samples and {len(count_table.columns)} features")
        
        try:
            # Create TreeSummarizedExperiment object
            r_counts = self._convert_to_r(count_table.T.astype(int))  # Features x samples
            r_metadata = self._convert_to_r(metadata)
            
            # Create TSE object
            r("""
            create_tse <- function(counts, metadata) {
                library(TreeSummarizedExperiment)
                tse <- TreeSummarizedExperiment(
                    assays = list(counts = as.matrix(counts)),
                    colData = metadata
                )
                return(tse)
            }
            """)
            
            r_tse = r['create_tse'](r_counts, r_metadata)
            
            # Run ANCOM-BC
            logger.info("Running ANCOM-BC analysis...")
            r_results = ancombc.ancombc(
                phyloseq=r_tse,
                formula=formula,
                alpha=alpha,
                **kwargs
            )
            
            # Extract results
            r_res_table = r_results.rx2('res')
            results_df = self._convert_from_r(r_res_table)
            
            # Standardize output format
            if 'lfc' in results_df.columns:
                results_df['log2fc'] = results_df['lfc']
            if 'p_val' in results_df.columns:
                results_df['pvalue'] = results_df['p_val']
            if 'q_val' in results_df.columns:
                results_df['padj'] = results_df['q_val']
            if 'W' in results_df.columns:
                results_df['statistic'] = results_df['W']
            
            # Add feature names
            results_df['feature'] = count_table.columns
            results_df = results_df.reset_index(drop=True)
            
            # Sort by p-value
            results_df = results_df.sort_values('pvalue')
            
            logger.info(f"ANCOM-BC analysis complete. Found {(results_df['padj'] < 0.05).sum()} significant features")
            
            return self.standardize_output(results_df)
            
        except Exception as e:
            logger.error(f"ANCOM-BC analysis failed: {e}")
            raise RuntimeError(f"ANCOM-BC analysis failed: {e}")
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get method parameters"""
        return {
            'group_column': {
                'type': str,
                'default': None,
                'description': 'Column name for grouping variable'
            },
            'formula': {
                'type': str,
                'default': None,
                'description': 'Formula for analysis (e.g., "~ group + batch")'
            },
            'alpha': {
                'type': float,
                'default': 0.05,
                'description': 'Significance level'
            }
        }
    
    def cite(self) -> str:
        """Return citation information"""
        return ("Lin, H. & Peddada, S.D. (2020). Analysis of compositions of microbiomes "
                "with bias correction. Nature Communications, 11(1), 3514.")


class DESeq2Method(RMethodBase):
    """DESeq2 implementation for differential abundance analysis"""
    
    def name(self) -> str:
        return "deseq2"
    
    def get_required_r_packages(self) -> list:
        return ["DESeq2"]
    
    def run(self, 
            count_table: pd.DataFrame, 
            metadata: pd.DataFrame, 
            formula: Optional[str] = None,
            group_column: Optional[str] = None,
            test: str = "Wald",
            fit_type: str = "parametric",
            **kwargs) -> pd.DataFrame:
        """
        Run DESeq2 analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        formula : str, optional
            Design formula
        group_column : str, optional
            Column name for grouping variable
        test : str
            Test for differential expression ('Wald' or 'LRT')
        fit_type : str
            Type of fitting ('parametric', 'local', 'mean')
        """
        
        self.validate_input(count_table, metadata)
        count_table, metadata = self._prepare_data_for_r(count_table, metadata)
        
        # Import DESeq2
        deseq2 = importr("DESeq2")
        
        # Determine grouping column and formula
        if group_column is None:
            group_column = metadata.columns[0]
            logger.info(f"Using '{group_column}' as grouping variable")
        
        if formula is None:
            formula = f"~ {group_column}"
            logger.info(f"Using design formula: {formula}")
        
        logger.info(f"Running DESeq2 with {len(count_table)} samples and {len(count_table.columns)} features")
        
        try:
            # Convert data to R
            r_counts = self._convert_to_r(count_table.T.astype(int))  # Features x samples
            r_metadata = self._convert_to_r(metadata)
            
            # Create DESeq2 dataset
            logger.info("Creating DESeqDataSet...")
            r_dds = deseq2.DESeqDataSetFromMatrix(
                countData=r_counts,
                colData=r_metadata,
                design=robjects.Formula(formula)
            )
            
            # Run DESeq2 analysis
            logger.info(f"Running DESeq2 analysis with {test} test...")
            r_dds = deseq2.DESeq(r_dds, test=test, fitType=fit_type)
            
            # Get results
            r_results = deseq2.results(r_dds)
            results_df = self._convert_from_r(r_results)
            
            # Standardize output format
            if 'log2FoldChange' in results_df.columns:
                results_df['log2fc'] = results_df['log2FoldChange']
            if 'pvalue' not in results_df.columns and 'pval' in results_df.columns:
                results_df['pvalue'] = results_df['pval']
            if 'padj' not in results_df.columns and 'FDR' in results_df.columns:
                results_df['padj'] = results_df['FDR']
            if 'stat' in results_df.columns:
                results_df['statistic'] = results_df['stat']
            
            # Add feature names and handle missing values
            results_df['feature'] = count_table.columns
            results_df = results_df.reset_index(drop=True)
            
            # Remove rows with NaN p-values
            valid_results = ~results_df['pvalue'].isna()
            if not valid_results.all():
                logger.warning(f"Removing {(~valid_results).sum()} features with missing p-values")
                results_df = results_df[valid_results]
            
            # Sort by p-value
            results_df = results_df.sort_values('pvalue')
            
            logger.info(f"DESeq2 analysis complete. Found {(results_df['padj'] < 0.05).sum()} significant features")
            
            return self.standardize_output(results_df)
            
        except Exception as e:
            logger.error(f"DESeq2 analysis failed: {e}")
            raise RuntimeError(f"DESeq2 analysis failed: {e}")
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get method parameters"""
        return {
            'group_column': {
                'type': str,
                'default': None,
                'description': 'Column name for grouping variable'
            },
            'formula': {
                'type': str,
                'default': None,
                'description': 'Design formula (e.g., "~ group + batch")'
            },
            'test': {
                'type': str,
                'default': 'Wald',
                'choices': ['Wald', 'LRT'],
                'description': 'Test for differential expression'
            },
            'fit_type': {
                'type': str,
                'default': 'parametric',
                'choices': ['parametric', 'local', 'mean'],
                'description': 'Type of fitting'
            }
        }
    
    def cite(self) -> str:
        """Return citation information"""
        return ("Love, M.I., Huber, W. & Anders, S. (2014). Moderated estimation of "
                "fold change and dispersion for RNA-seq data with DESeq2. "
                "Genome Biology, 15(12), 550.")


# Function to check R package availability
def check_r_package_availability():
    """Check which R packages are available"""
    if not R_AVAILABLE:
        return {}
    
    packages_to_check = {
        'ALDEx2': ALDEx2Method,
        'ANCOMBC': ANCOMBCMethod, 
        'DESeq2': DESeq2Method
    }
    
    available_methods = {}
    for package, method_class in packages_to_check.items():
        try:
            importr(package)
            available_methods[method_class().name()] = method_class
            logger.info(f"R package {package} is available")
        except RRuntimeError:
            logger.warning(f"R package {package} is not available")
    
    return available_methods