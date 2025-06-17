#!/usr/bin/env python3
"""
Template and examples for adding new differential abundance methods
"""

import pandas as pd
import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, Optional, Any
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import logging

logger = logging.getLogger(__name__)


class DAAMethodWrapper(ABC):
    """Abstract base class for all method wrappers"""
    
    @abstractmethod
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """Run differential abundance analysis"""
        pass
    
    @abstractmethod
    def name(self) -> str:
        """Return method name"""
        pass
    
    @abstractmethod
    def check_requirements(self) -> bool:
        """Check if method requirements are met"""
        pass
    
    def validate_input(self, count_table: pd.DataFrame, metadata: pd.DataFrame) -> None:
        """Common input validation"""
        if count_table.shape[0] != metadata.shape[0]:
            raise ValueError("Count table and metadata must have same number of samples")
        
        if not all(count_table.index == metadata.index):
            raise ValueError("Count table and metadata must have matching sample IDs")
        
        if count_table.isnull().any().any():
            raise ValueError("Count table contains null values")


class ALDEx2Wrapper(DAAMethodWrapper):
    """Wrapper for ALDEx2 method"""
    
    def __init__(self):
        self.r_aldex2 = None
        self._setup()
    
    def name(self) -> str:
        return "aldex2"
    
    def _setup(self):
        """Initialize R environment and load ALDEx2"""
        try:
            # Enable automatic pandas conversion
            pandas2ri.activate()
            
            # Import ALDEx2
            self.r_aldex2 = importr('ALDEx2')
            logger.info("ALDEx2 loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load ALDEx2: {e}")
            self.r_aldex2 = None
    
    def check_requirements(self) -> bool:
        """Check if ALDEx2 is available"""
        return self.r_aldex2 is not None
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """
        Run ALDEx2 analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Features x Samples count matrix (note: transposed)
        metadata : pd.DataFrame
            Sample metadata
        formula : str
            Not used for ALDEx2 (uses first metadata column)
        **kwargs : dict
            Additional parameters:
            - mc_samples: Number of Monte Carlo samples (default: 128)
            - test: Test type ('t', 'glm') (default: 't')
            - denom: Denominator ('all', 'iqlr', 'zero', 'lvha') (default: 'all')
        """
        
        self.validate_input(count_table.T, metadata)  # Transpose for validation
        
        if not self.check_requirements():
            raise RuntimeError("ALDEx2 not available. Please install in R.")
        
        # ALDEx2 expects features as rows
        if count_table.shape[0] == metadata.shape[0]:
            count_table = count_table.T
            logger.info("Transposed count table for ALDEx2")
        
        # Get parameters
        mc_samples = kwargs.get('mc_samples', 128)
        test = kwargs.get('test', 't')
        denom = kwargs.get('denom', 'all')
        
        # Get conditions from first metadata column
        conditions = metadata.iloc[:, 0].values
        
        try:
            # Convert to R objects
            r_counts = pandas2ri.py2rpy(count_table)
            r_conditions = ro.FactorVector(conditions)
            
            # Run ALDEx2
            logger.info(f"Running ALDEx2 with {mc_samples} MC samples...")
            
            # CLR transformation
            aldex_clr = self.r_aldex2.aldex_clr(
                r_counts, 
                r_conditions,
                mc_samples=mc_samples,
                denom=denom
            )
            
            # Statistical test
            if test == 't':
                aldex_res = self.r_aldex2.aldex_ttest(aldex_clr)
            else:
                aldex_res = self.r_aldex2.aldex_glm(aldex_clr, r_conditions)
            
            # Convert results back to pandas
            results = pandas2ri.rpy2py(aldex_res)
            results.index = count_table.index
            
            # Standardize output format
            output = pd.DataFrame({
                'feature': results.index,
                'pvalue': results['we.ep'] if 'we.ep' in results else results['glm.ep'],
                'padj': results['we.eBH'] if 'we.eBH' in results else results['glm.eBH'],
                'log2fc': results['diff.btw'] if 'diff.btw' in results else 0,
                'statistic': results['wi.ep'] if 'wi.ep' in results else 0
            })
            
            return output.sort_values('pvalue')
            
        except Exception as e:
            logger.error(f"ALDEx2 analysis failed: {e}")
            raise


class ANCOMBCWrapper(DAAMethodWrapper):
    """Wrapper for ANCOM-BC method"""
    
    def __init__(self):
        self.r_ancombc = None
        self._setup()
    
    def name(self) -> str:
        return "ancombc"
    
    def _setup(self):
        """Initialize R environment and load ANCOMBC"""
        try:
            pandas2ri.activate()
            self.r_ancombc = importr('ANCOMBC')
            logger.info("ANCOMBC loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load ANCOMBC: {e}")
            self.r_ancombc = None
    
    def check_requirements(self) -> bool:
        return self.r_ancombc is not None
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """Run ANCOM-BC analysis"""
        
        self.validate_input(count_table, metadata)
        
        if not self.check_requirements():
            raise RuntimeError("ANCOMBC not available. Please install in R.")
        
        # ANCOM-BC specific implementation
        # ... (similar pattern to ALDEx2)
        
        # Placeholder for now
        logger.warning("ANCOMBC wrapper not fully implemented yet")
        return pd.DataFrame()


class DESeq2Wrapper(DAAMethodWrapper):
    """Wrapper for DESeq2 method"""
    
    def __init__(self):
        self.r_deseq2 = None
        self._setup()
    
    def name(self) -> str:
        return "deseq2"
    
    def _setup(self):
        """Initialize R environment and load DESeq2"""
        try:
            pandas2ri.activate()
            self.r_deseq2 = importr('DESeq2')
            self.r_base = importr('base')
            self.r_stats = importr('stats')
            logger.info("DESeq2 loaded successfully")
        except Exception as e:
            logger.error(f"Failed to load DESeq2: {e}")
            self.r_deseq2 = None
    
    def check_requirements(self) -> bool:
        return self.r_deseq2 is not None
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """
        Run DESeq2 analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        formula : str
            Design formula (default: ~condition)
        **kwargs : dict
            Additional parameters:
            - test: Test type ('Wald', 'LRT') (default: 'Wald')
            - fitType: Fit type ('parametric', 'local', 'mean') (default: 'parametric')
        """
        
        self.validate_input(count_table, metadata)
        
        if not self.check_requirements():
            raise RuntimeError("DESeq2 not available. Please install in R.")
        
        # DESeq2 expects integer counts
        count_table = count_table.round().astype(int)
        
        # Default formula if not provided
        if formula is None:
            formula = f"~ {metadata.columns[0]}"
        
        # Get parameters
        test = kwargs.get('test', 'Wald')
        fitType = kwargs.get('fitType', 'parametric')
        
        try:
            # Convert to R objects
            r_counts = pandas2ri.py2rpy(count_table.T)  # Features as rows
            r_metadata = pandas2ri.py2rpy(metadata)
            
            # Create DESeqDataSet
            dds = self.r_deseq2.DESeqDataSetFromMatrix(
                countData=r_counts,
                colData=r_metadata,
                design=self.r_stats.formula(formula)
            )
            
            # Run DESeq2
            logger.info("Running DESeq2 analysis...")
            dds = self.r_deseq2.DESeq(dds, test=test, fitType=fitType, quiet=True)
            
            # Extract results
            res = self.r_deseq2.results(dds)
            results_df = pandas2ri.rpy2py(self.r_base.as_data_frame(res))
            results_df.index = count_table.columns
            
            # Standardize output
            output = pd.DataFrame({
                'feature': results_df.index,
                'pvalue': results_df['pvalue'],
                'padj': results_df['padj'],
                'log2fc': results_df['log2FoldChange'],
                'statistic': results_df['stat'],
                'baseMean': results_df['baseMean']
            })
            
            # Remove NA values
            output = output.dropna(subset=['pvalue'])
            
            return output.sort_values('pvalue')
            
        except Exception as e:
            logger.error(f"DESeq2 analysis failed: {e}")
            raise


class ZicoSeqWrapper(DAAMethodWrapper):
    """Wrapper for ZicoSeq method (Python implementation)"""
    
    def name(self) -> str:
        return "zicoseq"
    
    def check_requirements(self) -> bool:
        # Check if Python implementation is available
        try:
            # Would import actual ZicoSeq package here
            return True
        except:
            return False
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """
        Run ZicoSeq analysis
        
        This is a placeholder for the actual implementation
        ZicoSeq handles zero-inflation and compositional effects
        """
        
        self.validate_input(count_table, metadata)
        
        logger.warning("ZicoSeq wrapper not fully implemented yet")
        
        # Placeholder implementation
        # In reality, would call ZicoSeq functions here
        
        return pd.DataFrame({
            'feature': count_table.columns,
            'pvalue': np.random.uniform(0, 1, len(count_table.columns)),
            'padj': np.random.uniform(0, 1, len(count_table.columns)),
            'log2fc': np.random.normal(0, 1, len(count_table.columns))
        })


class CustomMethodWrapper(DAAMethodWrapper):
    """Template for adding your own custom method"""
    
    def __init__(self, method_name: str):
        self.method_name = method_name
        self.custom_function = None
    
    def name(self) -> str:
        return self.method_name
    
    def check_requirements(self) -> bool:
        """Check if your method's dependencies are available"""
        # Add your checks here
        return True
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """
        Implement your custom method here
        
        Required output format:
        - DataFrame with columns: feature, pvalue, padj, log2fc
        - Optional columns: statistic, any method-specific metrics
        """
        
        self.validate_input(count_table, metadata)
        
        # Your implementation here
        # Example structure:
        
        # 1. Prepare data
        # ...
        
        # 2. Run your statistical test
        # ...
        
        # 3. Format results
        results = pd.DataFrame({
            'feature': count_table.columns,
            'pvalue': [...],  # Your p-values
            'padj': [...],    # Adjusted p-values
            'log2fc': [...],  # Effect sizes
            # Add any additional columns
        })
        
        return results.sort_values('pvalue')


# Method registry for easy access
METHOD_REGISTRY = {
    'aldex2': ALDEx2Wrapper,
    'ancombc': ANCOMBCWrapper,
    'deseq2': DESeq2Wrapper,
    'zicoseq': ZicoSeqWrapper,
}


def register_method(name: str, wrapper_class: type) -> None:
    """Register a new method wrapper"""
    METHOD_REGISTRY[name] = wrapper_class
    logger.info(f"Registered new method: {name}")


def get_available_methods() -> Dict[str, bool]:
    """Check which methods are available"""
    available = {}
    for name, wrapper_class in METHOD_REGISTRY.items():
        try:
            wrapper = wrapper_class()
            available[name] = wrapper.check_requirements()
        except:
            available[name] = False
    return available


# Example of how to add a new method
if __name__ == "__main__":
    # Check available methods
    print("Available methods:")
    for method, available in get_available_methods().items():
        status = "✓" if available else "✗"
        print(f"  {status} {method}")
    
    # Example of registering a custom method
    class MyCustomMethod(DAAMethodWrapper):
        def name(self):
            return "my_method"
        
        def check_requirements(self):
            return True
        
        def run(self, count_table, metadata, formula=None, **kwargs):
            # Your implementation
            pass
    
    register_method("my_custom_method", MyCustomMethod)
