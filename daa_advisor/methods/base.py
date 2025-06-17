#!/usr/bin/env python3
"""
Base class for differential abundance methods
"""

import pandas as pd
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional


class DAAMethod(ABC):
    """Abstract base class for differential abundance methods"""
    
    @abstractmethod
    def run(self, 
            count_table: pd.DataFrame, 
            metadata: pd.DataFrame, 
            formula: Optional[str] = None,
            **kwargs) -> pd.DataFrame:
        """
        Run differential abundance analysis
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        formula : str, optional
            Formula for complex designs (e.g., "~ condition + batch")
        **kwargs : dict
            Method-specific parameters
            
        Returns:
        --------
        pd.DataFrame : Results with standardized columns:
            - feature: feature identifier
            - pvalue: raw p-value
            - padj: adjusted p-value (FDR)
            - log2fc: log2 fold change (when applicable)
            - statistic: test statistic
        """
        pass
    
    @abstractmethod
    def name(self) -> str:
        """Return method name"""
        pass
    
    def check_requirements(self, count_table: pd.DataFrame, metadata: pd.DataFrame) -> bool:
        """Check if data meets method requirements"""
        return True
    
    def validate_input(self, count_table: pd.DataFrame, metadata: pd.DataFrame) -> None:
        """Validate input data format"""
        # Basic validation
        if count_table.empty:
            raise ValueError("Count table is empty")
        
        if metadata.empty:
            raise ValueError("Metadata is empty")
        
        if not set(count_table.index).issubset(set(metadata.index)):
            raise ValueError("Sample IDs in count table must match metadata index")
        
        if (count_table < 0).any().any():
            raise ValueError("Count table contains negative values")
    
    def standardize_output(self, results: pd.DataFrame) -> pd.DataFrame:
        """Standardize output format across methods"""
        required_cols = ['feature', 'pvalue']
        
        for col in required_cols:
            if col not in results.columns:
                raise ValueError(f"Method output must contain '{col}' column")
        
        # Ensure standard columns exist
        if 'padj' not in results.columns:
            results['padj'] = results['pvalue']  # Fallback if no adjustment
        
        if 'log2fc' not in results.columns:
            results['log2fc'] = 0.0  # Fallback for methods without effect size
        
        if 'statistic' not in results.columns:
            results['statistic'] = 0.0  # Fallback for methods without test statistic
        
        return results[['feature', 'pvalue', 'padj', 'log2fc', 'statistic']].copy()
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get method-specific parameters and their descriptions"""
        return {}
    
    def cite(self) -> str:
        """Return citation information for the method"""
        return "Please cite the original method publication"