#!/usr/bin/env python3
"""
Wilcoxon rank-sum test implementation for DAAadvisor
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Optional
from scipy import stats
import logging

from .base import DAAMethod

logger = logging.getLogger(__name__)


class WilcoxonMethod(DAAMethod):
    """Wilcoxon rank-sum test for differential abundance"""
    
    def name(self) -> str:
        return "wilcoxon"
    
    def run(self, 
            count_table: pd.DataFrame, 
            metadata: pd.DataFrame, 
            formula: Optional[str] = None,
            group_column: Optional[str] = None,
            transform: str = 'none',
            **kwargs) -> pd.DataFrame:
        """
        Run Wilcoxon rank-sum test
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
        formula : str, optional
            Not used for Wilcoxon (simple two-group comparison only)
        group_column : str, optional
            Column name for grouping variable. If None, uses first column
        transform : str
            Data transformation: 'none', 'log', 'sqrt', 'clr'
        """
        
        self.validate_input(count_table, metadata)
        
        # Determine grouping column
        if group_column is None:
            group_column = metadata.columns[0]
            logger.info(f"Using '{group_column}' as grouping variable")
        
        if group_column not in metadata.columns:
            raise ValueError(f"Group column '{group_column}' not found in metadata")
        
        # Get groups
        groups = metadata[group_column].unique()
        groups = groups[~pd.isna(groups)]  # Remove NaN values
        
        if len(groups) != 2:
            raise ValueError(f"Wilcoxon test requires exactly 2 groups, found {len(groups)}: {groups}")
        
        # Apply transformation
        transformed_data = self._apply_transform(count_table, transform)
        
        # Get sample indices for each group
        group1_samples = metadata[metadata[group_column] == groups[0]].index
        group2_samples = metadata[metadata[group_column] == groups[1]].index
        
        # Ensure samples exist in count table
        group1_samples = group1_samples.intersection(transformed_data.index)
        group2_samples = group2_samples.intersection(transformed_data.index)
        
        if len(group1_samples) == 0 or len(group2_samples) == 0:
            raise ValueError("One or both groups have no samples in count table")
        
        logger.info(f"Comparing {len(group1_samples)} vs {len(group2_samples)} samples")
        logger.info(f"Groups: {groups[0]} vs {groups[1]}")
        
        results = []
        n_features = len(transformed_data.columns)
        
        for i, feature in enumerate(transformed_data.columns):
            if i % 100 == 0:
                logger.debug(f"Processing feature {i+1}/{n_features}")
            
            group1_values = transformed_data.loc[group1_samples, feature].values
            group2_values = transformed_data.loc[group2_samples, feature].values
            
            # Skip if all values are identical (no variation)
            if np.var(np.concatenate([group1_values, group2_values])) == 0:
                results.append({
                    'feature': feature,
                    'statistic': np.nan,
                    'pvalue': 1.0,
                    'mean_group1': np.mean(group1_values),
                    'mean_group2': np.mean(group2_values),
                    'log2fc': 0.0
                })
                continue
            
            # Perform Wilcoxon rank-sum test
            try:
                statistic, pvalue = stats.ranksums(group1_values, group2_values)
                
                # Handle perfect separation cases
                if np.isnan(pvalue):
                    pvalue = 1.0
                    statistic = 0.0
                
            except Exception as e:
                logger.warning(f"Error testing feature {feature}: {e}")
                statistic = np.nan
                pvalue = 1.0
            
            # Calculate effect size (log2 fold change)
            mean1 = np.mean(group1_values)
            mean2 = np.mean(group2_values)
            
            # For log2FC calculation, add pseudocount for zero values
            if transform == 'none':
                # Add pseudocount for zero values in original scale
                mean1_fc = mean1 + 1e-6
                mean2_fc = mean2 + 1e-6
                log2fc = np.log2(mean2_fc / mean1_fc) if mean1_fc > 0 else 0.0
            else:
                # For transformed data, calculate difference
                log2fc = mean2 - mean1
            
            results.append({
                'feature': feature,
                'statistic': statistic,
                'pvalue': pvalue,
                'mean_group1': mean1,
                'mean_group2': mean2,
                'log2fc': log2fc
            })
        
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        valid_pvalues = ~results_df['pvalue'].isna()
        
        if valid_pvalues.sum() > 0:
            _, padj_values, _, _ = multipletests(
                results_df.loc[valid_pvalues, 'pvalue'], 
                method='fdr_bh'
            )
            results_df.loc[valid_pvalues, 'padj'] = padj_values
            results_df.loc[~valid_pvalues, 'padj'] = np.nan
        else:
            results_df['padj'] = np.nan
        
        # Sort by p-value
        results_df = results_df.sort_values('pvalue')
        
        return self.standardize_output(results_df)
    
    def _apply_transform(self, count_table: pd.DataFrame, transform: str) -> pd.DataFrame:
        """Apply data transformation"""
        
        if transform == 'none':
            return count_table.copy()
        
        elif transform == 'log':
            # Log(x + 1) transformation
            return np.log(count_table + 1)
        
        elif transform == 'sqrt':
            # Square root transformation
            return np.sqrt(count_table)
        
        elif transform == 'clr':
            # Centered log-ratio transformation
            # Add pseudocount and apply CLR
            data_pseudo = count_table + 1e-6
            geometric_mean = np.exp(np.log(data_pseudo).mean(axis=1))
            clr_data = np.log(data_pseudo.div(geometric_mean, axis=0))
            return clr_data
        
        else:
            raise ValueError(f"Unknown transformation: {transform}")
    
    def check_requirements(self, count_table: pd.DataFrame, metadata: pd.DataFrame) -> bool:
        """Check if data meets Wilcoxon test requirements"""
        # Very minimal requirements - just need at least 3 samples per group
        if metadata.shape[1] == 0:
            return False
        
        group_sizes = metadata.iloc[:, 0].value_counts()
        return len(group_sizes) == 2 and all(size >= 3 for size in group_sizes)
    
    def get_parameters(self) -> Dict[str, Any]:
        """Get method parameters"""
        return {
            'group_column': {
                'type': str,
                'default': None,
                'description': 'Column name for grouping variable'
            },
            'transform': {
                'type': str,
                'default': 'none',
                'choices': ['none', 'log', 'sqrt', 'clr'],
                'description': 'Data transformation to apply'
            }
        }
    
    def cite(self) -> str:
        """Return citation information"""
        return ("Wilcoxon, F. (1945). Individual comparisons by ranking methods. "
                "Biometrics Bulletin, 1(6), 80-83.")