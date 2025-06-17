#!/usr/bin/env python3
"""
Core DAAadvisor functionality
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
import logging

from .profiler import DataProfiler
from .selector import MethodSelector
from .methods import MethodRegistry

logger = logging.getLogger(__name__)


class DifferentialAbundanceTool:
    """Main tool orchestrating the differential abundance analysis"""
    
    def __init__(self):
        self.profiler = DataProfiler()
        self.selector = MethodSelector()
        self.method_registry = MethodRegistry()
        self.results = {}
    
    def analyze(self, 
                count_table: pd.DataFrame, 
                metadata: pd.DataFrame,
                data_type: Optional[str] = None,
                use_consensus: bool = True,
                **kwargs) -> Dict:
        """
        Run differential abundance analysis with automatic method selection
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame  
            Sample metadata
        data_type : str, optional
            One of 'viral', 'asv', 'gene'
        use_consensus : bool
            Whether to run multiple methods for consensus
        
        Returns:
        --------
        dict : Results dictionary with recommendations and analyses
        """
        
        # Step 1: Profile the data
        logger.info("Profiling input data...")
        profile = self.profiler.profile_data(count_table, metadata, data_type)
        
        # Step 2: Get method recommendations
        logger.info("Selecting optimal methods...")
        recommendations = self.selector.recommend_methods(profile)
        
        logger.info(f"Recommended primary method: {recommendations.primary_method}")
        logger.info(f"Confidence: {recommendations.confidence:.2f}")
        
        # Step 3: Run analyses
        results = {
            'profile': profile,
            'recommendations': recommendations,
            'analyses': {}
        }
        
        # Run primary method
        primary_method = recommendations.primary_method
        if self.method_registry.has_method(primary_method):
            logger.info(f"Running {primary_method}...")
            method = self.method_registry.get_method(primary_method)
            results['analyses'][primary_method] = method.run(count_table, metadata, **kwargs)
        else:
            logger.warning(f"Method {primary_method} not available, using fallback")
            fallback_method = self.method_registry.get_fallback_method()
            results['analyses'][fallback_method.name()] = fallback_method.run(count_table, metadata, **kwargs)
        
        # Run secondary methods if consensus requested
        if use_consensus:
            for method_name in recommendations.secondary_methods:
                if self.method_registry.has_method(method_name):
                    logger.info(f"Running {method_name} for consensus...")
                    method = self.method_registry.get_method(method_name)
                    results['analyses'][method_name] = method.run(count_table, metadata, **kwargs)
        
        # Generate consensus results if multiple methods were run
        if len(results['analyses']) > 1:
            logger.info("Generating consensus results...")
            results['consensus'] = self._generate_consensus(results['analyses'])
        
        self.results = results
        return results
    
    def _generate_consensus(self, analyses: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Generate consensus results from multiple methods"""
        
        # Simple voting-based consensus
        all_features = set()
        for method_results in analyses.values():
            all_features.update(method_results['feature'].values)
        
        consensus_data = []
        for feature in all_features:
            feature_votes = 0
            feature_pvalues = []
            feature_log2fcs = []
            
            for method, results in analyses.items():
                feature_data = results[results['feature'] == feature]
                if not feature_data.empty:
                    padj = feature_data.get('padj', feature_data.get('qvalue', [np.nan])).iloc[0]
                    if padj < 0.05:
                        feature_votes += 1
                        feature_pvalues.append(feature_data['pvalue'].iloc[0])
                        if 'log2fc' in feature_data.columns:
                            feature_log2fcs.append(feature_data['log2fc'].iloc[0])
            
            consensus_data.append({
                'feature': feature,
                'n_significant': feature_votes,
                'n_methods': len(analyses),
                'mean_pvalue': np.mean(feature_pvalues) if feature_pvalues else np.nan,
                'mean_log2fc': np.mean(feature_log2fcs) if feature_log2fcs else np.nan,
                'consensus_significant': feature_votes >= len(analyses) / 2
            })
        
        return pd.DataFrame(consensus_data).sort_values('n_significant', ascending=False)
    
    def summarize_results(self) -> None:
        """Print summary of results"""
        if not self.results:
            logger.warning("No results to summarize. Run analyze() first.")
            return
        
        print("\n" + "="*60)
        print("DIFFERENTIAL ABUNDANCE ANALYSIS SUMMARY")
        print("="*60)
        
        profile = self.results['profile']
        print(f"\nData Profile:")
        print(f"  - Type: {profile.data_type}")
        print(f"  - Samples: {profile.n_samples}")
        print(f"  - Features: {profile.n_features}")
        print(f"  - Sparsity: {profile.sparsity:.1%}")
        print(f"  - Zero inflation: {profile.zero_inflation:.1%}")
        
        recommendations = self.results['recommendations']
        print(f"\nMethod Recommendations:")
        print(f"  - Primary: {recommendations.primary_method} (confidence: {recommendations.confidence:.2f})")
        print(f"  - Alternatives: {', '.join(recommendations.secondary_methods)}")
        
        print(f"\nAnalyses Run:")
        for method, results in self.results['analyses'].items():
            padj_col = 'padj' if 'padj' in results.columns else 'qvalue' if 'qvalue' in results.columns else None
            if padj_col:
                n_sig = (results[padj_col] < 0.05).sum()
            else:
                n_sig = (results['pvalue'] < 0.05).sum() if 'pvalue' in results.columns else 0
            print(f"  - {method}: {n_sig} significant features (FDR < 0.05)")
        
        if 'consensus' in self.results:
            consensus = self.results['consensus']
            n_consensus_sig = consensus['consensus_significant'].sum()
            print(f"\nConsensus Results:")
            print(f"  - {n_consensus_sig} features significant in majority of methods")
    
    def get_significant_features(self, alpha: float = 0.05, method: str = None) -> pd.DataFrame:
        """Get significant features from analysis results"""
        if not self.results:
            raise ValueError("No results available. Run analyze() first.")
        
        if method and method in self.results['analyses']:
            results = self.results['analyses'][method]
            padj_col = 'padj' if 'padj' in results.columns else 'qvalue' if 'qvalue' in results.columns else 'pvalue'
            return results[results[padj_col] < alpha]
        elif 'consensus' in self.results:
            return self.results['consensus'][self.results['consensus']['consensus_significant']]
        else:
            # Return results from primary method
            primary_method = list(self.results['analyses'].keys())[0]
            results = self.results['analyses'][primary_method]
            padj_col = 'padj' if 'padj' in results.columns else 'qvalue' if 'qvalue' in results.columns else 'pvalue'
            return results[results[padj_col] < alpha]