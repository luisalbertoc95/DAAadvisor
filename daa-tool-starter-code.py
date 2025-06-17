#!/usr/bin/env python3
"""
Differential Abundance Analysis Tool
Main implementation file demonstrating core architecture
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from abc import ABC, abstractmethod
import warnings
from scipy import stats
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class DataProfile:
    """Data characteristics profile"""
    data_type: str  # 'viral', 'asv', 'gene'
    n_samples: int
    n_features: int
    sparsity: float
    zero_inflation: float
    sequencing_depth_cv: float  # coefficient of variation
    group_sizes: Dict[str, int]
    compositional_bias: float
    metadata_factors: List[str]


@dataclass
class MethodRecommendation:
    """Method recommendation with reasoning"""
    primary_method: str
    secondary_methods: List[str]
    confidence: float
    reasoning: List[str]


class DataProfiler:
    """Analyzes input data characteristics"""
    
    def __init__(self):
        self.profile = None
    
    def profile_data(self, 
                    count_table: pd.DataFrame, 
                    metadata: pd.DataFrame,
                    data_type: Optional[str] = None) -> DataProfile:
        """
        Profile the input data to determine characteristics
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata with grouping variables
        data_type : str, optional
            One of 'viral', 'asv', 'gene'. If None, will attempt to detect
        """
        
        if data_type is None:
            data_type = self._detect_data_type(count_table)
        
        profile = DataProfile(
            data_type=data_type,
            n_samples=count_table.shape[0],
            n_features=count_table.shape[1],
            sparsity=self._calculate_sparsity(count_table),
            zero_inflation=self._assess_zero_inflation(count_table),
            sequencing_depth_cv=self._calculate_depth_variation(count_table),
            group_sizes=self._get_group_sizes(metadata),
            compositional_bias=self._detect_compositional_bias(count_table),
            metadata_factors=list(metadata.columns)
        )
        
        self.profile = profile
        logger.info(f"Data profiled: {profile.n_samples} samples, {profile.n_features} features")
        logger.info(f"Sparsity: {profile.sparsity:.2f}, Zero inflation: {profile.zero_inflation:.2f}")
        
        return profile
    
    def _detect_data_type(self, count_table: pd.DataFrame) -> str:
        """Attempt to detect data type from feature names and characteristics"""
        feature_names = count_table.columns
        
        # Check for common patterns
        if any('ASV' in str(f) or 'OTU' in str(f) for f in feature_names):
            return 'asv'
        elif any('gene' in str(f).lower() or 'KO' in str(f) for f in feature_names):
            return 'gene'
        elif self._calculate_sparsity(count_table) > 0.9:
            return 'viral'  # Very high sparsity often indicates viral
        else:
            logger.warning("Could not detect data type, defaulting to 'asv'")
            return 'asv'
    
    def _calculate_sparsity(self, count_table: pd.DataFrame) -> float:
        """Calculate proportion of zeros in the data"""
        return (count_table == 0).sum().sum() / (count_table.shape[0] * count_table.shape[1])
    
    def _assess_zero_inflation(self, count_table: pd.DataFrame) -> float:
        """Assess degree of zero inflation beyond expected from count distribution"""
        # Simplified metric: ratio of features with >50% zeros
        feature_zero_props = (count_table == 0).mean(axis=0)
        return (feature_zero_props > 0.5).mean()
    
    def _calculate_depth_variation(self, count_table: pd.DataFrame) -> float:
        """Calculate coefficient of variation in sequencing depth"""
        depths = count_table.sum(axis=1)
        return depths.std() / depths.mean()
    
    def _get_group_sizes(self, metadata: pd.DataFrame) -> Dict[str, int]:
        """Get sample sizes for each group in primary factor"""
        # Assuming first column is primary grouping factor
        if metadata.shape[1] > 0:
            primary_factor = metadata.iloc[:, 0]
            return primary_factor.value_counts().to_dict()
        return {}
    
    def _detect_compositional_bias(self, count_table: pd.DataFrame) -> float:
        """Detect potential compositional bias using correlation analysis"""
        # Simplified: check if high-abundance features negatively correlate with others
        if count_table.shape[1] < 10:
            return 0.0
        
        rel_abundance = count_table.div(count_table.sum(axis=1), axis=0)
        top_features = rel_abundance.mean().nlargest(5).index
        other_features = rel_abundance.columns.difference(top_features)
        
        correlations = []
        for top in top_features:
            for other in other_features[:10]:  # Sample for speed
                corr = rel_abundance[top].corr(rel_abundance[other])
                if not np.isnan(corr):
                    correlations.append(corr)
        
        if correlations:
            return abs(np.median(correlations))
        return 0.0


class MethodSelector:
    """Recommends best methods based on data profile"""
    
    def __init__(self):
        self.method_performance = self._load_method_performance()
    
    def _load_method_performance(self) -> Dict:
        """Load method performance characteristics from benchmarking studies"""
        return {
            'aldex2': {
                'good_for': ['asv', 'high_compositional_bias', 'medium_sample_size'],
                'bad_for': ['extreme_sparsity', 'very_small_sample_size'],
                'fdr_control': 'excellent',
                'power': 'moderate'
            },
            'ancom-bc': {
                'good_for': ['compositional_bias', 'asv', 'gene'],
                'bad_for': ['very_small_sample_size'],
                'fdr_control': 'good',
                'power': 'moderate'
            },
            'deseq2': {
                'good_for': ['gene', 'low_sparsity', 'small_sample_size'],
                'bad_for': ['high_compositional_bias', 'extreme_sparsity'],
                'fdr_control': 'moderate',
                'power': 'high'
            },
            'edger': {
                'good_for': ['gene', 'low_sparsity'],
                'bad_for': ['high_sparsity', 'compositional_bias'],
                'fdr_control': 'poor',
                'power': 'very_high'
            },
            'metagenomeseq': {
                'good_for': ['moderate_sparsity', 'gene', 'asv'],
                'bad_for': ['extreme_sparsity'],
                'fdr_control': 'good',
                'power': 'moderate'
            },
            'zicoseq': {
                'good_for': ['high_sparsity', 'zero_inflation', 'viral'],
                'bad_for': [],
                'fdr_control': 'excellent',
                'power': 'high'
            }
        }
    
    def recommend_methods(self, profile: DataProfile) -> MethodRecommendation:
        """Recommend methods based on data profile"""
        
        scores = {}
        reasoning = []
        
        # Score each method based on data characteristics
        for method, chars in self.method_performance.items():
            score = 0
            method_reasons = []
            
            # Data type compatibility
            if profile.data_type in chars['good_for']:
                score += 2
                method_reasons.append(f"Good for {profile.data_type} data")
            
            # Sample size considerations
            min_samples = min(profile.group_sizes.values()) if profile.group_sizes else profile.n_samples
            if min_samples < 10 and 'small_sample_size' in chars['good_for']:
                score += 1
                method_reasons.append("Handles small sample sizes well")
            elif min_samples < 10 and 'very_small_sample_size' in chars['bad_for']:
                score -= 2
                method_reasons.append("Poor with very small sample sizes")
            
            # Sparsity handling
            if profile.sparsity > 0.8:
                if 'high_sparsity' in chars['good_for']:
                    score += 2
                    method_reasons.append("Handles high sparsity")
                elif 'extreme_sparsity' in chars['bad_for']:
                    score -= 2
                    method_reasons.append("Poor with extreme sparsity")
            
            # Compositional bias
            if profile.compositional_bias > 0.3:
                if 'compositional_bias' in chars['good_for'] or 'high_compositional_bias' in chars['good_for']:
                    score += 2
                    method_reasons.append("Addresses compositional bias")
                elif 'compositional_bias' in chars['bad_for']:
                    score -= 1
            
            # Zero inflation
            if profile.zero_inflation > 0.5 and 'zero_inflation' in chars['good_for']:
                score += 1
                method_reasons.append("Handles zero inflation")
            
            scores[method] = (score, method_reasons)
        
        # Sort methods by score
        sorted_methods = sorted(scores.items(), key=lambda x: x[1][0], reverse=True)
        
        # Select primary and secondary methods
        primary_method = sorted_methods[0][0]
        secondary_methods = [m[0] for m in sorted_methods[1:4]]  # Top 3 alternatives
        
        # Generate reasoning
        reasoning.append(f"Primary method: {primary_method}")
        reasoning.extend(scores[primary_method][1])
        reasoning.append(f"Data characteristics: {profile.data_type} data, "
                        f"{profile.n_samples} samples, "
                        f"{profile.sparsity:.1%} sparsity")
        
        # Calculate confidence based on score differences
        if len(sorted_methods) > 1:
            score_diff = sorted_methods[0][1][0] - sorted_methods[1][1][0]
            confidence = min(0.95, 0.5 + score_diff * 0.1)
        else:
            confidence = 0.8
        
        return MethodRecommendation(
            primary_method=primary_method,
            secondary_methods=secondary_methods,
            confidence=confidence,
            reasoning=reasoning
        )


class DAAMethod(ABC):
    """Abstract base class for differential abundance methods"""
    
    @abstractmethod
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """Run differential abundance analysis"""
        pass
    
    @abstractmethod
    def name(self) -> str:
        """Return method name"""
        pass


class WilcoxonMethod(DAAMethod):
    """Simple Wilcoxon rank-sum test implementation"""
    
    def name(self) -> str:
        return "wilcoxon"
    
    def run(self, count_table: pd.DataFrame, metadata: pd.DataFrame, 
            formula: str = None, **kwargs) -> pd.DataFrame:
        """Run Wilcoxon test on each feature"""
        
        # Assume binary comparison using first metadata column
        group_col = metadata.columns[0]
        groups = metadata[group_col].unique()
        
        if len(groups) != 2:
            raise ValueError("Wilcoxon test requires exactly 2 groups")
        
        group1_samples = metadata[metadata[group_col] == groups[0]].index
        group2_samples = metadata[metadata[group_col] == groups[1]].index
        
        results = []
        for feature in count_table.columns:
            group1_values = count_table.loc[group1_samples, feature]
            group2_values = count_table.loc[group2_samples, feature]
            
            # Skip if all zeros
            if (group1_values == 0).all() and (group2_values == 0).all():
                results.append({
                    'feature': feature,
                    'statistic': np.nan,
                    'pvalue': np.nan,
                    'mean_group1': 0,
                    'mean_group2': 0,
                    'log2fc': 0
                })
                continue
            
            statistic, pvalue = stats.ranksums(group1_values, group2_values)
            
            # Calculate log2 fold change (with pseudocount)
            mean1 = group1_values.mean() + 1
            mean2 = group2_values.mean() + 1
            log2fc = np.log2(mean2 / mean1)
            
            results.append({
                'feature': feature,
                'statistic': statistic,
                'pvalue': pvalue,
                'mean_group1': group1_values.mean(),
                'mean_group2': group2_values.mean(),
                'log2fc': log2fc
            })
        
        results_df = pd.DataFrame(results)
        
        # Multiple testing correction
        from statsmodels.stats.multitest import multipletests
        if not results_df['pvalue'].isna().all():
            _, results_df['padj'], _, _ = multipletests(
                results_df['pvalue'].fillna(1), 
                method='fdr_bh'
            )
        else:
            results_df['padj'] = np.nan
        
        return results_df.sort_values('pvalue')


class DifferentialAbundanceTool:
    """Main tool orchestrating the analysis"""
    
    def __init__(self):
        self.profiler = DataProfiler()
        self.selector = MethodSelector()
        self.methods = {
            'wilcoxon': WilcoxonMethod(),
            # Add other method wrappers here
        }
        self.results = {}
    
    def analyze(self, 
                count_table: pd.DataFrame, 
                metadata: pd.DataFrame,
                data_type: Optional[str] = None,
                use_consensus: bool = True) -> Dict:
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
        profile = self.profiler.profile_data(count_table, metadata, data_type)
        
        # Step 2: Get method recommendations
        recommendations = self.selector.recommend_methods(profile)
        
        logger.info(f"Recommended primary method: {recommendations.primary_method}")
        logger.info(f"Confidence: {recommendations.confidence:.2f}")
        logger.info("Reasoning: " + "; ".join(recommendations.reasoning))
        
        # Step 3: Run analyses
        results = {
            'profile': profile,
            'recommendations': recommendations,
            'analyses': {}
        }
        
        # Run primary method
        primary_method = recommendations.primary_method
        if primary_method in self.methods:
            logger.info(f"Running {primary_method}...")
            results['analyses'][primary_method] = self.methods[primary_method].run(
                count_table, metadata
            )
        else:
            logger.warning(f"Method {primary_method} not yet implemented, falling back to Wilcoxon")
            results['analyses']['wilcoxon'] = self.methods['wilcoxon'].run(
                count_table, metadata
            )
        
        # Run secondary methods if consensus requested
        if use_consensus:
            for method in recommendations.secondary_methods:
                if method in self.methods:
                    logger.info(f"Running {method} for consensus...")
                    results['analyses'][method] = self.methods[method].run(
                        count_table, metadata
                    )
        
        # Generate consensus results if multiple methods were run
        if len(results['analyses']) > 1:
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
                if not feature_data.empty and feature_data['padj'].iloc[0] < 0.05:
                    feature_votes += 1
                    feature_pvalues.append(feature_data['pvalue'].iloc[0])
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
            n_sig = (results['padj'] < 0.05).sum() if 'padj' in results.columns else 0
            print(f"  - {method}: {n_sig} significant features (FDR < 0.05)")
        
        if 'consensus' in self.results:
            consensus = self.results['consensus']
            n_consensus_sig = consensus['consensus_significant'].sum()
            print(f"\nConsensus Results:")
            print(f"  - {n_consensus_sig} features significant in majority of methods")


# Example usage
if __name__ == "__main__":
    # Generate example data
    np.random.seed(42)
    
    # Simulate count data (100 samples x 500 features)
    n_samples = 100
    n_features = 500
    
    # Create sparse count matrix
    counts = np.random.negative_binomial(5, 0.3, size=(n_samples, n_features))
    counts[np.random.random((n_samples, n_features)) < 0.7] = 0  # Add zeros
    
    count_table = pd.DataFrame(
        counts,
        index=[f"Sample_{i}" for i in range(n_samples)],
        columns=[f"ASV_{i}" for i in range(n_features)]
    )
    
    # Create metadata with two groups
    metadata = pd.DataFrame({
        'condition': ['Control'] * 50 + ['Treatment'] * 50,
        'batch': np.random.choice(['A', 'B'], n_samples)
    }, index=count_table.index)
    
    # Run analysis
    tool = DifferentialAbundanceTool()
    results = tool.analyze(count_table, metadata, data_type='asv')
    tool.summarize_results()
