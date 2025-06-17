#!/usr/bin/env python3
"""
Data profiling module for DAAadvisor
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
from dataclasses import dataclass
import logging
from scipy import stats

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
    # Additional advanced metrics
    prevalence_distribution: Dict[str, float]
    feature_variance: float
    sample_diversity: Dict[str, float]
    batch_effects: Dict[str, float]


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
            metadata_factors=list(metadata.columns),
            prevalence_distribution=self._analyze_prevalence(count_table),
            feature_variance=self._calculate_feature_variance(count_table),
            sample_diversity=self._calculate_sample_diversity(count_table),
            batch_effects=self._detect_batch_effects(count_table, metadata)
        )
        
        self.profile = profile
        logger.info(f"Data profiled: {profile.n_samples} samples, {profile.n_features} features")
        logger.info(f"Sparsity: {profile.sparsity:.2f}, Zero inflation: {profile.zero_inflation:.2f}")
        logger.info(f"Data type: {profile.data_type}, Compositional bias: {profile.compositional_bias:.2f}")
        
        return profile
    
    def _detect_data_type(self, count_table: pd.DataFrame) -> str:
        """Attempt to detect data type from feature names and characteristics"""
        feature_names = count_table.columns
        
        # Check for common patterns
        if any('ASV' in str(f) or 'OTU' in str(f) for f in feature_names):
            return 'asv'
        elif any('gene' in str(f).lower() or 'KO' in str(f) or 'COG' in str(f) for f in feature_names):
            return 'gene'
        elif any('virus' in str(f).lower() or 'phage' in str(f).lower() for f in feature_names):
            return 'viral'
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
        # Calculate proportion of features with >50% zeros
        feature_zero_props = (count_table == 0).mean(axis=0)
        return (feature_zero_props > 0.5).mean()
    
    def _calculate_depth_variation(self, count_table: pd.DataFrame) -> float:
        """Calculate coefficient of variation in sequencing depth"""
        depths = count_table.sum(axis=1)
        if depths.mean() == 0:
            return 0.0
        return depths.std() / depths.mean()
    
    def _get_group_sizes(self, metadata: pd.DataFrame) -> Dict[str, int]:
        """Get sample sizes for each group in primary factor"""
        if metadata.shape[1] > 0:
            primary_factor = metadata.iloc[:, 0]
            return primary_factor.value_counts().to_dict()
        return {}
    
    def _detect_compositional_bias(self, count_table: pd.DataFrame) -> float:
        """Detect potential compositional bias using correlation analysis"""
        if count_table.shape[1] < 10:
            return 0.0
        
        # Convert to relative abundance
        rel_abundance = count_table.div(count_table.sum(axis=1), axis=0).fillna(0)
        
        # Find top abundant features
        mean_abundance = rel_abundance.mean()
        top_features = mean_abundance.nlargest(min(5, len(mean_abundance))).index
        other_features = rel_abundance.columns.difference(top_features)
        
        # Calculate correlations between top and other features
        correlations = []
        for top in top_features:
            for other in other_features[:min(20, len(other_features))]:  # Sample for speed
                corr = rel_abundance[top].corr(rel_abundance[other])
                if not np.isnan(corr):
                    correlations.append(corr)
        
        if correlations:
            # High negative correlation indicates compositional bias
            return abs(np.median([c for c in correlations if c < 0]))
        return 0.0
    
    def _analyze_prevalence(self, count_table: pd.DataFrame) -> Dict[str, float]:
        """Analyze feature prevalence distribution"""
        prevalence = (count_table > 0).mean(axis=0)
        
        return {
            'low_prevalence': (prevalence < 0.1).mean(),  # <10% samples
            'medium_prevalence': ((prevalence >= 0.1) & (prevalence < 0.5)).mean(),  # 10-50%
            'high_prevalence': (prevalence >= 0.5).mean(),  # >50%
            'core_features': (prevalence > 0.8).mean()  # >80% samples
        }
    
    def _calculate_feature_variance(self, count_table: pd.DataFrame) -> float:
        """Calculate mean coefficient of variation across features"""
        # Use relative abundance to normalize for sequencing depth
        rel_abundance = count_table.div(count_table.sum(axis=1), axis=0).fillna(0)
        
        feature_cvs = []
        for col in rel_abundance.columns:
            values = rel_abundance[col]
            if values.mean() > 0:
                cv = values.std() / values.mean()
                feature_cvs.append(cv)
        
        return np.mean(feature_cvs) if feature_cvs else 0.0
    
    def _calculate_sample_diversity(self, count_table: pd.DataFrame) -> Dict[str, float]:
        """Calculate alpha diversity metrics"""
        diversity_metrics = {}
        
        # Shannon diversity
        shannon_div = []
        for idx in count_table.index:
            counts = count_table.loc[idx]
            counts = counts[counts > 0]  # Remove zeros
            if len(counts) > 0:
                props = counts / counts.sum()
                shannon = -np.sum(props * np.log(props))
                shannon_div.append(shannon)
        
        diversity_metrics['shannon_mean'] = np.mean(shannon_div) if shannon_div else 0.0
        diversity_metrics['shannon_cv'] = (np.std(shannon_div) / np.mean(shannon_div)) if shannon_div and np.mean(shannon_div) > 0 else 0.0
        
        # Observed features (richness)
        richness = (count_table > 0).sum(axis=1)
        diversity_metrics['richness_mean'] = richness.mean()
        diversity_metrics['richness_cv'] = richness.std() / richness.mean() if richness.mean() > 0 else 0.0
        
        return diversity_metrics
    
    def _detect_batch_effects(self, count_table: pd.DataFrame, metadata: pd.DataFrame) -> Dict[str, float]:
        """Detect potential batch effects in the data"""
        batch_effects = {}
        
        # Look for batch-related columns
        batch_cols = [col for col in metadata.columns if 
                     any(term in col.lower() for term in ['batch', 'run', 'plate', 'lane', 'date'])]
        
        for batch_col in batch_cols:
            if batch_col in metadata.columns:
                batches = metadata[batch_col].unique()
                if len(batches) > 1:
                    # Calculate variance explained by batch effect
                    # Using a simplified approach with PC1
                    try:
                        from sklearn.decomposition import PCA
                        from sklearn.preprocessing import StandardScaler
                        
                        # Standardize data
                        scaler = StandardScaler()
                        scaled_data = scaler.fit_transform(count_table.fillna(0))
                        
                        # PCA
                        pca = PCA(n_components=1)
                        pc1 = pca.fit_transform(scaled_data).flatten()
                        
                        # Calculate correlation between PC1 and batch
                        batch_numeric = pd.Categorical(metadata[batch_col]).codes
                        correlation = np.corrcoef(pc1, batch_numeric)[0, 1]
                        batch_effects[batch_col] = abs(correlation) if not np.isnan(correlation) else 0.0
                        
                    except ImportError:
                        # Fallback to simple variance calculation
                        batch_means = []
                        for batch in batches:
                            batch_samples = metadata[metadata[batch_col] == batch].index
                            if len(batch_samples) > 0:
                                batch_mean = count_table.loc[batch_samples].mean().mean()
                                batch_means.append(batch_mean)
                        
                        if len(batch_means) > 1:
                            batch_cv = np.std(batch_means) / np.mean(batch_means) if np.mean(batch_means) > 0 else 0.0
                            batch_effects[batch_col] = min(batch_cv, 1.0)  # Cap at 1.0
        
        return batch_effects
    
    def print_profile_summary(self) -> None:
        """Print a detailed summary of the data profile"""
        if not self.profile:
            logger.warning("No profile available. Run profile_data() first.")
            return
        
        profile = self.profile
        
        print("\n" + "="*50)
        print("DATA PROFILE SUMMARY")
        print("="*50)
        
        print(f"\nBasic Characteristics:")
        print(f"  Data Type: {profile.data_type}")
        print(f"  Samples: {profile.n_samples}")
        print(f"  Features: {profile.n_features}")
        print(f"  Sparsity: {profile.sparsity:.1%}")
        print(f"  Zero Inflation: {profile.zero_inflation:.1%}")
        print(f"  Sequencing Depth CV: {profile.sequencing_depth_cv:.2f}")
        
        print(f"\nGroup Sizes:")
        for group, size in profile.group_sizes.items():
            print(f"  {group}: {size}")
        
        print(f"\nCompostional Analysis:")
        print(f"  Compositional Bias: {profile.compositional_bias:.2f}")
        
        print(f"\nPrevalence Distribution:")
        for category, prop in profile.prevalence_distribution.items():
            print(f"  {category.replace('_', ' ').title()}: {prop:.1%}")
        
        print(f"\nDiversity Metrics:")
        for metric, value in profile.sample_diversity.items():
            print(f"  {metric.replace('_', ' ').title()}: {value:.2f}")
        
        if profile.batch_effects:
            print(f"\nPotential Batch Effects:")
            for batch_col, effect_size in profile.batch_effects.items():
                print(f"  {batch_col}: {effect_size:.2f}")