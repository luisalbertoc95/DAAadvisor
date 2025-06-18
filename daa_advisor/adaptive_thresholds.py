#!/usr/bin/env python3
"""
Adaptive Thresholding Module Using Information Theory

This module provides adaptive significance thresholds based on information content
and data characteristics, improving upon standard fixed-threshold approaches.
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, Optional
import logging
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

class AdaptiveThresholdSelector:
    """
    Information theory-based adaptive threshold selection
    
    Methods:
    - Information content-weighted FDR correction
    - Data-adaptive significance levels
    - Effect size-informed thresholding
    - Entropy-based confidence scoring
    """
    
    def __init__(self, base_alpha: float = 0.05):
        self.base_alpha = base_alpha
        self.threshold_stats = {}
        
    def apply_adaptive_thresholds(self,
                                results: pd.DataFrame,
                                count_table: pd.DataFrame,
                                metadata: pd.DataFrame,
                                group_column: str,
                                method: str = 'information_weighted',
                                verbose: bool = True) -> pd.DataFrame:
        """
        Apply adaptive thresholds based on information theory
        
        Parameters:
        -----------
        results : pd.DataFrame
            Differential abundance results with p-values
        count_table : pd.DataFrame
            Original count data
        metadata : pd.DataFrame
            Sample metadata
        group_column : str
            Grouping variable column name
        method : str
            Adaptive method ('information_weighted', 'entropy_based', 'effect_size_weighted')
        verbose : bool
            Print threshold details
            
        Returns:
        --------
        pd.DataFrame
            Results with adaptive p-value corrections
        """
        
        if verbose:
            logger.info(f"🎯 Applying {method} adaptive thresholds")
            logger.info(f"📊 Base α = {self.base_alpha}")
        
        # Calculate information content for each feature
        info_weights = self._calculate_information_weights(
            results, count_table, metadata, group_column
        )
        
        # Apply selected adaptive method
        if method == 'information_weighted':
            corrected_results = self._information_weighted_fdr(
                results, info_weights, verbose
            )
        elif method == 'entropy_based':
            corrected_results = self._entropy_based_thresholds(
                results, count_table, info_weights, verbose
            )
        elif method == 'effect_size_weighted':
            corrected_results = self._effect_size_weighted_fdr(
                results, info_weights, verbose
            )
        else:
            raise ValueError(f"Unknown adaptive method: {method}")
        
        # Add confidence scores based on information content
        corrected_results = self._add_confidence_scores(
            corrected_results, info_weights, verbose
        )
        
        if verbose:
            n_significant = (corrected_results['padj_adaptive'] < self.base_alpha).sum()
            n_original = (corrected_results.get('padj', corrected_results.get('pvalue', [1]*len(corrected_results))) < self.base_alpha).sum()
            logger.info(f"📈 Significant features: {n_original} → {n_significant}")
        
        return corrected_results
    
    def _calculate_information_weights(self,
                                     results: pd.DataFrame,
                                     count_table: pd.DataFrame,
                                     metadata: pd.DataFrame,
                                     group_column: str) -> np.ndarray:
        """Calculate information content weights for each feature"""
        
        # Convert to compositions
        compositions = count_table.div(count_table.sum(axis=1), axis=0).fillna(0)
        groups = metadata[group_column].values
        unique_groups = np.unique(groups)
        
        if len(unique_groups) != 2:
            raise ValueError("Adaptive thresholds currently support only two-group comparisons")
        
        info_weights = []
        
        for feature in results['feature']:
            if feature not in compositions.columns:
                info_weights.append(0.0)
                continue
                
            feature_data = compositions[feature].values
            
            # Split by groups
            group1_data = feature_data[groups == unique_groups[0]]
            group2_data = feature_data[groups == unique_groups[1]]
            
            # Calculate multiple information measures
            weight = 0.0
            
            # 1. Jensen-Shannon divergence (between-group information)
            js_divergence = self._jensen_shannon_divergence(group1_data, group2_data)
            weight += js_divergence * 0.3
            
            # 2. Shannon entropy (within-group variability)
            entropy1 = self._shannon_entropy(group1_data)
            entropy2 = self._shannon_entropy(group2_data)
            avg_entropy = (entropy1 + entropy2) / 2
            weight += avg_entropy * 0.2
            
            # 3. Effect size (fold change)
            mean1 = np.mean(group1_data) + 1e-10
            mean2 = np.mean(group2_data) + 1e-10
            log_fold_change = abs(np.log2(mean2 / mean1))
            weight += log_fold_change * 0.2
            
            # 4. Prevalence information
            prevalence = np.mean(feature_data > 0)
            # Moderate prevalence (0.3-0.7) gets highest weight
            prevalence_weight = 1 - abs(prevalence - 0.5) * 2
            weight += prevalence_weight * 0.1
            
            # 5. Variance information
            total_variance = np.var(feature_data)
            weight += min(total_variance, 1.0) * 0.2  # Cap at 1.0
            
            info_weights.append(weight)
        
        info_weights = np.array(info_weights)
        
        # Normalize weights to [0.1, 2.0] range
        # Features with very low information get minimum weight
        # Features with high information get up to 2x weight
        if np.max(info_weights) > 0:
            info_weights = info_weights / np.max(info_weights)
        info_weights = 0.1 + 1.9 * info_weights
        
        return info_weights
    
    def _information_weighted_fdr(self,
                                results: pd.DataFrame,
                                info_weights: np.ndarray,
                                verbose: bool) -> pd.DataFrame:
        """Information content-weighted FDR correction"""
        
        # Standard Benjamini-Hochberg with information weighting
        pvalues = results['pvalue'].values
        n_tests = len(pvalues)
        
        # Sort by p-values
        sort_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sort_idx]
        sorted_weights = info_weights[sort_idx]
        
        # Apply weighted Benjamini-Hochberg procedure
        # Higher information content features get more lenient thresholds
        adjusted_pvals = np.ones(n_tests)
        
        for i in range(n_tests):
            # Adaptive alpha based on information content
            adaptive_alpha = self.base_alpha * sorted_weights[i]
            
            # Standard BH threshold with adaptive alpha
            bh_threshold = adaptive_alpha * (i + 1) / n_tests
            
            if sorted_pvals[i] <= bh_threshold:
                adjusted_pvals[i] = sorted_pvals[i] * n_tests / (i + 1)
            else:
                # Once we fail, all subsequent tests fail
                break
        
        # Ensure monotonicity
        for i in range(n_tests - 2, -1, -1):
            adjusted_pvals[i] = min(adjusted_pvals[i], adjusted_pvals[i + 1])
        
        # Map back to original order
        original_adjusted = np.ones(n_tests)
        original_adjusted[sort_idx] = adjusted_pvals
        
        # Cap at 1.0
        original_adjusted = np.minimum(original_adjusted, 1.0)
        
        results_copy = results.copy()
        results_copy['padj_adaptive'] = original_adjusted
        results_copy['info_weight'] = info_weights
        
        if verbose:
            logger.info(f"    Mean information weight: {np.mean(info_weights):.3f}")
            logger.info(f"    Weight range: [{np.min(info_weights):.3f}, {np.max(info_weights):.3f}]")
        
        return results_copy
    
    def _entropy_based_thresholds(self,
                                results: pd.DataFrame,
                                count_table: pd.DataFrame,
                                info_weights: np.ndarray,
                                verbose: bool) -> pd.DataFrame:
        """Entropy-based adaptive thresholds"""
        
        # Calculate dataset-level entropy characteristics
        compositions = count_table.div(count_table.sum(axis=1), axis=0).fillna(0)
        
        # Sample-wise entropy (diversity)
        sample_entropies = []
        for i in range(len(compositions)):
            sample_data = compositions.iloc[i].values
            entropy = self._shannon_entropy(sample_data)
            sample_entropies.append(entropy)
        
        dataset_entropy = np.mean(sample_entropies)
        entropy_cv = np.std(sample_entropies) / dataset_entropy if dataset_entropy > 0 else 0
        
        # Adaptive threshold based on dataset entropy characteristics
        # High entropy datasets get more lenient thresholds (more complex)
        # Low entropy datasets get stricter thresholds (simpler patterns)
        entropy_factor = 1.0 + 0.5 * (dataset_entropy / np.log(count_table.shape[1]))
        variability_factor = 1.0 - 0.3 * entropy_cv  # High variability = stricter
        
        adaptive_alpha = self.base_alpha * entropy_factor * variability_factor
        adaptive_alpha = np.clip(adaptive_alpha, 0.01, 0.2)  # Reasonable bounds
        
        # Apply standard FDR with adaptive alpha
        _, padj_adaptive, _, _ = multipletests(
            results['pvalue'], 
            alpha=adaptive_alpha, 
            method='fdr_bh'
        )
        
        results_copy = results.copy()
        results_copy['padj_adaptive'] = padj_adaptive
        results_copy['info_weight'] = info_weights
        results_copy['adaptive_alpha'] = adaptive_alpha
        
        if verbose:
            logger.info(f"    Dataset entropy: {dataset_entropy:.3f}")
            logger.info(f"    Entropy CV: {entropy_cv:.3f}")
            logger.info(f"    Adaptive α: {adaptive_alpha:.4f}")
        
        return results_copy
    
    def _effect_size_weighted_fdr(self,
                                results: pd.DataFrame,
                                info_weights: np.ndarray,
                                verbose: bool) -> pd.DataFrame:
        """Effect size-weighted FDR correction"""
        
        # Extract effect sizes (could be log fold change, test statistic, etc.)
        if 'log2fc' in results.columns:
            effect_sizes = np.abs(results['log2fc'].values)
        elif 'effect_size' in results.columns:
            effect_sizes = np.abs(results['effect_size'].values)
        elif 'statistic' in results.columns:
            effect_sizes = np.abs(results['statistic'].values)
        else:
            # Use information weights as proxy for effect size
            effect_sizes = info_weights
        
        # Normalize effect sizes
        if np.max(effect_sizes) > 0:
            effect_sizes = effect_sizes / np.max(effect_sizes)
        
        # Combine information weights and effect sizes
        combined_weights = 0.6 * info_weights + 0.4 * (effect_sizes + 0.1)
        
        # Apply weighted FDR
        pvalues = results['pvalue'].values
        n_tests = len(pvalues)
        
        # Sort by p-values
        sort_idx = np.argsort(pvalues)
        sorted_pvals = pvalues[sort_idx]
        sorted_weights = combined_weights[sort_idx]
        
        adjusted_pvals = np.ones(n_tests)
        
        for i in range(n_tests):
            # Weight-adjusted threshold
            weight_factor = sorted_weights[i] / np.mean(combined_weights)
            adaptive_threshold = self.base_alpha * weight_factor * (i + 1) / n_tests
            
            if sorted_pvals[i] <= adaptive_threshold:
                adjusted_pvals[i] = sorted_pvals[i] * n_tests / (i + 1)
            else:
                break
        
        # Ensure monotonicity and map back
        for i in range(n_tests - 2, -1, -1):
            adjusted_pvals[i] = min(adjusted_pvals[i], adjusted_pvals[i + 1])
        
        original_adjusted = np.ones(n_tests)
        original_adjusted[sort_idx] = adjusted_pvals
        original_adjusted = np.minimum(original_adjusted, 1.0)
        
        results_copy = results.copy()
        results_copy['padj_adaptive'] = original_adjusted
        results_copy['info_weight'] = info_weights
        results_copy['effect_weight'] = effect_sizes
        results_copy['combined_weight'] = combined_weights
        
        if verbose:
            logger.info(f"    Mean effect size: {np.mean(effect_sizes):.3f}")
            logger.info(f"    Mean combined weight: {np.mean(combined_weights):.3f}")
        
        return results_copy
    
    def _add_confidence_scores(self,
                             results: pd.DataFrame,
                             info_weights: np.ndarray,
                             verbose: bool) -> pd.DataFrame:
        """Add confidence scores based on information content"""
        
        # Confidence based on multiple factors
        pvalues = results['pvalue'].values
        padj = results['padj_adaptive'].values
        
        # P-value confidence (lower p-value = higher confidence)
        pval_confidence = -np.log10(pvalues + 1e-10)
        pval_confidence = pval_confidence / np.max(pval_confidence) if np.max(pval_confidence) > 0 else pval_confidence
        
        # Information content confidence
        info_confidence = (info_weights - np.min(info_weights)) / (np.max(info_weights) - np.min(info_weights))
        
        # Significance confidence (significant features get bonus)
        sig_confidence = (padj < self.base_alpha).astype(float) * 0.5
        
        # Combined confidence score
        confidence_scores = 0.5 * pval_confidence + 0.3 * info_confidence + 0.2 * sig_confidence
        
        results['confidence_score'] = confidence_scores
        results['high_confidence'] = confidence_scores > np.percentile(confidence_scores, 75)
        
        if verbose:
            n_high_conf = np.sum(results['high_confidence'])
            logger.info(f"    High confidence features: {n_high_conf}")
        
        return results
    
    def _shannon_entropy(self, data: np.ndarray) -> float:
        """Calculate Shannon entropy of data distribution"""
        # Remove zeros and normalize
        data_nonzero = data[data > 0]
        if len(data_nonzero) == 0:
            return 0.0
        
        # Normalize to probabilities
        probs = data_nonzero / np.sum(data_nonzero)
        
        # Calculate Shannon entropy
        entropy = -np.sum(probs * np.log2(probs))
        return entropy
    
    def _jensen_shannon_divergence(self, p: np.ndarray, q: np.ndarray) -> float:
        """Calculate Jensen-Shannon divergence"""
        # Normalize
        p = p / (np.sum(p) + 1e-10)
        q = q / (np.sum(q) + 1e-10)
        
        # Add small regularization
        p = p + 1e-10
        q = q + 1e-10
        
        # Re-normalize
        p = p / np.sum(p)
        q = q / np.sum(q)
        
        # Average distribution
        m = 0.5 * (p + q)
        
        # KL divergences
        kl_pm = np.sum(p * np.log(p / m))
        kl_qm = np.sum(q * np.log(q / m))
        
        # JS divergence
        js = 0.5 * kl_pm + 0.5 * kl_qm
        return js


def apply_adaptive_thresholds(results: pd.DataFrame,
                            count_table: pd.DataFrame,
                            metadata: pd.DataFrame,
                            group_column: str,
                            method: str = 'information_weighted',
                            base_alpha: float = 0.05,
                            verbose: bool = True) -> pd.DataFrame:
    """
    Convenience function for applying adaptive thresholds
    
    Parameters:
    -----------
    results : pd.DataFrame
        Differential abundance results
    count_table : pd.DataFrame
        Original count data
    metadata : pd.DataFrame
        Sample metadata
    group_column : str
        Grouping variable
    method : str
        Adaptive method ('information_weighted', 'entropy_based', 'effect_size_weighted')
    base_alpha : float
        Base significance threshold
    verbose : bool
        Print details
        
    Returns:
    --------
    pd.DataFrame
        Results with adaptive corrections
    """
    
    selector = AdaptiveThresholdSelector(base_alpha=base_alpha)
    return selector.apply_adaptive_thresholds(
        results=results,
        count_table=count_table,
        metadata=metadata,
        group_column=group_column,
        method=method,
        verbose=verbose
    )