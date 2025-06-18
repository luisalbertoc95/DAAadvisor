#!/usr/bin/env python3
"""
Information Theory Framework for Unified Differential Abundance Analysis

This module implements a unified information-theoretic approach to microbiome
differential abundance analysis, treating the problem as information extraction
from compositionally constrained data.

Based on principles from:
- Maximum entropy principle for compositional data
- Compositional Maximum Entropy (CME) 
- Information geometry on the simplex
- Transfer entropy for microbial interactions
"""

import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, Tuple, List
import logging
from scipy import stats
from scipy.optimize import minimize
from scipy.special import logsumexp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings

logger = logging.getLogger(__name__)

class CompositionInformationFramework:
    """
    Unified information theory framework for differential abundance analysis
    
    This framework treats microbiome data as messages from an information source
    with compositional constraints, using entropy-based measures to identify
    differentially abundant features while respecting the simplex geometry.
    """
    
    def __init__(self, regularization: float = 1e-6):
        """
        Initialize the information theory framework
        
        Parameters:
        -----------
        regularization : float
            Small value to add for numerical stability
        """
        self.regularization = regularization
        self.fitted = False
        
    def fit(self, count_table: pd.DataFrame, metadata: pd.DataFrame):
        """
        Fit the information theory model to the data
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata with grouping information
        """
        self.count_table = count_table.copy()
        self.metadata = metadata.copy()
        
        # Convert to relative abundances (compositional data)
        self.compositions = self._to_compositions(count_table)
        
        # Calculate information metrics
        self.entropy_metrics = self._calculate_entropy_metrics()
        self.interaction_network = self._infer_interaction_network()
        self.information_geometry = self._analyze_information_geometry()
        
        self.fitted = True
        return self
    
    def analyze_differential_information(self, 
                                       group_column: str,
                                       alpha: float = 0.05) -> pd.DataFrame:
        """
        Perform differential abundance analysis using information theory
        
        Parameters:
        -----------
        group_column : str
            Column name for grouping variable
        alpha : float
            Significance threshold
            
        Returns:
        --------
        pd.DataFrame
            Results with information-theoretic statistics
        """
        if not self.fitted:
            raise ValueError("Model must be fitted before analysis")
            
        groups = self.metadata[group_column].unique()
        if len(groups) != 2:
            raise ValueError(f"Only two-group comparisons supported, found {len(groups)} groups")
            
        results = []
        
        for feature_idx, feature_name in enumerate(self.count_table.columns):
            # Calculate information-theoretic measures
            result = self._analyze_feature_information(feature_idx, feature_name, group_column)
            results.append(result)
        
        results_df = pd.DataFrame(results)
        
        # Calculate FDR correction
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
        results_df['padj'] = padj
        
        # Sort by information divergence (effect size)
        results_df = results_df.sort_values('information_divergence', ascending=False)
        
        return results_df
    
    def _to_compositions(self, count_table: pd.DataFrame) -> pd.DataFrame:
        """Convert count data to compositional data (relative abundances)"""
        # Add pseudocount for zero handling
        counts_pseudo = count_table + self.regularization
        
        # Calculate relative abundances
        compositions = counts_pseudo.div(counts_pseudo.sum(axis=1), axis=0)
        
        return compositions
    
    def _calculate_entropy_metrics(self) -> Dict[str, np.ndarray]:
        """Calculate entropy-based metrics for each sample and feature"""
        compositions = self.compositions.values
        
        # Sample-wise entropy (diversity)
        sample_entropy = -np.sum(compositions * np.log(compositions + self.regularization), axis=1)
        
        # Feature-wise entropy (prevalence uncertainty)
        feature_entropy = -np.sum(compositions * np.log(compositions + self.regularization), axis=0)
        
        # Cross-entropy between samples
        n_samples = compositions.shape[0]
        cross_entropy_matrix = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(n_samples):
                if i != j:
                    # Cross-entropy H(P_i, P_j) = -sum(P_i * log(P_j))
                    cross_entropy_matrix[i, j] = -np.sum(
                        compositions[i] * np.log(compositions[j] + self.regularization)
                    )
        
        return {
            'sample_entropy': sample_entropy,
            'feature_entropy': feature_entropy,
            'cross_entropy_matrix': cross_entropy_matrix
        }
    
    def _infer_interaction_network(self) -> np.ndarray:
        """Infer feature interaction network using transfer entropy"""
        compositions = self.compositions.values
        n_features = compositions.shape[1]
        
        # Simplified transfer entropy calculation
        # In practice, this would use more sophisticated time-series methods
        interaction_matrix = np.zeros((n_features, n_features))
        
        for i in range(n_features):
            for j in range(n_features):
                if i != j:
                    # Calculate mutual information as proxy for transfer entropy
                    mi = self._calculate_mutual_information(
                        compositions[:, i], compositions[:, j]
                    )
                    interaction_matrix[i, j] = mi
        
        return interaction_matrix
    
    def _calculate_mutual_information(self, x: np.ndarray, y: np.ndarray) -> float:
        """Calculate mutual information between two variables"""
        # Discretize continuous variables for MI calculation
        n_bins = min(10, len(x) // 3)
        
        x_disc = pd.cut(x, bins=n_bins, labels=False)
        y_disc = pd.cut(y, bins=n_bins, labels=False)
        
        # Calculate mutual information
        contingency = pd.crosstab(x_disc, y_disc)
        contingency_norm = contingency / contingency.sum().sum()
        
        # Add small value to avoid log(0)
        contingency_norm = contingency_norm + self.regularization
        
        # MI = sum(p(x,y) * log(p(x,y) / (p(x) * p(y))))
        px = contingency_norm.sum(axis=1)
        py = contingency_norm.sum(axis=0)
        
        mi = 0
        for i in range(len(px)):
            for j in range(len(py)):
                if contingency_norm.iloc[i, j] > 0:
                    mi += contingency_norm.iloc[i, j] * np.log(
                        contingency_norm.iloc[i, j] / (px.iloc[i] * py.iloc[j])
                    )
        
        return mi
    
    def _analyze_information_geometry(self) -> Dict[str, Any]:
        """Analyze the information geometry of the compositional space"""
        compositions = self.compositions.values
        
        # Center log-ratio (CLR) transformation for Aitchison geometry
        clr_data = self._clr_transform(compositions)
        
        # Principal component analysis in CLR space
        pca = PCA()
        clr_pca = pca.fit_transform(clr_data)
        
        # Calculate information distance matrix using Jensen-Shannon divergence
        n_samples = compositions.shape[0]
        js_distance_matrix = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(i+1, n_samples):
                js_dist = self._jensen_shannon_divergence(compositions[i], compositions[j])
                js_distance_matrix[i, j] = js_dist
                js_distance_matrix[j, i] = js_dist
        
        return {
            'clr_data': clr_data,
            'pca_components': pca.components_,
            'pca_explained_variance': pca.explained_variance_ratio_,
            'js_distance_matrix': js_distance_matrix
        }
    
    def _clr_transform(self, compositions: np.ndarray) -> np.ndarray:
        """Centered log-ratio transformation"""
        # Add pseudocount for zero handling
        comp_pseudo = compositions + self.regularization
        
        # CLR transformation
        log_comp = np.log(comp_pseudo)
        geometric_mean = np.mean(log_comp, axis=1, keepdims=True)
        clr = log_comp - geometric_mean
        
        return clr
    
    def _jensen_shannon_divergence(self, p: np.ndarray, q: np.ndarray) -> float:
        """Calculate Jensen-Shannon divergence between two probability distributions"""
        # Normalize to ensure they sum to 1
        p = p / np.sum(p)
        q = q / np.sum(q)
        
        # Add small value for numerical stability
        p = p + self.regularization
        q = q + self.regularization
        
        # Average distribution
        m = 0.5 * (p + q)
        
        # Calculate KL divergences
        kl_pm = np.sum(p * np.log(p / m))
        kl_qm = np.sum(q * np.log(q / m))
        
        # Jensen-Shannon divergence
        js = 0.5 * kl_pm + 0.5 * kl_qm
        
        return js
    
    def _analyze_feature_information(self, 
                                   feature_idx: int, 
                                   feature_name: str,
                                   group_column: str) -> Dict[str, Any]:
        """Analyze information content of a single feature"""
        
        # Get feature data
        feature_data = self.compositions.iloc[:, feature_idx].values
        groups = self.metadata[group_column].values
        unique_groups = np.unique(groups)
        
        # Split by groups
        group1_data = feature_data[groups == unique_groups[0]]
        group2_data = feature_data[groups == unique_groups[1]]
        
        # Calculate group-specific entropy
        entropy1 = self._calculate_entropy(group1_data)
        entropy2 = self._calculate_entropy(group2_data)
        
        # Calculate information divergence (relative entropy/KL divergence)
        info_divergence = self._calculate_kl_divergence(group1_data, group2_data)
        
        # Calculate mutual information with group membership
        group_numeric = (groups == unique_groups[1]).astype(int)
        mutual_info = self._calculate_mutual_information(feature_data, group_numeric)
        
        # Statistical test using information-theoretic approach
        pvalue = self._information_statistical_test(group1_data, group2_data)
        
        # Effect size as information ratio
        effect_size = np.abs(np.log2(np.mean(group1_data + self.regularization) / 
                                   np.mean(group2_data + self.regularization)))
        
        return {
            'feature': feature_name,
            'entropy_group1': entropy1,
            'entropy_group2': entropy2,
            'information_divergence': info_divergence,
            'mutual_information': mutual_info,
            'effect_size': effect_size,
            'log2fc': effect_size * np.sign(np.mean(group1_data) - np.mean(group2_data)),
            'pvalue': pvalue,
            'statistic': info_divergence  # Use information divergence as test statistic
        }
    
    def _calculate_entropy(self, data: np.ndarray) -> float:
        """Calculate entropy of a data distribution
        
        For probability distributions (compositional data), calculate Shannon entropy directly.
        For count data, use histogram discretization.
        """
        # Check if data is already a probability distribution (sums to ~1)
        data_sum = np.sum(data)
        is_probability = np.abs(data_sum - 1.0) < 1e-6
        
        if is_probability:
            # Data is already a probability distribution - calculate Shannon entropy directly
            # Add small regularization to avoid log(0)
            data_reg = data + self.regularization
            data_reg = data_reg / np.sum(data_reg)  # Re-normalize after regularization
            
            # Calculate Shannon entropy: H(X) = -sum(p(x) * log2(p(x)))
            entropy = -np.sum(data_reg * np.log2(data_reg))
            return entropy
        
        else:
            # Data is count/continuous data - use histogram discretization
            n_bins = min(10, len(data) // 3) if len(data) > 10 else len(data)
            
            if n_bins <= 1:
                return 0.0
            
            counts, _ = np.histogram(data, bins=n_bins)
            # Remove zero counts
            counts = counts[counts > 0]
            
            if len(counts) == 0:
                return 0.0
            
            # Normalize to probabilities
            probabilities = counts / np.sum(counts)
            
            # Calculate entropy
            entropy = -np.sum(probabilities * np.log2(probabilities))
            
            return entropy
    
    def _calculate_kl_divergence(self, data1: np.ndarray, data2: np.ndarray) -> float:
        """Calculate KL divergence between two data distributions"""
        # Create histograms with same bins
        all_data = np.concatenate([data1, data2])
        n_bins = min(10, len(all_data) // 5)
        
        if n_bins <= 1:
            return 0.0
        
        bin_edges = np.histogram_bin_edges(all_data, bins=n_bins)
        
        counts1, _ = np.histogram(data1, bins=bin_edges)
        counts2, _ = np.histogram(data2, bins=bin_edges)
        
        # Convert to probabilities with pseudocounts
        p1 = (counts1 + self.regularization) / (np.sum(counts1) + n_bins * self.regularization)
        p2 = (counts2 + self.regularization) / (np.sum(counts2) + n_bins * self.regularization)
        
        # Calculate KL divergence
        kl_div = np.sum(p1 * np.log(p1 / p2))
        
        return kl_div
    
    def _information_statistical_test(self, group1: np.ndarray, group2: np.ndarray) -> float:
        """Statistical test based on information divergence"""
        
        # Use permutation test with information divergence as test statistic
        observed_divergence = self._calculate_kl_divergence(group1, group2)
        
        # Perform permutation test
        n_permutations = 1000
        all_data = np.concatenate([group1, group2])
        n1 = len(group1)
        
        permuted_divergences = []
        
        for _ in range(n_permutations):
            # Randomly permute group assignments
            np.random.shuffle(all_data)
            perm_group1 = all_data[:n1]
            perm_group2 = all_data[n1:]
            
            perm_divergence = self._calculate_kl_divergence(perm_group1, perm_group2)
            permuted_divergences.append(perm_divergence)
        
        # Calculate p-value
        permuted_divergences = np.array(permuted_divergences)
        pvalue = np.mean(permuted_divergences >= observed_divergence)
        
        # Ensure p-value is not exactly 0
        if pvalue == 0:
            pvalue = 1.0 / (n_permutations + 1)
        
        return pvalue


class MaximumEntropySelector:
    """
    Method selector based on maximum entropy principle
    
    This class uses information theory to select the most appropriate
    differential abundance method based on data characteristics.
    """
    
    def __init__(self):
        self.framework = CompositionInformationFramework()
        
    def select_method(self, 
                     count_table: pd.DataFrame, 
                     metadata: pd.DataFrame) -> Dict[str, Any]:
        """
        Select optimal method using maximum entropy principle
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Samples x Features count matrix
        metadata : pd.DataFrame
            Sample metadata
            
        Returns:
        --------
        Dict[str, Any]
            Method recommendation with information-theoretic justification
        """
        
        # Fit information framework
        self.framework.fit(count_table, metadata)
        
        # Calculate data characteristics from information perspective
        info_metrics = self._calculate_information_characteristics()
        
        # Score methods based on information theory principles
        method_scores = self._score_methods_information_theory(info_metrics)
        
        # Select method with maximum expected information preservation
        best_method = max(method_scores.items(), key=lambda x: x[1]['score'])
        
        return {
            'recommended_method': best_method[0],
            'score': best_method[1]['score'],
            'information_metrics': info_metrics,
            'method_scores': method_scores,
            'reasoning': best_method[1]['reasoning']
        }
    
    def _calculate_information_characteristics(self) -> Dict[str, float]:
        """Calculate information-theoretic characteristics of the data"""
        
        compositions = self.framework.compositions.values
        entropy_metrics = self.framework.entropy_metrics
        
        # Overall information content
        total_entropy = np.mean(entropy_metrics['sample_entropy'])
        
        # Information complexity (variability in entropy)
        entropy_cv = np.std(entropy_metrics['sample_entropy']) / np.mean(entropy_metrics['sample_entropy'])
        
        # Sparsity from information perspective
        # High sparsity = low information per feature
        feature_info_content = entropy_metrics['feature_entropy']
        sparsity_info = 1.0 - (np.mean(feature_info_content) / np.log(compositions.shape[0]))
        
        # Compositional constraint strength
        # Measured as deviation from uniform distribution
        uniform_entropy = np.log(compositions.shape[1])
        constraint_strength = 1.0 - (total_entropy / uniform_entropy)
        
        # Information interaction strength
        interaction_matrix = self.framework.interaction_network
        interaction_strength = np.mean(interaction_matrix[interaction_matrix > 0])
        
        return {
            'total_entropy': total_entropy,
            'entropy_complexity': entropy_cv,
            'sparsity_info': sparsity_info,
            'constraint_strength': constraint_strength,
            'interaction_strength': interaction_strength
        }
    
    def _score_methods_information_theory(self, info_metrics: Dict[str, float]) -> Dict[str, Dict[str, Any]]:
        """Score methods based on information theory principles"""
        
        scores = {}
        
        # ALDEx2 - Best for high compositional constraints and moderate complexity
        aldex2_score = (
            (1.0 - info_metrics['sparsity_info']) * 0.4 +  # Handles moderate sparsity well
            info_metrics['constraint_strength'] * 0.4 +     # Designed for compositional data
            (1.0 - info_metrics['entropy_complexity']) * 0.2  # Works with stable entropy
        )
        scores['aldex2'] = {
            'score': aldex2_score,
            'reasoning': 'Strong compositional awareness with CLR transformation preserves information geometry'
        }
        
        # DESeq2 - Best for high information content, low compositional constraints
        deseq2_score = (
            info_metrics['total_entropy'] * 0.4 +           # Needs high information content
            (1.0 - info_metrics['constraint_strength']) * 0.4 +  # Assumes count data properties
            (1.0 - info_metrics['sparsity_info']) * 0.2     # Works with dense data
        )
        scores['deseq2'] = {
            'score': deseq2_score,
            'reasoning': 'Maximum likelihood approach optimal for high-information, low-constraint scenarios'
        }
        
        # edgeR - Similar to DESeq2 but more robust to complexity
        edger_score = (
            info_metrics['total_entropy'] * 0.3 +
            (1.0 - info_metrics['constraint_strength']) * 0.3 +
            info_metrics['entropy_complexity'] * 0.2 +      # Handles complexity better
            (1.0 - info_metrics['sparsity_info']) * 0.2
        )
        scores['edger'] = {
            'score': edger_score,
            'reasoning': 'Robust negative binomial modeling handles information complexity well'
        }
        
        # Wilcoxon - Best for low information content, high uncertainty
        wilcoxon_score = (
            (1.0 - info_metrics['total_entropy']) * 0.4 +   # Non-parametric for low info
            info_metrics['sparsity_info'] * 0.3 +           # Robust to sparsity
            info_metrics['entropy_complexity'] * 0.3        # Handles complex distributions
        )
        scores['wilcoxon'] = {
            'score': wilcoxon_score,
            'reasoning': 'Non-parametric approach preserves rank information without distributional assumptions'
        }
        
        # ANCOM-BC - Balanced approach with bias correction
        ancombc_score = (
            info_metrics['constraint_strength'] * 0.3 +
            (1.0 - info_metrics['sparsity_info']) * 0.3 +
            info_metrics['interaction_strength'] * 0.2 +    # Accounts for feature interactions
            info_metrics['total_entropy'] * 0.2
        )
        scores['ancom-bc'] = {
            'score': ancombc_score,
            'reasoning': 'Bias correction maintains information integrity under compositional constraints'
        }
        
        return scores


def run_information_theory_analysis(count_table: pd.DataFrame,
                                   metadata: pd.DataFrame,
                                   group_column: str,
                                   alpha: float = 0.05) -> Dict[str, Any]:
    """
    Run complete information theory-based differential abundance analysis
    
    Parameters:
    -----------
    count_table : pd.DataFrame
        Samples x Features count matrix
    metadata : pd.DataFrame
        Sample metadata
    group_column : str
        Column name for grouping variable
    alpha : float
        Significance threshold
        
    Returns:
    --------
    Dict[str, Any]
        Complete analysis results including method selection and DA results
    """
    
    logger.info("Starting information theory-based differential abundance analysis...")
    
    # Initialize frameworks
    selector = MaximumEntropySelector()
    framework = CompositionInformationFramework()
    
    # Select optimal method using information theory
    logger.info("Selecting optimal method using maximum entropy principle...")
    method_selection = selector.select_method(count_table, metadata)
    
    # Run information theory-based differential abundance analysis
    logger.info("Running information theory differential abundance analysis...")
    framework.fit(count_table, metadata)
    da_results = framework.analyze_differential_information(group_column, alpha)
    
    # Calculate additional information metrics
    info_summary = {
        'total_samples': len(count_table),
        'total_features': len(count_table.columns),
        'information_entropy': np.mean(framework.entropy_metrics['sample_entropy']),
        'compositional_constraint': method_selection['information_metrics']['constraint_strength'],
        'sparsity_information': method_selection['information_metrics']['sparsity_info'],
        'significant_features': (da_results['padj'] < alpha).sum()
    }
    
    logger.info(f"Information theory analysis complete. Found {info_summary['significant_features']} significant features.")
    
    return {
        'method_selection': method_selection,
        'differential_abundance_results': da_results,
        'information_summary': info_summary,
        'framework': framework  # For further analysis
    }