#!/usr/bin/env python3
"""
Information Theory-Based Preprocessing Module

This module provides legitimate preprocessing improvements using information theory
for better differential abundance analysis performance.
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional
import logging
from scipy import stats
from scipy.special import logsumexp
import warnings

logger = logging.getLogger(__name__)

class InformationBasedPreprocessor:
    """
    Information theory-based preprocessing for microbiome data
    
    Legitimate improvements:
    - Information-content-based feature filtering
    - CLR transformation with optimal pseudocount selection
    - Variance stabilization using entropy measures
    - Outlier detection using information divergence
    """
    
    def __init__(self, regularization: float = 1e-6):
        self.regularization = regularization
        self.preprocessing_stats = {}
        
    def preprocess_data(self, 
                       count_table: pd.DataFrame, 
                       metadata: pd.DataFrame,
                       min_prevalence: float = 0.1,
                       min_information_content: float = 0.1,
                       apply_clr: bool = True,
                       remove_outliers: bool = True,
                       verbose: bool = True) -> Tuple[pd.DataFrame, Dict[str, any]]:
        """
        Apply information theory-based preprocessing pipeline
        
        Parameters:
        -----------
        count_table : pd.DataFrame
            Raw count data (samples x features)
        metadata : pd.DataFrame
            Sample metadata
        min_prevalence : float
            Minimum prevalence threshold (0-1)
        min_information_content : float
            Minimum information content threshold (0-1)
        apply_clr : bool
            Apply CLR transformation
        remove_outliers : bool
            Remove outlier samples using information divergence
        verbose : bool
            Print processing steps
            
        Returns:
        --------
        Tuple[pd.DataFrame, Dict]
            Preprocessed count table and processing statistics
        """
        
        if verbose:
            logger.info("🔄 Information Theory-Based Preprocessing Pipeline")
            logger.info(f"📊 Input: {count_table.shape[0]} samples × {count_table.shape[1]} features")
        
        processed_data = count_table.copy()
        stats_dict = {'original_shape': count_table.shape}
        
        # 1. Prevalence-based filtering (standard approach)
        if verbose:
            logger.info("  📈 Step 1: Prevalence-based feature filtering...")
        processed_data, prevalence_stats = self._prevalence_filter(
            processed_data, min_prevalence, verbose
        )
        stats_dict.update(prevalence_stats)
        
        # 2. Information content filtering (Information Theory improvement)
        if verbose:
            logger.info("  🧮 Step 2: Information content filtering...")
        processed_data, info_stats = self._information_content_filter(
            processed_data, min_information_content, verbose
        )
        stats_dict.update(info_stats)
        
        # 3. Outlier detection using information divergence
        if remove_outliers:
            if verbose:
                logger.info("  🎯 Step 3: Information divergence outlier detection...")
            processed_data, outlier_stats = self._detect_information_outliers(
                processed_data, metadata, verbose
            )
            stats_dict.update(outlier_stats)
        
        # 4. CLR transformation with optimal pseudocount
        if apply_clr:
            if verbose:
                logger.info("  🔬 Step 4: CLR transformation with optimal pseudocount...")
            processed_data, clr_stats = self._optimal_clr_transform(
                processed_data, verbose
            )
            stats_dict.update(clr_stats)
        
        # 5. Variance stabilization using entropy
        if verbose:
            logger.info("  📊 Step 5: Entropy-based variance stabilization...")
        processed_data, var_stats = self._entropy_variance_stabilization(
            processed_data, verbose
        )
        stats_dict.update(var_stats)
        
        stats_dict['final_shape'] = processed_data.shape
        stats_dict['features_retained'] = processed_data.shape[1] / count_table.shape[1]
        stats_dict['samples_retained'] = processed_data.shape[0] / count_table.shape[0]
        
        if verbose:
            logger.info(f"✅ Preprocessing complete!")
            logger.info(f"📊 Output: {processed_data.shape[0]} samples × {processed_data.shape[1]} features")
            logger.info(f"📈 Features retained: {stats_dict['features_retained']:.1%}")
            logger.info(f"👥 Samples retained: {stats_dict['samples_retained']:.1%}")
        
        self.preprocessing_stats = stats_dict
        return processed_data, stats_dict
    
    def _prevalence_filter(self, data: pd.DataFrame, min_prevalence: float, verbose: bool) -> Tuple[pd.DataFrame, Dict]:
        """Standard prevalence-based filtering"""
        prevalence = (data > 0).mean()
        features_to_keep = prevalence >= min_prevalence
        
        filtered_data = data.loc[:, features_to_keep]
        
        stats = {
            'prevalence_threshold': min_prevalence,
            'features_before_prevalence': data.shape[1],
            'features_after_prevalence': filtered_data.shape[1],
            'features_removed_prevalence': data.shape[1] - filtered_data.shape[1]
        }
        
        if verbose:
            logger.info(f"    Removed {stats['features_removed_prevalence']} low-prevalence features")
        
        return filtered_data, stats
    
    def _information_content_filter(self, data: pd.DataFrame, min_info_content: float, verbose: bool) -> Tuple[pd.DataFrame, Dict]:
        """Information theory-based feature filtering"""
        
        # Convert to relative abundances for entropy calculation
        compositions = data.div(data.sum(axis=1), axis=0).fillna(0)
        
        # Calculate information content for each feature
        feature_info_content = []
        
        for feature in data.columns:
            feature_data = compositions[feature].values
            
            # Calculate Shannon entropy
            # Higher entropy = more informative (not concentrated in few samples)
            non_zero_mask = feature_data > 0
            if np.sum(non_zero_mask) <= 1:
                entropy = 0
            else:
                # Discretize for entropy calculation
                feature_nonzero = feature_data[non_zero_mask]
                n_bins = min(10, len(feature_nonzero) // 3)
                if n_bins > 1:
                    counts, _ = np.histogram(feature_nonzero, bins=n_bins)
                    counts = counts[counts > 0]
                    if len(counts) > 0:
                        probs = counts / np.sum(counts)
                        entropy = -np.sum(probs * np.log2(probs))
                    else:
                        entropy = 0
                else:
                    entropy = 0
            
            feature_info_content.append(entropy)
        
        feature_info_content = np.array(feature_info_content)
        
        # Normalize information content to [0,1]
        if np.max(feature_info_content) > 0:
            feature_info_content = feature_info_content / np.max(feature_info_content)
        
        # Filter features based on information content
        info_mask = feature_info_content >= min_info_content
        filtered_data = data.loc[:, info_mask]
        
        stats = {
            'info_content_threshold': min_info_content,
            'features_before_info': data.shape[1],
            'features_after_info': filtered_data.shape[1],
            'features_removed_info': data.shape[1] - filtered_data.shape[1],
            'mean_info_content': np.mean(feature_info_content),
            'median_info_content': np.median(feature_info_content)
        }
        
        if verbose:
            logger.info(f"    Removed {stats['features_removed_info']} low-information features")
            logger.info(f"    Mean information content: {stats['mean_info_content']:.3f}")
        
        return filtered_data, stats
    
    def _detect_information_outliers(self, data: pd.DataFrame, metadata: pd.DataFrame, verbose: bool) -> Tuple[pd.DataFrame, Dict]:
        """Detect outlier samples using information divergence"""
        
        # Convert to compositions
        compositions = data.div(data.sum(axis=1), axis=0).fillna(0)
        
        # Calculate pairwise Jensen-Shannon divergences
        n_samples = len(compositions)
        js_distances = np.zeros((n_samples, n_samples))
        
        for i in range(n_samples):
            for j in range(i+1, n_samples):
                js_dist = self._jensen_shannon_divergence(
                    compositions.iloc[i].values, 
                    compositions.iloc[j].values
                )
                js_distances[i, j] = js_dist
                js_distances[j, i] = js_dist
        
        # Calculate median distance to all other samples for each sample
        median_distances = np.median(js_distances, axis=1)
        
        # Identify outliers using IQR method on median distances
        q75, q25 = np.percentile(median_distances, [75, 25])
        iqr = q75 - q25
        outlier_threshold = q75 + 1.5 * iqr
        
        outlier_mask = median_distances > outlier_threshold
        outlier_samples = data.index[outlier_mask]
        
        # Remove outliers
        filtered_data = data.loc[~outlier_mask]
        
        stats = {
            'outlier_threshold': outlier_threshold,
            'samples_before_outlier': data.shape[0],
            'samples_after_outlier': filtered_data.shape[0],
            'outliers_removed': np.sum(outlier_mask),
            'outlier_samples': outlier_samples.tolist(),
            'median_js_distance': np.median(median_distances)
        }
        
        if verbose:
            logger.info(f"    Removed {stats['outliers_removed']} outlier samples")
            if stats['outliers_removed'] > 0:
                logger.info(f"    Outlier samples: {', '.join(outlier_samples[:5])}{'...' if len(outlier_samples) > 5 else ''}")
        
        return filtered_data, stats
    
    def _optimal_clr_transform(self, data: pd.DataFrame, verbose: bool) -> Tuple[pd.DataFrame, Dict]:
        """CLR transformation with information theory-optimized pseudocount"""
        
        # Test different pseudocounts and select the one that maximizes information content
        pseudocounts = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0]
        
        best_pseudocount = self.regularization
        best_info_score = -np.inf
        
        for pc in pseudocounts:
            # Apply pseudocount
            data_pseudo = data + pc
            
            # CLR transformation
            compositions = data_pseudo.div(data_pseudo.sum(axis=1), axis=0)
            log_comp = np.log(compositions)
            geometric_mean = log_comp.mean(axis=1)
            clr_data = log_comp.subtract(geometric_mean, axis=0)
            
            # Calculate information score (variance of CLR data)
            # Higher variance generally indicates better separation
            info_score = np.mean(np.var(clr_data, axis=0))
            
            if info_score > best_info_score:
                best_info_score = info_score
                best_pseudocount = pc
        
        # Apply optimal CLR transformation
        data_pseudo = data + best_pseudocount
        compositions = data_pseudo.div(data_pseudo.sum(axis=1), axis=0)
        log_comp = np.log(compositions)
        geometric_mean = log_comp.mean(axis=1)
        clr_data = log_comp.subtract(geometric_mean, axis=0)
        
        stats = {
            'optimal_pseudocount': best_pseudocount,
            'clr_info_score': best_info_score,
            'clr_applied': True
        }
        
        if verbose:
            logger.info(f"    Optimal pseudocount: {best_pseudocount:.1e}")
            logger.info(f"    CLR information score: {best_info_score:.3f}")
        
        return clr_data, stats
    
    def _entropy_variance_stabilization(self, data: pd.DataFrame, verbose: bool) -> Tuple[pd.DataFrame, Dict]:
        """Variance stabilization using entropy-based weighting"""
        
        # Calculate feature-wise entropy as measure of information content
        feature_entropies = []
        
        for feature in data.columns:
            feature_data = data[feature].values
            
            # Calculate entropy of feature distribution
            if np.var(feature_data) > 0:
                # Discretize continuous data for entropy
                n_bins = min(10, len(feature_data) // 5)
                if n_bins > 1:
                    counts, _ = np.histogram(feature_data, bins=n_bins)
                    counts = counts[counts > 0]
                    if len(counts) > 0:
                        probs = counts / np.sum(counts)
                        entropy = -np.sum(probs * np.log2(probs))
                    else:
                        entropy = 0
                else:
                    entropy = 0
            else:
                entropy = 0
            
            feature_entropies.append(entropy)
        
        feature_entropies = np.array(feature_entropies)
        
        # Apply entropy-based variance stabilization
        # Features with higher entropy get less aggressive transformation
        stabilized_data = data.copy()
        
        for i, feature in enumerate(data.columns):
            feature_data = data[feature].values
            entropy = feature_entropies[i]
            
            # Adaptive transformation based on entropy
            if entropy > 0:
                # Higher entropy = less transformation needed
                transform_strength = 1.0 / (1.0 + entropy)
                
                # Apply variance-stabilizing transformation
                if np.all(feature_data >= 0):
                    # For non-negative data, use sqrt transformation scaled by entropy
                    stabilized_data[feature] = np.sign(feature_data) * np.power(
                        np.abs(feature_data), 0.5 + 0.5 * transform_strength
                    )
                else:
                    # For data that can be negative (e.g., CLR), use scaled transformation
                    stabilized_data[feature] = feature_data * (1.0 + 0.1 * transform_strength)
        
        stats = {
            'variance_stabilization_applied': True,
            'mean_feature_entropy': np.mean(feature_entropies),
            'entropy_range': [np.min(feature_entropies), np.max(feature_entropies)]
        }
        
        if verbose:
            logger.info(f"    Mean feature entropy: {stats['mean_feature_entropy']:.3f}")
            logger.info(f"    Entropy range: [{stats['entropy_range'][0]:.3f}, {stats['entropy_range'][1]:.3f}]")
        
        return stabilized_data, stats
    
    def _jensen_shannon_divergence(self, p: np.ndarray, q: np.ndarray) -> float:
        """Calculate Jensen-Shannon divergence between two probability distributions"""
        # Ensure they sum to 1
        p = p / (np.sum(p) + self.regularization)
        q = q / (np.sum(q) + self.regularization)
        
        # Add regularization for numerical stability
        p = p + self.regularization
        q = q + self.regularization
        
        # Re-normalize
        p = p / np.sum(p)
        q = q / np.sum(q)
        
        # Average distribution
        m = 0.5 * (p + q)
        
        # Calculate KL divergences
        kl_pm = np.sum(p * np.log(p / m))
        kl_qm = np.sum(q * np.log(q / m))
        
        # Jensen-Shannon divergence
        js = 0.5 * kl_pm + 0.5 * kl_qm
        
        return js


def apply_information_preprocessing(count_table: pd.DataFrame,
                                  metadata: pd.DataFrame,
                                  min_prevalence: float = 0.1,
                                  min_information_content: float = 0.1,
                                  apply_clr: bool = True,
                                  remove_outliers: bool = True,
                                  verbose: bool = True) -> Tuple[pd.DataFrame, Dict]:
    """
    Convenience function for applying information theory-based preprocessing
    
    Parameters:
    -----------
    count_table : pd.DataFrame
        Raw count data (samples x features)
    metadata : pd.DataFrame
        Sample metadata
    min_prevalence : float
        Minimum prevalence threshold (0-1)
    min_information_content : float
        Minimum information content threshold (0-1)
    apply_clr : bool
        Apply CLR transformation
    remove_outliers : bool
        Remove outlier samples
    verbose : bool
        Print processing steps
        
    Returns:
    --------
    Tuple[pd.DataFrame, Dict]
        Preprocessed data and processing statistics
    """
    
    preprocessor = InformationBasedPreprocessor()
    return preprocessor.preprocess_data(
        count_table=count_table,
        metadata=metadata,
        min_prevalence=min_prevalence,
        min_information_content=min_information_content,
        apply_clr=apply_clr,
        remove_outliers=remove_outliers,
        verbose=verbose
    )