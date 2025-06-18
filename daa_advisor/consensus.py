#!/usr/bin/env python3
"""
Advanced Consensus Analysis Module for DAAadvisor

This module provides sophisticated consensus analysis capabilities including:
- Advanced majority voting with method reliability weighting
- Inter-method agreement quantification
- Weighted confidence scoring based on method concordance
- Uncertainty quantification and confidence intervals
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import logging
from scipy import stats
from sklearn.metrics import cohen_kappa_score
import warnings

logger = logging.getLogger(__name__)

class AdvancedConsensusAnalyzer:
    """
    Advanced consensus analysis for differential abundance results
    
    Features:
    - Weighted majority voting based on method reliability
    - Inter-method agreement quantification (Cohen's kappa, correlation)
    - Confidence scoring based on method concordance
    - Uncertainty quantification with confidence intervals
    - Multiple voting strategies (simple, weighted, ranked)
    """
    
    def __init__(self, 
                 voting_strategy: str = 'weighted',
                 confidence_threshold: float = 0.7,
                 significance_threshold: float = 0.05):
        """
        Initialize advanced consensus analyzer
        
        Parameters:
        -----------
        voting_strategy : str
            Voting strategy ('simple', 'weighted', 'ranked')
        confidence_threshold : float
            Minimum confidence for high-confidence calls
        significance_threshold : float
            P-value threshold for significance
        """
        self.voting_strategy = voting_strategy
        self.confidence_threshold = confidence_threshold
        self.significance_threshold = significance_threshold
        
        # Method reliability weights (can be learned from benchmarking)
        self.method_weights = {
            'wilcoxon': 0.8,      # High reliability, conservative
            'deseq2': 0.9,        # High statistical power
            'edger': 0.85,        # Good balance
            'aldex2': 0.75,       # Compositional awareness
            'ancom-bc': 0.7,      # Very conservative
            'metagenomeseq': 0.8  # Good for zero-inflation
        }
        
    def generate_advanced_consensus(self, 
                                  analyses: Dict[str, pd.DataFrame],
                                  include_agreement_metrics: bool = True,
                                  include_confidence_scores: bool = True) -> Dict[str, Any]:
        """
        Generate advanced consensus analysis
        
        Parameters:
        -----------
        analyses : Dict[str, pd.DataFrame]
            Dictionary of method results
        include_agreement_metrics : bool
            Calculate inter-method agreement metrics
        include_confidence_scores : bool
            Calculate confidence scores
            
        Returns:
        --------
        Dict[str, Any]
            Advanced consensus results with metrics
        """
        
        logger.info(f"🤝 Generating advanced consensus from {len(analyses)} methods")
        logger.info(f"📊 Using {self.voting_strategy} voting strategy")
        
        # 1. Align features across all methods
        aligned_results = self._align_method_results(analyses)
        
        # 2. Generate consensus features
        if self.voting_strategy == 'simple':
            consensus_features = self._simple_majority_voting(aligned_results)
        elif self.voting_strategy == 'weighted':
            consensus_features = self._weighted_majority_voting(aligned_results)
        elif self.voting_strategy == 'ranked':
            consensus_features = self._ranked_voting(aligned_results)
        else:
            raise ValueError(f"Unknown voting strategy: {self.voting_strategy}")
        
        # 3. Calculate agreement metrics
        agreement_metrics = {}
        if include_agreement_metrics:
            agreement_metrics = self._calculate_agreement_metrics(aligned_results)
        
        # 4. Calculate confidence scores
        if include_confidence_scores:
            consensus_features = self._add_confidence_scores(consensus_features, aligned_results)
        
        # 5. Add uncertainty quantification
        consensus_features = self._add_uncertainty_quantification(consensus_features, aligned_results)
        
        # 6. Classify consensus strength
        consensus_features = self._classify_consensus_strength(consensus_features)
        
        logger.info(f"✅ Consensus analysis complete")
        logger.info(f"📈 {(consensus_features['consensus_significant']).sum()} consensus significant features")
        
        return {
            'consensus_features': consensus_features,
            'agreement_metrics': agreement_metrics,
            'method_weights': self.method_weights,
            'voting_strategy': self.voting_strategy,
            'n_methods': len(analyses),
            'method_names': list(analyses.keys())
        }
    
    def _align_method_results(self, analyses: Dict[str, pd.DataFrame]) -> pd.DataFrame:
        """Align results from different methods into a unified DataFrame"""
        
        # Get all unique features
        all_features = set()
        for results in analyses.values():
            all_features.update(results['feature'].values)
        
        all_features = sorted(list(all_features))
        
        # Create aligned results DataFrame
        aligned_data = []
        
        for feature in all_features:
            feature_data = {'feature': feature}
            
            # Extract data from each method
            for method_name, results in analyses.items():
                method_data = results[results['feature'] == feature]
                
                if not method_data.empty:
                    # Get significance status
                    padj_col = self._get_padj_column(method_data)
                    pvalue = method_data['pvalue'].iloc[0] if 'pvalue' in method_data.columns else np.nan
                    padj = method_data[padj_col].iloc[0] if padj_col else np.nan
                    is_significant = padj < self.significance_threshold if not pd.isna(padj) else False
                    
                    # Get effect size
                    effect_size = 0
                    if 'log2fc' in method_data.columns:
                        effect_size = method_data['log2fc'].iloc[0]
                    elif 'effect_size' in method_data.columns:
                        effect_size = method_data['effect_size'].iloc[0]
                    elif 'statistic' in method_data.columns:
                        effect_size = method_data['statistic'].iloc[0]
                    
                    feature_data.update({
                        f'{method_name}_pvalue': pvalue,
                        f'{method_name}_padj': padj,
                        f'{method_name}_significant': is_significant,
                        f'{method_name}_effect_size': effect_size
                    })
                else:
                    # Method did not test this feature
                    feature_data.update({
                        f'{method_name}_pvalue': np.nan,
                        f'{method_name}_padj': np.nan,
                        f'{method_name}_significant': False,
                        f'{method_name}_effect_size': np.nan
                    })
            
            aligned_data.append(feature_data)
        
        return pd.DataFrame(aligned_data)
    
    def _get_padj_column(self, results: pd.DataFrame) -> Optional[str]:
        """Get the adjusted p-value column name"""
        for col in ['padj', 'qvalue', 'adj_pvalue', 'padj_adaptive']:
            if col in results.columns:
                return col
        return None
    
    def _simple_majority_voting(self, aligned_results: pd.DataFrame) -> pd.DataFrame:
        """Simple majority voting (≥50% of methods)"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        consensus_data = []
        
        for _, row in aligned_results.iterrows():
            feature = row['feature']
            
            # Count significant calls
            significant_calls = [row[f'{method}_significant'] for method in method_names]
            n_significant = sum(significant_calls)
            n_methods = len(method_names)
            
            # Majority voting
            consensus_significant = n_significant >= (n_methods / 2)
            
            # Calculate mean statistics
            pvalues = [row[f'{method}_pvalue'] for method in method_names if not pd.isna(row[f'{method}_pvalue'])]
            effect_sizes = [row[f'{method}_effect_size'] for method in method_names if not pd.isna(row[f'{method}_effect_size'])]
            
            consensus_data.append({
                'feature': feature,
                'n_significant': n_significant,
                'n_methods': n_methods,
                'proportion_significant': n_significant / n_methods,
                'consensus_significant': consensus_significant,
                'mean_pvalue': np.mean(pvalues) if pvalues else np.nan,
                'median_pvalue': np.median(pvalues) if pvalues else np.nan,
                'mean_effect_size': np.mean(effect_sizes) if effect_sizes else np.nan,
                'median_effect_size': np.median(effect_sizes) if effect_sizes else np.nan
            })
        
        return pd.DataFrame(consensus_data).sort_values('n_significant', ascending=False)
    
    def _weighted_majority_voting(self, aligned_results: pd.DataFrame) -> pd.DataFrame:
        """Weighted majority voting based on method reliability"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        consensus_data = []
        
        for _, row in aligned_results.iterrows():
            feature = row['feature']
            
            # Calculate weighted votes
            weighted_votes = 0
            total_weight = 0
            
            for method in method_names:
                if row[f'{method}_significant']:
                    weight = self.method_weights.get(method, 0.5)  # Default weight for unknown methods
                    weighted_votes += weight
                total_weight += self.method_weights.get(method, 0.5)
            
            # Weighted consensus
            weighted_proportion = weighted_votes / total_weight if total_weight > 0 else 0
            consensus_significant = weighted_proportion >= 0.5
            
            # Count simple votes for comparison
            significant_calls = [row[f'{method}_significant'] for method in method_names]
            n_significant = sum(significant_calls)
            
            # Calculate weighted mean statistics
            weighted_pvalues = []
            weighted_effect_sizes = []
            
            for method in method_names:
                if not pd.isna(row[f'{method}_pvalue']):
                    weight = self.method_weights.get(method, 0.5)
                    weighted_pvalues.append((row[f'{method}_pvalue'], weight))
                if not pd.isna(row[f'{method}_effect_size']):
                    weight = self.method_weights.get(method, 0.5)
                    weighted_effect_sizes.append((row[f'{method}_effect_size'], weight))
            
            # Calculate weighted means
            if weighted_pvalues:
                weighted_mean_pvalue = sum(p * w for p, w in weighted_pvalues) / sum(w for _, w in weighted_pvalues)
            else:
                weighted_mean_pvalue = np.nan
            
            if weighted_effect_sizes:
                weighted_mean_effect = sum(e * w for e, w in weighted_effect_sizes) / sum(w for _, w in weighted_effect_sizes)
            else:
                weighted_mean_effect = np.nan
            
            consensus_data.append({
                'feature': feature,
                'n_significant': n_significant,
                'n_methods': len(method_names),
                'proportion_significant': n_significant / len(method_names),
                'weighted_votes': weighted_votes,
                'total_weight': total_weight,
                'weighted_proportion': weighted_proportion,
                'consensus_significant': consensus_significant,
                'weighted_mean_pvalue': weighted_mean_pvalue,
                'weighted_mean_effect_size': weighted_mean_effect
            })
        
        return pd.DataFrame(consensus_data).sort_values('weighted_votes', ascending=False)
    
    def _ranked_voting(self, aligned_results: pd.DataFrame) -> pd.DataFrame:
        """Ranked voting based on effect sizes and p-values"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        consensus_data = []
        
        for _, row in aligned_results.iterrows():
            feature = row['feature']
            
            # Collect method rankings
            method_scores = []
            
            for method in method_names:
                pvalue = row[f'{method}_pvalue']
                effect_size = abs(row[f'{method}_effect_size']) if not pd.isna(row[f'{method}_effect_size']) else 0
                
                if not pd.isna(pvalue):
                    # Combined score: -log10(p-value) + |effect_size|
                    score = -np.log10(pvalue + 1e-10) + effect_size
                    method_scores.append((method, score, row[f'{method}_significant']))
            
            # Sort by scores
            method_scores.sort(key=lambda x: x[1], reverse=True)
            
            # Ranked voting: top methods get more weight
            ranked_votes = 0
            for i, (method, score, is_significant) in enumerate(method_scores):
                if is_significant:
                    # Higher rank = more weight
                    rank_weight = len(method_scores) - i
                    ranked_votes += rank_weight
            
            max_possible_votes = sum(range(1, len(method_scores) + 1))
            ranked_proportion = ranked_votes / max_possible_votes if max_possible_votes > 0 else 0
            
            consensus_significant = ranked_proportion >= 0.5
            
            # Simple counts for comparison
            n_significant = sum(1 for _, _, is_sig in method_scores if is_sig)
            
            consensus_data.append({
                'feature': feature,
                'n_significant': n_significant,
                'n_methods': len(method_names),
                'ranked_votes': ranked_votes,
                'max_possible_votes': max_possible_votes,
                'ranked_proportion': ranked_proportion,
                'consensus_significant': consensus_significant,
                'top_method': method_scores[0][0] if method_scores else None,
                'top_score': method_scores[0][1] if method_scores else 0
            })
        
        return pd.DataFrame(consensus_data).sort_values('ranked_votes', ascending=False)
    
    def _calculate_agreement_metrics(self, aligned_results: pd.DataFrame) -> Dict[str, Any]:
        """Calculate inter-method agreement metrics"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        if len(method_names) < 2:
            return {'error': 'Need at least 2 methods for agreement calculation'}
        
        # Get significance calls matrix
        sig_matrix = aligned_results[[f'{method}_significant' for method in method_names]].values.astype(int)
        
        # 1. Pairwise Cohen's Kappa
        pairwise_kappa = {}
        kappa_values = []
        
        for i, method1 in enumerate(method_names):
            for j, method2 in enumerate(method_names):
                if i < j:  # Avoid duplicates
                    try:
                        kappa = cohen_kappa_score(sig_matrix[:, i], sig_matrix[:, j])
                        pairwise_kappa[f'{method1}_vs_{method2}'] = kappa
                        kappa_values.append(kappa)
                    except:
                        pairwise_kappa[f'{method1}_vs_{method2}'] = np.nan
        
        # 2. Overall agreement metrics
        mean_kappa = np.nanmean(kappa_values) if kappa_values else np.nan
        
        # 3. Proportion of features with unanimous agreement
        unanimous_agree = np.sum(np.all(sig_matrix == sig_matrix[:, 0:1], axis=1)) / len(sig_matrix)
        
        # 4. Proportion of features with majority agreement
        majority_agree = np.sum(np.sum(sig_matrix, axis=1) >= len(method_names)/2) / len(sig_matrix)
        
        # 5. Method-specific agreement (how often each method agrees with consensus)
        method_agreement = {}
        for i, method in enumerate(method_names):
            consensus_calls = np.sum(sig_matrix, axis=1) >= len(method_names)/2
            method_calls = sig_matrix[:, i].astype(bool)
            agreement = np.mean(consensus_calls == method_calls)
            method_agreement[method] = agreement
        
        # 6. Effect size correlations
        effect_size_correlations = {}
        effect_size_cols = [f'{method}_effect_size' for method in method_names]
        
        for i, method1 in enumerate(method_names):
            for j, method2 in enumerate(method_names):
                if i < j:
                    col1 = effect_size_cols[i]
                    col2 = effect_size_cols[j]
                    
                    # Get valid (non-NaN) pairs
                    valid_mask = ~(pd.isna(aligned_results[col1]) | pd.isna(aligned_results[col2]))
                    
                    if valid_mask.sum() > 5:  # Need at least 5 valid pairs
                        correlation = stats.spearmanr(
                            aligned_results.loc[valid_mask, col1],
                            aligned_results.loc[valid_mask, col2]
                        )[0]
                        effect_size_correlations[f'{method1}_vs_{method2}'] = correlation
        
        return {
            'pairwise_kappa': pairwise_kappa,
            'mean_kappa': mean_kappa,
            'kappa_interpretation': self._interpret_kappa(mean_kappa),
            'unanimous_agreement': unanimous_agree,
            'majority_agreement': majority_agree,
            'method_agreement': method_agreement,
            'effect_size_correlations': effect_size_correlations,
            'mean_effect_correlation': np.nanmean(list(effect_size_correlations.values())) if effect_size_correlations else np.nan
        }
    
    def _interpret_kappa(self, kappa: float) -> str:
        """Interpret Cohen's kappa value"""
        if pd.isna(kappa):
            return "Cannot calculate"
        elif kappa < 0:
            return "Poor agreement"
        elif kappa < 0.20:
            return "Slight agreement"
        elif kappa < 0.40:
            return "Fair agreement"
        elif kappa < 0.60:
            return "Moderate agreement"
        elif kappa < 0.80:
            return "Substantial agreement"
        else:
            return "Almost perfect agreement"
    
    def _add_confidence_scores(self, consensus_features: pd.DataFrame, aligned_results: pd.DataFrame) -> pd.DataFrame:
        """Add confidence scores based on method agreement"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        confidence_scores = []
        
        for _, row in consensus_features.iterrows():
            feature = row['feature']
            
            # Get feature data from aligned results
            feature_data = aligned_results[aligned_results['feature'] == feature].iloc[0]
            
            # Confidence factors
            confidence_factors = []
            
            # 1. Proportion of methods agreeing
            if 'proportion_significant' in row:
                agreement_factor = min(row['proportion_significant'], 1 - row['proportion_significant']) * 2
                confidence_factors.append(agreement_factor)
            
            # 2. P-value consistency (lower variance = higher confidence)
            pvalues = [feature_data[f'{method}_pvalue'] for method in method_names if not pd.isna(feature_data[f'{method}_pvalue'])]
            if len(pvalues) > 1:
                pvalue_cv = np.std(pvalues) / np.mean(pvalues) if np.mean(pvalues) > 0 else 1
                pvalue_consistency = 1 / (1 + pvalue_cv)  # Higher consistency = lower CV
                confidence_factors.append(pvalue_consistency)
            
            # 3. Effect size consistency
            effect_sizes = [feature_data[f'{method}_effect_size'] for method in method_names if not pd.isna(feature_data[f'{method}_effect_size'])]
            if len(effect_sizes) > 1:
                effect_cv = np.std(effect_sizes) / max(np.mean(np.abs(effect_sizes)), 0.1)
                effect_consistency = 1 / (1 + effect_cv)
                confidence_factors.append(effect_consistency)
            
            # 4. Method reliability weighting
            if row['consensus_significant']:
                significant_methods = [method for method in method_names if feature_data[f'{method}_significant']]
                avg_reliability = np.mean([self.method_weights.get(method, 0.5) for method in significant_methods])
                confidence_factors.append(avg_reliability)
            
            # Combined confidence score
            confidence_score = np.mean(confidence_factors) if confidence_factors else 0.5
            confidence_scores.append(confidence_score)
        
        consensus_features['confidence_score'] = confidence_scores
        consensus_features['high_confidence'] = consensus_features['confidence_score'] >= self.confidence_threshold
        
        return consensus_features
    
    def _add_uncertainty_quantification(self, consensus_features: pd.DataFrame, aligned_results: pd.DataFrame) -> pd.DataFrame:
        """Add uncertainty quantification measures"""
        
        method_names = [col.replace('_significant', '') for col in aligned_results.columns if col.endswith('_significant')]
        
        uncertainty_measures = []
        
        for _, row in consensus_features.iterrows():
            feature = row['feature']
            feature_data = aligned_results[aligned_results['feature'] == feature].iloc[0]
            
            # P-value uncertainty (confidence interval)
            pvalues = [feature_data[f'{method}_pvalue'] for method in method_names if not pd.isna(feature_data[f'{method}_pvalue'])]
            
            if len(pvalues) > 1:
                pvalue_ci = np.percentile(pvalues, [25, 75])
                pvalue_uncertainty = pvalue_ci[1] - pvalue_ci[0]
            else:
                pvalue_uncertainty = np.nan
            
            # Effect size uncertainty
            effect_sizes = [feature_data[f'{method}_effect_size'] for method in method_names if not pd.isna(feature_data[f'{method}_effect_size'])]
            
            if len(effect_sizes) > 1:
                effect_ci = np.percentile(effect_sizes, [25, 75])
                effect_uncertainty = effect_ci[1] - effect_ci[0]
            else:
                effect_uncertainty = np.nan
            
            uncertainty_measures.append({
                'pvalue_uncertainty': pvalue_uncertainty,
                'effect_uncertainty': effect_uncertainty
            })
        
        uncertainty_df = pd.DataFrame(uncertainty_measures)
        consensus_features = pd.concat([consensus_features, uncertainty_df], axis=1)
        
        return consensus_features
    
    def _classify_consensus_strength(self, consensus_features: pd.DataFrame) -> pd.DataFrame:
        """Classify consensus strength categories"""
        
        def classify_strength(row):
            if not row['consensus_significant']:
                return 'Not Significant'
            
            confidence = row.get('confidence_score', 0.5)
            proportion = row.get('proportion_significant', 0)
            
            if confidence >= 0.8 and proportion >= 0.8:
                return 'Strong Consensus'
            elif confidence >= 0.6 and proportion >= 0.6:
                return 'Moderate Consensus'
            elif confidence >= 0.4 or proportion >= 0.5:
                return 'Weak Consensus'
            else:
                return 'Conflicting Evidence'
        
        consensus_features['consensus_strength'] = consensus_features.apply(classify_strength, axis=1)
        
        return consensus_features


def generate_advanced_consensus(analyses: Dict[str, pd.DataFrame],
                              voting_strategy: str = 'weighted',
                              confidence_threshold: float = 0.7,
                              significance_threshold: float = 0.05) -> Dict[str, Any]:
    """
    Convenience function for generating advanced consensus analysis
    
    Parameters:
    -----------
    analyses : Dict[str, pd.DataFrame]
        Dictionary of method results
    voting_strategy : str
        Voting strategy ('simple', 'weighted', 'ranked')
    confidence_threshold : float
        Minimum confidence for high-confidence calls
    significance_threshold : float
        P-value threshold for significance
        
    Returns:
    --------
    Dict[str, Any]
        Advanced consensus results
    """
    
    analyzer = AdvancedConsensusAnalyzer(
        voting_strategy=voting_strategy,
        confidence_threshold=confidence_threshold,
        significance_threshold=significance_threshold
    )
    
    return analyzer.generate_advanced_consensus(
        analyses=analyses,
        include_agreement_metrics=True,
        include_confidence_scores=True
    )