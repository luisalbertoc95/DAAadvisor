#!/usr/bin/env python3
"""
Method selection module for DAAadvisor
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple
from dataclasses import dataclass
import logging

from .profiler import DataProfile

logger = logging.getLogger(__name__)


@dataclass
class MethodRecommendation:
    """Method recommendation with reasoning"""
    primary_method: str
    secondary_methods: List[str]
    confidence: float
    reasoning: List[str]
    scores: Dict[str, float]


class MethodSelector:
    """Recommends best methods based on data profile"""
    
    def __init__(self):
        self.method_performance = self._load_method_performance()
        self.decision_rules = self._load_decision_rules()
    
    def _load_method_performance(self) -> Dict:
        """Load method performance characteristics from benchmarking studies"""
        return {
            'aldex2': {
                'good_for': ['asv', 'high_compositional_bias', 'medium_sample_size'],
                'bad_for': ['extreme_sparsity', 'very_small_sample_size'],
                'fdr_control': 'excellent',
                'power': 'moderate',
                'min_samples': 8,
                'max_sparsity': 0.95,
                'handles_compositionality': True,
                'robust_to_outliers': True
            },
            'ancom-bc': {
                'good_for': ['compositional_bias', 'asv', 'gene', 'balanced_design'],
                'bad_for': ['very_small_sample_size', 'extreme_unbalanced'],
                'fdr_control': 'good',
                'power': 'moderate',
                'min_samples': 10,
                'max_sparsity': 0.90,
                'handles_compositionality': True,
                'robust_to_outliers': False
            },
            'deseq2': {
                'good_for': ['gene', 'low_sparsity', 'small_sample_size', 'complex_design'],
                'bad_for': ['high_compositional_bias', 'extreme_sparsity'],
                'fdr_control': 'moderate',
                'power': 'high',
                'min_samples': 6,
                'max_sparsity': 0.70,
                'handles_compositionality': False,
                'robust_to_outliers': False
            },
            'edger': {
                'good_for': ['gene', 'low_sparsity', 'large_sample_size'],
                'bad_for': ['high_sparsity', 'compositional_bias', 'small_sample_size'],
                'fdr_control': 'poor',
                'power': 'very_high',
                'min_samples': 10,
                'max_sparsity': 0.60,
                'handles_compositionality': False,
                'robust_to_outliers': False
            },
            'metagenomeseq': {
                'good_for': ['moderate_sparsity', 'gene', 'asv', 'uneven_depth'],
                'bad_for': ['extreme_sparsity', 'very_small_sample_size'],
                'fdr_control': 'good',
                'power': 'moderate',
                'min_samples': 8,
                'max_sparsity': 0.85,
                'handles_compositionality': True,
                'robust_to_outliers': True
            },
            'zicoseq': {
                'good_for': ['high_sparsity', 'zero_inflation', 'viral', 'small_sample_size'],
                'bad_for': ['low_sparsity'],
                'fdr_control': 'excellent',
                'power': 'high',
                'min_samples': 5,
                'max_sparsity': 0.99,
                'handles_compositionality': True,
                'robust_to_outliers': True
            },
            'linda': {
                'good_for': ['asv', 'gene', 'complex_design', 'confounders'],
                'bad_for': ['extreme_sparsity', 'very_small_sample_size'],
                'fdr_control': 'good',
                'power': 'moderate',
                'min_samples': 10,
                'max_sparsity': 0.80,
                'handles_compositionality': True,
                'robust_to_outliers': True
            },
            'maaslin3': {
                'good_for': ['asv', 'gene', 'complex_design', 'mixed_effects'],
                'bad_for': ['extreme_sparsity', 'small_sample_size'],
                'fdr_control': 'good',
                'power': 'moderate',
                'min_samples': 15,
                'max_sparsity': 0.75,
                'handles_compositionality': True,
                'robust_to_outliers': True
            },
            'wilcoxon': {
                'good_for': ['small_sample_size', 'non_parametric', 'robust'],
                'bad_for': ['complex_design', 'continuous_variables'],
                'fdr_control': 'moderate',
                'power': 'low',
                'min_samples': 3,
                'max_sparsity': 0.99,
                'handles_compositionality': False,
                'robust_to_outliers': True
            }
        }
    
    def _load_decision_rules(self) -> Dict:
        """Load decision tree rules for method selection"""inactive
        return {
            'viral': {
                'high_sparsity': 'zicoseq',
                'medium_sparsity_small_n': 'ancom-bc',
                'medium_sparsity_large_n': 'aldex2',
                'fallback': 'wilcoxon'
            },
            'asv': {
                'high_compositional_bias': 'aldex2',
                'medium_compositional_bias_balanced': 'ancom-bc',
                'low_sparsity': 'linda',
                'small_sample_size': 'aldex2',
                'fallback': 'aldex2'
            },
            'gene': {
                'low_sparsity_large_n': 'deseq2',
                'low_sparsity_small_n': 'metagenomeseq',
                'high_sparsity': 'zicoseq',
                'complex_design': 'linda',
                'fallback': 'deseq2'
            }
        }
    
    def recommend_methods(self, profile: DataProfile) -> MethodRecommendation:
        """Recommend methods based on data profile"""
        
        scores = {}
        reasoning = []
        
        # Score each method based on data characteristics
        for method, chars in self.method_performance.items():
            score, method_reasons = self._score_method(method, chars, profile)
            scores[method] = score
            
            if score > 0:  # Only add reasoning for viable methods
                reasoning.extend([f"{method}: {reason}" for reason in method_reasons])
        
        # Sort methods by score
        sorted_methods = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        
        # Filter out methods with negative scores
        viable_methods = [(method, score) for method, score in sorted_methods if score > 0]
        
        if not viable_methods:
            # Fallback to most robust method
            primary_method = 'wilcoxon'
            secondary_methods = ['aldex2']
            confidence = 0.3
            reasoning.append("No methods scored positively, using robust fallback")
        else:
            # Select primary and secondary methods
            primary_method = viable_methods[0][0]
            secondary_methods = [m[0] for m in viable_methods[1:4]]  # Top 3 alternatives
            
            # Calculate confidence based on score differences and absolute score
            top_score = viable_methods[0][1]
            second_score = viable_methods[1][1] if len(viable_methods) > 1 else 0
            
            score_gap = top_score - second_score
            confidence = min(0.95, 0.4 + (top_score * 0.1) + (score_gap * 0.05))
        
        # Add data-specific reasoning
        reasoning.insert(0, f"Data type: {profile.data_type}")
        reasoning.insert(1, f"Samples: {profile.n_samples}, Features: {profile.n_features}")
        reasoning.insert(2, f"Sparsity: {profile.sparsity:.1%}, Zero inflation: {profile.zero_inflation:.1%}")
        
        if profile.compositional_bias > 0.3:
            reasoning.append("High compositional bias detected")
        
        if profile.batch_effects:
            max_effect = max(profile.batch_effects.values())
            if max_effect > 0.3:
                reasoning.append(f"Potential batch effects detected (max: {max_effect:.2f})")
        
        return MethodRecommendation(
            primary_method=primary_method,
            secondary_methods=secondary_methods,
            confidence=confidence,
            reasoning=reasoning,
            scores=scores
        )
    
    def _score_method(self, method: str, characteristics: Dict, profile: DataProfile) -> Tuple[float, List[str]]:
        """Score a method based on data profile"""
        score = 0
        reasons = []
        
        # Basic compatibility checks
        min_samples = min(profile.group_sizes.values()) if profile.group_sizes else profile.n_samples
        
        # Sample size check
        if min_samples < characteristics.get('min_samples', 5):
            score -= 3
            reasons.append(f"insufficient samples (need ≥{characteristics['min_samples']})")
            return score, reasons  # Early return for critical failure
        else:
            score += 1
            reasons.append("adequate sample size")
        
        # Sparsity handling
        max_sparsity = characteristics.get('max_sparsity', 0.95)
        if profile.sparsity > max_sparsity:
            score -= 2
            reasons.append(f"sparsity too high ({profile.sparsity:.1%} > {max_sparsity:.1%})")
        elif profile.sparsity > 0.8 and 'high_sparsity' in characteristics['good_for']:
            score += 2
            reasons.append("handles high sparsity well")
        elif profile.sparsity < 0.5 and 'low_sparsity' in characteristics['good_for']:
            score += 2
            reasons.append("optimized for low sparsity")
        
        # Data type compatibility
        if profile.data_type in characteristics['good_for']:
            score += 3
            reasons.append(f"designed for {profile.data_type} data")
        elif profile.data_type in characteristics.get('bad_for', []):
            score -= 1
            reasons.append(f"not optimal for {profile.data_type} data")
        
        # Compositional bias handling
        if profile.compositional_bias > 0.3:
            if characteristics.get('handles_compositionality', False):
                score += 2
                reasons.append("addresses compositional bias")
            else:
                score -= 1
                reasons.append("may be affected by compositional bias")
        
        # Zero inflation
        if profile.zero_inflation > 0.6 and 'zero_inflation' in characteristics['good_for']:
            score += 1
            reasons.append("handles zero inflation")
        
        # Sample balance
        if profile.group_sizes:
            group_sizes = list(profile.group_sizes.values())
            balance_ratio = min(group_sizes) / max(group_sizes) if max(group_sizes) > 0 else 0
            
            if balance_ratio < 0.3:  # Highly unbalanced
                if 'extreme_unbalanced' in characteristics.get('bad_for', []):
                    score -= 1
                    reasons.append("sensitive to unbalanced design")
                elif characteristics.get('robust_to_outliers', False):
                    score += 1
                    reasons.append("robust to unbalanced design")
        
        # Batch effects
        if profile.batch_effects:
            max_batch_effect = max(profile.batch_effects.values())
            if max_batch_effect > 0.3:
                if 'confounders' in characteristics['good_for'] or 'complex_design' in characteristics['good_for']:
                    score += 1
                    reasons.append("can handle batch effects")
                else:
                    score -= 0.5
                    reasons.append("may be affected by batch effects")
        
        # Sequencing depth variation
        if profile.sequencing_depth_cv > 1.0:  # High variation
            if 'uneven_depth' in characteristics['good_for']:
                score += 1
                reasons.append("handles uneven sequencing depth")
        
        # Sample size categories
        if min_samples < 10:
            if 'small_sample_size' in characteristics['good_for']:
                score += 1
                reasons.append("good for small sample sizes")
            elif 'very_small_sample_size' in characteristics.get('bad_for', []):
                score -= 1
                reasons.append("prefers larger sample sizes")
        elif min_samples > 50:
            if 'large_sample_size' in characteristics['good_for']:
                score += 1
                reasons.append("leverages large sample size")
        
        # Method-specific bonuses
        fdr_quality = characteristics.get('fdr_control', 'moderate')
        if fdr_quality == 'excellent':
            score += 1
            reasons.append("excellent FDR control")
        elif fdr_quality == 'poor':
            score -= 0.5
            reasons.append("poor FDR control")
        
        return score, reasons
    
    def explain_recommendation(self, recommendation: MethodRecommendation) -> None:
        """Print detailed explanation of the recommendation"""
        
        print("\n" + "="*60)
        print("METHOD RECOMMENDATION EXPLANATION")
        print("="*60)
        
        print(f"\nPrimary Recommendation: {recommendation.primary_method.upper()}")
        print(f"Confidence: {recommendation.confidence:.1%}")
        
        print(f"\nAlternative Methods:")
        for i, method in enumerate(recommendation.secondary_methods, 1):
            score = recommendation.scores.get(method, 0)
            print(f"  {i}. {method} (score: {score:.1f})")
        
        print(f"\nReasoning:")
        for reason in recommendation.reasoning:
            print(f"  • {reason}")
        
        print(f"\nMethod Scores:")
        sorted_scores = sorted(recommendation.scores.items(), key=lambda x: x[1], reverse=True)
        for method, score in sorted_scores:
            if score > 0:
                print(f"  {method}: {score:.1f}")
            else:
                print(f"  {method}: {score:.1f} (not recommended)")
    
    def get_method_info(self, method: str) -> Dict:
        """Get detailed information about a specific method"""
        if method not in self.method_performance:
            return {}
        
        info = self.method_performance[method]
        return {
            'name': method,
            'good_for': info['good_for'],
            'bad_for': info.get('bad_for', []),
            'fdr_control': info['fdr_control'],
            'power': info['power'],
            'min_samples': info['min_samples'],
            'max_sparsity': info['max_sparsity'],
            'handles_compositionality': info['handles_compositionality'],
            'robust_to_outliers': info['robust_to_outliers']
        }