#!/usr/bin/env python3
"""
Visualization module for DAAadvisor
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Dict, List, Optional, Tuple
import logging
from pathlib import Path

from .profiler import DataProfile

logger = logging.getLogger(__name__)


class DAAVisualizer:
    """Comprehensive visualization for differential abundance analysis"""
    
    def __init__(self, style: str = 'seaborn-v0_8', figsize: Tuple[int, int] = (10, 8)):
        """
        Initialize visualizer
        
        Parameters:
        -----------
        style : str
            Matplotlib style to use
        figsize : tuple
            Default figure size
        """
        self.style = style
        self.figsize = figsize
        
        # Set matplotlib style
        try:
            plt.style.use(style)
        except:
            plt.style.use('default')
            logger.warning(f"Style {style} not available, using default")
        
        # Set seaborn palette
        sns.set_palette("husl")
    
    def plot_data_characteristics(self, 
                                profile: DataProfile, 
                                save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot comprehensive data characteristics
        
        Parameters:
        -----------
        profile : DataProfile
            Data profile object
        save_path : str, optional
            Path to save the plot
            
        Returns:
        --------
        matplotlib.Figure
        """
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'Data Characteristics Summary - {profile.data_type.upper()} Data', 
                     fontsize=16, fontweight='bold')
        
        # Basic statistics
        ax = axes[0, 0]
        stats_data = [
            profile.n_samples,
            profile.n_features,
            profile.sparsity * 100,
            profile.zero_inflation * 100,
            profile.compositional_bias * 100
        ]
        stats_labels = ['Samples', 'Features', 'Sparsity%', 'Zero Inflation%', 'Comp. Bias%']
        
        bars = ax.bar(stats_labels, stats_data)
        ax.set_title('Basic Statistics')
        ax.set_ylabel('Count / Percentage')
        
        # Color bars differently
        colors = ['skyblue', 'lightgreen', 'orange', 'lightcoral', 'plum']
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        
        # Add value labels on bars
        for bar, value in zip(bars, stats_data):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(stats_data)*0.01,
                   f'{value:.1f}' if value > 1 else f'{value:.3f}',
                   ha='center', va='bottom')
        
        # Group sizes
        ax = axes[0, 1]
        if profile.group_sizes:
            groups = list(profile.group_sizes.keys())
            sizes = list(profile.group_sizes.values())
            
            pie = ax.pie(sizes, labels=groups, autopct='%1.1f%%', startangle=90)
            ax.set_title('Group Distribution')
        else:
            ax.text(0.5, 0.5, 'No group information', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Group Distribution')
        
        # Prevalence distribution
        ax = axes[0, 2]
        if hasattr(profile, 'prevalence_distribution') and profile.prevalence_distribution:
            prev_categories = list(profile.prevalence_distribution.keys())
            prev_values = [v * 100 for v in profile.prevalence_distribution.values()]
            
            bars = ax.bar(prev_categories, prev_values)
            ax.set_title('Feature Prevalence Distribution')
            ax.set_ylabel('Percentage of Features')
            ax.tick_params(axis='x', rotation=45)
            
            # Color bars
            colors = sns.color_palette("viridis", len(bars))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
        else:
            ax.text(0.5, 0.5, 'Prevalence data not available', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Feature Prevalence Distribution')
        
        # Sample diversity
        ax = axes[1, 0]
        if hasattr(profile, 'sample_diversity') and profile.sample_diversity:
            diversity_metrics = ['shannon_mean', 'shannon_cv', 'richness_mean', 'richness_cv']
            diversity_values = [profile.sample_diversity.get(metric, 0) for metric in diversity_metrics]
            diversity_labels = ['Shannon\nMean', 'Shannon\nCV', 'Richness\nMean', 'Richness\nCV']
            
            bars = ax.bar(diversity_labels, diversity_values)
            ax.set_title('Sample Diversity Metrics')
            ax.set_ylabel('Value')
            
            for bar, value in zip(bars, diversity_values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + max(diversity_values)*0.01,
                       f'{value:.2f}', ha='center', va='bottom')
        else:
            ax.text(0.5, 0.5, 'Diversity data not available', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Sample Diversity Metrics')
        
        # Batch effects
        ax = axes[1, 1]
        if hasattr(profile, 'batch_effects') and profile.batch_effects:
            batch_factors = list(profile.batch_effects.keys())
            batch_effects = list(profile.batch_effects.values())
            
            bars = ax.bar(batch_factors, batch_effects)
            ax.set_title('Potential Batch Effects')
            ax.set_ylabel('Effect Size')
            ax.tick_params(axis='x', rotation=45)
            
            # Add threshold line
            ax.axhline(y=0.3, color='red', linestyle='--', alpha=0.7, label='Concern threshold')
            ax.legend()
        else:
            ax.text(0.5, 0.5, 'No batch effects detected', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Potential Batch Effects')
        
        # Data quality summary
        ax = axes[1, 2]
        quality_scores = self._calculate_quality_scores(profile)
        quality_labels = list(quality_scores.keys())
        quality_values = list(quality_scores.values())
        
        bars = ax.barh(quality_labels, quality_values)
        ax.set_title('Data Quality Assessment')
        ax.set_xlabel('Quality Score (0-1)')
        ax.set_xlim(0, 1)
        
        # Color bars by quality
        for bar, value in zip(bars, quality_values):
            if value >= 0.8:
                bar.set_color('green')
            elif value >= 0.6:
                bar.set_color('orange')
            else:
                bar.set_color('red')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Data characteristics plot saved to {save_path}")
        
        return fig
    
    def plot_volcano(self, 
                    results: pd.DataFrame, 
                    alpha: float = 0.05,
                    fc_threshold: float = 1.0,
                    title: str = "Volcano Plot",
                    save_path: Optional[str] = None) -> plt.Figure:
        """
        Create volcano plot of differential abundance results
        
        Parameters:
        -----------
        results : pd.DataFrame
            Results from differential abundance analysis
        alpha : float
            Significance threshold
        fc_threshold : float
            Fold change threshold (log2)
        title : str
            Plot title
        save_path : str, optional
            Path to save the plot
        """
        
        fig, ax = plt.subplots(figsize=self.figsize)
        
        # Prepare data
        log2fc = results['log2fc'].fillna(0)
        pvalues = results['pvalue'].fillna(1)
        padj = results.get('padj', pvalues)
        
        # Calculate -log10(p-value)
        neg_log_p = -np.log10(pvalues.replace(0, 1e-300))  # Avoid log(0)
        
        # Classify points
        significant = padj < alpha
        high_fc = np.abs(log2fc) > fc_threshold
        
        # Different point types
        ns_points = ~significant & ~high_fc  # Not significant, low FC
        fc_only = ~significant & high_fc     # High FC but not significant
        sig_only = significant & ~high_fc    # Significant but low FC
        sig_fc = significant & high_fc       # Significant and high FC
        
        # Plot points
        if ns_points.any():
            ax.scatter(log2fc[ns_points], neg_log_p[ns_points], 
                      c='lightgray', alpha=0.6, s=20, label='Not significant')
        
        if fc_only.any():
            ax.scatter(log2fc[fc_only], neg_log_p[fc_only], 
                      c='blue', alpha=0.7, s=30, label=f'|log2FC| > {fc_threshold}')
        
        if sig_only.any():
            ax.scatter(log2fc[sig_only], neg_log_p[sig_only], 
                      c='orange', alpha=0.7, s=30, label=f'FDR < {alpha}')
        
        if sig_fc.any():
            ax.scatter(log2fc[sig_fc], neg_log_p[sig_fc], 
                      c='red', alpha=0.8, s=40, label=f'Significant & |log2FC| > {fc_threshold}')
        
        # Add threshold lines
        ax.axhline(y=-np.log10(alpha), color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=fc_threshold, color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=-fc_threshold, color='red', linestyle='--', alpha=0.7)
        
        # Labels and formatting
        ax.set_xlabel('log2(Fold Change)')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add text with counts
        n_sig = significant.sum()
        n_total = len(results)
        ax.text(0.02, 0.98, f'Significant: {n_sig}/{n_total} ({n_sig/n_total*100:.1f}%)',
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Volcano plot saved to {save_path}")
        
        return fig
    
    def plot_method_comparison(self, 
                              method_results: Dict[str, pd.DataFrame],
                              alpha: float = 0.05,
                              save_path: Optional[str] = None) -> plt.Figure:
        """
        Compare results from multiple methods
        
        Parameters:
        -----------
        method_results : dict
            Dictionary of method results
        alpha : float
            Significance threshold
        save_path : str, optional
            Path to save the plot
        """
        
        n_methods = len(method_results)
        if n_methods < 2:
            raise ValueError("Need at least 2 methods for comparison")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Method Comparison', fontsize=16, fontweight='bold')
        
        # Prepare data
        method_names = list(method_results.keys())
        method_data = {}
        
        for method, results in method_results.items():
            padj_col = 'padj' if 'padj' in results.columns else 'qvalue' if 'qvalue' in results.columns else 'pvalue'
            significant = results[results[padj_col] < alpha]
            method_data[method] = {
                'n_significant': len(significant),
                'significant_features': set(significant['feature']),
                'pvalues': results['pvalue'].fillna(1),
                'log2fc': results.get('log2fc', pd.Series([0] * len(results)))
            }
        
        # Number of significant features
        ax = axes[0, 0]
        n_sig = [method_data[method]['n_significant'] for method in method_names]
        bars = ax.bar(method_names, n_sig)
        ax.set_title('Number of Significant Features')
        ax.set_ylabel('Count')
        ax.tick_params(axis='x', rotation=45)
        
        # Add value labels
        for bar, value in zip(bars, n_sig):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(n_sig)*0.01,
                   f'{value}', ha='center', va='bottom')
        
        # Overlap analysis (Venn-like)
        ax = axes[0, 1]
        if n_methods == 2:
            # Simple two-set overlap
            set1, set2 = [method_data[method]['significant_features'] for method in method_names]
            overlap = len(set1 & set2)
            only1 = len(set1 - set2)
            only2 = len(set2 - set1)
            
            categories = [f'Only {method_names[0]}', 'Overlap', f'Only {method_names[1]}']
            values = [only1, overlap, only2]
            
            bars = ax.bar(categories, values)
            ax.set_title('Feature Overlap')
            ax.set_ylabel('Number of Features')
            ax.tick_params(axis='x', rotation=45)
        else:
            # Multi-method overlap (simplified)
            all_features = set()
            for method_data_dict in method_data.values():
                all_features.update(method_data_dict['significant_features'])
            
            overlap_counts = []
            for feature in all_features:
                count = sum(1 for method_data_dict in method_data.values() 
                           if feature in method_data_dict['significant_features'])
                overlap_counts.append(count)
            
            bins = range(1, n_methods + 2)
            ax.hist(overlap_counts, bins=bins, alpha=0.7, edgecolor='black')
            ax.set_title('Feature Detection Frequency')
            ax.set_xlabel('Number of Methods Detecting Feature')
            ax.set_ylabel('Number of Features')
        
        # P-value distribution comparison
        ax = axes[1, 0]
        for method in method_names:
            pvals = method_data[method]['pvalues']
            ax.hist(pvals, bins=50, alpha=0.6, label=method, density=True)
        
        ax.set_xlabel('P-value')
        ax.set_ylabel('Density')
        ax.set_title('P-value Distributions')
        ax.legend()
        ax.axvline(x=alpha, color='red', linestyle='--', alpha=0.7)
        
        # Effect size comparison (if available)
        ax = axes[1, 1]
        log2fc_data = []
        labels = []
        
        for method in method_names:
            fc_values = method_data[method]['log2fc']
            if fc_values.std() > 0:  # Only if there's variation
                log2fc_data.append(fc_values)
                labels.append(method)
        
        if log2fc_data:
            ax.boxplot(log2fc_data, labels=labels)
            ax.set_title('Log2 Fold Change Distributions')
            ax.set_ylabel('log2(Fold Change)')
            ax.tick_params(axis='x', rotation=45)
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.7)
        else:
            ax.text(0.5, 0.5, 'No fold change data available', 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Log2 Fold Change Distributions')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Method comparison plot saved to {save_path}")
        
        return fig
    
    def create_interactive_dashboard(self, 
                                   analysis_results: Dict,
                                   save_path: Optional[str] = None) -> go.Figure:
        """
        Create interactive dashboard with plotly
        
        Parameters:
        -----------
        analysis_results : dict
            Complete analysis results from DifferentialAbundanceTool
        save_path : str, optional
            Path to save HTML file
        """
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                'Data Profile Summary',
                'Method Recommendations', 
                'Volcano Plot',
                'Results Overview'
            ),
            specs=[
                [{"type": "bar"}, {"type": "scatter"}],
                [{"type": "scatter"}, {"type": "bar"}]
            ]
        )
        
        profile = analysis_results['profile']
        recommendations = analysis_results['recommendations']
        analyses = analysis_results['analyses']
        
        # Data profile summary
        profile_stats = [
            profile.n_samples,
            profile.n_features, 
            profile.sparsity * 100,
            profile.zero_inflation * 100
        ]
        profile_labels = ['Samples', 'Features', 'Sparsity %', 'Zero Inflation %']
        
        fig.add_trace(
            go.Bar(
                x=profile_labels,
                y=profile_stats,
                name='Data Profile',
                showlegend=False,
                text=[f'{v:.1f}' for v in profile_stats],
                textposition='auto'
            ),
            row=1, col=1
        )
        
        # Method recommendations
        scores = recommendations.scores
        method_names = list(scores.keys())
        score_values = list(scores.values())
        
        colors = ['red' if name == recommendations.primary_method else 'blue' 
                 for name in method_names]
        
        fig.add_trace(
            go.Scatter(
                x=method_names,
                y=score_values,
                mode='markers',
                marker=dict(
                    size=15,
                    color=colors,
                    line=dict(width=2, color='black')
                ),
                name='Method Scores',
                showlegend=False,
                text=[f'Score: {s:.1f}' for s in score_values],
                hovertemplate='Method: %{x}<br>Score: %{y:.1f}<extra></extra>'
            ),
            row=1, col=2
        )
        
        # Volcano plot (using first analysis result)
        if analyses:
            first_method = list(analyses.keys())[0]
            results = analyses[first_method]
            
            log2fc = results.get('log2fc', pd.Series([0] * len(results)))
            pvalues = results['pvalue'].fillna(1)
            padj = results.get('padj', pvalues)
            
            neg_log_p = -np.log10(pvalues.replace(0, 1e-300))
            
            # Color by significance
            colors = ['red' if p < 0.05 else 'blue' for p in padj]
            
            fig.add_trace(
                go.Scatter(
                    x=log2fc,
                    y=neg_log_p,
                    mode='markers',
                    marker=dict(
                        color=colors,
                        size=6,
                        opacity=0.7
                    ),
                    name=f'{first_method} Results',
                    showlegend=False,
                    text=results['feature'],
                    hovertemplate='Feature: %{text}<br>' +
                                  'log2FC: %{x:.3f}<br>' +
                                  '-log10(p): %{y:.3f}<extra></extra>'
                ),
                row=2, col=1
            )
            
            # Add significance threshold line
            fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="red", 
                         row=2, col=1, annotation_text="p=0.05")
        
        # Results overview
        method_sig_counts = []
        method_labels = []
        
        for method, results in analyses.items():
            padj_col = 'padj' if 'padj' in results.columns else 'pvalue'
            n_sig = (results[padj_col] < 0.05).sum()
            method_sig_counts.append(n_sig)
            method_labels.append(method)
        
        fig.add_trace(
            go.Bar(
                x=method_labels,
                y=method_sig_counts,
                name='Significant Features',
                showlegend=False,
                text=[f'{v}' for v in method_sig_counts],
                textposition='auto'
            ),
            row=2, col=2
        )
        
        # Update layout
        fig.update_layout(
            title_text="DAAadvisor Analysis Dashboard",
            height=800,
            showlegend=False
        )
        
        # Update axis labels
        fig.update_xaxes(title_text="Metric", row=1, col=1)
        fig.update_yaxes(title_text="Value", row=1, col=1)
        
        fig.update_xaxes(title_text="Method", row=1, col=2)
        fig.update_yaxes(title_text="Score", row=1, col=2)
        
        fig.update_xaxes(title_text="log2(Fold Change)", row=2, col=1)
        fig.update_yaxes(title_text="-log10(p-value)", row=2, col=1)
        
        fig.update_xaxes(title_text="Method", row=2, col=2)
        fig.update_yaxes(title_text="Significant Features", row=2, col=2)
        
        if save_path:
            fig.write_html(save_path)
            logger.info(f"Interactive dashboard saved to {save_path}")
        
        return fig
    
    def _calculate_quality_scores(self, profile: DataProfile) -> Dict[str, float]:
        """Calculate data quality scores"""
        
        scores = {}
        
        # Sample size score
        if profile.n_samples >= 50:
            scores['Sample Size'] = 1.0
        elif profile.n_samples >= 20:
            scores['Sample Size'] = 0.8
        elif profile.n_samples >= 10:
            scores['Sample Size'] = 0.6
        else:
            scores['Sample Size'] = 0.3
        
        # Balance score
        if profile.group_sizes:
            sizes = list(profile.group_sizes.values())
            if len(sizes) >= 2:
                balance_ratio = min(sizes) / max(sizes)
                scores['Group Balance'] = balance_ratio
            else:
                scores['Group Balance'] = 0.5
        else:
            scores['Group Balance'] = 0.0
        
        # Sparsity score (lower is better for most methods)
        if profile.sparsity < 0.3:
            scores['Sparsity'] = 1.0
        elif profile.sparsity < 0.6:
            scores['Sparsity'] = 0.8
        elif profile.sparsity < 0.8:
            scores['Sparsity'] = 0.6
        else:
            scores['Sparsity'] = 0.4
        
        # Sequencing depth consistency
        if profile.sequencing_depth_cv < 0.3:
            scores['Depth Consistency'] = 1.0
        elif profile.sequencing_depth_cv < 0.5:
            scores['Depth Consistency'] = 0.8
        elif profile.sequencing_depth_cv < 1.0:
            scores['Depth Consistency'] = 0.6
        else:
            scores['Depth Consistency'] = 0.4
        
        return scores


def create_comprehensive_report(analysis_results: Dict, 
                              output_dir: str = "analysis_report") -> None:
    """
    Create comprehensive analysis report with all visualizations
    
    Parameters:
    -----------
    analysis_results : dict
        Complete analysis results from DifferentialAbundanceTool
    output_dir : str
        Directory to save report
    """
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    visualizer = DAAVisualizer()
    
    # Data characteristics plot
    profile_fig = visualizer.plot_data_characteristics(
        analysis_results['profile'],
        save_path=output_path / 'data_characteristics.png'
    )
    plt.close(profile_fig)
    
    # Method comparison (if multiple methods)
    if len(analysis_results['analyses']) > 1:
        comparison_fig = visualizer.plot_method_comparison(
            analysis_results['analyses'],
            save_path=output_path / 'method_comparison.png'
        )
        plt.close(comparison_fig)
    
    # Volcano plots for each method
    for method_name, results in analysis_results['analyses'].items():
        volcano_fig = visualizer.plot_volcano(
            results,
            title=f"Volcano Plot - {method_name.upper()}",
            save_path=output_path / f'volcano_{method_name}.png'
        )
        plt.close(volcano_fig)
    
    # Interactive dashboard
    dashboard = visualizer.create_interactive_dashboard(
        analysis_results,
        save_path=output_path / 'interactive_dashboard.html'
    )
    
    logger.info(f"Comprehensive report created in {output_dir}")