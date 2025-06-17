#!/usr/bin/env python3
"""
Visualization module for Differential Abundance Analysis Tool
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots


class DAAVisualizer:
    """Visualization tools for differential abundance analysis"""
    
    def __init__(self, style='seaborn'):
        plt.style.use(style)
        self.colors = sns.color_palette("husl", 8)
    
    def plot_data_characteristics(self, profile, save_path: Optional[str] = None) -> plt.Figure:
        """Create a comprehensive view of data characteristics"""
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle(f'Data Profile Summary - {profile.data_type.upper()} Data', fontsize=16)
        
        # 1. Sparsity heatmap
        ax = axes[0, 0]
        sparsity_data = pd.DataFrame({
            'Metric': ['Overall\nSparsity', 'Zero\nInflation', 'Depth\nVariation'],
            'Value': [profile.sparsity, profile.zero_inflation, profile.sequencing_depth_cv]
        })
        sns.barplot(data=sparsity_data, x='Metric', y='Value', ax=ax, palette='viridis')
        ax.set_title('Data Sparsity Metrics')
        ax.set_ylim(0, 1)
        
        # 2. Sample size distribution
        ax = axes[0, 1]
        if profile.group_sizes:
            groups = list(profile.group_sizes.keys())
            sizes = list(profile.group_sizes.values())
            ax.bar(groups, sizes, color=self.colors[:len(groups)])
            ax.set_title('Sample Sizes by Group')
            ax.set_ylabel('Number of Samples')
            
            # Add text labels
            for i, (g, s) in enumerate(zip(groups, sizes)):
                ax.text(i, s + 0.5, str(s), ha='center')
        
        # 3. Feature distribution info
        ax = axes[0, 2]
        info_text = f"""Data Summary:
        
Total Samples: {profile.n_samples}
Total Features: {profile.n_features}
Data Type: {profile.data_type}

Compositional Bias: {profile.compositional_bias:.3f}
Metadata Factors: {len(profile.metadata_factors)}
"""
        ax.text(0.1, 0.5, info_text, transform=ax.transAxes, 
                fontsize=12, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        ax.axis('off')
        
        # 4. Method recommendation gauge
        ax = axes[1, 0]
        ax.axis('off')
        
        # 5. Sparsity pattern visualization
        ax = axes[1, 1]
        # Create synthetic sparsity pattern for visualization
        pattern = np.random.random((20, 30))
        pattern[pattern < profile.sparsity] = 0
        pattern[pattern > 0] = 1
        
        im = ax.imshow(pattern, cmap='Greys', aspect='auto')
        ax.set_title(f'Sparsity Pattern Example\n({profile.sparsity:.1%} zeros)')
        ax.set_xlabel('Features')
        ax.set_ylabel('Samples')
        
        # 6. Recommendations text
        ax = axes[1, 2]
        ax.axis('off')
        rec_text = """Recommended Methods:
        
Based on data characteristics,
methods addressing compositional
bias and sparsity are preferred.

See full recommendations in
the analysis output.
"""
        ax.text(0.1, 0.5, rec_text, transform=ax.transAxes,
                fontsize=11, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_method_comparison(self, results: Dict, save_path: Optional[str] = None) -> plt.Figure:
        """Compare results across different methods"""
        
        analyses = results.get('analyses', {})
        if len(analyses) < 2:
            print("Need at least 2 methods for comparison")
            return None
        
        # Prepare comparison data
        method_summary = []
        for method, df in analyses.items():
            if 'padj' in df.columns:
                n_sig = (df['padj'] < 0.05).sum()
                n_total = len(df)
                mean_effect = df['log2fc'].abs().mean() if 'log2fc' in df.columns else 0
                
                method_summary.append({
                    'Method': method,
                    'Significant Features': n_sig,
                    'Proportion Significant': n_sig / n_total,
                    'Mean |log2FC|': mean_effect
                })
        
        summary_df = pd.DataFrame(method_summary)
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        
        # 1. Number of significant features
        ax = axes[0]
        sns.barplot(data=summary_df, x='Method', y='Significant Features', ax=ax)
        ax.set_title('Number of Significant Features (FDR < 0.05)')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        
        # 2. Proportion significant
        ax = axes[1]
        sns.barplot(data=summary_df, x='Method', y='Proportion Significant', ax=ax)
        ax.set_title('Proportion of Features Significant')
        ax.set_ylim(0, max(0.1, summary_df['Proportion Significant'].max() * 1.1))
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        
        # 3. Effect size distribution
        ax = axes[2]
        for method, df in analyses.items():
            if 'log2fc' in df.columns and 'padj' in df.columns:
                sig_features = df[df['padj'] < 0.05]
                if not sig_features.empty:
                    ax.hist(sig_features['log2fc'], alpha=0.5, label=method, bins=20)
        
        ax.set_xlabel('log2 Fold Change')
        ax.set_ylabel('Count')
        ax.set_title('Effect Size Distribution (Significant Features)')
        ax.legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_volcano(self, results_df: pd.DataFrame, method_name: str = "", 
                     save_path: Optional[str] = None) -> plt.Figure:
        """Create volcano plot for a single method's results"""
        
        if 'log2fc' not in results_df.columns or 'padj' not in results_df.columns:
            print("Results must contain 'log2fc' and 'padj' columns")
            return None
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Prepare data
        results_df = results_df.copy()
        results_df['-log10(padj)'] = -np.log10(results_df['padj'].clip(lower=1e-300))
        
        # Define significance thresholds
        padj_threshold = 0.05
        fc_threshold = 1.0  # log2FC threshold
        
        # Color points by significance
        colors = []
        for _, row in results_df.iterrows():
            if row['padj'] < padj_threshold:
                if abs(row['log2fc']) > fc_threshold:
                    colors.append('red')  # Significant and high fold change
                else:
                    colors.append('orange')  # Significant but low fold change
            else:
                colors.append('gray')  # Not significant
        
        # Create scatter plot
        scatter = ax.scatter(results_df['log2fc'], results_df['-log10(padj)'], 
                           c=colors, alpha=0.6, s=20)
        
        # Add threshold lines
        ax.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=fc_threshold, color='black', linestyle='--', alpha=0.5)
        ax.axvline(x=-fc_threshold, color='black', linestyle='--', alpha=0.5)
        
        # Labels and title
        ax.set_xlabel('log2 Fold Change')
        ax.set_ylabel('-log10(adjusted p-value)')
        ax.set_title(f'Volcano Plot - {method_name}' if method_name else 'Volcano Plot')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', label='Significant & |log2FC| > 1'),
            Patch(facecolor='orange', label='Significant & |log2FC| ≤ 1'),
            Patch(facecolor='gray', label='Not significant')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        # Add annotations for top features
        top_features = results_df.nsmallest(5, 'padj')
        for _, row in top_features.iterrows():
            if row['padj'] < padj_threshold:
                ax.annotate(row['feature'], 
                          xy=(row['log2fc'], row['-log10(padj)']),
                          xytext=(5, 5), textcoords='offset points',
                          fontsize=8, alpha=0.7)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def plot_consensus_heatmap(self, results: Dict, top_n: int = 50,
                              save_path: Optional[str] = None) -> plt.Figure:
        """Create heatmap showing agreement between methods"""
        
        if 'consensus' not in results or len(results['analyses']) < 2:
            print("Consensus results with multiple methods required")
            return None
        
        consensus = results['consensus']
        analyses = results['analyses']
        
        # Get top features by consensus
        top_features = consensus.nlargest(top_n, 'n_significant')['feature'].values
        
        # Create agreement matrix
        methods = list(analyses.keys())
        agreement_matrix = pd.DataFrame(0, index=top_features, columns=methods)
        
        for method, df in analyses.items():
            for feature in top_features:
                feature_data = df[df['feature'] == feature]
                if not feature_data.empty and feature_data['padj'].iloc[0] < 0.05:
                    agreement_matrix.loc[feature, method] = 1
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(8, 12))
        
        sns.heatmap(agreement_matrix, cmap='RdYlBu_r', cbar_kws={'label': 'Significant'},
                   ax=ax, linewidths=0.5, linecolor='gray')
        
        ax.set_title(f'Method Agreement - Top {top_n} Features by Consensus')
        ax.set_xlabel('Method')
        ax.set_ylabel('Feature')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig
    
    def create_interactive_dashboard(self, results: Dict) -> go.Figure:
        """Create interactive Plotly dashboard"""
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Method Performance', 'P-value Distribution',
                          'Effect Size by Method', 'Feature Rankings'),
            specs=[[{'type': 'bar'}, {'type': 'histogram'}],
                   [{'type': 'box'}, {'type': 'scatter'}]]
        )
        
        analyses = results.get('analyses', {})
        colors = px.colors.qualitative.Set3[:len(analyses)]
        
        # 1. Method performance (number of significant features)
        methods = []
        n_sig = []
        for method, df in analyses.items():
            if 'padj' in df.columns:
                methods.append(method)
                n_sig.append((df['padj'] < 0.05).sum())
        
        fig.add_trace(
            go.Bar(x=methods, y=n_sig, name='Significant Features',
                  marker_color=colors[:len(methods)]),
            row=1, col=1
        )
        
        # 2. P-value distribution
        for i, (method, df) in enumerate(analyses.items()):
            if 'pvalue' in df.columns:
                fig.add_trace(
                    go.Histogram(x=df['pvalue'], name=method, opacity=0.5,
                               marker_color=colors[i % len(colors)]),
                    row=1, col=2
                )
        
        # 3. Effect size distribution
        for i, (method, df) in enumerate(analyses.items()):
            if 'log2fc' in df.columns and 'padj' in df.columns:
                sig_features = df[df['padj'] < 0.05]
                if not sig_features.empty:
                    fig.add_trace(
                        go.Box(y=sig_features['log2fc'], name=method,
                             marker_color=colors[i % len(colors)]),
                        row=2, col=1
                    )
        
        # 4. Feature rankings scatter
        if 'consensus' in results:
            consensus = results['consensus'].head(20)
            fig.add_trace(
                go.Scatter(
                    x=consensus['mean_log2fc'],
                    y=consensus['n_significant'],
                    mode='markers+text',
                    text=consensus['feature'],
                    textposition='top center',
                    marker=dict(
                        size=10,
                        color=consensus['n_significant'],
                        colorscale='Viridis',
                        showscale=True
                    ),
                    name='Features'
                ),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            title_text="Differential Abundance Analysis Dashboard",
            showlegend=True,
            height=800
        )
        
        fig.update_xaxes(title_text="Method", row=1, col=1)
        fig.update_xaxes(title_text="P-value", row=1, col=2)
        fig.update_xaxes(title_text="Mean log2FC", row=2, col=2)
        
        fig.update_yaxes(title_text="Count", row=1, col=1)
        fig.update_yaxes(title_text="Count", row=1, col=2)
        fig.update_yaxes(title_text="log2FC", row=2, col=1)
        fig.update_yaxes(title_text="# Methods Significant", row=2, col=2)
        
        return fig
    
    def plot_abundance_patterns(self, count_table: pd.DataFrame, 
                              metadata: pd.DataFrame,
                              top_features: List[str],
                              save_path: Optional[str] = None) -> plt.Figure:
        """Plot abundance patterns for top differential features"""
        
        n_features = min(len(top_features), 9)  # Max 9 subplots
        fig, axes = plt.subplots(3, 3, figsize=(15, 12))
        axes = axes.flatten()
        
        group_col = metadata.columns[0]
        
        for i, feature in enumerate(top_features[:n_features]):
            ax = axes[i]
            
            # Prepare data
            plot_data = pd.DataFrame({
                'Abundance': count_table[feature],
                'Group': metadata[group_col]
            })
            
            # Create violin plot with points
            sns.violinplot(data=plot_data, x='Group', y='Abundance', ax=ax)
            sns.stripplot(data=plot_data, x='Group', y='Abundance', 
                         color='black', alpha=0.5, size=3, ax=ax)
            
            ax.set_title(feature)
            ax.set_xlabel('')
            
            # Add log scale if appropriate
            if plot_data['Abundance'].max() > 100:
                ax.set_yscale('log')
                ax.set_ylabel('Abundance (log scale)')
        
        # Hide unused subplots
        for i in range(n_features, 9):
            axes[i].axis('off')
        
        plt.suptitle('Abundance Patterns - Top Differential Features', fontsize=16)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        
        return fig


# Example usage with the main tool
if __name__ == "__main__":
    # This would typically be called after running the main analysis
    # Here we create mock results for demonstration
    
    # Create mock profile
    from collections import namedtuple
    DataProfile = namedtuple('DataProfile', 
                           ['data_type', 'n_samples', 'n_features', 'sparsity',
                            'zero_inflation', 'sequencing_depth_cv', 'group_sizes',
                            'compositional_bias', 'metadata_factors'])
    
    mock_profile = DataProfile(
        data_type='asv',
        n_samples=100,
        n_features=500,
        sparsity=0.75,
        zero_inflation=0.6,
        sequencing_depth_cv=0.3,
        group_sizes={'Control': 50, 'Treatment': 50},
        compositional_bias=0.4,
        metadata_factors=['condition', 'batch']
    )
    
    # Create visualizer
    viz = DAAVisualizer()
    
    # Plot data characteristics
    fig = viz.plot_data_characteristics(mock_profile)
    plt.show()
    
    # Create mock results for method comparison
    mock_results = {
        'analyses': {
            'aldex2': pd.DataFrame({
                'feature': [f'ASV_{i}' for i in range(100)],
                'pvalue': np.random.beta(0.5, 10, 100),
                'padj': np.random.beta(0.5, 10, 100),
                'log2fc': np.random.normal(0, 1.5, 100)
            }),
            'deseq2': pd.DataFrame({
                'feature': [f'ASV_{i}' for i in range(100)],
                'pvalue': np.random.beta(0.3, 10, 100),
                'padj': np.random.beta(0.3, 10, 100),
                'log2fc': np.random.normal(0, 2, 100)
            })
        }
    }
    
    # Plot method comparison
    fig = viz.plot_method_comparison(mock_results)
    plt.show()
    
    # Create volcano plot
    fig = viz.plot_volcano(mock_results['analyses']['aldex2'], 'ALDEx2')
    plt.show()
