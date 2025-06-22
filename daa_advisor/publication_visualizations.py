#!/usr/bin/env python3
"""
Publication-Quality Visualizations for DAAadvisor Benchmarking

This module creates publication-ready figures and tables suitable for 
high-impact journals including Nature Methods, Bioinformatics, and Microbiome.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Any
import logging
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.figure_factory as ff
from scipy import stats
import warnings

logger = logging.getLogger(__name__)

class PublicationVisualizer:
    """
    Create publication-quality visualizations for benchmarking results
    
    Features:
    - High-resolution publication figures (300+ DPI)
    - Journal-style formatting and color schemes
    - Statistical significance annotations
    - Multi-panel composite figures
    - Interactive supplementary materials
    """
    
    def __init__(self, output_dir: str = "publication_figures", dpi: int = 300):
        """
        Initialize publication visualizer
        
        Parameters:
        -----------
        output_dir : str
            Directory for saving figures
        dpi : int
            Resolution for publication figures
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.dpi = dpi
        
        # Set publication-quality style
        self._setup_publication_style()
        
        logger.info(f"📊 Publication visualizer initialized: {output_dir}")
    
    def _setup_publication_style(self):
        """Setup matplotlib for publication-quality figures"""
        
        # Publication style settings
        plt.style.use('default')
        
        # Set font properties
        plt.rcParams.update({
            'font.size': 12,
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'Liberation Serif', 'serif'],
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'figure.titlesize': 18,
            'figure.dpi': self.dpi,
            'savefig.dpi': self.dpi,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            'axes.linewidth': 1.2,
            'axes.edgecolor': 'black',
            'axes.facecolor': 'white',
            'figure.facecolor': 'white'
        })
        
        # Define color palette
        self.method_colors = {
            'wilcoxon': '#1f77b4',      # Blue
            'deseq2': '#ff7f0e',        # Orange  
            'edger': '#2ca02c',         # Green
            'aldex2': '#d62728',        # Red
            'ancom-bc': '#9467bd',      # Purple
            'metagenomeseq': '#8c564b', # Brown
            'consensus': '#e377c2'      # Pink
        }
        
        self.dataset_colors = {
            'IBD': '#e41a1c',
            'CRC': '#377eb8', 
            'T2D': '#4daf4a',
            'antibiotic': '#984ea3',
            'diet': '#ff7f00',
            'controlled': '#808080'
        }
    
    def create_main_performance_figure(self, 
                                     summary_df: pd.DataFrame,
                                     save_path: str = None) -> plt.Figure:
        """
        Create main performance comparison figure (Figure 1)
        
        Parameters:
        -----------
        summary_df : pd.DataFrame
            Summary results from benchmarking
        save_path : str, optional
            Path to save figure
            
        Returns:
        --------
        plt.Figure
            Main performance figure
        """
        
        logger.info("📊 Creating main performance figure...")
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('DAAadvisor Performance Across Real-World Microbiome Datasets', 
                    fontsize=20, fontweight='bold', y=0.95)
        
        # Extract numeric values from confidence interval strings
        def extract_mean(value_str):
            try:
                return float(value_str.split(' ± ')[0])
            except:
                return np.nan
        
        def extract_std(value_str):
            try:
                return float(value_str.split(' ± ')[1])
            except:
                return np.nan
        
        # Prepare data
        plot_data = summary_df.copy()
        
        metrics = ['F1_Score', 'Sensitivity', 'Specificity', 'AUROC', 'AUPRC', 'Precision']
        metric_labels = ['F1 Score', 'Sensitivity', 'Specificity', 'AUROC', 'AUPRC', 'Precision']
        
        for i, (metric, label) in enumerate(zip(metrics, metric_labels)):
            ax = axes[i // 3, i % 3]
            
            # Extract means and stds
            means = plot_data[metric].apply(extract_mean)
            stds = plot_data[metric].apply(extract_std)
            
            # Group by method
            methods = plot_data['Method'].unique()
            
            x_pos = np.arange(len(methods))
            
            method_means = []
            method_stds = []
            
            for method in methods:
                method_data = plot_data[plot_data['Method'] == method]
                method_mean = means[method_data.index].mean()
                method_std = stds[method_data.index].mean()
                method_means.append(method_mean)
                method_stds.append(method_std)
            
            # Create bar plot
            bars = ax.bar(x_pos, method_means, yerr=method_stds, 
                         capsize=5, capthick=2, 
                         color=[self.method_colors.get(m, '#808080') for m in methods],
                         alpha=0.8, edgecolor='black', linewidth=1)
            
            ax.set_xlabel('Method', fontweight='bold')
            ax.set_ylabel(label, fontweight='bold')
            ax.set_title(f'{label} Performance', fontweight='bold')
            ax.set_xticks(x_pos)
            ax.set_xticklabels(methods, rotation=45, ha='right')
            ax.grid(axis='y', alpha=0.3)
            
            # Add value labels on bars
            for bar, mean_val in zip(bars, method_means):
                if not np.isnan(mean_val):
                    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                           f'{mean_val:.3f}', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # Save figure
        if save_path is None:
            save_path = self.output_dir / "main_performance_figure.png"
        
        plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"✅ Main performance figure saved: {save_path}")
        
        return fig
    
    def create_dataset_comparison_figure(self,
                                       all_results: Dict,
                                       save_path: str = None) -> plt.Figure:
        """
        Create dataset-specific performance comparison (Figure 2)
        
        Parameters:
        -----------
        all_results : Dict
            Complete benchmark results
        save_path : str, optional
            Path to save figure
            
        Returns:
        --------
        plt.Figure
            Dataset comparison figure
        """
        
        logger.info("📊 Creating dataset comparison figure...")
        
        # Organize data by dataset type
        dataset_types = {
            'Disease States': [k for k in all_results.keys() if any(disease in k.lower() for disease in ['ibd', 'crc', 't2d'])],
            'Antibiotic Studies': [k for k in all_results.keys() if 'antibiotic' in k.lower()],
            'Controlled Experiments': [k for k in all_results.keys() if 'controlled' in k.lower()],
            'Validation Studies': [k for k in all_results.keys() if any(val in k.lower() for val in ['obesity', 'diet'])]
        }
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Performance Across Different Study Types', 
                    fontsize=20, fontweight='bold', y=0.95)
        
        for idx, (study_type, datasets) in enumerate(dataset_types.items()):
            if not datasets:
                continue
                
            ax = axes[idx // 2, idx % 2]
            
            # Collect F1 scores for this study type
            method_f1_scores = {}
            
            for dataset_name in datasets:
                if dataset_name in all_results:
                    for method, stats in all_results[dataset_name].items():
                        if 'f1_score_mean' in stats:
                            if method not in method_f1_scores:
                                method_f1_scores[method] = []
                            method_f1_scores[method].append(stats['f1_score_mean'])
            
            # Create violin plot
            data_for_violin = []
            labels_for_violin = []
            
            for method, scores in method_f1_scores.items():
                data_for_violin.extend(scores)
                labels_for_violin.extend([method] * len(scores))
            
            if data_for_violin:
                violin_data = pd.DataFrame({
                    'F1_Score': data_for_violin,
                    'Method': labels_for_violin
                })
                
                sns.violinplot(data=violin_data, x='Method', y='F1_Score', ax=ax,
                             palette=[self.method_colors.get(m, '#808080') for m in violin_data['Method'].unique()])
                
                ax.set_title(f'{study_type}', fontweight='bold', fontsize=14)
                ax.set_xlabel('Method', fontweight='bold')
                ax.set_ylabel('F1 Score', fontweight='bold')
                ax.tick_params(axis='x', rotation=45)
                ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        # Save figure
        if save_path is None:
            save_path = self.output_dir / "dataset_comparison_figure.png"
        
        plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"✅ Dataset comparison figure saved: {save_path}")
        
        return fig
    
    def create_statistical_significance_figure(self,
                                             all_results: Dict,
                                             save_path: str = None) -> plt.Figure:
        """
        Create statistical significance comparison figure
        
        Parameters:
        -----------
        all_results : Dict
            Complete benchmark results
        save_path : str, optional
            Path to save figure
            
        Returns:
        --------
        plt.Figure
            Statistical significance figure
        """
        
        logger.info("📊 Creating statistical significance figure...")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle('Statistical Rigor and Confidence Intervals', 
                    fontsize=18, fontweight='bold')
        
        # Collect data for confidence intervals
        methods_data = {}
        
        for dataset_name, dataset_results in all_results.items():
            for method, stats in dataset_results.items():
                if 'f1_score_mean' in stats and 'f1_score_ci_lower' in stats:
                    if method not in methods_data:
                        methods_data[method] = {'means': [], 'ci_lower': [], 'ci_upper': []}
                    
                    methods_data[method]['means'].append(stats['f1_score_mean'])
                    methods_data[method]['ci_lower'].append(stats['f1_score_ci_lower'])
                    methods_data[method]['ci_upper'].append(stats['f1_score_ci_upper'])
        
        # Plot 1: Confidence intervals
        methods = list(methods_data.keys())
        x_pos = np.arange(len(methods))
        
        for i, method in enumerate(methods):
            means = np.array(methods_data[method]['means'])
            ci_lower = np.array(methods_data[method]['ci_lower'])
            ci_upper = np.array(methods_data[method]['ci_upper'])
            
            # Plot individual points
            y_positions = np.full(len(means), i) + np.random.normal(0, 0.05, len(means))
            ax1.scatter(means, y_positions, 
                       color=self.method_colors.get(method, '#808080'),
                       alpha=0.6, s=50)
            
            # Plot confidence intervals
            for j, (mean, lower, upper) in enumerate(zip(means, ci_lower, ci_upper)):
                ax1.plot([lower, upper], [y_positions[j], y_positions[j]], 
                        color=self.method_colors.get(method, '#808080'),
                        alpha=0.8, linewidth=2)
        
        ax1.set_yticks(x_pos)
        ax1.set_yticklabels(methods)
        ax1.set_xlabel('F1 Score', fontweight='bold')
        ax1.set_ylabel('Method', fontweight='bold')
        ax1.set_title('95% Confidence Intervals', fontweight='bold')
        ax1.grid(axis='x', alpha=0.3)
        
        # Plot 2: Bootstrap success rates
        bootstrap_success = {}
        for dataset_name, dataset_results in all_results.items():
            for method, stats in dataset_results.items():
                if 'n_bootstrap_success' in stats:
                    if method not in bootstrap_success:
                        bootstrap_success[method] = []
                    bootstrap_success[method].append(stats['n_bootstrap_success'])
        
        if bootstrap_success:
            methods = list(bootstrap_success.keys())
            success_rates = [np.mean(bootstrap_success[method]) for method in methods]
            
            bars = ax2.bar(methods, success_rates,
                          color=[self.method_colors.get(m, '#808080') for m in methods],
                          alpha=0.8, edgecolor='black', linewidth=1)
            
            ax2.set_xlabel('Method', fontweight='bold')
            ax2.set_ylabel('Average Bootstrap Success Rate', fontweight='bold')
            ax2.set_title('Bootstrap Reliability', fontweight='bold')
            ax2.tick_params(axis='x', rotation=45)
            ax2.grid(axis='y', alpha=0.3)
            
            # Add value labels
            for bar, rate in zip(bars, success_rates):
                ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                        f'{rate:.1f}', ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        
        # Save figure
        if save_path is None:
            save_path = self.output_dir / "statistical_significance_figure.png"
        
        plt.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
        logger.info(f"✅ Statistical significance figure saved: {save_path}")
        
        return fig
    
    def create_interactive_dashboard(self,
                                   all_results: Dict,
                                   summary_df: pd.DataFrame,
                                   save_path: str = None) -> str:
        """
        Create interactive HTML dashboard for supplementary materials
        
        Parameters:
        -----------
        all_results : Dict
            Complete benchmark results
        summary_df : pd.DataFrame
            Summary results table
        save_path : str, optional
            Path to save HTML file
            
        Returns:
        --------
        str
            Path to saved HTML file
        """
        
        logger.info("🌐 Creating interactive dashboard...")
        
        # Create subplots
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=('Method Performance Overview', 'Dataset-Specific Results',
                          'Bootstrap Confidence Intervals', 'ROC/PR Curves',
                          'Effect Size Accuracy', 'Computational Efficiency'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Extract data for plotting
        methods = summary_df['Method'].unique()
        datasets = summary_df['Dataset'].unique()
        
        # Plot 1: Overall performance radar chart
        metrics = ['F1_Score', 'Sensitivity', 'Specificity', 'AUROC']
        
        for method in methods:
            method_data = summary_df[summary_df['Method'] == method]
            
            values = []
            for metric in metrics:
                try:
                    mean_val = method_data[metric].apply(lambda x: float(x.split(' ± ')[0])).mean()
                    values.append(mean_val)
                except:
                    values.append(0)
            
            # Close the radar chart
            values.append(values[0])
            theta = metrics + [metrics[0]]
            
            fig.add_trace(
                go.Scatterpolar(
                    r=values,
                    theta=theta,
                    fill='toself',
                    name=method,
                    line_color=self.method_colors.get(method, '#808080')
                ),
                row=1, col=1
            )
        
        # Plot 2: Heatmap of method x dataset performance
        heatmap_data = []
        heatmap_datasets = []
        heatmap_methods = []
        
        for dataset in datasets:
            dataset_data = summary_df[summary_df['Dataset'] == dataset]
            for method in methods:
                method_data = dataset_data[dataset_data['Method'] == method]
                if not method_data.empty:
                    try:
                        f1_score = float(method_data['F1_Score'].iloc[0].split(' ± ')[0])
                    except:
                        f1_score = 0
                    heatmap_data.append(f1_score)
                    heatmap_datasets.append(dataset)
                    heatmap_methods.append(method)
        
        if heatmap_data:
            fig.add_trace(
                go.Heatmap(
                    z=heatmap_data,
                    x=heatmap_methods,
                    y=heatmap_datasets,
                    colorscale='Viridis',
                    name='F1 Score Heatmap'
                ),
                row=1, col=2
            )
        
        # Update layout
        fig.update_layout(
            title_text="DAAadvisor Comprehensive Benchmarking Dashboard",
            title_x=0.5,
            title_font_size=20,
            height=1200,
            showlegend=True
        )
        
        # Save dashboard
        if save_path is None:
            save_path = self.output_dir / "interactive_dashboard.html"
        
        fig.write_html(str(save_path))
        logger.info(f"✅ Interactive dashboard saved: {save_path}")
        
        return str(save_path)
    
    def create_publication_table(self,
                               summary_df: pd.DataFrame,
                               save_path: str = None) -> pd.DataFrame:
        """
        Create publication-ready summary table
        
        Parameters:
        -----------
        summary_df : pd.DataFrame
            Summary results from benchmarking
        save_path : str, optional
            Path to save table
            
        Returns:
        --------
        pd.DataFrame
            Formatted publication table
        """
        
        logger.info("📋 Creating publication table...")
        
        # Group by method and calculate overall statistics
        pub_table_data = []
        
        methods = summary_df['Method'].unique()
        
        for method in methods:
            method_data = summary_df[summary_df['Method'] == method]
            
            # Extract numeric values
            def extract_mean(col):
                try:
                    return method_data[col].apply(lambda x: float(x.split(' ± ')[0])).mean()
                except:
                    return np.nan
            
            def extract_std(col):
                try:
                    return method_data[col].apply(lambda x: float(x.split(' ± ')[1])).mean()
                except:
                    return np.nan
            
            row = {
                'Method': method.upper(),
                'N_Datasets': len(method_data),
                'F1_Score': f"{extract_mean('F1_Score'):.3f} ± {extract_std('F1_Score'):.3f}",
                'Sensitivity': f"{extract_mean('Sensitivity'):.3f} ± {extract_std('Sensitivity'):.3f}",
                'Specificity': f"{extract_mean('Specificity'):.3f} ± {extract_std('Specificity'):.3f}",
                'Precision': f"{extract_mean('Precision'):.3f} ± {extract_std('Precision'):.3f}",
                'AUROC': f"{extract_mean('AUROC'):.3f} ± {extract_std('AUROC'):.3f}",
                'AUPRC': f"{extract_mean('AUPRC'):.3f} ± {extract_std('AUPRC'):.3f}",
                'FDR': f"{extract_mean('FDR'):.3f} ± {extract_std('FDR'):.3f}",
                'Bootstrap_Success': f"{method_data['Bootstrap_N'].mean():.0f}/{method_data['Bootstrap_N'].max():.0f}"
            }
            pub_table_data.append(row)
        
        pub_table = pd.DataFrame(pub_table_data)
        
        # Sort by F1 score (descending)
        pub_table['_f1_sort'] = pub_table['F1_Score'].apply(lambda x: float(x.split(' ± ')[0]))
        pub_table = pub_table.sort_values('_f1_sort', ascending=False).drop('_f1_sort', axis=1)
        
        # Save table
        if save_path is None:
            save_path = self.output_dir / "publication_table.csv"
        
        pub_table.to_csv(save_path, index=False)
        
        # Create LaTeX version
        latex_path = save_path.with_suffix('.tex')
        latex_table = pub_table.to_latex(index=False, float_format="%.3f", 
                                       caption="DAAadvisor Performance Across Real-World Microbiome Datasets",
                                       label="tab:daaadvisor_performance")
        
        with open(latex_path, 'w') as f:
            f.write(latex_table)
        
        logger.info(f"✅ Publication table saved: {save_path}")
        logger.info(f"✅ LaTeX table saved: {latex_path}")
        
        return pub_table
    
    def generate_all_publication_figures(self,
                                       all_results: Dict,
                                       summary_df: pd.DataFrame) -> Dict[str, str]:
        """
        Generate complete set of publication figures
        
        Parameters:
        -----------
        all_results : Dict
            Complete benchmark results
        summary_df : pd.DataFrame
            Summary results table
            
        Returns:
        --------
        Dict[str, str]
            Dictionary of figure names and paths
        """
        
        logger.info("🎨 Generating complete publication figure set...")
        
        figure_paths = {}
        
        # Main performance figure
        main_fig = self.create_main_performance_figure(summary_df)
        figure_paths['main_performance'] = str(self.output_dir / "Figure1_main_performance.png")
        
        # Dataset comparison figure
        dataset_fig = self.create_dataset_comparison_figure(all_results)
        figure_paths['dataset_comparison'] = str(self.output_dir / "Figure2_dataset_comparison.png")
        
        # Statistical significance figure
        stats_fig = self.create_statistical_significance_figure(all_results)
        figure_paths['statistical_significance'] = str(self.output_dir / "Figure3_statistical_significance.png")
        
        # Interactive dashboard
        dashboard_path = self.create_interactive_dashboard(all_results, summary_df)
        figure_paths['interactive_dashboard'] = dashboard_path
        
        # Publication table
        pub_table = self.create_publication_table(summary_df)
        figure_paths['publication_table'] = str(self.output_dir / "Table1_publication_table.csv")
        
        # Close all figures to save memory
        plt.close('all')
        
        logger.info(f"✅ Generated {len(figure_paths)} publication figures")
        
        return figure_paths


def create_publication_figures(benchmark_results: Dict,
                             output_dir: str = "publication_figures") -> Dict[str, str]:
    """
    Convenience function to create all publication figures
    
    Parameters:
    -----------
    benchmark_results : Dict
        Results from run_publication_benchmark
    output_dir : str
        Directory for saving figures
        
    Returns:
    --------
    Dict[str, str]
        Dictionary of figure paths
    """
    
    visualizer = PublicationVisualizer(output_dir=output_dir)
    
    return visualizer.generate_all_publication_figures(
        benchmark_results['results'],
        benchmark_results['summary']
    )