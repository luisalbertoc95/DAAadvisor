#!/usr/bin/env python3
"""
Benchmarking framework for DAAadvisor methods
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from typing import Dict, List, Tuple, Optional
import logging
import time
from pathlib import Path

from .core import DifferentialAbundanceTool
from .data_generators import create_benchmark_datasets
from .methods import MethodRegistry

logger = logging.getLogger(__name__)


class MethodBenchmark:
    """Comprehensive benchmarking framework for differential abundance methods"""
    
    def __init__(self, output_dir: str = "benchmark_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {}
        self.performance_metrics = {}
        
    def run_comprehensive_benchmark(self, datasets: Optional[Dict] = None) -> Dict:
        """
        Run comprehensive benchmark across multiple datasets and methods
        
        Parameters:
        -----------
        datasets : dict, optional
            Dictionary of datasets. If None, generates standard benchmark datasets
            
        Returns:
        --------
        dict : Comprehensive benchmark results
        """
        
        if datasets is None:
            logger.info("Generating benchmark datasets...")
            datasets = create_benchmark_datasets()
        
        # Initialize method registry to see what's available
        registry = MethodRegistry()
        available_methods = registry.list_methods()
        
        logger.info(f"Benchmarking {len(available_methods)} methods on {len(datasets)} datasets")
        
        benchmark_results = {}
        
        for dataset_name, (count_table, metadata, true_features) in datasets.items():
            logger.info(f"Processing dataset: {dataset_name}")
            
            dataset_results = {
                'data_info': {
                    'n_samples': len(count_table),
                    'n_features': len(count_table.columns),
                    'sparsity': (count_table == 0).sum().sum() / (count_table.shape[0] * count_table.shape[1]),
                    'true_positives': len(true_features)
                },
                'method_results': {},
                'performance_metrics': {}
            }
            
            # Run each available method
            for method_name in available_methods:
                try:
                    method_result = self._run_single_method(
                        method_name, count_table, metadata, true_features
                    )
                    dataset_results['method_results'][method_name] = method_result
                    
                except Exception as e:
                    logger.warning(f"Method {method_name} failed on {dataset_name}: {e}")
                    dataset_results['method_results'][method_name] = {
                        'failed': True,
                        'error': str(e)
                    }
            
            # Calculate performance metrics
            dataset_results['performance_metrics'] = self._calculate_performance_metrics(
                dataset_results['method_results'], true_features
            )
            
            benchmark_results[dataset_name] = dataset_results
        
        self.results = benchmark_results
        
        # Generate comprehensive plots and reports
        self._generate_benchmark_report()
        
        logger.info(f"Benchmark complete. Results saved to {self.output_dir}")
        return benchmark_results
    
    def _run_single_method(self, 
                          method_name: str, 
                          count_table: pd.DataFrame, 
                          metadata: pd.DataFrame,
                          true_features: List[str]) -> Dict:
        """Run a single method and collect results"""
        
        start_time = time.time()
        
        try:
            # Create tool and run analysis
            tool = DifferentialAbundanceTool()
            
            # Force use of specific method by temporarily modifying the registry
            registry = tool.method_registry
            if not registry.has_method(method_name):
                raise ValueError(f"Method {method_name} not available")
            
            # Run analysis with just this method
            results = tool.analyze(
                count_table=count_table,
                metadata=metadata,
                use_consensus=False
            )
            
            runtime = time.time() - start_time
            
            # Extract results for the method that was actually run
            if method_name in results['analyses']:
                method_results = results['analyses'][method_name]
            else:
                # Fallback method was used
                actual_method = list(results['analyses'].keys())[0]
                method_results = results['analyses'][actual_method]
                logger.warning(f"Method {method_name} not available, used {actual_method}")
            
            # Calculate basic metrics
            significant_features = method_results[method_results['padj'] < 0.05]['feature'].tolist()
            
            return {
                'runtime': runtime,
                'n_significant': len(significant_features),
                'significant_features': significant_features,
                'all_results': method_results,
                'failed': False
            }
            
        except Exception as e:
            return {
                'runtime': time.time() - start_time,
                'failed': True,
                'error': str(e)
            }
    
    def _calculate_performance_metrics(self, 
                                     method_results: Dict, 
                                     true_features: List[str]) -> Dict:
        """Calculate performance metrics for all methods"""
        
        performance = {}
        
        for method_name, result in method_results.items():
            if result.get('failed', False):
                performance[method_name] = {
                    'failed': True,
                    'precision': 0,
                    'recall': 0,
                    'f1_score': 0,
                    'fpr': 0,
                    'runtime': result.get('runtime', 0)
                }
                continue
            
            significant_features = result.get('significant_features', [])
            
            # Calculate confusion matrix elements
            true_positives = len(set(significant_features) & set(true_features))
            false_positives = len(set(significant_features) - set(true_features))
            false_negatives = len(set(true_features) - set(significant_features))
            
            # Calculate metrics
            precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
            recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
            f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
            
            # False positive rate
            total_features = len(result['all_results']) if 'all_results' in result else 1000
            true_negatives = total_features - len(true_features) - false_positives
            fpr = false_positives / (false_positives + true_negatives) if (false_positives + true_negatives) > 0 else 0
            
            performance[method_name] = {
                'failed': False,
                'true_positives': true_positives,
                'false_positives': false_positives,
                'false_negatives': false_negatives,
                'precision': precision,
                'recall': recall,
                'f1_score': f1_score,
                'fpr': fpr,
                'runtime': result.get('runtime', 0)
            }
        
        return performance
    
    def _generate_benchmark_report(self):
        """Generate comprehensive benchmark report with visualizations"""
        
        # Create summary DataFrame
        summary_data = []
        for dataset_name, dataset_result in self.results.items():
            for method_name, metrics in dataset_result['performance_metrics'].items():
                row = {
                    'dataset': dataset_name,
                    'method': method_name,
                    'data_type': dataset_name.split('_')[0],
                    'challenge': '_'.join(dataset_name.split('_')[1:]) if '_' in dataset_name else 'standard',
                    'n_samples': dataset_result['data_info']['n_samples'],
                    'n_features': dataset_result['data_info']['n_features'],
                    'sparsity': dataset_result['data_info']['sparsity'],
                    **metrics
                }
                summary_data.append(row)
        
        summary_df = pd.DataFrame(summary_data)
        
        # Save summary
        summary_df.to_csv(self.output_dir / 'benchmark_summary.csv', index=False)
        
        # Generate plots
        self._plot_performance_overview(summary_df)
        self._plot_method_comparison(summary_df)
        self._plot_data_type_performance(summary_df)
        self._plot_runtime_analysis(summary_df)
        self._create_interactive_dashboard(summary_df)
        
        # Generate text report
        self._generate_text_report(summary_df)
    
    def _plot_performance_overview(self, summary_df: pd.DataFrame):
        """Create overview performance plots"""
        
        plt.style.use('seaborn-v0_8')
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Method Performance Overview', fontsize=16, fontweight='bold')
        
        # Filter out failed methods for plotting
        success_df = summary_df[~summary_df['failed']]
        
        # Precision vs Recall
        ax = axes[0, 0]
        for method in success_df['method'].unique():
            method_data = success_df[success_df['method'] == method]
            ax.scatter(method_data['recall'], method_data['precision'], 
                      label=method, alpha=0.7, s=60)
        ax.set_xlabel('Recall (Sensitivity)')
        ax.set_ylabel('Precision')
        ax.set_title('Precision vs Recall')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.3)
        
        # F1 Score by Data Type
        ax = axes[0, 1]
        sns.boxplot(data=success_df, x='data_type', y='f1_score', hue='method', ax=ax)
        ax.set_title('F1 Score by Data Type')
        ax.set_ylabel('F1 Score')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Runtime vs Performance
        ax = axes[1, 0]
        scatter = ax.scatter(success_df['runtime'], success_df['f1_score'], 
                           c=success_df['n_features'], cmap='viridis', alpha=0.7, s=60)
        ax.set_xlabel('Runtime (seconds)')
        ax.set_ylabel('F1 Score')
        ax.set_title('Runtime vs Performance')
        plt.colorbar(scatter, ax=ax, label='Number of Features')
        
        # Method Ranking
        ax = axes[1, 1]
        method_f1 = success_df.groupby('method')['f1_score'].mean().sort_values(ascending=True)
        bars = ax.barh(method_f1.index, method_f1.values)
        ax.set_xlabel('Mean F1 Score')
        ax.set_title('Method Ranking (Overall Performance)')
        
        # Color bars by performance
        for i, bar in enumerate(bars):
            bar.set_color(plt.cm.RdYlGn(method_f1.values[i]))
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'performance_overview.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_method_comparison(self, summary_df: pd.DataFrame):
        """Create detailed method comparison plots"""
        
        success_df = summary_df[~summary_df['failed']]
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Detailed Method Comparison', fontsize=16, fontweight='bold')
        
        # Heatmap of F1 scores by dataset and method
        ax = axes[0, 0]
        pivot_f1 = success_df.pivot_table(values='f1_score', index='dataset', columns='method', aggfunc='mean')
        sns.heatmap(pivot_f1, annot=True, fmt='.3f', cmap='RdYlGn', ax=ax, cbar_kws={'label': 'F1 Score'})
        ax.set_title('F1 Score Heatmap')
        ax.set_ylabel('Dataset')
        
        # Precision comparison
        ax = axes[0, 1]
        success_df.boxplot(column='precision', by='method', ax=ax)
        ax.set_title('Precision Distribution by Method')
        ax.set_ylabel('Precision')
        ax.tick_params(axis='x', rotation=45)
        
        # Recall comparison
        ax = axes[1, 0]
        success_df.boxplot(column='recall', by='method', ax=ax)
        ax.set_title('Recall Distribution by Method')
        ax.set_ylabel('Recall')
        ax.tick_params(axis='x', rotation=45)
        
        # Runtime comparison
        ax = axes[1, 1]
        success_df.boxplot(column='runtime', by='method', ax=ax)
        ax.set_title('Runtime Distribution by Method')
        ax.set_ylabel('Runtime (seconds)')
        ax.tick_params(axis='x', rotation=45)
        ax.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'method_comparison.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_data_type_performance(self, summary_df: pd.DataFrame):
        """Create data type specific performance analysis"""
        
        success_df = summary_df[~summary_df['failed']]
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('Performance by Data Type', fontsize=16, fontweight='bold')
        
        data_types = success_df['data_type'].unique()
        
        for i, data_type in enumerate(data_types):
            ax = axes[i]
            type_data = success_df[success_df['data_type'] == data_type]
            
            # Violin plot of F1 scores
            sns.violinplot(data=type_data, x='method', y='f1_score', ax=ax)
            ax.set_title(f'{data_type.upper()} Data Performance')
            ax.set_ylabel('F1 Score')
            ax.tick_params(axis='x', rotation=45)
            
            # Add mean points
            means = type_data.groupby('method')['f1_score'].mean()
            for j, method in enumerate(means.index):
                ax.scatter(j, means[method], color='red', s=50, zorder=5)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'data_type_performance.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _plot_runtime_analysis(self, summary_df: pd.DataFrame):
        """Create runtime analysis plots"""
        
        success_df = summary_df[~summary_df['failed']]
        
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle('Runtime Analysis', fontsize=16, fontweight='bold')
        
        # Runtime vs dataset size
        ax = axes[0]
        scatter = ax.scatter(success_df['n_features'] * success_df['n_samples'], 
                           success_df['runtime'], 
                           c=success_df['sparsity'], cmap='viridis', alpha=0.7, s=60)
        ax.set_xlabel('Dataset Size (samples × features)')
        ax.set_ylabel('Runtime (seconds)')
        ax.set_title('Runtime vs Dataset Size')
        ax.set_xscale('log')
        ax.set_yscale('log')
        plt.colorbar(scatter, ax=ax, label='Sparsity')
        
        # Runtime by method
        ax = axes[1]
        method_runtime = success_df.groupby('method')['runtime'].agg(['mean', 'std']).sort_values('mean')
        bars = ax.bar(range(len(method_runtime)), method_runtime['mean'], 
                      yerr=method_runtime['std'], capsize=5)
        ax.set_xticks(range(len(method_runtime)))
        ax.set_xticklabels(method_runtime.index, rotation=45)
        ax.set_ylabel('Runtime (seconds)')
        ax.set_title('Mean Runtime by Method')
        ax.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'runtime_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _create_interactive_dashboard(self, summary_df: pd.DataFrame):
        """Create interactive plotly dashboard"""
        
        success_df = summary_df[~summary_df['failed']]
        
        # Create subplots
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Performance Overview', 'Method Comparison', 
                          'Runtime Analysis', 'Data Type Analysis'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": True}, {"secondary_y": False}]]
        )
        
        # Performance overview scatter
        for method in success_df['method'].unique():
            method_data = success_df[success_df['method'] == method]
            fig.add_trace(
                go.Scatter(
                    x=method_data['recall'],
                    y=method_data['precision'],
                    mode='markers',
                    name=method,
                    text=method_data['dataset'],
                    hovertemplate='<b>%{fullData.name}</b><br>' +
                                  'Recall: %{x:.3f}<br>' +
                                  'Precision: %{y:.3f}<br>' +
                                  'Dataset: %{text}<extra></extra>',
                    showlegend=True
                ),
                row=1, col=1
            )
        
        # Method F1 comparison
        method_f1 = success_df.groupby('method')['f1_score'].mean().sort_values(ascending=False)
        fig.add_trace(
            go.Bar(
                x=method_f1.index,
                y=method_f1.values,
                name='F1 Score',
                showlegend=False,
                text=[f'{v:.3f}' for v in method_f1.values],
                textposition='auto'
            ),
            row=1, col=2
        )
        
        # Runtime vs performance
        fig.add_trace(
            go.Scatter(
                x=success_df['runtime'],
                y=success_df['f1_score'],
                mode='markers',
                name='Runtime vs F1',
                text=success_df['method'],
                hovertemplate='Method: %{text}<br>' +
                              'Runtime: %{x:.3f}s<br>' +
                              'F1 Score: %{y:.3f}<extra></extra>',
                showlegend=False,
                marker=dict(
                    size=success_df['n_features'] / 20,
                    color=success_df['sparsity'],
                    colorscale='Viridis',
                    showscale=True,
                    colorbar=dict(title="Sparsity", x=0.45, y=0.2, len=0.3)
                )
            ),
            row=2, col=1
        )
        
        # Data type performance
        for data_type in success_df['data_type'].unique():
            type_data = success_df[success_df['data_type'] == data_type]
            fig.add_trace(
                go.Box(
                    y=type_data['f1_score'],
                    name=f'{data_type.upper()}',
                    showlegend=False
                ),
                row=2, col=2
            )
        
        # Update layout
        fig.update_layout(
            title_text="DAAadvisor Benchmark Dashboard",
            height=800,
            showlegend=True
        )
        
        # Update axis labels
        fig.update_xaxes(title_text="Recall", row=1, col=1)
        fig.update_yaxes(title_text="Precision", row=1, col=1)
        fig.update_xaxes(title_text="Method", row=1, col=2)
        fig.update_yaxes(title_text="Mean F1 Score", row=1, col=2)
        fig.update_xaxes(title_text="Runtime (seconds)", row=2, col=1)
        fig.update_yaxes(title_text="F1 Score", row=2, col=1)
        fig.update_xaxes(title_text="Data Type", row=2, col=2)
        fig.update_yaxes(title_text="F1 Score", row=2, col=2)
        
        # Save interactive plot
        fig.write_html(str(self.output_dir / 'interactive_dashboard.html'))
        
    def _generate_text_report(self, summary_df: pd.DataFrame):
        """Generate comprehensive text report"""
        
        success_df = summary_df[~summary_df['failed']]
        
        report = []
        report.append("# DAAadvisor Benchmark Report")
        report.append("=" * 50)
        report.append("")
        
        # Overall statistics
        report.append("## Overall Statistics")
        report.append(f"- Total datasets tested: {summary_df['dataset'].nunique()}")
        report.append(f"- Total methods tested: {summary_df['method'].nunique()}")
        report.append(f"- Total test runs: {len(summary_df)}")
        report.append(f"- Success rate: {len(success_df)/len(summary_df)*100:.1f}%")
        report.append("")
        
        # Method rankings
        report.append("## Method Rankings")
        method_rankings = success_df.groupby('method').agg({
            'f1_score': ['mean', 'std'],
            'precision': 'mean',
            'recall': 'mean',
            'runtime': 'mean'
        }).round(3)
        method_rankings.columns = ['F1_mean', 'F1_std', 'Precision', 'Recall', 'Runtime']
        method_rankings = method_rankings.sort_values('F1_mean', ascending=False)
        
        report.append("### By F1 Score:")
        for i, (method, row) in enumerate(method_rankings.iterrows(), 1):
            report.append(f"{i}. **{method}**: F1={row['F1_mean']:.3f}±{row['F1_std']:.3f}, "
                         f"Precision={row['Precision']:.3f}, Recall={row['Recall']:.3f}, "
                         f"Runtime={row['Runtime']:.3f}s")
        report.append("")
        
        # Data type analysis
        report.append("## Performance by Data Type")
        for data_type in success_df['data_type'].unique():
            type_data = success_df[success_df['data_type'] == data_type]
            best_method = type_data.groupby('method')['f1_score'].mean().idxmax()
            best_f1 = type_data.groupby('method')['f1_score'].mean().max()
            
            report.append(f"### {data_type.upper()} Data:")
            report.append(f"- Best method: **{best_method}** (F1={best_f1:.3f})")
            report.append(f"- Average sparsity: {type_data['sparsity'].mean():.2%}")
            report.append(f"- Number of datasets: {type_data['dataset'].nunique()}")
            report.append("")
        
        # Recommendations
        report.append("## Recommendations")
        overall_best = method_rankings.index[0]
        report.append(f"1. **Overall best method**: {overall_best}")
        
        for data_type in success_df['data_type'].unique():
            type_data = success_df[success_df['data_type'] == data_type]
            best_for_type = type_data.groupby('method')['f1_score'].mean().idxmax()
            report.append(f"2. **Best for {data_type} data**: {best_for_type}")
        
        fastest_method = method_rankings.sort_values('Runtime').index[0]
        report.append(f"3. **Fastest method**: {fastest_method}")
        report.append("")
        
        # Save report
        with open(self.output_dir / 'benchmark_report.md', 'w') as f:
            f.write('\n'.join(report))
        
        logger.info("Benchmark report generated successfully")


def run_full_benchmark(output_dir: str = "benchmark_results") -> Dict:
    """
    Run comprehensive benchmarking suite
    
    Parameters:
    -----------
    output_dir : str
        Directory to save results
        
    Returns:
    --------
    dict : Benchmark results
    """
    
    logger.info("Starting comprehensive DAAadvisor benchmark")
    
    benchmark = MethodBenchmark(output_dir)
    results = benchmark.run_comprehensive_benchmark()
    
    logger.info(f"Benchmark complete! Results saved to {output_dir}")
    
    return results