#!/usr/bin/env python3
"""
Command-line interface for DAAadvisor
"""

import click
import pandas as pd
import numpy as np
import logging
import sys
from pathlib import Path

from .core import DifferentialAbundanceTool
from .profiler import DataProfiler
from .selector import MethodSelector

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
def main(verbose):
    """DAAadvisor: Intelligent differential abundance analysis"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)


@main.command()
@click.argument('count_table', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--data-type', type=click.Choice(['asv', 'gene', 'viral']), 
              help='Data type (auto-detected if not specified)')
@click.option('--output', '-o', default='daa_results', 
              help='Output directory')
@click.option('--consensus', is_flag=True, default=True,
              help='Run consensus analysis with multiple methods')
@click.option('--alpha', default=0.05, type=float,
              help='Significance threshold')
def analyze(count_table, metadata, data_type, output, consensus, alpha):
    """
    Run differential abundance analysis
    
    COUNT_TABLE: Path to count matrix CSV (samples x features)
    METADATA: Path to metadata CSV
    """
    
    click.echo(f"🧬 DAAadvisor Analysis")
    click.echo(f"Count table: {count_table}")
    click.echo(f"Metadata: {metadata}")
    
    try:
        # Load data
        click.echo("📊 Loading data...")
        counts = pd.read_csv(count_table, index_col=0)
        meta = pd.read_csv(metadata, index_col=0)
        
        click.echo(f"Loaded {counts.shape[0]} samples x {counts.shape[1]} features")
        
        # Run analysis
        click.echo("🔬 Running analysis...")
        tool = DifferentialAbundanceTool()
        results = tool.analyze(
            count_table=counts,
            metadata=meta,
            data_type=data_type,
            use_consensus=consensus
        )
        
        # Create output directory
        output_path = Path(output)
        output_path.mkdir(exist_ok=True)
        
        # Save results
        click.echo(f"💾 Saving results to {output_path}")
        
        # Save primary results
        primary_method = list(results['analyses'].keys())[0]
        primary_results = results['analyses'][primary_method]
        primary_results.to_csv(output_path / f'{primary_method}_results.csv', index=False)
        
        # Save significant features
        significant = tool.get_significant_features(alpha=alpha)
        significant.to_csv(output_path / 'significant_features.csv', index=False)
        
        # Save consensus if available
        if 'consensus' in results:
            results['consensus'].to_csv(output_path / 'consensus_results.csv', index=False)
        
        # Print summary
        click.echo("\n" + "="*50)
        tool.summarize_results()
        
        click.echo(f"\n✅ Analysis complete! Results saved to: {output_path}")
        click.echo(f"📈 Found {len(significant)} significant features (α = {alpha})")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument('count_table', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--data-type', type=click.Choice(['asv', 'gene', 'viral']),
              help='Data type (auto-detected if not specified)')
def profile(count_table, metadata, data_type):
    """
    Profile data characteristics only
    
    COUNT_TABLE: Path to count matrix CSV (samples x features)
    METADATA: Path to metadata CSV
    """
    
    click.echo("📊 Data Profiling")
    
    try:
        # Load data
        counts = pd.read_csv(count_table, index_col=0)
        meta = pd.read_csv(metadata, index_col=0)
        
        # Profile data
        profiler = DataProfiler()
        profile = profiler.profile_data(counts, meta, data_type)
        profiler.print_profile_summary()
        
        # Get recommendations
        selector = MethodSelector()
        recommendations = selector.recommend_methods(profile)
        selector.explain_recommendation(recommendations)
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@main.command()
def methods():
    """List available methods and their characteristics"""
    
    click.echo("🔬 Available Methods")
    click.echo("=" * 50)
    
    from .methods import MethodRegistry
    registry = MethodRegistry()
    
    available_methods = registry.list_methods()
    for method_name in available_methods:
        method_info = registry.get_method_info(method_name)
        click.echo(f"\n{method_name.upper()}")
        click.echo(f"  Citation: {method_info['citation']}")
        if method_info['parameters']:
            click.echo(f"  Parameters: {', '.join(method_info['parameters'])}")
    
    click.echo(f"\n📋 Total available methods: {len(available_methods)}")


@main.command()
@click.option('--samples', default=50, help='Number of samples')
@click.option('--features', default=200, help='Number of features')
@click.option('--sparsity', default=0.7, help='Sparsity level (0-1)')
@click.option('--output', default='example_data', help='Output directory')
@click.option('--data-type', type=click.Choice(['asv', 'gene', 'viral']), default='asv', help='Type of data to generate')
def generate_example(samples, features, sparsity, output, data_type):
    """Generate example data for testing"""
    
    from .data_generators import MicrobiomeDataGenerator
    
    click.echo(f"🎲 Generating {data_type} example data")
    click.echo(f"Samples: {samples}, Features: {features}, Sparsity: {sparsity}")
    
    generator = MicrobiomeDataGenerator(random_seed=42)
    
    if data_type == 'asv':
        count_table, metadata, diff_features = generator.generate_asv_data(
            n_samples=samples, n_features=features, sparsity=sparsity
        )
    elif data_type == 'gene':
        count_table, metadata, diff_features = generator.generate_gene_data(
            n_samples=samples, n_features=features, sparsity=sparsity
        )
    else:  # viral
        count_table, metadata, diff_features = generator.generate_viral_data(
            n_samples=samples, n_features=features, sparsity=sparsity
        )
    
    # Create output directory
    output_path = Path(output)
    output_path.mkdir(exist_ok=True)
    
    # Save files
    count_table.to_csv(output_path / 'counts.csv')
    metadata.to_csv(output_path / 'metadata.csv')
    
    # Save ground truth
    with open(output_path / 'differential_features.txt', 'w') as f:
        f.write('\n'.join(diff_features))
    
    click.echo(f"✅ Example data saved to: {output_path}")
    click.echo(f"  - counts.csv: {count_table.shape}")
    click.echo(f"  - metadata.csv: {metadata.shape}")
    click.echo(f"  - differential_features.txt: {len(diff_features)} features")
    click.echo(f"\nTo analyze: daaadvisor analyze {output_path}/counts.csv {output_path}/metadata.csv")


@main.command()
@click.argument('count_table', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--output', '-o', default='visualization_report', help='Output directory')
@click.option('--data-type', type=click.Choice(['asv', 'gene', 'viral']), help='Data type')
@click.option('--interactive', is_flag=True, help='Create interactive dashboard')
def visualize(count_table, metadata, output, data_type, interactive):
    """
    Create comprehensive visualizations for your data
    
    COUNT_TABLE: Path to count matrix CSV
    METADATA: Path to metadata CSV
    """
    
    click.echo("📊 Creating visualizations...")
    
    try:
        # Load data and run analysis
        counts = pd.read_csv(count_table, index_col=0)
        meta = pd.read_csv(metadata, index_col=0)
        
        tool = DifferentialAbundanceTool()
        results = tool.analyze(
            count_table=counts,
            metadata=meta,
            data_type=data_type,
            use_consensus=True
        )
        
        # Create visualizations
        from .visualization import create_comprehensive_report
        create_comprehensive_report(results, output)
        
        click.echo(f"✅ Visualizations created in: {output}")
        click.echo(f"  - Data characteristics plot")
        click.echo(f"  - Volcano plots")
        if len(results['analyses']) > 1:
            click.echo(f"  - Method comparison plots")
        if interactive:
            click.echo(f"  - Interactive dashboard: {output}/interactive_dashboard.html")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--output', '-o', default='benchmark_results', help='Output directory')
@click.option('--quick', is_flag=True, help='Run quick benchmark (fewer datasets)')
@click.option('--data-types', default='asv,gene,viral', help='Comma-separated data types to test')
def benchmark(output, quick, data_types):
    """
    Run comprehensive method benchmarking
    """
    
    click.echo("🏁 Starting comprehensive benchmarking...")
    
    try:
        from .benchmarking import run_full_benchmark
        from .data_generators import create_benchmark_datasets, MicrobiomeDataGenerator
        
        # Create datasets
        if quick:
            # Quick benchmark with smaller datasets
            generator = MicrobiomeDataGenerator()
            datasets = {}
            
            for data_type in data_types.split(','):
                if data_type == 'asv':
                    datasets[f'{data_type}_quick'] = generator.generate_asv_data(
                        n_samples=30, n_features=50, n_differential=5
                    )
                elif data_type == 'gene':
                    datasets[f'{data_type}_quick'] = generator.generate_gene_data(
                        n_samples=30, n_features=50, n_differential=5
                    )
                elif data_type == 'viral':
                    datasets[f'{data_type}_quick'] = generator.generate_viral_data(
                        n_samples=30, n_features=50, n_differential=5
                    )
        else:
            # Full benchmark
            datasets = None  # Will generate all benchmark datasets
        
        # Run benchmark
        results = run_full_benchmark(output)
        
        click.echo(f"✅ Benchmark complete! Results in: {output}")
        click.echo(f"  - Summary CSV")
        click.echo(f"  - Performance plots")
        click.echo(f"  - Interactive dashboard")
        click.echo(f"  - Detailed report")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.argument('count_table', type=click.Path(exists=True))
@click.argument('metadata', type=click.Path(exists=True))
@click.option('--group-column', default=None, help='Column for grouping (auto-detect if not specified)')
@click.option('--output', '-o', default='info_theory_results', help='Output directory')
@click.option('--alpha', default=0.05, type=float, help='Significance threshold')
def info_theory(count_table, metadata, group_column, output, alpha):
    """
    Run information theory-based differential abundance analysis
    
    COUNT_TABLE: Path to count matrix CSV
    METADATA: Path to metadata CSV
    """
    
    click.echo("🧮 Information Theory Analysis")
    
    try:
        # Load data
        counts = pd.read_csv(count_table, index_col=0)
        meta = pd.read_csv(metadata, index_col=0)
        
        # Auto-detect group column if not specified
        if group_column is None:
            group_column = meta.columns[0]
            click.echo(f"Using '{group_column}' as grouping variable")
        
        # Run information theory analysis
        from .information_theory import run_information_theory_analysis
        
        click.echo("Running information theory framework...")
        results = run_information_theory_analysis(
            count_table=counts,
            metadata=meta,
            group_column=group_column,
            alpha=alpha
        )
        
        # Create output directory
        output_path = Path(output)
        output_path.mkdir(exist_ok=True)
        
        # Save results
        da_results = results['differential_abundance_results']
        da_results.to_csv(output_path / 'info_theory_results.csv', index=False)
        
        # Save method selection info
        method_info = results['method_selection']
        with open(output_path / 'method_selection.txt', 'w') as f:
            f.write(f"Recommended method: {method_info['recommended_method']}\n")
            f.write(f"Confidence score: {method_info['score']:.3f}\n")
            f.write(f"Reasoning: {method_info['reasoning']}\n\n")
            
            f.write("Information Metrics:\n")
            for key, value in method_info['information_metrics'].items():
                f.write(f"  {key}: {value:.3f}\n")
        
        # Print summary
        info_summary = results['information_summary']
        significant = da_results[da_results['padj'] < alpha]
        
        click.echo(f"\n✅ Information theory analysis complete!")
        click.echo(f"📊 Analyzed {info_summary['total_features']} features")
        click.echo(f"📈 Found {len(significant)} significant features (α = {alpha})")
        click.echo(f"🧮 Mean information entropy: {info_summary['information_entropy']:.3f}")
        click.echo(f"🔬 Compositional constraint: {info_summary['compositional_constraint']:.3f}")
        click.echo(f"🎯 Recommended method: {method_info['recommended_method']}")
        click.echo(f"📊 Method confidence: {method_info['score']:.3f}")
        
        if len(significant) > 0:
            click.echo(f"\n🔍 Top 5 features by information divergence:")
            for _, row in significant.head().iterrows():
                click.echo(f"  {row['feature']}: divergence={row['information_divergence']:.3f}, p={row['pvalue']:.2e}")
        
        click.echo(f"\n📁 Results saved to: {output_path}")
        
    except Exception as e:
        click.echo(f"❌ Error: {e}", err=True)
        sys.exit(1)


@main.command()
@click.option('--output', '-o', default='test_results', help='Output directory')
def test(output):
    """Run comprehensive test suite"""
    
    click.echo("🧪 Running comprehensive tests...")
    
    try:
        # Import and run tests
        from tests.test_comprehensive import run_comprehensive_tests
        
        success = run_comprehensive_tests()
        
        if success:
            click.echo("✅ All tests passed!")
        else:
            click.echo("❌ Some tests failed!")
            sys.exit(1)
            
    except ImportError:
        click.echo("❌ Test module not found. Run from package root directory.")
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ Error running tests: {e}", err=True)
        sys.exit(1)


if __name__ == '__main__':
    main()