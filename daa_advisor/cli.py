#!/usr/bin/env python3
"""
Command-line interface for DAAadvisor
"""

import click
import pandas as pd
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
def generate_example(samples, features, sparsity, output):
    """Generate example data for testing"""
    
    import numpy as np
    
    click.echo(f"🎲 Generating example data")
    click.echo(f"Samples: {samples}, Features: {features}, Sparsity: {sparsity}")
    
    # Set seed for reproducibility
    np.random.seed(42)
    
    # Generate count data with specified sparsity
    counts = np.random.negative_binomial(5, 0.3, size=(samples, features))
    
    # Add zeros for sparsity
    zero_mask = np.random.random((samples, features)) < sparsity
    counts[zero_mask] = 0
    
    # Create count table
    count_df = pd.DataFrame(
        counts,
        index=[f"Sample_{i+1}" for i in range(samples)],
        columns=[f"ASV_{i+1}" for i in range(features)]
    )
    
    # Create metadata
    metadata_df = pd.DataFrame({
        'condition': ['Control'] * (samples//2) + ['Treatment'] * (samples//2),
        'batch': np.random.choice(['A', 'B'], samples),
        'age': np.random.randint(20, 70, samples)
    }, index=count_df.index)
    
    # Create output directory
    output_path = Path(output)
    output_path.mkdir(exist_ok=True)
    
    # Save files
    count_df.to_csv(output_path / 'counts.csv')
    metadata_df.to_csv(output_path / 'metadata.csv')
    
    click.echo(f"✅ Example data saved to: {output_path}")
    click.echo(f"  - counts.csv: {count_df.shape}")
    click.echo(f"  - metadata.csv: {metadata_df.shape}")
    click.echo(f"\nTo analyze: daaadvisor analyze {output_path}/counts.csv {output_path}/metadata.csv")


if __name__ == '__main__':
    main()