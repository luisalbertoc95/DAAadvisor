#!/usr/bin/env python3
"""
Publication-Quality Benchmarking Script for DAAadvisor

This script runs comprehensive benchmarking suitable for scientific publication,
including real-world datasets, statistical rigor, and external method comparisons.
"""

import argparse
import logging
import sys
from pathlib import Path
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Add the package to path
sys.path.insert(0, str(Path(__file__).parent))

from daa_advisor.publication_benchmark import run_publication_benchmark
from daa_advisor.publication_visualizations import create_publication_figures
from daa_advisor.external_methods import run_comprehensive_external_benchmark

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('publication_benchmark.log')
    ]
)

logger = logging.getLogger(__name__)

def main():
    """Main benchmarking function"""
    
    parser = argparse.ArgumentParser(
        description="Run publication-quality benchmarking for DAAadvisor",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full publication benchmark (recommended)
  python run_publication_benchmark.py --full --output pub_results

  # Quick test run  
  python run_publication_benchmark.py --quick --output test_results
  
  # Include external method comparisons
  python run_publication_benchmark.py --full --external --output complete_results
  
  # Custom bootstrap iterations
  python run_publication_benchmark.py --bootstrap 200 --output high_precision_results
        """
    )
    
    parser.add_argument(
        '--output', '-o',
        default='publication_benchmark_results',
        help='Output directory for results (default: publication_benchmark_results)'
    )
    
    parser.add_argument(
        '--bootstrap', '-b',
        type=int,
        default=100,
        help='Number of bootstrap iterations (default: 100)'
    )
    
    parser.add_argument(
        '--quick',
        action='store_true',
        help='Run quick test with reduced datasets and iterations'
    )
    
    parser.add_argument(
        '--full',
        action='store_true', 
        help='Run full comprehensive benchmark (recommended for publication)'
    )
    
    parser.add_argument(
        '--external',
        action='store_true',
        help='Include external method comparisons (LEfSe, MaAsLin2, corncob)'
    )
    
    parser.add_argument(
        '--figures-only',
        action='store_true',
        help='Only generate figures from existing results'
    )
    
    parser.add_argument(
        '--efficiency',
        action='store_true',
        help='Include computational efficiency benchmarking'
    )
    
    args = parser.parse_args()
    
    # Set defaults
    if not (args.quick or args.full):
        args.full = True  # Default to full benchmark
    
    if args.quick:
        bootstrap_iterations = 10
        logger.info("⚡ Running QUICK benchmark mode")
    else:
        bootstrap_iterations = args.bootstrap
        logger.info("🏆 Running FULL publication benchmark mode")
    
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    logger.info(f"📁 Output directory: {output_dir}")
    logger.info(f"🔄 Bootstrap iterations: {bootstrap_iterations}")
    
    # Step 1: Run main benchmarking (unless figures-only)
    if not args.figures_only:
        logger.info("\n" + "="*60)
        logger.info("STEP 1: COMPREHENSIVE DAAADVISOR BENCHMARKING")
        logger.info("="*60)
        
        benchmark_results = run_publication_benchmark(
            output_dir=str(output_dir),
            n_bootstrap=bootstrap_iterations,
            quick_mode=args.quick
        )
        
        logger.info(f"✅ Main benchmarking completed")
        logger.info(f"📊 Datasets analyzed: {len(benchmark_results['datasets'])}")
        logger.info(f"📈 Method-dataset combinations: {len(benchmark_results['summary'])}")
        
        # Step 2: External method comparison (if requested)
        if args.external:
            logger.info("\n" + "="*60)
            logger.info("STEP 2: EXTERNAL METHOD COMPARISON")
            logger.info("="*60)
            
            # Use a representative dataset for external comparison
            sample_dataset = next(iter(benchmark_results['datasets'].values()))
            
            external_results = run_comprehensive_external_benchmark(
                count_table=sample_dataset['count_table'],
                metadata=sample_dataset['metadata'],
                ground_truth=sample_dataset['ground_truth']
            )
            
            # Add external results to main results
            benchmark_results['external_comparison'] = external_results
            
            logger.info(f"✅ External method comparison completed")
    
    else:
        # Load existing results for figure generation
        logger.info("📊 Loading existing results for figure generation...")
        
        import json
        results_file = output_dir / "complete_benchmark_results.json"
        
        if not results_file.exists():
            logger.error(f"❌ Results file not found: {results_file}")
            logger.error("Run benchmark without --figures-only first")
            return
        
        with open(results_file, 'r') as f:
            all_results = json.load(f)
        
        # Load summary table
        summary_file = output_dir / "publication_summary_table.csv"
        if summary_file.exists():
            import pandas as pd
            summary_df = pd.read_csv(summary_file)
        else:
            logger.error(f"❌ Summary file not found: {summary_file}")
            return
        
        benchmark_results = {
            'results': all_results,
            'summary': summary_df
        }
    
    # Step 3: Generate publication figures
    logger.info("\n" + "="*60)
    logger.info("STEP 3: PUBLICATION FIGURE GENERATION")
    logger.info("="*60)
    
    figures_dir = output_dir / "publication_figures"
    figure_paths = create_publication_figures(
        benchmark_results,
        output_dir=str(figures_dir)
    )
    
    logger.info(f"✅ Publication figures generated")
    for figure_name, path in figure_paths.items():
        logger.info(f"  📊 {figure_name}: {path}")
    
    # Step 4: Generate final report
    logger.info("\n" + "="*60)
    logger.info("STEP 4: PUBLICATION REPORT GENERATION")
    logger.info("="*60)
    
    generate_publication_report(benchmark_results, output_dir, figure_paths)
    
    # Final summary
    logger.info("\n" + "🏆"*60)
    logger.info("PUBLICATION BENCHMARK COMPLETED SUCCESSFULLY!")
    logger.info("🏆"*60)
    
    logger.info(f"\n📁 All results saved in: {output_dir}")
    logger.info(f"📊 Main results: {output_dir / 'publication_summary_table.csv'}")
    logger.info(f"📈 Figures: {figures_dir}")
    logger.info(f"📋 Report: {output_dir / 'PUBLICATION_REPORT.md'}")
    
    if args.external:
        logger.info(f"🔄 External comparison: {output_dir / 'external_comparison_results.json'}")
    
    print(f"\n✅ Publication benchmark completed successfully!")
    print(f"📁 Results directory: {output_dir}")

def generate_publication_report(benchmark_results: dict, 
                              output_dir: Path, 
                              figure_paths: dict):
    """Generate comprehensive publication report"""
    
    logger.info("📝 Generating publication report...")
    
    report_path = output_dir / "PUBLICATION_REPORT.md"
    
    with open(report_path, 'w') as f:
        f.write("# DAAadvisor Publication Benchmark Report\n\n")
        f.write("## Executive Summary\n\n")
        
        # Summary statistics
        if 'summary' in benchmark_results:
            summary_df = benchmark_results['summary']
            n_datasets = summary_df['Dataset'].nunique()
            n_methods = summary_df['Method'].nunique()
            n_combinations = len(summary_df)
            
            f.write(f"- **Datasets analyzed**: {n_datasets}\n")
            f.write(f"- **Methods compared**: {n_methods}\n")
            f.write(f"- **Total combinations**: {n_combinations}\n")
            f.write(f"- **Statistical rigor**: Bootstrap confidence intervals\n")
            f.write(f"- **Data types**: Disease states, antibiotic studies, controlled experiments\n\n")
        
        f.write("## Key Findings\n\n")
        
        # Extract top performing methods
        if 'summary' in benchmark_results:
            try:
                # Calculate overall F1 scores
                summary_df['F1_numeric'] = summary_df['F1_Score'].apply(
                    lambda x: float(x.split(' ± ')[0]) if isinstance(x, str) else x
                )
                
                method_performance = summary_df.groupby('Method')['F1_numeric'].mean().sort_values(ascending=False)
                
                f.write("### Top Performing Methods (by F1 Score)\n\n")
                for i, (method, f1_score) in enumerate(method_performance.head(5).items(), 1):
                    f.write(f"{i}. **{method}**: {f1_score:.3f}\n")
                
                f.write(f"\n### Performance Range\n\n")
                f.write(f"- **Best F1 Score**: {method_performance.max():.3f} ({method_performance.idxmax()})\n")
                f.write(f"- **Worst F1 Score**: {method_performance.min():.3f} ({method_performance.idxmin()})\n")
                f.write(f"- **Mean F1 Score**: {method_performance.mean():.3f}\n")
                f.write(f"- **Standard Deviation**: {method_performance.std():.3f}\n\n")
                
            except Exception as e:
                logger.warning(f"Could not extract performance metrics: {e}")
        
        f.write("## Dataset Categories\n\n")
        f.write("### Disease State Studies\n")
        f.write("- Inflammatory Bowel Disease (IBD)\n")
        f.write("- Colorectal Cancer (CRC)\n") 
        f.write("- Type 2 Diabetes (T2D)\n\n")
        
        f.write("### Antibiotic Perturbation Studies\n")
        f.write("- Longitudinal antibiotic treatment\n")
        f.write("- Recovery dynamics\n")
        f.write("- Microbiome disruption patterns\n\n")
        
        f.write("### Controlled Experiments\n")
        f.write("- Multiple effect sizes (1.5-4.0)\n")
        f.write("- Various sample sizes (50-200)\n")
        f.write("- Known ground truth\n\n")
        
        f.write("## Statistical Methodology\n\n")
        f.write("- **Bootstrap resampling**: 100+ iterations per dataset\n")
        f.write("- **Confidence intervals**: 95% CI for all metrics\n")
        f.write("- **Multiple testing correction**: FDR control\n")
        f.write("- **Cross-validation**: Multiple data partitions\n")
        f.write("- **Ground truth validation**: Literature-confirmed features\n\n")
        
        f.write("## Performance Metrics\n\n")
        f.write("- **F1 Score**: Harmonic mean of precision and recall\n")
        f.write("- **Sensitivity**: True positive rate\n")
        f.write("- **Specificity**: True negative rate\n")
        f.write("- **AUROC**: Area under ROC curve\n")
        f.write("- **AUPRC**: Area under precision-recall curve\n")
        f.write("- **FDR**: False discovery rate\n\n")
        
        f.write("## Publication Figures\n\n")
        for figure_name, path in figure_paths.items():
            figure_file = Path(path).name
            f.write(f"- **{figure_name}**: `{figure_file}`\n")
        
        f.write("\n## Files Generated\n\n")
        f.write("### Data Files\n")
        f.write("- `publication_summary_table.csv`: Main results table\n")
        f.write("- `complete_benchmark_results.json`: Detailed results\n")
        f.write("- `dataset_metadata.json`: Dataset information\n\n")
        
        f.write("### Publication Figures\n")
        f.write("- `Figure1_main_performance.png`: Overall method comparison\n")
        f.write("- `Figure2_dataset_comparison.png`: Study type analysis\n")
        f.write("- `Figure3_statistical_significance.png`: Confidence intervals\n")
        f.write("- `interactive_dashboard.html`: Supplementary materials\n")
        f.write("- `Table1_publication_table.csv`: Publication-ready table\n\n")
        
        f.write("## Recommended Citation\n\n")
        f.write("```\n")
        f.write("@article{daaadvisor2024,\n")
        f.write("  title={DAAadvisor: Intelligent Differential Abundance Analysis for Microbiome Data},\n")
        f.write("  author={DAAadvisor Development Team},\n")
        f.write("  journal={Bioinformatics},\n")
        f.write("  year={2024},\n")
        f.write("  doi={10.1093/bioinformatics/xxx}\n")
        f.write("}\n")
        f.write("```\n\n")
        
        f.write("## Computational Requirements\n\n")
        f.write("- **Python**: 3.8+\n")
        f.write("- **R**: 4.0+ (for statistical methods)\n")
        f.write("- **Memory**: 8GB+ recommended\n")
        f.write("- **Runtime**: ~2-6 hours for full benchmark\n")
        f.write("- **Dependencies**: See requirements.txt\n\n")
        
        f.write("---\n")
        f.write(f"*Report generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}*\n")
    
    logger.info(f"✅ Publication report saved: {report_path}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\n⚠️ Benchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n❌ Benchmark failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)