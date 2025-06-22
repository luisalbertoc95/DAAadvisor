#!/usr/bin/env python3
"""
Quick test of the publication benchmarking framework
"""

import pandas as pd
import numpy as np
from daa_advisor import (
    run_publication_benchmark, 
    create_publication_figures,
    MicrobiomeDataGenerator
)
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def test_publication_benchmark():
    """Test the publication benchmarking framework"""
    
    print("🧪 Testing Publication Benchmarking Framework")
    print("=" * 60)
    
    # Test with quick mode for demonstration
    print("⚡ Running quick benchmark test...")
    
    try:
        # Run publication benchmark in quick mode
        results = run_publication_benchmark(
            output_dir="test_publication_results",
            n_bootstrap=5,  # Very quick for testing
            quick_mode=True
        )
        
        print(f"\n✅ Benchmark completed successfully!")
        print(f"📊 Datasets analyzed: {len(results['datasets'])}")
        print(f"📈 Summary table rows: {len(results['summary'])}")
        
        # Show sample results
        if not results['summary'].empty:
            print(f"\n📋 Sample Results:")
            print(results['summary'].head())
            
            # Test figure generation
            print(f"\n🎨 Testing figure generation...")
            figure_paths = create_publication_figures(
                results,
                output_dir="test_publication_figures"
            )
            
            print(f"✅ Generated {len(figure_paths)} publication figures:")
            for name, path in figure_paths.items():
                print(f"  📊 {name}: {path}")
        
        print(f"\n🏆 Publication benchmark test completed successfully!")
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def demo_comprehensive_metrics():
    """Demonstrate comprehensive metrics calculation"""
    
    print("\n📊 Demonstrating Comprehensive Metrics")
    print("-" * 40)
    
    from daa_advisor.publication_benchmark import PublicationBenchmark
    
    # Create sample data
    generator = MicrobiomeDataGenerator()
    count_table, metadata, ground_truth = generator.generate_asv_data(
        n_samples=40,
        n_features=100,
        n_differential=20
    )
    
    # Create mock results
    mock_results = {}
    
    # Wilcoxon results
    features = list(count_table.columns)
    wilcoxon_results = pd.DataFrame({
        'feature': features,
        'pvalue': np.random.uniform(0.001, 0.5, len(features)),
        'padj': np.random.uniform(0.01, 0.3, len(features)),
        'log2fc': np.random.normal(0, 1.5, len(features))
    })
    
    # Make some features clearly significant  
    significant_mask = wilcoxon_results['feature'].isin(ground_truth[:10])
    wilcoxon_results.loc[significant_mask, 'pvalue'] = np.random.uniform(0.0001, 0.01, significant_mask.sum())
    wilcoxon_results.loc[significant_mask, 'padj'] = np.random.uniform(0.001, 0.04, significant_mask.sum())
    
    mock_results['wilcoxon'] = wilcoxon_results
    
    # Calculate comprehensive metrics
    benchmark = PublicationBenchmark()
    metrics = benchmark.calculate_comprehensive_metrics(
        mock_results,
        ground_truth,
        features
    )
    
    print(f"🎯 Comprehensive Metrics for Mock Data:")
    for method, method_metrics in metrics.items():
        print(f"\n📈 {method.upper()}:")
        for metric, value in method_metrics.items():
            if isinstance(value, float):
                print(f"  • {metric}: {value:.3f}")
            else:
                print(f"  • {metric}: {value}")
    
    return metrics

if __name__ == "__main__":
    print("🏆 DAAadvisor Publication Benchmarking Test Suite")
    print("=" * 60)
    
    # Test 1: Basic benchmark
    success1 = test_publication_benchmark()
    
    # Test 2: Metrics demonstration  
    demo_comprehensive_metrics()
    
    if success1:
        print(f"\n🎉 All tests passed! Publication benchmarking ready for use.")
        print(f"💡 To run full benchmark: python run_publication_benchmark.py --full")
    else:
        print(f"\n⚠️  Some tests failed. Check the error messages above.")