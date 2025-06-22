#!/usr/bin/env python3
"""
Minimal test of publication benchmarking core functionality
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add package to path
sys.path.insert(0, str(Path(__file__).parent))

def test_minimal_benchmark():
    """Test core benchmark functionality without full execution"""
    
    print("🧪 Testing Core Publication Benchmark Functionality")
    print("=" * 55)
    
    try:
        # Test imports
        print("📦 Testing imports...")
        from daa_advisor.publication_benchmark import PublicationBenchmark
        from daa_advisor.publication_visualizations import PublicationVisualizer
        from daa_advisor.external_methods import ExternalMethodBenchmark
        print("✅ All imports successful")
        
        # Test benchmark initialization
        print("\n🏗️ Testing benchmark initialization...")
        benchmark = PublicationBenchmark(
            output_dir="test_minimal_results",
            n_bootstrap=5,
            n_cv_folds=3
        )
        print("✅ Benchmark initialized successfully")
        
        # Test metrics calculation
        print("\n📊 Testing metrics calculation...")
        
        # Create sample data
        features = [f"Feature_{i}" for i in range(20)]
        ground_truth = features[:5]  # First 5 are true positives
        
        # Mock method results
        mock_results = {
            'test_method': pd.DataFrame({
                'feature': features,
                'pvalue': np.concatenate([
                    np.random.uniform(0.001, 0.01, 5),    # Significant (true positives)
                    np.random.uniform(0.1, 0.5, 15)       # Non-significant
                ]),
                'padj': np.concatenate([
                    np.random.uniform(0.001, 0.04, 5),    # Significant (true positives)
                    np.random.uniform(0.1, 0.5, 15)       # Non-significant
                ]),
                'log2fc': np.random.normal(0, 1.5, 20)
            })
        }
        
        # Calculate metrics
        metrics = benchmark.calculate_comprehensive_metrics(
            mock_results, ground_truth, features
        )
        
        print("✅ Metrics calculation successful")
        print(f"📈 Sample metrics for test_method:")
        for metric, value in metrics['test_method'].items():
            if isinstance(value, float):
                print(f"   • {metric}: {value:.3f}")
        
        # Test visualizer initialization
        print("\n🎨 Testing visualizer...")
        visualizer = PublicationVisualizer(output_dir="test_minimal_figures")
        print("✅ Visualizer initialized successfully")
        
        # Test external methods
        print("\n🔄 Testing external methods framework...")
        external_benchmark = ExternalMethodBenchmark()
        print(f"✅ External methods available: {list(external_benchmark.available_methods.keys())}")
        
        print("\n🎉 All core functionality tests passed!")
        print("💡 Ready for full publication benchmarking")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_dataset_generation():
    """Test realistic dataset generation"""
    
    print("\n📊 Testing Dataset Generation")
    print("-" * 30)
    
    try:
        from daa_advisor.publication_benchmark import PublicationBenchmark
        
        benchmark = PublicationBenchmark(output_dir="test_datasets")
        
        # Test individual dataset type generation
        print("🦠 Testing disease datasets...")
        disease_datasets = benchmark._get_disease_datasets()
        print(f"✅ Generated {len(disease_datasets)} disease datasets")
        
        print("💊 Testing antibiotic datasets...")
        antibiotic_datasets = benchmark._get_antibiotic_datasets()
        print(f"✅ Generated {len(antibiotic_datasets)} antibiotic datasets")
        
        print("🎯 Testing controlled datasets...")
        controlled_datasets = benchmark._generate_controlled_datasets()
        print(f"✅ Generated {len(controlled_datasets)} controlled datasets")
        
        # Show sample dataset info
        sample_dataset = next(iter(disease_datasets.values()))
        print(f"\n📋 Sample dataset characteristics:")
        print(f"   • Count table shape: {sample_dataset['count_table'].shape}")
        print(f"   • Metadata columns: {list(sample_dataset['metadata'].columns)}")
        print(f"   • Ground truth features: {len(sample_dataset['ground_truth'])}")
        
        return True
        
    except Exception as e:
        print(f"❌ Dataset generation test failed: {e}")
        return False

if __name__ == "__main__":
    print("🏆 DAAadvisor Publication Benchmark - Minimal Test Suite")
    print("=" * 60)
    
    # Run tests
    test1_passed = test_minimal_benchmark()
    test2_passed = test_dataset_generation()
    
    print("\n" + "="*60)
    if test1_passed and test2_passed:
        print("🎉 ALL TESTS PASSED - Publication benchmark ready!")
        print("📋 To run full benchmark:")
        print("   python run_publication_benchmark.py --quick")
        print("   python run_publication_benchmark.py --full")
    else:
        print("⚠️  Some tests failed - check errors above")
    print("="*60)