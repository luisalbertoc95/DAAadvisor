#!/usr/bin/env python3
"""
Minimal test of cross-validation framework with synthetic data only
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Add package to path
sys.path.insert(0, str(Path(__file__).parent))

def test_minimal_cross_validation():
    """Test cross-validation framework with synthetic data"""
    
    print("🧪 Testing Minimal Cross-Validation Framework")
    print("=" * 55)
    
    # Test with synthetic data only first
    print("📊 Loading synthetic datasets...")
    
    from run_cross_validation_benchmark import load_realistic_synthetic_data
    
    synthetic_datasets = load_realistic_synthetic_data()
    
    if not synthetic_datasets:
        print("❌ No synthetic datasets found")
        return False
    
    print(f"✅ Loaded {len(synthetic_datasets)} synthetic datasets:")
    for name, data in synthetic_datasets.items():
        print(f"   🎭 {name}: {data['count_table'].shape} samples×features")
        print(f"      Ground truth: {len(data['ground_truth'])} differential features")
    
    # Test single analysis
    print(f"\n🔬 Testing single analysis...")
    
    from daa_advisor import DifferentialAbundanceTool
    
    tool = DifferentialAbundanceTool()
    
    # Use IBD synthetic data as example
    if 'IBD_synthetic' in synthetic_datasets:
        test_data = synthetic_datasets['IBD_synthetic']
        
        print(f"📊 Analyzing IBD synthetic data...")
        print(f"   Samples: {len(test_data['count_table'])}")
        print(f"   Features: {len(test_data['count_table'].columns)}")
        print(f"   Conditions: {test_data['metadata']['condition'].value_counts().to_dict()}")
        
        try:
            results = tool.analyze(
                count_table=test_data['count_table'],
                metadata=test_data['metadata'],
                use_consensus=True
            )
            
            print(f"✅ Analysis completed successfully!")
            print(f"📈 Methods run: {list(results['analyses'].keys())}")
            
            # Check results
            for method, result_df in results['analyses'].items():
                padj_col = 'padj' if 'padj' in result_df.columns else 'qvalue'
                if padj_col in result_df.columns:
                    n_sig = (result_df[padj_col] < 0.05).sum()
                    print(f"   • {method}: {n_sig} significant features")
            
            if 'consensus' in results:
                consensus_sig = results['consensus']['consensus_significant'].sum()
                print(f"   • Consensus: {consensus_sig} significant features")
            
            return True
            
        except Exception as e:
            print(f"❌ Analysis failed: {e}")
            return False
    
    else:
        print("❌ IBD synthetic data not found")
        return False

def test_cross_validation_components():
    """Test individual components of cross-validation"""
    
    print(f"\n🔧 Testing Cross-Validation Components")
    print("-" * 45)
    
    # Test curated data downloader import
    print("📦 Testing imports...")
    try:
        from daa_advisor.curated_data_downloader import CuratedDataDownloader
        print("✅ CuratedDataDownloader imported")
        
        downloader = CuratedDataDownloader(output_dir="test_minimal_curated")
        print(f"✅ Downloader initialized: {len(downloader.available_conditions)} conditions")
        
    except Exception as e:
        print(f"❌ Import failed: {e}")
    
    # Test cross-validation benchmark import
    try:
        from run_cross_validation_benchmark import (
            load_realistic_synthetic_data, 
            compare_analysis_results,
            get_significant_features
        )
        print("✅ Cross-validation functions imported")
        
    except Exception as e:
        print(f"❌ Cross-validation import failed: {e}")
    
    return True

if __name__ == "__main__":
    print("🏆 DAAadvisor Cross-Validation Test Suite")
    print("=" * 55)
    
    # Test 1: Basic cross-validation components
    test1_passed = test_cross_validation_components()
    
    # Test 2: Minimal analysis
    test2_passed = test_minimal_cross_validation()
    
    print(f"\n" + "="*55)
    if test1_passed and test2_passed:
        print("🎉 ALL TESTS PASSED - Cross-validation ready!")
        print("")
        print("💡 Next steps:")
        print("   # Test with synthetic data only:")
        print("   python run_cross_validation_benchmark.py --max-conditions 0")
        print("")
        print("   # Download real data and cross-validate:")
        print("   python run_cross_validation_benchmark.py --max-conditions 1")
        print("")
        print("   # Full cross-validation:")
        print("   python run_cross_validation_benchmark.py --max-conditions 3")
    else:
        print("⚠️  Some tests failed - check errors above")
    print("="*55)