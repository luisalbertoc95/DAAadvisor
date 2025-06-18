#!/usr/bin/env python3
"""
Test complete workflow with all 6 methods and consensus analysis
"""

import sys
sys.path.insert(0, '.')

import pandas as pd
import numpy as np
from daa_advisor import DifferentialAbundanceTool
from daa_advisor.methods.registry import MethodRegistry
import time

def test_complete_workflow():
    """Test all 6 methods with consensus analysis"""
    
    print("🚀 TESTING COMPLETE WORKFLOW WITH ALL 6 METHODS + CONSENSUS")
    print("=" * 70)
    
    # Create comprehensive test data
    np.random.seed(42)
    n_samples, n_features = 40, 80
    n_differential = 10
    
    # Generate realistic microbiome-like data
    baseline = np.random.negative_binomial(25, 0.25, (n_samples, n_features)) + 3
    
    # Add strong differential signal for first features
    treatment_mask = np.array([False] * 20 + [True] * 20)
    fold_changes = np.random.uniform(3.5, 5.5, n_differential)
    
    for i in range(n_differential):
        baseline[treatment_mask, i] = baseline[treatment_mask, i] * fold_changes[i]
    
    # Create DataFrames
    count_table = pd.DataFrame(
        baseline,
        columns=[f'ASV_{i:03d}' for i in range(n_features)],
        index=[f'Sample_{i:03d}' for i in range(n_samples)]
    )
    
    metadata = pd.DataFrame({
        'condition': ['Control'] * 20 + ['Treatment'] * 20,
        'batch': ['A', 'B'] * 20,
        'patient_id': [f'P{i//2}' for i in range(n_samples)]
    }, index=count_table.index)
    
    print(f"📊 Dataset: {count_table.shape[0]} samples × {count_table.shape[1]} features")
    print(f"🎯 True differential features: {n_differential}")
    print(f"📈 Sparsity: {((count_table == 0).sum().sum() / (count_table.shape[0] * count_table.shape[1])):.1%}")
    print(f"📋 Groups: {metadata['condition'].value_counts().to_dict()}")
    print()
    
    # Test individual methods first
    print("🔬 TESTING ALL 6 METHODS INDIVIDUALLY:")
    print("-" * 50)
    
    registry = MethodRegistry()
    all_methods = registry.list_methods()
    
    method_results = {}
    true_differential = [f'ASV_{i:03d}' for i in range(n_differential)]
    
    for method_name in all_methods:
        print(f"Testing {method_name.upper()}...", end=" ")
        try:
            start_time = time.time()
            method = registry.get_method(method_name)
            result = method.run(count_table, metadata, group_column='condition')
            runtime = time.time() - start_time
            
            # Calculate performance
            significant = result[result['padj'] < 0.05]
            sig_features = set(significant['feature'].tolist())
            true_set = set(true_differential)
            
            tp = len(true_set & sig_features)
            fp = len(sig_features - true_set)
            fn = len(true_set - sig_features)
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            method_results[method_name] = {
                'success': True,
                'significant': len(significant),
                'tp': tp,
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'runtime': runtime
            }
            
            print(f"✅ SUCCESS ({len(significant)} sig, F1={f1:.3f}, {runtime:.2f}s)")
            
        except Exception as e:
            method_results[method_name] = {'success': False, 'error': str(e)[:60]}
            print(f"❌ FAILED: {str(e)[:60]}...")
    
    print()
    print("📊 INDIVIDUAL METHOD PERFORMANCE:")
    print("-" * 50)
    working_methods = []
    for method, result in method_results.items():
        if result['success']:
            working_methods.append(method)
            print(f"✅ {method.upper():12}: {result['significant']:2d} sig, "
                  f"P={result['precision']:.3f}, R={result['recall']:.3f}, "
                  f"F1={result['f1']:.3f}, {result['runtime']:.2f}s")
        else:
            print(f"❌ {method.upper():12}: FAILED")
    
    print()
    print(f"🎯 SUMMARY: {len(working_methods)}/6 methods working ({len(working_methods)/6*100:.1f}%)")
    
    if len(working_methods) >= 2:
        print()
        print("🤝 TESTING CONSENSUS ANALYSIS:")
        print("-" * 40)
        
        try:
            # Initialize tool for consensus analysis
            tool = DifferentialAbundanceTool()
            
            print("Running complete analysis with consensus...")
            start_time = time.time()
            results = tool.analyze(
                count_table=count_table,
                metadata=metadata,
                data_type='asv',
                use_consensus=True,
                group_column='condition'
            )
            consensus_runtime = time.time() - start_time
            
            print("✅ Consensus analysis completed!")
            print()
            
            # Display results
            print("📋 CONSENSUS RESULTS:")
            consensus_features = tool.get_significant_features(alpha=0.05)
            print(f"Consensus significant features: {len(consensus_features)}")
            
            # Calculate consensus performance
            if hasattr(consensus_features, 'feature'):
                consensus_set = set(consensus_features['feature'])
            elif isinstance(consensus_features, list):
                consensus_set = set(consensus_features)
            else:
                consensus_set = set()
            
            true_set = set(true_differential)
            tp = len(true_set & consensus_set)
            fp = len(consensus_set - true_set)
            fn = len(true_set - consensus_set)
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            print(f"Consensus Performance:")
            print(f"  True Positives: {tp}/{len(true_differential)}")
            print(f"  Precision: {precision:.3f}")
            print(f"  Recall: {recall:.3f}")
            print(f"  F1 Score: {f1:.3f}")
            print(f"  Runtime: {consensus_runtime:.2f}s")
            
            # Show detailed summary from tool
            print()
            print("📊 DETAILED ANALYSIS SUMMARY:")
            tool.summarize_results()
            
        except Exception as e:
            print(f"❌ Consensus analysis failed: {e}")
    
    print()
    print("🎉 COMPLETE WORKFLOW TEST FINISHED!")
    print("=" * 70)
    
    return method_results, len(working_methods)

if __name__ == "__main__":
    results, working_count = test_complete_workflow()
    print(f"\n🏆 FINAL RESULT: {working_count}/6 methods functional!")