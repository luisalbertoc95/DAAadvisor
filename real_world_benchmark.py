#!/usr/bin/env python3
"""
Real-world comprehensive benchmark with increased sample sizes and dataset diversity
"""

import sys
sys.path.insert(0, '.')

import pandas as pd
import numpy as np
from daa_advisor import DifferentialAbundanceTool
from daa_advisor.methods.registry import MethodRegistry
from daa_advisor.data_generators import MicrobiomeDataGenerator
import time
import os
from pathlib import Path

def create_real_world_datasets():
    """Create diverse real-world-like datasets with varying characteristics"""
    
    datasets = {}
    generator = MicrobiomeDataGenerator()
    
    # 1. Large ASV dataset (typical 16S study)
    print("Generating large ASV dataset...")
    count_table, metadata, true_features = generator.generate_asv_data(
        n_samples=200,
        n_features=500,
        n_differential=50,
        sparsity=0.85,
        effect_size=2.5
    )
    datasets['large_asv'] = {
        'data': (count_table, metadata),
        'description': 'Large 16S rRNA study (200 samples, 500 ASVs)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    # 2. Deep gene functional dataset
    print("Generating deep gene functional dataset...")
    count_table, metadata, true_features = generator.generate_gene_data(
        n_samples=150,
        n_features=1000,
        n_differential=80,
        sparsity=0.45,
        effect_size=1.8
    )
    datasets['deep_gene'] = {
        'data': (count_table, metadata),
        'description': 'Deep functional analysis (150 samples, 1000 genes)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    # 3. Multi-batch ASV study
    print("Generating multi-batch ASV study...")
    count_table, metadata, true_features = generator.generate_asv_data(
        n_samples=120,
        n_features=400,
        n_differential=40,
        sparsity=0.80,
        effect_size=3.0
    )
    datasets['multi_batch'] = {
        'data': (count_table, metadata),
        'description': 'Multi-batch ASV study (120 samples, 400 ASVs)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    # 4. Sparse viral dataset
    print("Generating sparse viral dataset...")
    count_table, metadata, true_features = generator.generate_viral_data(
        n_samples=80,
        n_features=300,
        n_differential=25,
        sparsity=0.95,
        effect_size=4.0
    )
    datasets['sparse_viral'] = {
        'data': (count_table, metadata),
        'description': 'Sparse viral analysis (80 samples, 300 viruses)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    # 5. Longitudinal-like dataset (unbalanced)
    print("Generating longitudinal-like dataset...")
    count_table, metadata, true_features = generator.generate_asv_data(
        n_samples=180,
        n_features=350,
        n_differential=35,
        sparsity=0.75,
        effect_size=2.2
    )
    datasets['longitudinal'] = {
        'data': (count_table, metadata),
        'description': 'Longitudinal-like study (180 samples, 350 ASVs)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    # 6. High-diversity gene dataset
    print("Generating high-diversity gene dataset...")
    count_table, metadata, true_features = generator.generate_gene_data(
        n_samples=100,
        n_features=800,
        n_differential=60,
        sparsity=0.50,
        effect_size=1.5
    )
    datasets['high_diversity'] = {
        'data': (count_table, metadata),
        'description': 'High-diversity gene analysis (100 samples, 800 genes)',
        'true_diff': len(true_features),
        'true_features': true_features
    }
    
    return datasets

def benchmark_single_dataset(dataset_name, dataset_info, output_dir):
    """Benchmark all methods on a single dataset"""
    
    print(f"\n🔬 Benchmarking: {dataset_info['description']}")
    print("-" * 60)
    
    count_table, metadata = dataset_info['data']
    true_features = dataset_info['true_features']
    true_diff = len(true_features)
    
    # Dataset statistics
    sparsity = (count_table == 0).sum().sum() / (count_table.shape[0] * count_table.shape[1])
    print(f"📊 Dataset stats: {count_table.shape[0]} samples × {count_table.shape[1]} features")
    print(f"📈 Sparsity: {sparsity:.1%}, True differential: {true_diff}")
    
    # Test individual methods
    registry = MethodRegistry()
    all_methods = registry.list_methods()
    
    results = {}
    
    for method_name in all_methods:
        print(f"  Testing {method_name}...", end=" ")
        try:
            start_time = time.time()
            method = registry.get_method(method_name)
            result = method.run(count_table, metadata, group_column='condition')
            runtime = time.time() - start_time
            
            # Calculate performance
            significant = result[result['padj'] < 0.05]
            sig_features = set(significant['feature'].tolist())
            true_set = set(true_features)
            
            tp = len(true_set & sig_features)
            fp = len(sig_features - true_set)
            fn = len(true_set - sig_features)
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            results[method_name] = {
                'success': True,
                'n_significant': len(significant),
                'tp': tp,
                'fp': fp,
                'fn': fn,
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'runtime': runtime
            }
            
            print(f"✅ {len(significant)} sig, F1={f1:.3f}")
            
        except Exception as e:
            results[method_name] = {
                'success': False,
                'error': str(e)[:80],
                'runtime': 0
            }
            print(f"❌ FAILED")
    
    # Test consensus analysis
    print("  Testing consensus...", end=" ")
    try:
        start_time = time.time()
        tool = DifferentialAbundanceTool()
        consensus_results = tool.analyze(
            count_table=count_table,
            metadata=metadata,
            data_type=dataset_info['description'].split()[0].lower(),
            use_consensus=True,
            group_column='condition'
        )
        consensus_runtime = time.time() - start_time
        
        consensus_features = tool.get_significant_features(alpha=0.05)
        
        # Calculate consensus performance
        if hasattr(consensus_features, 'feature'):
            consensus_set = set(consensus_features['feature'])
        elif isinstance(consensus_features, list):
            consensus_set = set(consensus_features)
        else:
            consensus_set = set()
        
        true_set = set(true_features)
        tp = len(true_set & consensus_set)
        fp = len(consensus_set - true_set)
        fn = len(true_set - consensus_set)
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        
        results['consensus'] = {
            'success': True,
            'n_significant': len(consensus_set),
            'tp': tp,
            'fp': fp,
            'fn': fn,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'runtime': consensus_runtime
        }
        
        print(f"✅ {len(consensus_set)} sig, F1={f1:.3f}")
        
    except Exception as e:
        results['consensus'] = {
            'success': False,
            'error': str(e)[:80],
            'runtime': 0
        }
        print(f"❌ FAILED")
    
    # Save detailed results
    results_df = pd.DataFrame(results).T
    results_df.to_csv(f"{output_dir}/{dataset_name}_detailed_results.csv")
    
    return results

def run_comprehensive_benchmark():
    """Run comprehensive real-world benchmark"""
    
    print("🌍 REAL-WORLD COMPREHENSIVE BENCHMARK")
    print("=" * 60)
    print("Testing all 6 methods with diverse, large-scale datasets")
    print()
    
    # Create output directory
    output_dir = "real_world_benchmark_results"
    Path(output_dir).mkdir(exist_ok=True)
    
    # Generate datasets
    print("📊 Generating real-world datasets...")
    datasets = create_real_world_datasets()
    
    # Benchmark each dataset
    all_results = {}
    
    for dataset_name, dataset_info in datasets.items():
        results = benchmark_single_dataset(dataset_name, dataset_info, output_dir)
        all_results[dataset_name] = results
    
    # Generate summary report
    print("\n📋 COMPREHENSIVE BENCHMARK SUMMARY")
    print("=" * 60)
    
    # Calculate overall performance
    method_performance = {}
    all_methods = ['wilcoxon', 'aldex2', 'ancom-bc', 'deseq2', 'edger', 'metagenomeseq', 'consensus']
    
    for method in all_methods:
        successes = 0
        total_f1 = 0
        total_runtime = 0
        count = 0
        
        for dataset_results in all_results.values():
            if method in dataset_results and dataset_results[method]['success']:
                successes += 1
                total_f1 += dataset_results[method]['f1']
                total_runtime += dataset_results[method]['runtime']
                count += 1
        
        if count > 0:
            method_performance[method] = {
                'success_rate': successes / len(datasets),
                'avg_f1': total_f1 / count,
                'avg_runtime': total_runtime / count,
                'datasets_tested': count
            }
    
    # Display summary table
    print("\n📊 METHOD PERFORMANCE SUMMARY:")
    print("-" * 80)
    print(f"{'Method':<15} {'Success Rate':<12} {'Avg F1':<8} {'Avg Runtime':<12} {'Datasets'}")
    print("-" * 80)
    
    for method, perf in method_performance.items():
        success_pct = perf['success_rate'] * 100
        print(f"{method.upper():<15} {success_pct:>8.1f}%    {perf['avg_f1']:>6.3f}   {perf['avg_runtime']:>8.2f}s    {perf['datasets_tested']}/6")
    
    # Save summary
    summary_df = pd.DataFrame(method_performance).T
    summary_df.to_csv(f"{output_dir}/benchmark_summary.csv")
    
    # Generate detailed dataset-by-dataset results
    print(f"\n📈 DATASET-BY-DATASET RESULTS:")
    print("-" * 80)
    
    for dataset_name, dataset_info in datasets.items():
        print(f"\n🔬 {dataset_info['description']}:")
        results = all_results[dataset_name]
        
        working_methods = [m for m, r in results.items() if r['success']]
        print(f"  Working methods: {len(working_methods)}/7")
        
        if working_methods:
            best_method = max(working_methods, key=lambda m: results[m]['f1'])
            best_f1 = results[best_method]['f1']
            print(f"  Best performer: {best_method.upper()} (F1={best_f1:.3f})")
            
            if 'consensus' in working_methods:
                consensus_f1 = results['consensus']['f1']
                improvement = consensus_f1 - best_f1
                print(f"  Consensus F1: {consensus_f1:.3f} ({improvement:+.3f} vs best individual)")
    
    print(f"\n✅ Benchmark complete! Results saved to: {output_dir}/")
    print(f"📁 Generated files:")
    print(f"  - benchmark_summary.csv (overall performance)")
    for dataset_name in datasets.keys():
        print(f"  - {dataset_name}_detailed_results.csv")
    
    return all_results, method_performance

if __name__ == "__main__":
    results, performance = run_comprehensive_benchmark()
    
    print(f"\n🎯 FINAL ASSESSMENT:")
    print(f"All 6 statistical methods tested across 6 diverse real-world datasets")
    total_tests = len(results) * 7  # 6 datasets × 7 methods
    successful_tests = sum(sum(1 for r in dataset_results.values() if r['success']) for dataset_results in results.values())
    print(f"Overall success rate: {successful_tests}/{total_tests} ({successful_tests/total_tests*100:.1f}%)")