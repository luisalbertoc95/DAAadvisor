#!/usr/bin/env python3
"""
Basic usage example for DAAadvisor
"""

import numpy as np
import pandas as pd
import sys
import os

# Add the parent directory to the path so we can import daa_advisor
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from daa_advisor import DifferentialAbundanceTool, DataProfiler, MethodSelector

def generate_example_data():
    """Generate example microbiome data"""
    print("🎲 Generating example data...")
    
    # Set seed for reproducibility
    np.random.seed(42)
    
    # Generate count data (100 samples x 200 ASVs)
    n_samples = 100
    n_features = 200
    
    # Create sparse count matrix with realistic microbiome characteristics
    counts = np.random.negative_binomial(3, 0.3, size=(n_samples, n_features))
    
    # Add high sparsity (typical for microbiome data)
    zero_mask = np.random.random((n_samples, n_features)) < 0.75
    counts[zero_mask] = 0
    
    # Create some truly differential features
    differential_features = np.random.choice(n_features, 20, replace=False)
    treatment_samples = np.arange(50, 100)  # Second half are treatment
    
    for feature in differential_features:
        # Increase abundance in treatment group
        multiplier = np.random.uniform(2, 5)
        counts[treatment_samples, feature] = (counts[treatment_samples, feature] * multiplier).astype(int)
    
    # Create DataFrame
    count_table = pd.DataFrame(
        counts,
        index=[f"Sample_{i+1}" for i in range(n_samples)],
        columns=[f"ASV_{i+1}" for i in range(n_features)]
    )
    
    # Create metadata
    metadata = pd.DataFrame({
        'condition': ['Control'] * 50 + ['Treatment'] * 50,
        'batch': np.random.choice(['A', 'B', 'C'], n_samples),
        'age': np.random.randint(20, 70, n_samples),
        'bmi': np.random.normal(25, 5, n_samples)
    }, index=count_table.index)
    
    print(f"✅ Generated {n_samples} samples x {n_features} features")
    print(f"   Sparsity: {(counts == 0).mean():.1%}")
    print(f"   Differential features added: {len(differential_features)}")
    
    return count_table, metadata

def main():
    """Main example workflow"""
    
    print("🧬 DAAadvisor Basic Usage Example")
    print("=" * 50)
    
    # Generate example data
    count_table, metadata = generate_example_data()
    
    # Step 1: Profile the data
    print("\n📊 Step 1: Data Profiling")
    print("-" * 30)
    
    profiler = DataProfiler()
    profile = profiler.profile_data(count_table, metadata, data_type='asv')
    profiler.print_profile_summary()
    
    # Step 2: Get method recommendations
    print("\n🎯 Step 2: Method Recommendations")
    print("-" * 35)
    
    selector = MethodSelector()
    recommendations = selector.recommend_methods(profile)
    selector.explain_recommendation(recommendations)
    
    # Step 3: Run differential abundance analysis
    print("\n🔬 Step 3: Differential Abundance Analysis")
    print("-" * 45)
    
    tool = DifferentialAbundanceTool()
    results = tool.analyze(
        count_table=count_table,
        metadata=metadata,
        data_type='asv',
        use_consensus=True  # Run multiple methods for consensus
    )
    
    # Display results summary
    tool.summarize_results()
    
    # Step 4: Examine significant features
    print("\n📈 Step 4: Significant Features")
    print("-" * 32)
    
    significant_features = tool.get_significant_features(alpha=0.05)
    print(f"Found {len(significant_features)} significant features at FDR < 0.05")
    
    if len(significant_features) > 0:
        print("\nTop 5 most significant features:")
        top_features = significant_features.head()
        for idx, row in top_features.iterrows():
            print(f"  {row['feature']}: p={row['pvalue']:.2e}, padj={row['padj']:.2e}, log2FC={row['log2fc']:.2f}")
    
    # Step 5: Method comparison
    if 'consensus' in results:
        print("\n🤝 Step 5: Consensus Analysis")
        print("-" * 28)
        
        consensus = results['consensus']
        consensus_sig = consensus[consensus['consensus_significant']]
        print(f"Consensus significant features: {len(consensus_sig)}")
        
        if len(consensus_sig) > 0:
            print("\nFeatures significant across multiple methods:")
            for idx, row in consensus_sig.head().iterrows():
                print(f"  {row['feature']}: {row['n_significant']}/{row['n_methods']} methods")
    
    print("\n✅ Analysis complete!")
    print("\n💡 Next steps:")
    print("  - Save results: results['analyses']['primary_method'].to_csv('results.csv')")
    print("  - Visualize: Create volcano plots, heatmaps, etc.")
    print("  - Validate: Check significant features in independent dataset")

if __name__ == "__main__":
    main()