#!/usr/bin/env python3
"""
Advanced DAAadvisor Usage Examples

This script demonstrates:
1. R method integration (ALDEx2, ANCOM-BC, DESeq2, edgeR, metagenomeSeq)
2. Information theory framework for unified analysis
3. Method comparison and selection
"""

import pandas as pd
import numpy as np
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Import DAAadvisor components
from daa_advisor import (
    DifferentialAbundanceTool,
    MicrobiomeDataGenerator, 
    CompositionInformationFramework,
    MaximumEntropySelector,
    run_information_theory_analysis
)

def demonstrate_r_methods():
    """Demonstrate R method integration"""
    print("\n" + "="*60)
    print("R METHODS INTEGRATION DEMONSTRATION")
    print("="*60)
    
    # Generate example data
    generator = MicrobiomeDataGenerator(random_seed=42)
    count_table, metadata, true_features = generator.generate_gene_data(
        n_samples=40, n_features=100, n_differential=10, sparsity=0.4
    )
    
    print(f"Generated dataset: {count_table.shape[0]} samples × {count_table.shape[1]} features")
    print(f"True differential features: {len(true_features)}")
    print(f"Groups: {metadata['condition'].value_counts().to_dict()}")
    
    # Initialize tool
    tool = DifferentialAbundanceTool()
    
    # Test each R method (if available)
    r_methods = ['aldex2', 'ancom-bc', 'deseq2', 'edger', 'metagenomeseq']
    
    results = {}
    for method in r_methods:
        try:
            print(f"\n--- Testing {method.upper()} ---")
            
            # Run analysis with specific method
            result = tool.analyze(
                count_table=count_table,
                metadata=metadata,
                data_type='gene',
                use_consensus=False,
                methods=[method]  # Force specific method
            )
            
            if method in result['analyses']:
                method_results = result['analyses'][method]
                n_significant = (method_results['padj'] < 0.05).sum()
                print(f"✅ {method.upper()}: {n_significant} significant features found")
                results[method] = method_results
            else:
                print(f"❌ {method.upper()}: Method not available (R package missing)")
                
        except Exception as e:
            print(f"❌ {method.upper()}: Failed - {str(e)}")
    
    # Compare results if multiple methods worked
    if len(results) > 1:
        print(f"\n--- METHOD COMPARISON ---")
        print("Method\t\tSignificant\tMean |log2FC|\tMean p-value")
        print("-" * 55)
        
        for method, result in results.items():
            significant = result[result['padj'] < 0.05]
            n_sig = len(significant)
            mean_fc = np.mean(np.abs(result['log2fc']))
            mean_pval = np.mean(result['pvalue'])
            
            print(f"{method:<12}\t{n_sig:<8}\t{mean_fc:.3f}\t\t{mean_pval:.2e}")
    
    return results

def demonstrate_information_theory():
    """Demonstrate information theory framework"""
    print("\n" + "="*60)
    print("INFORMATION THEORY FRAMEWORK DEMONSTRATION")
    print("="*60)
    
    # Generate example data with different characteristics
    generator = MicrobiomeDataGenerator(random_seed=123)
    
    datasets = {
        'ASV (high sparsity)': generator.generate_asv_data(
            n_samples=50, n_features=80, sparsity=0.85, n_differential=8
        ),
        'Gene (moderate sparsity)': generator.generate_gene_data(
            n_samples=50, n_features=80, sparsity=0.4, n_differential=8  
        ),
        'Viral (very high sparsity)': generator.generate_viral_data(
            n_samples=50, n_features=80, sparsity=0.95, n_differential=8
        )
    }
    
    print("\n--- Information Theory Method Selection ---")
    
    selector = MaximumEntropySelector()
    
    for dataset_name, (count_table, metadata, true_features) in datasets.items():
        print(f"\n{dataset_name}:")
        print(f"  Data shape: {count_table.shape}")
        print(f"  Sparsity: {(count_table == 0).sum().sum() / (count_table.shape[0] * count_table.shape[1]):.1%}")
        
        # Get method recommendation using information theory
        recommendation = selector.select_method(count_table, metadata)
        
        print(f"  Recommended method: {recommendation['recommended_method']}")
        print(f"  Confidence score: {recommendation['score']:.3f}")
        print(f"  Reasoning: {recommendation['reasoning']}")
        
        # Show information metrics
        metrics = recommendation['information_metrics']
        print(f"  Information entropy: {metrics['total_entropy']:.3f}")
        print(f"  Compositional constraint: {metrics['constraint_strength']:.3f}")
        print(f"  Sparsity (info): {metrics['sparsity_info']:.3f}")
    
    print("\n--- Full Information Theory Analysis ---")
    
    # Run complete information theory analysis on one dataset
    count_table, metadata, true_features = datasets['Gene (moderate sparsity)']
    
    info_results = run_information_theory_analysis(
        count_table=count_table,
        metadata=metadata,
        group_column='condition',
        alpha=0.05
    )
    
    # Display results
    da_results = info_results['differential_abundance_results']
    info_summary = info_results['information_summary']
    
    print(f"\nInformation Theory Analysis Results:")
    print(f"  Total features analyzed: {info_summary['total_features']}")
    print(f"  Significant features: {info_summary['significant_features']}")
    print(f"  Overall information entropy: {info_summary['information_entropy']:.3f}")
    print(f"  Compositional constraint strength: {info_summary['compositional_constraint']:.3f}")
    
    print(f"\nTop 5 features by information divergence:")
    top_features = da_results.head()
    for _, row in top_features.iterrows():
        print(f"  {row['feature']}: divergence={row['information_divergence']:.3f}, "
              f"p={row['pvalue']:.2e}, log2FC={row['log2fc']:.2f}")
    
    return info_results

def demonstrate_unified_analysis():
    """Demonstrate unified analysis combining traditional and information theory approaches"""
    print("\n" + "="*60)
    print("UNIFIED ANALYSIS DEMONSTRATION")
    print("="*60)
    
    # Generate complex dataset
    generator = MicrobiomeDataGenerator(random_seed=456)
    count_table, metadata, true_features = generator.generate_asv_data(
        n_samples=60, n_features=120, n_differential=12, sparsity=0.7
    )
    
    print(f"Dataset: {count_table.shape[0]} samples × {count_table.shape[1]} features")
    print(f"True differential features: {len(true_features)}")
    
    # 1. Traditional DAAadvisor analysis
    print("\n--- Traditional Analysis ---")
    tool = DifferentialAbundanceTool()
    traditional_results = tool.analyze(
        count_table=count_table,
        metadata=metadata,
        data_type='asv',
        use_consensus=True
    )
    
    tool.summarize_results()
    traditional_significant = tool.get_significant_features(alpha=0.05)
    
    # 2. Information theory analysis
    print("\n--- Information Theory Analysis ---")
    info_results = run_information_theory_analysis(
        count_table=count_table,
        metadata=metadata,
        group_column='condition',
        alpha=0.05
    )
    
    info_significant = info_results['differential_abundance_results']
    info_significant = info_significant[info_significant['padj'] < 0.05]
    
    # 3. Compare approaches
    print("\n--- Comparison of Approaches ---")
    
    traditional_features = set(traditional_significant['feature'])
    info_features = set(info_significant['feature'])
    true_feature_set = set(true_features)
    
    # Calculate overlaps
    traditional_overlap = len(traditional_features & true_feature_set)
    info_overlap = len(info_features & true_feature_set)
    method_overlap = len(traditional_features & info_features)
    
    print(f"Traditional method:")
    print(f"  Significant features: {len(traditional_features)}")
    print(f"  True positives: {traditional_overlap}/{len(true_feature_set)} ({traditional_overlap/len(true_feature_set)*100:.1f}%)")
    
    print(f"Information theory method:")
    print(f"  Significant features: {len(info_features)}")
    print(f"  True positives: {info_overlap}/{len(true_feature_set)} ({info_overlap/len(true_feature_set)*100:.1f}%)")
    
    print(f"Method agreement:")
    print(f"  Common significant features: {method_overlap}")
    print(f"  Agreement rate: {method_overlap/max(len(traditional_features), len(info_features))*100:.1f}%")
    
    # Show unique findings from each approach
    traditional_only = traditional_features - info_features
    info_only = info_features - traditional_features
    
    if traditional_only:
        print(f"\nFeatures unique to traditional analysis: {len(traditional_only)}")
        for feature in list(traditional_only)[:3]:  # Show first 3
            print(f"  {feature}")
    
    if info_only:
        print(f"\nFeatures unique to information theory: {len(info_only)}")
        for feature in list(info_only)[:3]:  # Show first 3
            print(f"  {feature}")
    
    return {
        'traditional': traditional_results,
        'information_theory': info_results,
        'comparison': {
            'traditional_features': traditional_features,
            'info_features': info_features,
            'true_features': true_feature_set
        }
    }

def main():
    """Run all demonstrations"""
    print("DAAadvisor Advanced Analysis Demonstrations")
    print("==========================================")
    
    try:
        # 1. Test R methods integration
        r_results = demonstrate_r_methods()
        
        # 2. Test information theory framework  
        info_results = demonstrate_information_theory()
        
        # 3. Test unified analysis
        unified_results = demonstrate_unified_analysis()
        
        print("\n" + "="*60)
        print("ALL DEMONSTRATIONS COMPLETED SUCCESSFULLY!")
        print("="*60)
        
        print("\nKey Insights:")
        print("• R methods provide complementary approaches with different assumptions")
        print("• Information theory framework offers principled method selection")
        print("• Unified analysis can improve feature discovery robustness")
        print("• Different methods may find different but valid biological signals")
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Some features may require additional dependencies (rpy2, R packages)")
        
    except Exception as e:
        print(f"❌ Error during demonstration: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()