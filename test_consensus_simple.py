#!/usr/bin/env python3
"""
Simple test for advanced consensus analysis integration
"""

import numpy as np
import pandas as pd
from daa_advisor import DifferentialAbundanceTool, MicrobiomeDataGenerator

def test_consensus_simple():
    """Test the advanced consensus with methods that work reliably"""
    
    print("🧪 Testing Advanced Consensus Analysis (Simple)")
    print("=" * 60)
    
    # Generate smaller, cleaner test data
    generator = MicrobiomeDataGenerator()
    count_table, metadata, differential_features = generator.generate_asv_data(
        n_samples=40,
        n_features=50,
        effect_size=3.0,
        n_differential=10,
        sparsity=0.6  # Less sparsity for more reliable results
    )
    
    print(f"📊 Generated test data: {count_table.shape} samples x features")
    print(f"🎯 {metadata['condition'].value_counts().to_dict()} samples per group")
    
    # Manual analysis with reliable methods to test consensus
    from daa_advisor.methods import MethodRegistry
    registry = MethodRegistry()
    
    # Run only reliable methods
    reliable_methods = ['wilcoxon', 'deseq2', 'edger']
    analyses = {}
    
    for method_name in reliable_methods:
        if registry.has_method(method_name):
            try:
                print(f"\n🔬 Running {method_name}...")
                method = registry.get_method(method_name)
                result = method.run(count_table, metadata)
                analyses[method_name] = result
                
                # Count significant results
                padj_col = 'padj' if 'padj' in result.columns else 'qvalue'
                n_sig = (result[padj_col] < 0.05).sum() if padj_col in result.columns else 0
                print(f"   ✅ {method_name}: {n_sig} significant features")
                
            except Exception as e:
                print(f"   ❌ {method_name} failed: {e}")
    
    # Test advanced consensus directly
    if len(analyses) >= 2:
        print(f"\n🤝 Testing Advanced Consensus with {len(analyses)} methods...")
        
        from daa_advisor.consensus import AdvancedConsensusAnalyzer
        
        # Test different voting strategies
        for strategy in ['simple', 'weighted', 'ranked']:
            print(f"\n📊 Testing {strategy} voting strategy:")
            
            analyzer = AdvancedConsensusAnalyzer(voting_strategy=strategy)
            consensus_results = analyzer.generate_advanced_consensus(analyses)
            
            consensus = consensus_results['consensus_features']
            metrics = consensus_results['agreement_metrics']
            
            print(f"   • Total features: {len(consensus)}")
            print(f"   • Consensus significant: {consensus['consensus_significant'].sum()}")
            
            if 'confidence_score' in consensus.columns:
                print(f"   • High confidence calls: {consensus['high_confidence'].sum()}")
                print(f"   • Mean confidence: {consensus['confidence_score'].mean():.3f}")
            
            if 'consensus_strength' in consensus.columns:
                strength_dist = consensus['consensus_strength'].value_counts()
                print(f"   • Consensus strength: {dict(strength_dist)}")
            
            # Agreement metrics
            if 'mean_kappa' in metrics:
                print(f"   • Mean Cohen's κ: {metrics['mean_kappa']:.3f} ({metrics['kappa_interpretation']})")
            
            if 'majority_agreement' in metrics:
                print(f"   • Majority agreement: {metrics['majority_agreement']:.1%}")
    
    else:
        print(f"❌ Need at least 2 methods for consensus, got {len(analyses)}")
    
    print(f"\n✅ Advanced consensus test completed!")
    return analyses

if __name__ == "__main__":
    test_consensus_simple()