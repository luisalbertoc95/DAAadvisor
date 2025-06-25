#!/usr/bin/env python3
"""
Test script for advanced consensus analysis integration
"""

import numpy as np
import pandas as pd
from daa_advisor import DifferentialAbundanceTool, MicrobiomeDataGenerator

def test_advanced_consensus():
    """Test the integrated advanced consensus analysis"""
    
    print("🧪 Testing Advanced Consensus Analysis Integration")
    print("=" * 60)
    
    # Generate test data
    generator = MicrobiomeDataGenerator()
    count_table, metadata, differential_features = generator.generate_asv_data(
        n_samples=50,
        n_features=100,
        effect_size=2.0,
        n_differential=20
    )
    
    print(f"📊 Generated test data: {count_table.shape} samples x features")
    print(f"🎯 {metadata['condition'].value_counts().to_dict()} samples per group")
    print(f"🔬 {len(differential_features)} differential features generated")
    
    # Run analysis with consensus
    tool = DifferentialAbundanceTool()
    results = tool.analyze(
        count_table=count_table,
        metadata=metadata,
        data_type='asv',
        use_consensus=True
    )
    
    print(f"\n🔬 Methods run: {list(results['analyses'].keys())}")
    
    # Test consensus results
    if 'consensus' in results:
        consensus = results['consensus']
        print(f"\n🤝 Advanced Consensus Results:")
        print(f"   • Total features analyzed: {len(consensus)}")
        print(f"   • Consensus significant: {consensus['consensus_significant'].sum()}")
        
        if 'confidence_score' in consensus.columns:
            print(f"   • High confidence calls: {consensus['high_confidence'].sum()}")
            print(f"   • Mean confidence score: {consensus['confidence_score'].mean():.3f}")
        
        if 'consensus_strength' in consensus.columns:
            strength_dist = consensus['consensus_strength'].value_counts()
            print(f"   • Consensus strength distribution:")
            for strength, count in strength_dist.items():
                print(f"     - {strength}: {count}")
    
    # Test agreement metrics
    if 'consensus_metrics' in results:
        metrics = results['consensus_metrics']
        print(f"\n📈 Agreement Metrics:")
        
        if 'mean_kappa' in metrics:
            print(f"   • Mean Cohen's κ: {metrics['mean_kappa']:.3f} ({metrics['kappa_interpretation']})")
        
        if 'majority_agreement' in metrics:
            print(f"   • Majority agreement: {metrics['majority_agreement']:.1%}")
        
        if 'unanimous_agreement' in metrics:
            print(f"   • Unanimous agreement: {metrics['unanimous_agreement']:.1%}")
        
        if 'method_agreement' in metrics:
            print(f"   • Method-specific agreement:")
            for method, agreement in metrics['method_agreement'].items():
                print(f"     - {method}: {agreement:.3f}")
    
    # Display summary
    print(f"\n📋 Analysis Summary:")
    tool.summarize_results()
    
    print(f"\n✅ Advanced consensus analysis integration test completed!")
    return results

if __name__ == "__main__":
    test_advanced_consensus()