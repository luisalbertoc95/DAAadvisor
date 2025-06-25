#!/usr/bin/env python3
"""
Demo test with artificial results to showcase advanced consensus features
"""

import numpy as np
import pandas as pd
from daa_advisor.consensus import AdvancedConsensusAnalyzer

def create_demo_results():
    """Create artificial method results for demonstration"""
    
    features = [f"Feature_{i+1}" for i in range(20)]
    
    # Create realistic method results with some overlap
    wilcoxon_results = pd.DataFrame({
        'feature': features,
        'pvalue': np.random.uniform(0.001, 0.3, 20),
        'padj': np.random.uniform(0.01, 0.2, 20),
        'log2fc': np.random.normal(0, 2, 20)
    })
    
    deseq2_results = pd.DataFrame({
        'feature': features,
        'pvalue': np.random.uniform(0.001, 0.4, 20),
        'padj': np.random.uniform(0.01, 0.3, 20),
        'log2fc': np.random.normal(0, 1.5, 20)
    })
    
    edger_results = pd.DataFrame({
        'feature': features,
        'pvalue': np.random.uniform(0.001, 0.5, 20),
        'padj': np.random.uniform(0.01, 0.4, 20),
        'log2fc': np.random.normal(0, 1.8, 20)
    })
    
    # Make some features clearly significant across methods
    significant_features = ['Feature_1', 'Feature_3', 'Feature_7', 'Feature_12']
    
    for df in [wilcoxon_results, deseq2_results, edger_results]:
        mask = df['feature'].isin(significant_features)
        df.loc[mask, 'pvalue'] = np.random.uniform(0.0001, 0.01, mask.sum())
        df.loc[mask, 'padj'] = np.random.uniform(0.001, 0.04, mask.sum())
        df.loc[mask, 'log2fc'] = np.random.uniform(1.5, 3.0, mask.sum())
    
    return {
        'wilcoxon': wilcoxon_results,
        'deseq2': deseq2_results,
        'edger': edger_results
    }

def demo_advanced_consensus():
    """Demonstrate advanced consensus analysis features"""
    
    print("🧪 Advanced Consensus Analysis Demo")
    print("=" * 50)
    
    # Create demo results
    analyses = create_demo_results()
    
    print(f"📊 Demo data: {len(analyses)} methods, {len(analyses['wilcoxon'])} features")
    
    # Show method-specific significant features
    for method, results in analyses.items():
        n_sig = (results['padj'] < 0.05).sum()
        print(f"   • {method}: {n_sig} significant features")
    
    # Test all voting strategies
    strategies = ['simple', 'weighted', 'ranked']
    
    for strategy in strategies:
        print(f"\n🗳️  {strategy.upper()} VOTING STRATEGY")
        print("-" * 40)
        
        analyzer = AdvancedConsensusAnalyzer(
            voting_strategy=strategy,
            confidence_threshold=0.7
        )
        
        results = analyzer.generate_advanced_consensus(analyses)
        
        consensus = results['consensus_features']
        metrics = results['agreement_metrics']
        
        # Consensus results
        print(f"📈 Consensus Results:")
        print(f"   • Total features analyzed: {len(consensus)}")
        print(f"   • Consensus significant: {consensus['consensus_significant'].sum()}")
        
        if 'high_confidence' in consensus.columns:
            print(f"   • High confidence calls: {consensus['high_confidence'].sum()}")
        
        if 'confidence_score' in consensus.columns:
            print(f"   • Mean confidence score: {consensus['confidence_score'].mean():.3f}")
        
        # Show top consensus features
        top_features = consensus.nlargest(5, 'n_significant')[['feature', 'n_significant', 'consensus_significant']]
        if 'confidence_score' in consensus.columns:
            top_features = consensus.nlargest(5, 'confidence_score')[['feature', 'n_significant', 'consensus_significant', 'confidence_score']]
        
        print(f"   • Top features:")
        for _, row in top_features.iterrows():
            if 'confidence_score' in row:
                print(f"     - {row['feature']}: {row['n_significant']}/{len(analyses)} methods, conf={row['confidence_score']:.3f}")
            else:
                print(f"     - {row['feature']}: {row['n_significant']}/{len(analyses)} methods")
        
        # Agreement metrics
        print(f"🤝 Agreement Metrics:")
        if 'mean_kappa' in metrics and not pd.isna(metrics['mean_kappa']):
            print(f"   • Mean Cohen's κ: {metrics['mean_kappa']:.3f} ({metrics['kappa_interpretation']})")
        
        if 'majority_agreement' in metrics:
            print(f"   • Majority agreement: {metrics['majority_agreement']:.1%}")
        
        if 'unanimous_agreement' in metrics:
            print(f"   • Unanimous agreement: {metrics['unanimous_agreement']:.1%}")
        
        # Consensus strength distribution
        if 'consensus_strength' in consensus.columns:
            strength_dist = consensus['consensus_strength'].value_counts()
            print(f"   • Consensus strength:")
            for strength, count in strength_dist.items():
                print(f"     - {strength}: {count}")
    
    print(f"\n✅ Advanced consensus demo completed!")
    
    return results

if __name__ == "__main__":
    demo_advanced_consensus()