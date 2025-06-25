#!/usr/bin/env python3
"""
Create DAAadvisor Workflow Diagram using Graphviz
Professional, publication-quality diagrams
"""

try:
    import graphviz
    GRAPHVIZ_AVAILABLE = True
except ImportError:
    GRAPHVIZ_AVAILABLE = False

def create_core_workflow():
    """Create core analysis workflow using Graphviz"""
    
    if not GRAPHVIZ_AVAILABLE:
        print("❌ Graphviz not available. Install with: pip install graphviz")
        return None
    
    # Create directed graph
    dot = graphviz.Digraph(comment='DAAadvisor Core Workflow')
    dot.attr(rankdir='TB', size='10,12', dpi='300')
    dot.attr('node', shape='box', style='rounded,filled', fontname='Arial', fontsize='11')
    dot.attr('edge', fontname='Arial', fontsize='9', color='#2E3440')
    
    # Define modern color scheme
    colors = {
        'input': '#88C0D0',
        'profiling': '#D08770', 
        'selection': '#A3BE8C',
        'analysis': '#5E81AC',
        'consensus': '#B48EAD'
    }
    
    # Add nodes
    dot.node('input', 
             '📊 Data Input\\n\\nCount Table: samples × features matrix\\nMetadata: experimental conditions\\nSupports: ASV/16S, gene, viral data', 
             fillcolor=colors['input'], fontcolor='#2E3440')
    
    dot.node('profiling', 
             '🔍 Comprehensive Data Profiling\\n\\nStatistical characterization:\\n• Zero-inflation: quantify sparsity (0-95%)\\n• Distribution analysis: mean, variance, skewness\\n• Library size normalization assessment\\n• Compositional bias detection\\n• Data type classification (ASV/gene/viral)', 
             fillcolor=colors['profiling'], fontcolor='#2E3440')
    
    dot.node('selection', 
             '🧮 Information-Theoretic Method Selection\\n\\nMathematical optimization framework:\\n• Maximum entropy principle: Method* = argmax H(X|θ)\\n• Jensen-Shannon divergence: between-group analysis\\n• CLR transformation compatibility assessment\\n• Method suitability scoring (0-100%)\\n• Uncertainty quantification', 
             fillcolor=colors['selection'], fontcolor='#2E3440')
    
    dot.node('analysis', 
             '⚡ Parallel Multi-Method Analysis\\n\\nSimultaneous statistical testing:\\n• Wilcoxon: rank-based non-parametric\\n• ALDEx2: CLR + Monte Carlo (n=128)\\n• DESeq2: negative binomial GLM\\n• edgeR: TMM normalization + QLF\\n• metagenomeSeq: zero-inflated log-normal\\n• ANCOM-BC: bias-corrected compositional', 
             fillcolor=colors['analysis'], fontcolor='white')
    
    dot.node('consensus', 
             '🤝 Advanced Consensus Integration\\n\\nMulti-method result synthesis:\\n• Voting strategies: simple, weighted, ranked\\n• Cohen kappa: inter-method agreement\\n• Effect size concordance analysis\\n• P-value consistency assessment\\n• Confidence intervals: bootstrap (n=50-100)\\n• Final recommendation scoring', 
             fillcolor=colors['consensus'], fontcolor='#2E3440')
    
    # Add edges with labels
    dot.edge('input', 'profiling', label='Raw microbiome\\ncount matrix')
    dot.edge('input', 'selection', label='Data characteristics\\n& constraints')
    dot.edge('profiling', 'analysis', label='Statistical profile\\n& recommendations')
    dot.edge('selection', 'consensus', label='Optimal method\\nselection & weights')
    dot.edge('analysis', 'consensus', label='Individual method\\nresults & p-values')
    
    return dot

def create_validation_workflow():
    """Create validation workflow using Graphviz"""
    
    if not GRAPHVIZ_AVAILABLE:
        return None
    
    # Create directed graph
    dot = graphviz.Digraph(comment='DAAadvisor Validation Workflow')
    dot.attr(rankdir='TB', size='10,12', dpi='300')
    dot.attr('node', shape='box', style='rounded,filled', fontname='Arial', fontsize='11')
    dot.attr('edge', fontname='Arial', fontsize='9', color='#2E3440')
    
    # Define colors
    colors = {
        'input': '#ECEFF4',
        'real_data': '#88C0D0',
        'synthetic': '#D08770',
        'cross_val': '#A3BE8C',
        'validation': '#5E81AC',
        'reporting': '#B48EAD',
        'output': '#8FBCBB'
    }
    
    # Add nodes
    dot.node('input', 
             '⬆️ From Core Analysis\\n\\nStatistical results: p-values, effect sizes\\nMethod performance: success rates, runtimes\\nConsensus calls: agreement scores', 
             fillcolor=colors['input'], fontcolor='#2E3440')
    
    dot.node('real_data', 
             '🧬 Real-World Data Integration\\n\\nData source: curatedMetagenomicData (Bioconductor)\\n• IBD study: 1,627 samples (1,201 IBD + 426 controls)\\n• Ground truth: published biomarkers (Faecalibacterium, Escherichia)\\n• Data format: TreeSummarizedExperiment\\n• Quality control: sample alignment validation\\n• Metadata: disease status, demographics', 
             fillcolor=colors['real_data'], fontcolor='#2E3440')
    
    dot.node('synthetic', 
             '🎭 Controlled Data Simulation\\n\\nRealistic microbiome simulation:\\n• Literature effect sizes: IBD (log2FC: 0.5-2.0)\\n• Compositional structure: Dirichlet-multinomial\\n• Known differential features: 50-200 taxa\\n• Sparsity patterns: 70-90% zeros\\n• Sample sizes: 20-200 per group\\n• Controlled confounders: batch, age, gender', 
             fillcolor=colors['synthetic'], fontcolor='#2E3440')
    
    dot.node('cross_val', 
             '🔄 Cross-Validation Framework\\n\\nMethod performance evaluation:\\n• Real vs synthetic correlation: Pearson r\\n• Bootstrap iterations: 50-100 replicates\\n• Performance metrics: sensitivity, specificity, F1\\n• Ground truth recovery: known feature detection\\n• Robustness: consistency across datasets\\n• Statistical testing: paired t-tests, Wilcoxon', 
             fillcolor=colors['cross_val'], fontcolor='#2E3440')
    
    dot.node('validation', 
             '🏆 Publication-Quality Validation\\n\\nScientific rigor assessment:\\n• Statistical confidence: 95% bootstrap CI\\n• Literature confirmation: biomarker validation\\n• Performance metrics: AUROC, AUPRC, MCC\\n• Effect size validation: Cohen d, log2FC\\n• FDR control: Benjamini-Hochberg correction\\n• Reproducibility: seed-controlled analysis', 
             fillcolor=colors['validation'], fontcolor='white')
    
    dot.node('reporting', 
             '📊 Comprehensive Reporting\\n\\nPublication-ready documentation:\\n• Interactive HTML dashboards (plotly)\\n• High-resolution figures (300 DPI, PNG/SVG)\\n• Statistical tables: CSV/Excel export\\n• Method comparison matrices\\n• Bootstrap result summaries\\n• Reproducible analysis scripts', 
             fillcolor=colors['reporting'], fontcolor='#2E3440')
    
    dot.node('output', 
             '🎯 Validated Scientific Results\\n\\nPeer-review ready output:\\n• Cross-validated performance (real + synthetic)\\n• Literature-confirmed biomarker detection\\n• Statistical significance: p < 0.05, FDR < 0.1\\n• Effect size confidence intervals\\n• Method reliability assessment\\n• Complete reproducibility documentation', 
             fillcolor=colors['output'], fontcolor='#2E3440')
    
    # Add edges
    dot.edge('input', 'real_data', label='Core analysis results\\nfor validation testing')
    dot.edge('input', 'synthetic', label='Method performance\\nbaseline metrics')
    dot.edge('real_data', 'cross_val', label='Real dataset + literature\\nbiomarker ground truth')
    dot.edge('synthetic', 'cross_val', label='Synthetic dataset +\\nknown differential features')
    dot.edge('cross_val', 'validation', label='Cross-validation metrics\\n& performance correlation')
    dot.edge('cross_val', 'reporting', label='Bootstrap results &\\nstatistical summaries')
    dot.edge('validation', 'output', label='Publication-quality\\nvalidation results')
    dot.edge('reporting', 'output', label='Complete documentation\\n& visualizations')
    
    return dot

def main():
    """Create and save workflow diagrams"""
    
    if not GRAPHVIZ_AVAILABLE:
        print("❌ Graphviz not available.")
        print("💡 Install with: pip install graphviz")
        print("💡 Also need system Graphviz: brew install graphviz (macOS) or apt-get install graphviz (Ubuntu)")
        return
    
    try:
        # Create core workflow
        core_dot = create_core_workflow()
        if core_dot:
            core_dot.render('daaadvisor_graphviz_core_workflow', format='png', cleanup=True)
            core_dot.render('daaadvisor_graphviz_core_workflow', format='svg', cleanup=True)
            print("✅ Core workflow created: daaadvisor_graphviz_core_workflow.png/.svg")
        
        # Create validation workflow
        validation_dot = create_validation_workflow()
        if validation_dot:
            validation_dot.render('daaadvisor_graphviz_validation_workflow', format='png', cleanup=True)
            validation_dot.render('daaadvisor_graphviz_validation_workflow', format='svg', cleanup=True)
            print("✅ Validation workflow created: daaadvisor_graphviz_validation_workflow.png/.svg")
            
    except Exception as e:
        print(f"❌ Error creating diagrams: {e}")

if __name__ == "__main__":
    main()