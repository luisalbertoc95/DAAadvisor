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
             '📊 Data Input\\n\\nCount Table + Metadata\\nYour microbiome dataset', 
             fillcolor=colors['input'], fontcolor='#2E3440')
    
    dot.node('profiling', 
             '🔍 Smart Data Profiling\\n\\nAutomatic characterization\\n• Sparsity assessment\\n• Data type detection\\n• Quality metrics', 
             fillcolor=colors['profiling'], fontcolor='#2E3440')
    
    dot.node('selection', 
             '🧠 AI-Powered Method Selection\\n\\nInformation theory optimization\\n• Maximum entropy principle\\n• Divergence analysis\\n• Confidence scoring', 
             fillcolor=colors['selection'], fontcolor='#2E3440')
    
    dot.node('analysis', 
             '⚡ Parallel Statistical Testing\\n\\nMulti-method analysis\\n• Wilcoxon (non-parametric)\\n• ALDEx2 (compositional)\\n• DESeq2/edgeR (count-based)', 
             fillcolor=colors['analysis'], fontcolor='white')
    
    dot.node('consensus', 
             '🤝 Intelligent Result Integration\\n\\nAdvanced consensus\\n• Voting strategies\\n• Agreement metrics\\n• Confidence scoring', 
             fillcolor=colors['consensus'], fontcolor='#2E3440')
    
    # Add edges with labels
    dot.edge('input', 'profiling', label='Raw\\nMicrobiome\\nData')
    dot.edge('input', 'selection', label='Dataset\\nCharacteristics')
    dot.edge('profiling', 'analysis', label='Data\\nProfile')
    dot.edge('selection', 'consensus', label='Method\\nRecommendation')
    dot.edge('analysis', 'consensus', label='Statistical\\nResults')
    
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
             '⬆️ From Core Analysis\\n\\nStatistical results\\nMethod performance', 
             fillcolor=colors['input'], fontcolor='#2E3440')
    
    dot.node('real_data', 
             '🧬 Real-World Data\\n\\ncuratedMetagenomicData\\n• IBD: 1,627 samples\\n• Literature biomarkers\\n• Quality control', 
             fillcolor=colors['real_data'], fontcolor='#2E3440')
    
    dot.node('synthetic', 
             '🎭 Controlled Simulation\\n\\nRealistic synthetic data\\n• Literature-based effects\\n• Known ground truth\\n• Controlled conditions', 
             fillcolor=colors['synthetic'], fontcolor='#2E3440')
    
    dot.node('cross_val', 
             '🔄 Cross-Validation Engine\\n\\nPerformance comparison\\n• Real vs synthetic\\n• Bootstrap testing\\n• Robustness assessment', 
             fillcolor=colors['cross_val'], fontcolor='#2E3440')
    
    dot.node('validation', 
             '🏆 Publication Validation\\n\\nScientific rigor\\n• Statistical confidence\\n• Literature confirmation\\n• Comprehensive metrics', 
             fillcolor=colors['validation'], fontcolor='white')
    
    dot.node('reporting', 
             '📊 Comprehensive Reporting\\n\\nPublication-ready output\\n• Interactive dashboards\\n• Journal figures\\n• Statistical summaries', 
             fillcolor=colors['reporting'], fontcolor='#2E3440')
    
    dot.node('output', 
             '🎯 Publication-Ready Results\\n\\nValidated scientific output\\n• Real-world tested\\n• Literature confirmed\\n• Statistical confidence', 
             fillcolor=colors['output'], fontcolor='#2E3440')
    
    # Add edges
    dot.edge('input', 'real_data', label='Core Results\\nfor Testing')
    dot.edge('input', 'synthetic', label='Performance\\nBaseline')
    dot.edge('real_data', 'cross_val', label='Real Dataset\\n+ Biomarkers')
    dot.edge('synthetic', 'cross_val', label='Synthetic Dataset\\n+ Known Truth')
    dot.edge('cross_val', 'validation', label='Performance\\nMetrics')
    dot.edge('cross_val', 'reporting', label='Validation\\nResults')
    dot.edge('validation', 'output', label='Scientific\\nValidation')
    dot.edge('reporting', 'output', label='Reports +\\nVisualizations')
    
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