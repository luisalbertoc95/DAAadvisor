#!/usr/bin/env python3
"""
Create DAAadvisor Workflow Diagram
Comprehensive visual representation of the 8-step framework
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
import numpy as np

# Set up the figure with high resolution for GitHub
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'

fig, ax = plt.subplots(1, 1, figsize=(16, 20))
ax.set_xlim(0, 10)
ax.set_ylim(0, 24)
ax.axis('off')

# Color scheme
colors = {
    'input': '#E3F2FD',      # Light blue
    'profiling': '#F3E5F5',   # Light purple
    'selection': '#E8F5E8',   # Light green
    'analysis': '#FFF3E0',    # Light orange
    'consensus': '#FCE4EC',   # Light pink
    'real_data': '#E0F2F1',   # Light teal
    'cross_val': '#FFF8E1',   # Light yellow
    'validation': '#FFEBEE',  # Light red
    'output': '#F1F8E9'       # Light lime
}

def create_step_box(x, y, width, height, title, details, color, step_num):
    """Create a detailed step box with title and bullet points"""
    
    # Main box
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.1",
        facecolor=color,
        edgecolor='black',
        linewidth=1.5
    )
    ax.add_patch(box)
    
    # Step number circle
    circle = plt.Circle((x + 0.3, y + height - 0.3), 0.2, 
                       facecolor='white', edgecolor='black', linewidth=2)
    ax.add_patch(circle)
    ax.text(x + 0.3, y + height - 0.3, str(step_num), 
            ha='center', va='center', fontweight='bold', fontsize=12)
    
    # Title
    ax.text(x + 0.7, y + height - 0.3, title, 
            ha='left', va='center', fontweight='bold', fontsize=12)
    
    # Details
    detail_y = y + height - 0.8
    for detail in details:
        ax.text(x + 0.2, detail_y, f"• {detail}", 
                ha='left', va='top', fontsize=9, wrap=True)
        detail_y -= 0.4

def create_arrow(start_x, start_y, end_x, end_y, style='->'):
    """Create connection arrow between steps"""
    arrow = ConnectionPatch(
        (start_x, start_y), (end_x, end_y), 
        "data", "data",
        arrowstyle=style, shrinkA=5, shrinkB=5,
        mutation_scale=20, fc="black", linewidth=2
    )
    ax.add_patch(arrow)

# Title
ax.text(5, 23, 'DAAadvisor: Complete Methodology Workflow', 
        ha='center', va='center', fontsize=20, fontweight='bold')
ax.text(5, 22.5, '8-Step Framework with Cross-Validation & Real Data Integration', 
        ha='center', va='center', fontsize=14, style='italic')

# Step 1: Data Input & Profiling
create_step_box(
    0.5, 20, 4, 2.5,
    "📊 Data Assessment & Profiling",
    [
        "Sparsity Analysis: Zero-inflation (83%)",
        "Data Type: ASV/Gene/Viral detection",
        "Compositional Bias: Library size variation",
        "Information Preprocessing: Adaptive thresholds"
    ],
    colors['profiling'], 1
)

# Step 2: Information-Theoretic Selection
create_step_box(
    5.5, 20, 4, 2.5,
    "🧮 Information-Theoretic Selection",
    [
        "Maximum Entropy: Method* = argmax H(X|θ)",
        "Jensen-Shannon Divergence: Between-group",
        "CLR Transformation: log(x/g(x))",
        "Confidence Scoring: Selection certainty"
    ],
    colors['selection'], 2
)

# Step 3: Multi-Method Analysis
create_step_box(
    0.5, 17, 4, 2.5,
    "🔬 Multi-Method Statistical Analysis",
    [
        "Wilcoxon: Non-parametric (100% success)",
        "ALDEx2: CLR + Monte Carlo sampling",
        "DESeq2: Negative binomial modeling",
        "edgeR/metagenomeSeq: TMM/Zero-inflated"
    ],
    colors['analysis'], 3
)

# Step 4: Advanced Consensus
create_step_box(
    5.5, 17, 4, 2.5,
    "🤝 Advanced Consensus Analysis",
    [
        "Sophisticated Voting: Weighted reliability",
        "Cohen's Kappa: κ = 0.436 agreement",
        "Confidence Scoring: Method concordance",
        "Consensus Strength: Strong/Moderate/Weak"
    ],
    colors['consensus'], 4
)

# Step 5: Real Data Integration
create_step_box(
    0.5, 14, 4, 2.5,
    "🧬 Real Data Integration",
    [
        "curatedMetagenomicData: R/Bioconductor",
        "IBD Dataset: 1,627 samples validated",
        "Ground Truth: Literature biomarkers",
        "Quality Control: Sample alignment"
    ],
    colors['real_data'], 5
)

# Step 6: Cross-Validation Framework
create_step_box(
    5.5, 14, 4, 2.5,
    "🔄 Cross-Validation Framework",
    [
        "Real vs Synthetic: Performance correlation",
        "Bootstrap: 50-100 iterations",
        "Ground Truth Recovery: Method validation",
        "Robustness: Consistency evaluation"
    ],
    colors['cross_val'], 6
)

# Step 7: Publication Validation
create_step_box(
    0.5, 11, 4, 2.5,
    "🏆 Publication-Quality Validation",
    [
        "Statistical Rigor: Confidence intervals",
        "Literature Confirmation: Published biomarkers",
        "Comprehensive Metrics: AUROC, F1, AUPRC",
        "Real-World: 1,627 IBD samples tested"
    ],
    colors['validation'], 7
)

# Step 8: Results & Reporting
create_step_box(
    5.5, 11, 4, 2.5,
    "📈 Results & Comprehensive Reporting",
    [
        "Interactive Dashboards: HTML reports",
        "Cross-Validation: Real vs synthetic plots",
        "Publication Figures: Journal-ready",
        "Bootstrap Results: Statistical significance"
    ],
    colors['output'], 8
)

# Add connecting arrows
# Vertical flow
create_arrow(2.5, 20, 2.5, 19.5)  # 1 to 3
create_arrow(7.5, 20, 7.5, 19.5)  # 2 to 4
create_arrow(2.5, 17, 2.5, 16.5)  # 3 to 5
create_arrow(7.5, 17, 7.5, 16.5)  # 4 to 6
create_arrow(2.5, 14, 2.5, 13.5)  # 5 to 7
create_arrow(7.5, 14, 7.5, 13.5)  # 6 to 8

# Horizontal connections
create_arrow(4.5, 21.2, 5.5, 21.2)  # 1 to 2
create_arrow(4.5, 18.2, 5.5, 18.2)  # 3 to 4
create_arrow(4.5, 15.2, 5.5, 15.2)  # 5 to 6
create_arrow(4.5, 12.2, 5.5, 12.2)  # 7 to 8

# Cross-validation connections (showing integration)
create_arrow(2.5, 16.5, 7, 15.5, style='->')  # Real data to cross-val
create_arrow(7, 15.5, 2.5, 13.5, style='->')  # Cross-val to validation

# Add key achievements box
achievements_box = FancyBboxPatch(
    (1, 8), 8, 2.5,
    boxstyle="round,pad=0.2",
    facecolor='#F5F5F5',
    edgecolor='black',
    linewidth=2
)
ax.add_patch(achievements_box)

ax.text(5, 9.8, '🎉 Key Achievements & Validation Results', 
        ha='center', va='center', fontsize=14, fontweight='bold')

achievements_text = [
    "✅ 100% Complete Framework: All 8 steps implemented and validated",
    "✅ Real Data Integration: 1,627 IBD samples (1,201 IBD + 426 controls) from curatedMetagenomicData",
    "✅ Cross-Validation Pipeline: Real vs synthetic data comparison with bootstrap validation",
    "✅ Literature Validation: Known biomarkers confirmed (Faecalibacterium, Escherichia)",
    "✅ Publication Ready: Statistical rigor with confidence intervals and comprehensive reporting",
    "✅ Method Performance: 6/6 statistical methods functional with consensus analysis"
]

y_pos = 9.4
for achievement in achievements_text:
    ax.text(1.2, y_pos, achievement, ha='left', va='center', fontsize=10)
    y_pos -= 0.3

# Add technical specs box
specs_box = FancyBboxPatch(
    (1, 5), 8, 2.5,
    boxstyle="round,pad=0.2",
    facecolor='#E8F5E8',
    edgecolor='black',
    linewidth=2
)
ax.add_patch(specs_box)

ax.text(5, 6.8, '🔬 Technical Specifications', 
        ha='center', va='center', fontsize=14, fontweight='bold')

specs_text = [
    "📊 Data Types: ASV/16S, Gene/Functional, Viral microbiome data",
    "🧮 Methods: Wilcoxon, ALDEx2, DESeq2, edgeR, metagenomeSeq, ANCOM-BC",
    "📈 Validation: Bootstrap confidence intervals (50-100 iterations)",
    "🧬 Real Data: curatedMetagenomicData (Bioconductor) with automated R integration",
    "🎯 Consensus: Cohen's kappa agreement (κ = 0.436) with weighted voting",
    "📋 Output: Interactive HTML dashboards + publication-ready figures"
]

y_pos = 6.4
for spec in specs_text:
    ax.text(1.2, y_pos, spec, ha='left', va='center', fontsize=10)
    y_pos -= 0.3

# Add usage examples
usage_box = FancyBboxPatch(
    (1, 2), 8, 2.5,
    boxstyle="round,pad=0.2",
    facecolor='#FFF3E0',
    edgecolor='black',
    linewidth=2
)
ax.add_patch(usage_box)

ax.text(5, 3.8, '🚀 Usage Examples', 
        ha='center', va='center', fontsize=14, fontweight='bold')

usage_text = [
    "# Basic analysis with intelligent method selection",
    "from daa_advisor import DifferentialAbundanceTool",
    "tool = DifferentialAbundanceTool()",
    "results = tool.analyze(count_table, metadata, use_consensus=True)",
    "",
    "# Cross-validation with real data",
    "python run_cross_validation_benchmark.py --max-conditions 1"
]

y_pos = 3.4
for usage in usage_text:
    if usage.startswith('#'):
        ax.text(1.2, y_pos, usage, ha='left', va='center', fontsize=10, 
                style='italic', color='green')
    elif usage.startswith('python'):
        ax.text(1.2, y_pos, usage, ha='left', va='center', fontsize=10, 
                family='monospace', color='blue')
    elif usage == "":
        y_pos -= 0.1
        continue
    else:
        ax.text(1.2, y_pos, usage, ha='left', va='center', fontsize=10, 
                family='monospace')
    y_pos -= 0.25

# Add footer
ax.text(5, 0.5, 'DAAadvisor: Intelligent Differential Abundance Analysis for Microbiome Data', 
        ha='center', va='center', fontsize=12, fontweight='bold')
ax.text(5, 0.1, 'Complete framework with cross-validation and real data integration', 
        ha='center', va='center', fontsize=10, style='italic')

plt.tight_layout()
plt.savefig('daaadvisor_workflow_diagram.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Workflow diagram created: daaadvisor_workflow_diagram.png")
print("📊 High-resolution PNG ready for GitHub")
print("🎯 Comprehensive 8-step framework visualized")