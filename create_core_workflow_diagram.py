#!/usr/bin/env python3
"""
Create DAAadvisor Core Workflow Diagram
Main 8-step analysis framework
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

fig, ax = plt.subplots(1, 1, figsize=(16, 18))
ax.set_xlim(0, 12)
ax.set_ylim(0, 20)
ax.axis('off')

# Color scheme
colors = {
    'input': '#E3F2FD',      # Light blue
    'profiling': '#F3E5F5',   # Light purple
    'selection': '#E8F5E8',   # Light green
    'analysis': '#FFF3E0',    # Light orange
    'consensus': '#FCE4EC',   # Light pink
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
        detail_y -= 0.35

def create_arrow(start_x, start_y, end_x, end_y, style='->', label='', label_pos='mid'):
    """Create connection arrow between steps with optional label"""
    arrow = ConnectionPatch(
        (start_x, start_y), (end_x, end_y), 
        "data", "data",
        arrowstyle=style, shrinkA=5, shrinkB=5,
        mutation_scale=20, fc="black", linewidth=2
    )
    ax.add_patch(arrow)
    
    # Add label if provided
    if label:
        if label_pos == 'mid':
            label_x = (start_x + end_x) / 2
            label_y = (start_y + end_y) / 2
        elif label_pos == 'start':
            label_x = start_x + 0.2
            label_y = start_y + 0.1
        else:  # 'end'
            label_x = end_x - 0.2
            label_y = end_y - 0.1
            
        ax.text(label_x, label_y, label, ha='center', va='center', 
                fontsize=8, bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))

# Title
ax.text(6, 19, 'DAAadvisor: Core Analysis Framework', 
        ha='center', va='center', fontsize=18, fontweight='bold')
ax.text(6, 18.5, 'Intelligent Method Selection & Multi-Method Analysis Pipeline', 
        ha='center', va='center', fontsize=12, style='italic')

# Step 1: Data Input & Profiling
create_step_box(
    0.5, 15.5, 4.5, 2.5,
    "📊 Data Assessment & Profiling",
    [
        "Sparsity analysis & zero-inflation detection",
        "Data type classification (ASV/Gene/Viral)", 
        "Compositional bias assessment",
        "Library size & distribution profiling"
    ],
    colors['profiling'], 1
)

# Step 2: Information-Theoretic Selection
create_step_box(
    7, 15.5, 4.5, 2.5,
    "🧮 Information-Theoretic Selection",
    [
        "Maximum entropy principle optimization",
        "Jensen-Shannon divergence calculation",
        "CLR transformation assessment",
        "Method confidence scoring"
    ],
    colors['selection'], 2
)

# Step 3: Multi-Method Analysis
create_step_box(
    0.5, 12, 4.5, 2.5,
    "🔬 Multi-Method Statistical Analysis",
    [
        "Wilcoxon: Non-parametric testing",
        "ALDEx2: CLR + Monte Carlo sampling",
        "DESeq2: Negative binomial modeling", 
        "edgeR/metagenomeSeq: TMM/Zero-inflated"
    ],
    colors['analysis'], 3
)

# Step 4: Advanced Consensus
create_step_box(
    7, 12, 4.5, 2.5,
    "🤝 Advanced Consensus Analysis",
    [
        "Sophisticated voting strategies",
        "Cohen's kappa agreement (κ = 0.436)",
        "Method concordance scoring",
        "Consensus strength classification"
    ],
    colors['consensus'], 4
)

# Input Data Box
input_box = FancyBboxPatch(
    (3, 17.7), 6, 0.8,
    boxstyle="round,pad=0.1",
    facecolor=colors['input'],
    edgecolor='black',
    linewidth=2
)
ax.add_patch(input_box)
ax.text(6, 18.1, '📋 Input: Count Table + Metadata', 
        ha='center', va='center', fontweight='bold', fontsize=12)

# Results Box  
results_box = FancyBboxPatch(
    (2.5, 8.5), 7, 2.5,
    boxstyle="round,pad=0.1",
    facecolor=colors['output'],
    edgecolor='black',
    linewidth=2
)
ax.add_patch(results_box)

ax.text(6, 10.5, '📈 Core Analysis Results', 
        ha='center', va='center', fontweight='bold', fontsize=14)

results_details = [
    "• Individual method results: p-values, effect sizes, significant features",
    "• Method-specific performance metrics & diagnostic information", 
    "• Consensus calls with confidence scores (if consensus enabled)",
    "• Statistical summary & method recommendations for final use"
]

y_pos = 10.1
for detail in results_details:
    ax.text(2.7, y_pos, detail, ha='left', va='center', fontsize=10)
    y_pos -= 0.3

# Add connecting arrows with labels
# Input to processing
create_arrow(5, 17.7, 2.8, 18, label='Count Matrix +\nMetadata', label_pos='mid')  # Input to Step 1
create_arrow(7, 17.7, 9.2, 18, label='Data\nCharacteristics', label_pos='mid')  # Input to Step 2

# Between steps
create_arrow(5, 16.7, 7, 16.7, label='Data Profile +\nRecommendations', label_pos='mid')  # 1 to 2
create_arrow(2.75, 15.5, 2.75, 14.5, label='Method\nSelection', label_pos='mid')    # 1 to 3
create_arrow(9.25, 15.5, 9.25, 14.5, label='Selected Methods +\nConfidence', label_pos='mid')    # 2 to 4
create_arrow(5, 13.2, 7, 13.2, label='Individual Method\nResults', label_pos='mid')  # 3 to 4

# To results
create_arrow(2.75, 12, 4.5, 11, label='Method-Specific\nResults', label_pos='mid')    # 3 to results
create_arrow(9.25, 12, 7.5, 11, label='Consensus Results +\nConfidence Scores', label_pos='mid')    # 4 to results

# Add method icons
method_y = 6.5
methods = [
    "🔬 Wilcoxon", "🧪 ALDEx2", "📊 DESeq2", 
    "📈 edgeR", "🔍 metagenomeSeq", "⚖️ ANCOM-BC"
]

ax.text(6, 7.5, 'Integrated Statistical Methods', 
        ha='center', va='center', fontweight='bold', fontsize=12)

for i, method in enumerate(methods):
    x_pos = 2 + (i % 3) * 2.5
    y_pos = method_y - (i // 3) * 0.5
    
    method_box = FancyBboxPatch(
        (x_pos - 0.6, y_pos - 0.15), 1.2, 0.3,
        boxstyle="round,pad=0.05",
        facecolor='white',
        edgecolor='gray',
        linewidth=1
    )
    ax.add_patch(method_box)
    ax.text(x_pos, y_pos, method, ha='center', va='center', fontsize=9)

# Add R Integration note
ax.text(6, 5, '🔗 Complete R Integration via rpy2', 
        ha='center', va='center', fontsize=11, style='italic', color='blue')

# Add footer
ax.text(6, 0.5, 'DAAadvisor Core Framework: Intelligent Method Selection & Analysis', 
        ha='center', va='center', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('daaadvisor_core_workflow.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Core workflow diagram created: daaadvisor_core_workflow.png")