#!/usr/bin/env python3
"""
Create Enhanced DAAadvisor Core Workflow Diagram
Sophisticated layout with clear explanations of core characteristics
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch, Rectangle
import numpy as np

# Set up the figure with high resolution
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'sans-serif'

fig, ax = plt.subplots(1, 1, figsize=(18, 20))
ax.set_xlim(0, 14)
ax.set_ylim(0, 22)
ax.axis('off')

# Modern color scheme
colors = {
    'header': '#2E3440',
    'input': '#88C0D0',
    'profiling': '#D08770', 
    'selection': '#A3BE8C',
    'analysis': '#5E81AC',
    'consensus': '#B48EAD',
    'output': '#8FBCBB',
    'accent': '#BF616A',
    'text': '#2E3440',
    'light_bg': '#ECEFF4'
}

def create_modern_box(x, y, width, height, title, description, details, color, icon=''):
    """Create a modern, sophisticated step box"""
    
    # Main box with shadow effect
    shadow = FancyBboxPatch(
        (x + 0.1, y - 0.1), width, height,
        boxstyle="round,pad=0.15",
        facecolor='#D8DEE9',
        alpha=0.3,
        zorder=1
    )
    ax.add_patch(shadow)
    
    # Main box
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.15",
        facecolor=color,
        edgecolor=colors['text'],
        linewidth=2,
        zorder=2
    )
    ax.add_patch(box)
    
    # Header section
    header_box = FancyBboxPatch(
        (x + 0.1, y + height - 0.8), width - 0.2, 0.6,
        boxstyle="round,pad=0.05",
        facecolor=colors['header'],
        alpha=0.1,
        zorder=3
    )
    ax.add_patch(header_box)
    
    # Title with icon
    ax.text(x + 0.3, y + height - 0.5, f"{icon} {title}", 
            ha='left', va='center', fontweight='bold', fontsize=12, 
            color=colors['text'])
    
    # Description
    ax.text(x + 0.3, y + height - 1.1, description, 
            ha='left', va='top', fontsize=10, style='italic',
            color=colors['text'], wrap=True)
    
    # Details
    detail_y = y + height - 1.5
    for detail in details:
        # Add bullet point
        ax.text(x + 0.3, detail_y, "●", ha='left', va='center', 
                fontsize=8, color=colors['accent'])
        ax.text(x + 0.5, detail_y, detail, ha='left', va='center', 
                fontsize=9, color=colors['text'], wrap=True)
        detail_y -= 0.35

def create_modern_arrow(start_x, start_y, end_x, end_y, label, curve=0):
    """Create a modern arrow with label"""
    if curve != 0:
        # Create curved arrow
        from matplotlib.patches import FancyArrowPatch
        arrow = FancyArrowPatch(
            (start_x, start_y), (end_x, end_y),
            arrowstyle='->', mutation_scale=20,
            color=colors['text'], linewidth=2.5,
            connectionstyle=f"arc3,rad={curve}",
            zorder=4
        )
    else:
        # Straight arrow
        arrow = ConnectionPatch(
            (start_x, start_y), (end_x, end_y), 
            "data", "data",
            arrowstyle='->', shrinkA=5, shrinkB=5,
            mutation_scale=20, fc=colors['text'], linewidth=2.5,
            color=colors['text'], zorder=4
        )
    ax.add_patch(arrow)
    
    # Add label with background
    if curve != 0:
        label_x = (start_x + end_x) / 2 + curve
        label_y = (start_y + end_y) / 2 + abs(curve) * 0.5
    else:
        label_x = (start_x + end_x) / 2
        label_y = (start_y + end_y) / 2
    
    ax.text(label_x, label_y, label, ha='center', va='center', 
            fontsize=8, fontweight='bold',
            bbox=dict(boxstyle="round,pad=0.4", facecolor='white', 
                     edgecolor=colors['text'], alpha=0.9),
            zorder=5)

# Title section
title_box = FancyBboxPatch(
    (1, 20), 12, 1.5,
    boxstyle="round,pad=0.2",
    facecolor=colors['header'],
    edgecolor='none',
    zorder=2
)
ax.add_patch(title_box)

ax.text(7, 20.75, 'DAAadvisor: Intelligent Core Analysis Framework', 
        ha='center', va='center', fontsize=16, fontweight='bold', color='white')
ax.text(7, 20.25, 'From Raw Data to Publication-Ready Results', 
        ha='center', va='center', fontsize=11, style='italic', color='white')

# Input section
create_modern_box(
    5.5, 18, 3, 1.5,
    "Data Input",
    "Your microbiome dataset",
    [
        "Count table (samples × features)",
        "Sample metadata with conditions"
    ],
    colors['input'], "📊"
)

# Step 1: Data Profiling
create_modern_box(
    0.5, 15.5, 5.5, 2.8,
    "Smart Data Profiling",
    "Automatically characterizes your dataset to understand its structure",
    [
        "Sparsity Assessment: How many zeros? (typical: 70-90%)",
        "Data Type Detection: ASV/16S, genes, or viral sequences?", 
        "Compositional Analysis: Library size variation patterns",
        "Quality Metrics: Distribution shapes and outlier detection"
    ],
    colors['profiling'], "🔍"
)

# Step 2: Intelligent Method Selection  
create_modern_box(
    8, 15.5, 5.5, 2.8,
    "AI-Powered Method Selection",
    "Uses information theory to choose the best statistical approach",
    [
        "Maximum Entropy: Finds optimal method given data constraints",
        "Divergence Analysis: Measures group separation potential",
        "Method Confidence: Quantifies selection certainty (0-100%)",
        "Recommendation Engine: Suggests primary + backup methods"
    ],
    colors['selection'], "🧠"
)

# Step 3: Multi-Method Analysis
create_modern_box(
    0.5, 12, 5.5, 2.8,
    "Parallel Statistical Testing",
    "Runs multiple proven methods simultaneously for robust results",
    [
        "Non-parametric: Wilcoxon (works on any distribution)",
        "Compositional: ALDEx2 (handles microbiome ratios correctly)",
        "Count-based: DESeq2, edgeR (for abundance data)",
        "Zero-inflation: metagenomeSeq (handles sparse data)"
    ],
    colors['analysis'], "⚡"
)

# Step 4: Smart Consensus
create_modern_box(
    8, 12, 5.5, 2.8,
    "Intelligent Result Integration",
    "Combines multiple methods using advanced agreement algorithms",
    [
        "Voting Strategies: Simple majority, weighted by reliability",
        "Agreement Metrics: Cohen's kappa inter-method concordance",
        "Confidence Scoring: How sure are we? (High/Medium/Low)",
        "Conflict Resolution: Handles disagreements systematically"
    ],
    colors['consensus'], "🤝"
)

# Results section - removed as requested

# Add sophisticated arrows with labels
create_modern_arrow(6.5, 18, 3.5, 18.3, "Raw\nMicrobiome\nData")
create_modern_arrow(7.5, 18, 10.5, 18.3, "Dataset\nCharacteristics")

create_modern_arrow(3.25, 15.5, 3.25, 14.8, "Data\nProfile")
create_modern_arrow(10.75, 15.5, 10.75, 14.8, "Method\nRecommendation")

create_modern_arrow(6, 16.9, 8, 16.9, "Profiling\nResults", curve=0.3)
create_modern_arrow(6, 13.4, 8, 13.4, "Statistical\nResults", curve=0.3)

# Key characteristics panel - moved up to fill space
char_box = FancyBboxPatch(
    (1, 8), 12, 2.5,
    boxstyle="round,pad=0.2",
    facecolor=colors['light_bg'],
    edgecolor=colors['text'],
    linewidth=2,
    zorder=2
)
ax.add_patch(char_box)

ax.text(7, 10.2, '💡 Why This Framework is Powerful', 
        ha='center', va='center', fontweight='bold', fontsize=13, color=colors['text'])

characteristics = [
    "🎯 Automatic Intelligence: No guessing which method to use - AI chooses based on your data",
    "🔬 Scientific Rigor: Multiple methods + consensus = more reliable than any single test",
    "⚡ Handles All Data Types: Works with ASV, gene, viral, sparse, or dense microbiome data",
    "📊 Publication Ready: Bootstrap confidence intervals + comprehensive reporting included"
]

y_pos = 9.8
for char in characteristics:
    ax.text(1.5, y_pos, char, ha='left', va='center', fontsize=11, 
            color=colors['text'], fontweight='500')
    y_pos -= 0.35

# Footer with integration info - moved up to fill space
footer_box = FancyBboxPatch(
    (2, 6), 10, 1.5,
    boxstyle="round,pad=0.15",
    facecolor='white',
    edgecolor=colors['accent'],
    linewidth=2,
    zorder=2
)
ax.add_patch(footer_box)

ax.text(7, 7.1, '🔗 Complete R Integration', 
        ha='center', va='center', fontweight='bold', fontsize=12, color=colors['accent'])
ax.text(7, 6.6, 'Seamlessly integrates 6 statistical methods via rpy2', 
        ha='center', va='center', fontsize=11, style='italic', color=colors['text'])
ax.text(7, 6.3, 'Wilcoxon • ALDEx2 • DESeq2 • edgeR • metagenomeSeq • ANCOM-BC', 
        ha='center', va='center', fontsize=10, color=colors['text'])

# Usage example
ax.text(7, 5.3, '🚀 Simple Usage', 
        ha='center', va='center', fontweight='bold', fontsize=12, color=colors['text'])
ax.text(7, 4.9, 'tool = DifferentialAbundanceTool()', 
        ha='center', va='center', fontsize=10, family='monospace', 
        color=colors['text'])
ax.text(7, 4.6, 'results = tool.analyze(count_table, metadata, use_consensus=True)', 
        ha='center', va='center', fontsize=10, family='monospace', 
        color=colors['text'])

# Final footer
ax.text(7, 3.7, 'DAAadvisor: Making Microbiome Analysis Intelligent and Reproducible', 
        ha='center', va='center', fontsize=12, fontweight='bold', color=colors['header'])

plt.tight_layout()
plt.savefig('daaadvisor_enhanced_core_workflow.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Enhanced core workflow diagram created: daaadvisor_enhanced_core_workflow.png")
print("🎨 Modern, sophisticated layout with clear explanations")
print("💡 Intuitive descriptions of core characteristics")