#!/usr/bin/env python3
"""
Create Enhanced DAAadvisor Validation Workflow Diagram
Sophisticated layout with clear explanations of validation characteristics
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, ConnectionPatch, Rectangle
import numpy as np

# Set up the figure with high resolution
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 9
plt.rcParams['font.family'] = 'sans-serif'

fig, ax = plt.subplots(1, 1, figsize=(18, 20))
ax.set_xlim(0, 14)
ax.set_ylim(0, 22)
ax.axis('off')

# Modern color scheme (matching core workflow)
colors = {
    'header': '#2E3440',
    'real_data': '#88C0D0',
    'synthetic': '#D08770', 
    'cross_val': '#A3BE8C',
    'validation': '#5E81AC',
    'reporting': '#B48EAD',
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

ax.text(7, 20.75, 'DAAadvisor: Scientific Validation Framework', 
        ha='center', va='center', fontsize=16, fontweight='bold', color='white')
ax.text(7, 20.25, 'Ensuring Your Results are Publication-Ready and Reproducible', 
        ha='center', va='center', fontsize=11, style='italic', color='white')

# Input from core analysis
input_box = FancyBboxPatch(
    (5.5, 18), 3, 1.2,
    boxstyle="round,pad=0.1",
    facecolor=colors['light_bg'],
    edgecolor=colors['text'],
    linewidth=2,
    zorder=2
)
ax.add_patch(input_box)
ax.text(7, 18.6, '⬆️ From Core Analysis', 
        ha='center', va='center', fontweight='bold', fontsize=11, color=colors['text'])
ax.text(7, 18.3, 'Statistical results + method performance', 
        ha='center', va='center', fontsize=9, style='italic', color=colors['text'])

# Real Data Integration
create_modern_box(
    0.5, 15, 5.5, 2.8,
    "Real-World Data Integration",
    "Downloads and processes actual microbiome studies for validation",
    [
        "Live Data Source: curatedMetagenomicData (Bioconductor)",
        "Proven Datasets: IBD, Cancer, Diabetes studies (1,000+ samples)",
        "Ground Truth: Known biomarkers from published literature",
        "Quality Control: Standardized formatting and sample validation"
    ],
    colors['real_data'], "🧬"
)

# Synthetic Data Generation
create_modern_box(
    8, 15, 5.5, 2.8,
    "Controlled Data Simulation",
    "Creates realistic synthetic datasets with known ground truth",
    [
        "Literature-Based: Effect sizes from real studies (IBD, CRC)",
        "Realistic Patterns: Mimics true microbiome characteristics",
        "Known Truth: We know exactly which features should be significant",
        "Controlled Testing: Perfect for method performance evaluation"
    ],
    colors['synthetic'], "🎭"
)

# Cross-Validation Engine
create_modern_box(
    3, 11.5, 8, 2.8,
    "Cross-Validation Engine",
    "Compares method performance on real vs synthetic data",
    [
        "Performance Correlation: Do methods agree on real and synthetic data?",
        "Bootstrap Testing: 50-100 iterations for statistical confidence",
        "Ground Truth Recovery: How well do methods find known answers?",
        "Robustness Assessment: Are results consistent across datasets?"
    ],
    colors['cross_val'], "🔄"
)

# Publication Validation
create_modern_box(
    0.5, 8, 5.5, 2.8,
    "Publication-Quality Validation",
    "Ensures results meet scientific publication standards",
    [
        "Statistical Rigor: Bootstrap confidence intervals (95% CI)",
        "Literature Confirmation: Validates against known biomarkers",
        "Comprehensive Metrics: AUROC, AUPRC, F1, sensitivity, specificity",
        "Peer Review Ready: Meets journal requirements for rigor"
    ],
    colors['validation'], "🏆"
)

# Comprehensive Reporting
create_modern_box(
    8, 8, 5.5, 2.8,
    "Comprehensive Reporting",
    "Generates publication-ready reports and visualizations",
    [
        "Interactive Dashboards: HTML reports with drill-down capability",
        "Journal Figures: High-resolution plots ready for submission",
        "Method Comparisons: Side-by-side performance analysis",
        "Statistical Summaries: Complete results with uncertainty quantification"
    ],
    colors['reporting'], "📊"
)

# Final Results
results_box = FancyBboxPatch(
    (2, 4.5), 10, 2.8,
    boxstyle="round,pad=0.2",
    facecolor=colors['output'],
    edgecolor=colors['text'],
    linewidth=3,
    zorder=2
)
ax.add_patch(results_box)

results_shadow = FancyBboxPatch(
    (2.1, 4.4), 10, 2.8,
    boxstyle="round,pad=0.2",
    facecolor='#D8DEE9',
    alpha=0.3,
    zorder=1
)
ax.add_patch(results_shadow)

ax.text(7, 6.8, '🎯 Publication-Ready Scientific Output', 
        ha='center', va='center', fontweight='bold', fontsize=14, color=colors['text'])

results_content = [
    "📈 Validated Performance: Real-world testing with 1,627+ samples",
    "🔬 Literature Confirmed: Results validated against known biomarkers", 
    "📊 Statistical Confidence: Bootstrap intervals and rigorous testing",
    "📋 Complete Documentation: Ready for peer review and publication"
]

y_pos = 6.3
for item in results_content:
    ax.text(2.5, y_pos, item, ha='left', va='center', fontsize=10, 
            color=colors['text'], fontweight='500')
    y_pos -= 0.35

# Add sophisticated arrows
create_modern_arrow(6.5, 18, 3.5, 17.8, "Core Results\nfor Testing")
create_modern_arrow(7.5, 18, 10.5, 17.8, "Performance\nBaseline")

create_modern_arrow(3.25, 15, 5.5, 14.3, "Real Data +\nBiomarkers", curve=0.3)
create_modern_arrow(10.75, 15, 8.5, 14.3, "Synthetic Data +\nKnown Truth", curve=-0.3)

create_modern_arrow(5, 11.5, 3.5, 10.8, "Cross-Val\nResults")
create_modern_arrow(9, 11.5, 10.5, 10.8, "Performance\nMetrics")

create_modern_arrow(3.25, 8, 5, 7.3, "Scientific\nValidation")
create_modern_arrow(10.75, 8, 9, 7.3, "Reports +\nVisualizations")

# Validation highlights
highlight_box = FancyBboxPatch(
    (1, 1.5), 12, 2.5,
    boxstyle="round,pad=0.2",
    facecolor=colors['light_bg'],
    edgecolor=colors['text'],
    linewidth=2,
    zorder=2
)
ax.add_patch(highlight_box)

ax.text(7, 3.7, '💡 Why This Validation Matters', 
        ha='center', va='center', fontweight='bold', fontsize=13, color=colors['text'])

validation_points = [
    "🎯 Confidence in Results: Know your findings are robust, not just lucky",
    "🔬 Publication Ready: Meets peer review standards with proper statistical testing",
    "📊 Real-World Validated: Tested on actual microbiome studies, not just simulations",
    "🏆 Scientific Rigor: Bootstrap confidence intervals + literature confirmation"
]

y_pos = 3.3
for point in validation_points:
    ax.text(1.5, y_pos, point, ha='left', va='center', fontsize=10, 
            color=colors['text'], fontweight='500')
    y_pos -= 0.35

# Footer
ax.text(7, 0.5, 'DAAadvisor: Bringing Scientific Rigor to Microbiome Analysis', 
        ha='center', va='center', fontsize=12, fontweight='bold', color=colors['header'])

plt.tight_layout()
plt.savefig('daaadvisor_enhanced_validation_workflow.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Enhanced validation workflow diagram created: daaadvisor_enhanced_validation_workflow.png")
print("🎨 Modern, sophisticated layout with clear explanations")
print("💡 Intuitive descriptions of validation characteristics")