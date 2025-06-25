#!/usr/bin/env python3
"""
Create DAAadvisor Validation & Testing Workflow Diagram
Cross-validation, real data integration, and reporting pipeline
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

fig, ax = plt.subplots(1, 1, figsize=(14, 16))
ax.set_xlim(0, 10)
ax.set_ylim(0, 18)
ax.axis('off')

# Color scheme
colors = {
    'real_data': '#E0F2F1',   # Light teal
    'synthetic': '#E8EAF6',   # Light indigo
    'cross_val': '#FFF8E1',   # Light yellow
    'validation': '#FFEBEE',  # Light red
    'reporting': '#F1F8E9',   # Light lime
    'results': '#FCE4EC'      # Light pink
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

def create_arrow(start_x, start_y, end_x, end_y, style='->', color='black'):
    """Create connection arrow between steps"""
    arrow = ConnectionPatch(
        (start_x, start_y), (end_x, end_y), 
        "data", "data",
        arrowstyle=style, shrinkA=5, shrinkB=5,
        mutation_scale=20, fc=color, linewidth=2, color=color
    )
    ax.add_patch(arrow)

# Title
ax.text(5, 17, 'DAAadvisor: Validation & Testing Framework', 
        ha='center', va='center', fontsize=18, fontweight='bold')
ax.text(5, 16.5, 'Cross-Validation, Real Data Integration & Publication-Quality Reporting', 
        ha='center', va='center', fontsize=12, style='italic')

# Input from Core Analysis
input_box = FancyBboxPatch(
    (2, 15.2), 6, 0.8,
    boxstyle="round,pad=0.1",
    facecolor='#F5F5F5',
    edgecolor='black',
    linewidth=2
)
ax.add_patch(input_box)
ax.text(5, 15.6, '⬆️ From Core Analysis: Method Results & Consensus', 
        ha='center', va='center', fontweight='bold', fontsize=12)

# Step 5: Real Data Integration
create_step_box(
    0.5, 12.5, 4, 2.2,
    "🧬 Real Data Integration",
    [
        "curatedMetagenomicData download",
        "IBD: 1,627 samples (1,201 IBD + 426 controls)",
        "Literature ground truth biomarkers",
        "Data standardization & quality control"
    ],
    colors['real_data'], 5
)

# Synthetic Data Generation
create_step_box(
    5.5, 12.5, 4, 2.2,
    "🎭 Synthetic Data Generation",
    [
        "Literature-based realistic simulations",
        "IBD/CRC/Antibiotic effect patterns",
        "Known differential features (ground truth)",
        "Controlled experimental conditions"
    ],
    colors['synthetic'], "S"
)

# Step 6: Cross-Validation Framework
create_step_box(
    2.25, 9.8, 5.5, 2.2,
    "🔄 Cross-Validation Framework",
    [
        "Real vs synthetic performance comparison",
        "Bootstrap validation (50-100 iterations)",
        "Ground truth recovery assessment",
        "Method robustness & consistency evaluation"
    ],
    colors['cross_val'], 6
)

# Step 7: Publication Validation
create_step_box(
    0.5, 7.1, 4, 2.2,
    "🏆 Publication-Quality Validation",
    [
        "Statistical rigor with confidence intervals",
        "Literature biomarker confirmation",
        "Comprehensive metrics (AUROC, F1, AUPRC)",
        "Real-world testing validation"
    ],
    colors['validation'], 7
)

# Step 8: Results & Reporting
create_step_box(
    5.5, 7.1, 4, 2.2,
    "📈 Results & Comprehensive Reporting",
    [
        "Interactive HTML dashboards",
        "Cross-validation reports & plots",
        "Publication-ready figures",
        "Bootstrap statistical summaries"
    ],
    colors['reporting'], 8
)

# Final Output Box
output_box = FancyBboxPatch(
    (1.5, 4.5), 7, 2,
    boxstyle="round,pad=0.15",
    facecolor=colors['results'],
    edgecolor='black',
    linewidth=2
)
ax.add_patch(output_box)

ax.text(5, 6, '🎯 Publication-Ready Results', 
        ha='center', va='center', fontweight='bold', fontsize=14)

output_details = [
    "• Cross-validated performance metrics with confidence intervals",
    "• Real vs synthetic data comparison reports",
    "• Literature-confirmed biomarker validation",
    "• Interactive dashboards & journal-ready visualizations"
]

y_pos = 5.6
for detail in output_details:
    ax.text(1.7, y_pos, detail, ha='left', va='center', fontsize=10)
    y_pos -= 0.25

# Add connecting arrows
# From input to data sources
create_arrow(4, 15.2, 2.5, 14.7)  # Input to real data
create_arrow(6, 15.2, 7.5, 14.7)  # Input to synthetic

# Data sources to cross-validation
create_arrow(2.5, 12.5, 4, 12)    # Real data to cross-val
create_arrow(7.5, 12.5, 6, 12)    # Synthetic to cross-val

# Cross-validation to outputs
create_arrow(3.5, 9.8, 2.5, 9.3)  # Cross-val to validation
create_arrow(6.5, 9.8, 7.5, 9.3)  # Cross-val to reporting

# To final results
create_arrow(2.5, 7.1, 4, 6.5)    # Validation to results
create_arrow(7.5, 7.1, 6, 6.5)    # Reporting to results

# Add validation datasets
datasets_y = 2.5
ax.text(5, 3.5, 'Validated Real-World Datasets', 
        ha='center', va='center', fontweight='bold', fontsize=12)

datasets = [
    "🦠 IBD (1,627 samples)", "🔬 Colorectal Cancer", "💊 Type 2 Diabetes",
    "⚖️ Obesity Studies", "🫀 Liver Cirrhosis", "💉 Antibiotic Treatment"
]

for i, dataset in enumerate(datasets):
    x_pos = 1 + (i % 3) * 2.5
    y_pos = datasets_y - (i // 3) * 0.4
    
    dataset_box = FancyBboxPatch(
        (x_pos - 0.7, y_pos - 0.15), 1.4, 0.3,
        boxstyle="round,pad=0.05",
        facecolor='white',
        edgecolor='teal',
        linewidth=1
    )
    ax.add_patch(dataset_box)
    ax.text(x_pos, y_pos, dataset, ha='center', va='center', fontsize=9)

# Add performance highlights
ax.text(5, 1.2, '🏆 Key Achievements: 100% Success Rate • 6/6 Methods Functional • Real Data Validated', 
        ha='center', va='center', fontsize=11, fontweight='bold', color='green')

# Add footer
ax.text(5, 0.5, 'DAAadvisor Validation Framework: Ensuring Scientific Rigor & Reproducibility', 
        ha='center', va='center', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('daaadvisor_validation_workflow.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.close()

print("✅ Validation workflow diagram created: daaadvisor_validation_workflow.png")