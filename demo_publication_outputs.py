#!/usr/bin/env python3
"""
Demo script to show publication output structure without full benchmark
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

def create_demo_outputs():
    """Create demo publication outputs to show structure"""
    
    output_dir = Path("demo_publication_outputs")
    output_dir.mkdir(exist_ok=True)
    
    print("📁 Creating demo publication outputs...")
    
    # 1. Main Summary Table
    print("📊 Creating publication_summary_table.csv...")
    summary_data = []
    
    datasets = ["IBD_pediatric", "CRC_fecal", "antibiotic_longitudinal", "controlled_es2.0_n100"]
    methods = ["wilcoxon", "deseq2", "edger", "aldex2"]
    
    for dataset in datasets:
        for method in methods:
            # Simulate realistic performance metrics
            f1 = np.random.uniform(0.1, 0.8)
            sensitivity = np.random.uniform(0.2, 0.9)
            specificity = np.random.uniform(0.7, 0.95)
            
            summary_data.append({
                'Dataset': dataset,
                'Study_Type': 'Disease State' if 'IBD' in dataset or 'CRC' in dataset else 
                            'Antibiotic' if 'antibiotic' in dataset else 'Controlled',
                'Data_Type': '16S V4' if 'IBD' in dataset else 'Shotgun' if 'CRC' in dataset else '16S V3-V4',
                'Sample_Size': np.random.choice([50, 80, 100, 150]),
                'Method': method,
                'F1_Score': f"{f1:.3f} ± {f1*0.1:.3f}",
                'F1_CI': f"[{f1-0.05:.3f}, {f1+0.05:.3f}]",
                'Sensitivity': f"{sensitivity:.3f} ± {sensitivity*0.1:.3f}",
                'Specificity': f"{specificity:.3f} ± {specificity*0.1:.3f}",
                'Precision': f"{f1*0.9:.3f} ± {f1*0.09:.3f}",
                'AUROC': f"{np.random.uniform(0.6, 0.9):.3f} ± {0.05:.3f}",
                'AUPRC': f"{np.random.uniform(0.4, 0.8):.3f} ± {0.05:.3f}",
                'FDR': f"{np.random.uniform(0.1, 0.4):.3f} ± {0.03:.3f}",
                'Bootstrap_N': np.random.choice([95, 98, 100])
            })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / "publication_summary_table.csv", index=False)
    
    # 2. Publication-ready LaTeX table
    print("📋 Creating Table1_publication_table.csv...")
    pub_table_data = []
    for method in methods:
        method_data = summary_df[summary_df['Method'] == method]
        
        pub_table_data.append({
            'Method': method.upper(),
            'N_Datasets': len(method_data),
            'F1_Score': f"{method_data['F1_Score'].apply(lambda x: float(x.split(' ± ')[0])).mean():.3f} ± {0.05:.3f}",
            'Sensitivity': f"{method_data['Sensitivity'].apply(lambda x: float(x.split(' ± ')[0])).mean():.3f} ± {0.05:.3f}",
            'Specificity': f"{method_data['Specificity'].apply(lambda x: float(x.split(' ± ')[0])).mean():.3f} ± {0.05:.3f}",
            'AUROC': f"{method_data['AUROC'].apply(lambda x: float(x.split(' ± ')[0])).mean():.3f} ± {0.05:.3f}",
            'Bootstrap_Success': f"{method_data['Bootstrap_N'].mean():.0f}/100"
        })
    
    pub_table = pd.DataFrame(pub_table_data)
    pub_table.to_csv(output_dir / "Table1_publication_table.csv", index=False)
    
    # Create LaTeX version
    latex_table = pub_table.to_latex(
        index=False, 
        float_format="%.3f",
        caption="DAAadvisor Performance Across Real-World Microbiome Datasets",
        label="tab:daaadvisor_performance"
    )
    with open(output_dir / "Table1_publication_table.tex", 'w') as f:
        f.write(latex_table)
    
    # 3. Dataset metadata
    print("📋 Creating dataset_metadata.json...")
    metadata = {
        "IBD_pediatric": {
            "study": "Pediatric IBD microbiome",
            "disease": "Inflammatory Bowel Disease", 
            "tissue": "Fecal",
            "sequencing": "16S V4",
            "sample_size": 100,
            "expected_differential": 50
        },
        "CRC_fecal": {
            "study": "Colorectal cancer microbiome",
            "disease": "Colorectal Cancer",
            "tissue": "Fecal", 
            "sequencing": "16S V3-V4",
            "sample_size": 150,
            "expected_differential": 75
        },
        "antibiotic_longitudinal": {
            "study": "Antibiotic perturbation longitudinal",
            "treatment": "Broad-spectrum antibiotics",
            "tissue": "Fecal",
            "sequencing": "16S V4", 
            "sample_size": 120,
            "expected_differential": 200
        }
    }
    
    with open(output_dir / "dataset_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    # 4. Complete benchmark results (simplified)
    print("📋 Creating complete_benchmark_results.json...")
    complete_results = {
        "IBD_pediatric": {
            "wilcoxon": {
                "f1_score_mean": 0.456,
                "f1_score_std": 0.045,
                "f1_score_ci_lower": 0.411,
                "f1_score_ci_upper": 0.501,
                "sensitivity_mean": 0.678,
                "specificity_mean": 0.823,
                "n_bootstrap_success": 98
            },
            "deseq2": {
                "f1_score_mean": 0.534,
                "f1_score_std": 0.053,
                "f1_score_ci_lower": 0.481,
                "f1_score_ci_upper": 0.587,
                "sensitivity_mean": 0.734,
                "specificity_mean": 0.856,
                "n_bootstrap_success": 95
            }
        }
    }
    
    with open(output_dir / "complete_benchmark_results.json", 'w') as f:
        json.dump(complete_results, f, indent=2)
    
    # 5. Create placeholder figures directory
    figures_dir = output_dir / "publication_figures"
    figures_dir.mkdir(exist_ok=True)
    
    # Create figure placeholder files
    figure_names = [
        "Figure1_main_performance.png",
        "Figure2_dataset_comparison.png", 
        "Figure3_statistical_significance.png",
        "interactive_dashboard.html"
    ]
    
    for fig_name in figure_names:
        placeholder_file = figures_dir / fig_name
        with open(placeholder_file, 'w') as f:
            if fig_name.endswith('.html'):
                f.write(f"<html><body><h1>Interactive Dashboard</h1><p>Publication-quality interactive benchmark results would be here.</p></body></html>")
            else:
                f.write(f"# Placeholder for {fig_name}\n# This would contain a high-resolution publication figure")
    
    # 6. Publication report
    print("📋 Creating PUBLICATION_REPORT.md...")
    report_content = f"""# DAAadvisor Publication Benchmark Report

## Executive Summary

- **Datasets analyzed**: {len(datasets)}
- **Methods compared**: {len(methods)}  
- **Total combinations**: {len(summary_df)}
- **Statistical rigor**: Bootstrap confidence intervals
- **Data types**: Disease states, antibiotic studies, controlled experiments

## Key Findings

### Top Performing Methods (by F1 Score)

1. **DESEQ2**: 0.534
2. **WILCOXON**: 0.456  
3. **EDGER**: 0.445
4. **ALDEX2**: 0.423

### Performance Range

- **Best F1 Score**: 0.534 (DESEQ2)
- **Worst F1 Score**: 0.423 (ALDEX2)
- **Mean F1 Score**: 0.465
- **Standard Deviation**: 0.047

## Dataset Categories

### Disease State Studies
- Inflammatory Bowel Disease (IBD)
- Colorectal Cancer (CRC)

### Antibiotic Perturbation Studies  
- Longitudinal antibiotic treatment
- Recovery dynamics
- Microbiome disruption patterns

### Controlled Experiments
- Multiple effect sizes (1.5-4.0)
- Various sample sizes (50-200)
- Known ground truth

## Files Generated

### Data Files
- `publication_summary_table.csv`: Main results table
- `complete_benchmark_results.json`: Detailed results
- `dataset_metadata.json`: Dataset information

### Publication Figures
- `Figure1_main_performance.png`: Overall method comparison
- `Figure2_dataset_comparison.png`: Study type analysis
- `Figure3_statistical_significance.png`: Confidence intervals
- `interactive_dashboard.html`: Supplementary materials
- `Table1_publication_table.csv`: Publication-ready table

---
*Report generated on {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""
    
    with open(output_dir / "PUBLICATION_REPORT.md", 'w') as f:
        f.write(report_content)
    
    print(f"\n✅ Demo publication outputs created in: {output_dir}")
    print(f"📊 Files generated:")
    
    for file_path in sorted(output_dir.rglob("*")):
        if file_path.is_file():
            print(f"   📄 {file_path.relative_to(output_dir)}")
    
    return output_dir

if __name__ == "__main__":
    print("🏆 DAAadvisor Publication Output Demo")
    print("=" * 50)
    
    output_dir = create_demo_outputs()
    
    print(f"\n💡 This shows the structure of outputs from:")
    print(f"   python run_publication_benchmark.py --full")
    print(f"")
    print(f"🎯 For actual results, run the full benchmark!")