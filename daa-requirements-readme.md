# Differential Abundance Analysis Tool

An intelligent tool for microbiome differential abundance analysis that automatically selects the best statistical method based on your data characteristics.

## Features

- **Automatic Method Selection**: Analyzes your data characteristics and recommends the most appropriate differential abundance method
- **Multi-Method Integration**: Supports ALDEx2, ANCOM-BC, DESeq2, edgeR, metagenomeSeq, ZicoSeq, and more
- **Data Type Awareness**: Optimized for different data types (16S/ASV, metagenomics/genes, viral)
- **Consensus Analysis**: Runs multiple methods and provides consensus results
- **Comprehensive Visualization**: Data profiling, method comparisons, volcano plots, and interactive dashboards
- **Detailed Reporting**: Explanations for method choices and comprehensive result summaries

## Installation

### Requirements

Create a conda environment with the necessary dependencies:

```bash
# Create conda environment
conda create -n daa-tool python=3.9
conda activate daa-tool

# Install Python packages
pip install -r requirements.txt

# Install R and R packages
conda install -c conda-forge r-base=4.3
```

### requirements.txt

```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=0.24.0
statsmodels>=0.12.0
matplotlib>=3.4.0
seaborn>=0.11.0
plotly>=5.0.0
biom-format>=2.1.10
rpy2>=3.5.0
click>=8.0.0
tqdm>=4.62.0
pytest>=6.2.0
jupyter>=1.0.0
```

### R Package Installation

```R
# In R or RStudio
install.packages("BiocManager")
BiocManager::install(c(
    "ALDEx2",
    "ANCOMBC", 
    "DESeq2",
    "edgeR",
    "metagenomeSeq"
))

# Additional packages
install.packages(c("tidyverse", "vegan"))
```

## Quick Start

```python
from daa_tool import DifferentialAbundanceTool
import pandas as pd

# Load your data
count_table = pd.read_csv("counts.csv", index_col=0)  # Samples x Features
metadata = pd.read_csv("metadata.csv", index_col=0)   # Sample metadata

# Run analysis
tool = DifferentialAbundanceTool()
results = tool.analyze(
    count_table=count_table,
    metadata=metadata,
    data_type='asv',  # or 'gene' or 'viral'
    use_consensus=True
)

# Summarize results
tool.summarize_results()

# Visualize
from daa_tool.visualization import DAAVisualizer
viz = DAAVisualizer()
viz.plot_data_characteristics(results['profile'])
viz.plot_method_comparison(results)
```

## Data Format

### Count Table (CSV)
```
        ASV_1   ASV_2   ASV_3   ...
Sample1   10      0      45    ...
Sample2    0     23       0    ...
Sample3   15      5      12    ...
```

### Metadata (CSV)
```
        condition   batch   age
Sample1   Control     A     25
Sample2 Treatment     A     30
Sample3   Control     B     28
```

## Method Selection Logic

The tool automatically selects methods based on:

1. **Data Type**
   - ASV/16S: Prefers compositionally-aware methods (ALDEx2, ANCOM-BC)
   - Genes: Can use RNA-seq methods (DESeq2, edgeR)
   - Viral: Handles extreme sparsity (ZicoSeq, ANCOM-BC)

2. **Sample Size**
   - Small (<20): Methods robust to small samples
   - Large (>50): Methods with better power

3. **Data Characteristics**
   - High sparsity: Zero-inflated methods
   - Compositional bias: Methods with bias correction
   - Uneven sequencing depth: Proper normalization

## Understanding Results

### Profile Summary
- **Sparsity**: Proportion of zeros in your data
- **Zero Inflation**: Features with excessive zeros
- **Compositional Bias**: Correlation between abundant features

### Method Recommendations
- **Primary Method**: Best suited for your data
- **Confidence**: How certain the recommendation is
- **Reasoning**: Why this method was chosen

### Analysis Results
- **Individual Methods**: Results from each method
- **Consensus**: Features significant across multiple methods
- **Visualizations**: Volcano plots, method comparisons

## Advanced Usage

### Custom Method Parameters

```python
# Run with custom parameters
results = tool.analyze(
    count_table=count_table,
    metadata=metadata,
    method_params={
        'deseq2': {'test': 'Wald', 'fitType': 'local'},
        'aldex2': {'mc.samples': 128, 'denom': 'all'}
    }
)
```

### Batch Processing

```python
# Process multiple datasets
datasets = ['study1.csv', 'study2.csv', 'study3.csv']
all_results = {}

for dataset in datasets:
    counts = pd.read_csv(dataset, index_col=0)
    results = tool.analyze(counts, metadata)
    all_results[dataset] = results
```

### Export Results

```python
# Save results
results['analyses']['primary'].to_csv('primary_results.csv')
results['consensus'].to_csv('consensus_results.csv')

# Save visualizations
viz.plot_volcano(results['analyses']['primary'], save_path='volcano.png')
viz.create_interactive_dashboard(results).write_html('dashboard.html')
```

## Troubleshooting

### Common Issues

1. **R packages not found**
   - Ensure R is in your PATH
   - Check R package installation with `library(package_name)` in R

2. **Memory errors with large datasets**
   - Filter low-prevalence features
   - Use chunked processing
   - Increase available memory

3. **Method-specific errors**
   - Check input data format
   - Ensure no missing values in metadata
   - Verify group sizes meet method requirements

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
# Clone repository
git clone https://github.com/yourusername/daa-tool.git
cd daa-tool

# Install in development mode
pip install -e .

# Run tests
pytest tests/
```

## Citation

If you use this tool in your research, please cite:

```
@software{daa_tool,
  title = {Differential Abundance Analysis Tool},
  author = {Your Name},
  year = {2024},
  url = {https://github.com/yourusername/daa-tool}
}
```

Also cite the individual methods used:
- ALDEx2: Fernandes et al. (2014)
- ANCOM-BC: Lin & Peddada (2020)
- DESeq2: Love et al. (2014)
- etc.

## License

MIT License - see LICENSE file for details

## Support

- Documentation: [Read the Docs](https://daa-tool.readthedocs.io)
- Issues: [GitHub Issues](https://github.com/yourusername/daa-tool/issues)
- Discussions: [GitHub Discussions](https://github.com/yourusername/daa-tool/discussions)

## Acknowledgments

This tool integrates methods developed by many research groups. We thank all the original authors for their contributions to the field.