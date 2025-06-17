<p align="center">
  <img src="daaadvisor_logo.png" alt="DAAadvisor Logo" width="400"/>
</p>

# DAAadvisor

**Differential Abundance Analysis Advisor for Microbiome Data**

An intelligent tool that automatically selects the best statistical method for your microbiome differential abundance analysis based on data characteristics.

## Features

- **🧠 Intelligent Method Selection**: Automatically recommends the best statistical method based on your data characteristics
- **📊 Comprehensive Data Profiling**: Analyzes sparsity, zero-inflation, compositional bias, and other key metrics
- **🔄 Multi-Method Support**: Integrates ALDEx2, ANCOM-BC, DESeq2, edgeR, metagenomeSeq, ZicoSeq, and more
- **🎯 Consensus Analysis**: Combines results from multiple methods for robust findings
- **📈 Rich Visualizations**: Generate comprehensive plots and interactive dashboards
- **⚡ Easy to Use**: Simple Python API with sensible defaults

## Quick Start

```python
import pandas as pd
from daa_advisor import DifferentialAbundanceTool

# Load your data
count_table = pd.read_csv("counts.csv", index_col=0)  # Samples x Features
metadata = pd.read_csv("metadata.csv", index_col=0)   # Sample metadata

# Run analysis with automatic method selection
tool = DifferentialAbundanceTool()
results = tool.analyze(
    count_table=count_table,
    metadata=metadata,
    data_type='asv',  # or 'gene' or 'viral'
    use_consensus=True
)

# View results
tool.summarize_results()

# Get significant features
significant = tool.get_significant_features(alpha=0.05)
print(f"Found {len(significant)} significant features")
```

## Installation

### Basic Installation

```bash
# Install from PyPI (when available)
pip install daaadvisor

# Or install from source
git clone https://github.com/yourusername/daaadvisor.git
cd daaadvisor
pip install -e .
```

### With R Methods Support

To use R-based methods (ALDEx2, ANCOM-BC, DESeq2, etc.), you'll need R and specific packages:

```bash
# Install R dependencies
conda install -c conda-forge r-base=4.3

# Install Python package with R support
pip install daaadvisor[r]
```

Then in R:
```r
# Install required R packages
install.packages("BiocManager")
BiocManager::install(c("ALDEx2", "ANCOMBC", "DESeq2", "edgeR", "metagenomeSeq"))
```

## How It Works

### 1. Data Profiling
DAAadvisor analyzes your data characteristics:
- **Sparsity**: Proportion of zeros
- **Zero inflation**: Features with excessive zeros  
- **Compositional bias**: Correlation between abundant features
- **Sample size**: Number of samples per group
- **Data type**: ASV/16S, genes, or viral

### 2. Method Selection
Based on benchmarking studies, the tool scores each method:
- **ALDEx2**: Best for compositional data with moderate sparsity
- **ANCOM-BC**: Good for compositional bias correction
- **DESeq2**: Powerful for gene data with low sparsity
- **ZicoSeq**: Excels with high sparsity and zero inflation
- **And more...**

### 3. Consensus Analysis
Optionally runs multiple methods and combines results for robust findings.

## Supported Methods

| Method | Best For | Handles Compositionality | Min Samples |
|--------|----------|-------------------------|-------------|
| ALDEx2 | ASV data, compositional bias | ✅ | 8 |
| ANCOM-BC | ASV/gene, balanced design | ✅ | 10 |
| DESeq2 | Gene data, complex designs | ❌ | 6 |
| edgeR | Gene data, large samples | ❌ | 10 |
| metagenomeSeq | Moderate sparsity | ✅ | 8 |
| ZicoSeq | High sparsity, viral data | ✅ | 5 |
| Wilcoxon | Small samples, robust | ❌ | 3 |

## Examples

### Basic Analysis
```python
from daa_advisor import DifferentialAbundanceTool

tool = DifferentialAbundanceTool()
results = tool.analyze(count_table, metadata)
tool.summarize_results()
```

### Advanced Usage
```python
# Custom parameters
results = tool.analyze(
    count_table=counts,
    metadata=meta,
    data_type='gene',
    use_consensus=True,
    # Method-specific parameters
    method_params={
        'deseq2': {'test': 'Wald', 'fitType': 'local'},
        'aldex2': {'mc.samples': 128}
    }
)

# Detailed method information
from daa_advisor import MethodSelector
selector = MethodSelector()
print(selector.get_method_info('aldex2'))
```

### Data Profiling Only
```python
from daa_advisor import DataProfiler

profiler = DataProfiler()
profile = profiler.profile_data(count_table, metadata)
profiler.print_profile_summary()
```

## Data Format

### Count Table (CSV)
Samples as rows, features as columns:
```
        ASV_1   ASV_2   ASV_3
Sample1   10      0      45
Sample2    0     23       0  
Sample3   15      5      12
```

### Metadata (CSV)
```
        condition   batch   age
Sample1   Control     A     25
Sample2 Treatment     A     30
Sample3   Control     B     28
```

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md).

## Citation

If you use DAAadvisor in your research, please cite:

```
@software{daaadvisor,
  title = {DAAadvisor: Intelligent Differential Abundance Analysis for Microbiome Data},
  author = {DAAadvisor Team},
  year = {2024},
  url = {https://github.com/yourusername/daaadvisor}
}
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- 📖 **Documentation**: [Read the Docs](https://daaadvisor.readthedocs.io)
- 🐛 **Issues**: [GitHub Issues](https://github.com/yourusername/daaadvisor/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/yourusername/daaadvisor/discussions)

---

**DAAadvisor** - Making differential abundance analysis intelligent and reproducible! 🧬✨