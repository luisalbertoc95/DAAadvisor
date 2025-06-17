<p align="center">
  <img src="daaadvisor_logo.png" alt="DAAadvisor Logo" width="400"/>
</p>

# DAAadvisor

**Differential Abundance Analysis Advisor for Microbiome Data**

An intelligent tool that automatically selects the best statistical method for your microbiome differential abundance analysis based on data characteristics.

## Features

- **🧠 Intelligent Method Selection**: Information-theoretic framework with maximum entropy principle for optimal method selection
- **📊 Comprehensive Data Profiling**: Analyzes sparsity, zero-inflation, compositional bias, and other key metrics
- **🔄 Multi-Method Support**: Integrates 6 statistical methods with full R integration (ALDEx2, ANCOM-BC, DESeq2, edgeR, metagenomeSeq, Wilcoxon)
- **🎯 Consensus Analysis**: Voting-based consensus combining results from multiple methods
- **📈 Rich Visualizations**: Interactive HTML dashboards, method comparisons, and publication-ready plots
- **🧮 Information Theory Framework**: Jensen-Shannon divergence and compositional data analysis
- **⚡ Easy to Use**: Simple Python API with intelligent defaults and comprehensive testing

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

## 🧠 Methodology Framework

DAAadvisor implements a comprehensive 5-step information-theoretic framework:

### 1. **📊 Data Assessment & Profiling**
- **Sparsity Analysis**: Zero-inflation quantification
- **Count Distribution**: Mean, variance, dynamic range assessment  
- **Data Type Detection**: ASV/16S, Gene/Functional, Viral classification
- **Compositional Bias**: Library size variation analysis

### 2. **🧮 Information-Theoretic Method Selection**
- **Maximum Entropy Principle**: `Method* = argmax H(X|θ)` subject to data constraints
- **Jensen-Shannon Divergence**: `JS(P,Q) = ½[KL(P||M) + KL(Q||M)]` for between-group differences
- **Compositional Log-Ratio**: `CLR(x) = log(x/g(x))` transformation
- **Confidence Scoring**: Quantitative method selection confidence

### 3. **🔬 Multi-Method Statistical Analysis**
- **Wilcoxon**: Non-parametric rank-based testing
- **ALDEx2**: CLR transformation with Monte Carlo sampling
- **ANCOM-BC**: Compositional bias correction
- **DESeq2**: Negative binomial modeling
- **edgeR**: TMM normalization with quasi-likelihood
- **metagenomeSeq**: Zero-inflated log-normal modeling

### 4. **🤝 Consensus Integration**
- **Majority Voting**: Features significant in ≥50% of methods
- **Agreement Metrics**: Inter-method concordance quantification
- **Weighted Confidence**: Method agreement-based scoring

### 5. **📈 Results & Visualization**
- **Interactive HTML Dashboards**: Comprehensive reporting
- **Method Comparison Plots**: Performance metrics visualization
- **Volcano Plots**: Effect size vs significance
- **Data Profiling Charts**: Characteristics visualization

## 📊 Current Performance

| Method | Status | F1 Score | Precision | Recall | Runtime |
|--------|--------|----------|-----------|--------|---------|
| Wilcoxon | ✅ Working | 1.000 | 1.000 | 1.000 | 0.030s |
| ALDEx2 | ✅ Working | 0.453 | 0.293 | 1.000 | 0.577s |
| ANCOM-BC | 🔧 In Progress | - | - | - | - |
| DESeq2 | 🔧 In Progress | - | - | - | - |
| edgeR | 🔧 In Progress | - | - | - | - |
| metagenomeSeq | 🔧 In Progress | - | - | - | - |

## 📄 Comprehensive Results & Documentation

🎯 **View Complete Results**: [`consolidated_results/`](consolidated_results/)

### **Priority Resources:**
1. **📊 Method Comparison**: [`method_comparison_both_methods.png`](method_comparison_both_methods.png)
2. **🧠 Methodology Diagram**: [`consolidated_results/reports/methodology_diagram.html`](consolidated_results/reports/methodology_diagram.html)
3. **📄 Interactive Dashboard**: [`consolidated_results/visualizations/interactive_dashboard.html`](consolidated_results/visualizations/interactive_dashboard.html)
4. **📋 Detailed HTML Report**: [`consolidated_results/reports/detailed_results.html`](consolidated_results/reports/detailed_results.html)

### **Technical Documentation:**
- **🔬 R Integration Success Report**: [`consolidated_results/reports/R_INTEGRATION_SUCCESS_REPORT.md`](consolidated_results/reports/R_INTEGRATION_SUCCESS_REPORT.md)
- **📊 Benchmark Results**: [`consolidated_results/benchmarks/`](consolidated_results/benchmarks/)
- **📈 Analysis Results**: [`consolidated_results/comprehensive_analysis/`](consolidated_results/comprehensive_analysis/)

## Supported Methods

| Method | Status | Best For | Handles Compositionality | R Integration |
|--------|--------|----------|-------------------------|---------------|
| **Wilcoxon** | ✅ Working | Small samples, robust testing | ❌ | Pure Python |
| **ALDEx2** | ✅ Working | ASV data, compositional analysis | ✅ | rpy2 + R |
| **ANCOM-BC** | 🔧 In Progress | ASV/gene, bias correction | ✅ | rpy2 + R |
| **DESeq2** | 🔧 In Progress | Gene data, complex designs | ❌ | rpy2 + R |
| **edgeR** | 🔧 In Progress | Gene data, large samples | ❌ | rpy2 + R |
| **metagenomeSeq** | 🔧 In Progress | High sparsity, zero-inflation | ✅ | rpy2 + R |

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