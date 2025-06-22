# Real Microbiome Data for Publication Benchmarking

## 🎯 Overview

For publication-quality benchmarking, DAAadvisor supports both **realistic simulated data** and **real published datasets**. This guide shows you how to obtain and use real microbiome data.

## 📊 Available Real Datasets

### 🦠 Disease State Studies

| Dataset | Study | Samples | Data Type | Ground Truth | Download |
|---------|-------|---------|-----------|--------------|----------|
| **IBD Multi-omics** | Franzosa et al. 2019 | 132 | Metagenomics | IBD vs Healthy | [IBDMDB.org](https://ibdmdb.org/downloads) |
| **CRC Meta-analysis** | Wirbel et al. 2019 | 768 | Metagenomics | Cancer vs Healthy | [Zenodo](https://zenodo.org/record/3250873) |
| **Type 2 Diabetes** | Qin et al. 2012 | 344 | Metagenomics | T2D vs Healthy | [NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA179209) |

### 💊 Antibiotic Studies

| Dataset | Study | Samples | Data Type | Ground Truth | Download |
|---------|-------|---------|-----------|--------------|----------|
| **Antibiotic Perturbation** | Dethlefsen & Relman 2008 | 50 | 16S | Pre/During/Post | [NCBI SRA](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA20573) |
| **Recovery Dynamics** | Raymond et al. 2016 | 82 | 16S | Baseline/Post-treatment | [Qiita](https://qiita.ucsd.edu/study/description/10349) |

### 🌍 Large Collections

| Dataset | Study | Samples | Data Type | Ground Truth | Download |
|---------|-------|---------|-----------|--------------|----------|
| **American Gut** | McDonald et al. 2018 | 11,842 | 16S | Healthy variation | [Qiita](https://qiita.ucsd.edu/) |
| **curatedMetagenomicData** | Pasolli et al. 2017 | 5,716 | Metagenomics | Multiple studies | [Bioconductor](https://bioconductor.org/packages/curatedMetagenomicData/) |

## 🚀 Quick Start: Realistic Data

### Option 1: Generate Realistic Demo Data (Immediate)

```bash
# Generate realistic data based on published studies
python download_real_data.py

# This creates:
# realistic_demo_data/
# ├── ibd_count_table.csv          # IBD study (132 samples)
# ├── ibd_metadata.csv
# ├── crc_count_table.csv          # CRC study (200 samples) 
# ├── crc_metadata.csv
# ├── antibiotic_count_table.csv   # Antibiotic study (60 samples)
# └── antibiotic_metadata.csv

# Run benchmark with realistic data
python run_publication_benchmark.py --full --output realistic_benchmark_results
```

### Option 2: Download Real Published Data

## 📥 Download Methods

### 🔧 Method 1: Direct Download (Easiest)

```bash
# Install required tools
conda install wget curl

# IBD data from IBDMDB
mkdir -p real_data/ibd_franzosa_2019
cd real_data/ibd_franzosa_2019

wget "https://ibdmdb.org/downloads/HMP2/metadata/hmp2_metadata.csv"
wget "https://ibdmdb.org/downloads/HMP2/mgx/taxonomic_profiles.tsv.gz"

# CRC data from Zenodo  
mkdir -p ../crc_wirbel_2019
cd ../crc_wirbel_2019

wget "https://zenodo.org/record/3250873/files/species_abundance.tsv"
wget "https://zenodo.org/record/3250873/files/metadata.tsv"
```

### 🧬 Method 2: SRA Toolkit (For NCBI Data)

```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Download antibiotic study data
mkdir -p real_data/antibiotic_dethlefsen_2008

# Get sample accession list from NCBI BioProject PRJNA20573
esearch -db sra -query "PRJNA20573" | efetch -format runinfo > runinfo.csv

# Download first few samples
fastq-dump --split-files --gzip SRR013294 SRR013295 SRR013296
```

### 🌐 Method 3: Qiita Database

```bash
# 1. Register at https://qiita.ucsd.edu/
# 2. Search for American Gut Project (Study 10317)
# 3. Download processed OTU table and metadata
# 4. Convert BIOM to TSV format

biom convert -i otu_table.biom -o otu_table.tsv --to-tsv
```

### 🧬 Method 4: curatedMetagenomicData (R Package)

```r
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("curatedMetagenomicData")

library(curatedMetagenomicData)

# Download IBD data
ibd_data <- curatedMetagenomicData("IBD", dataType="relative_abundance", rownames="short")

# Save for DAAadvisor
write.csv(assay(ibd_data), "ibd_abundance.csv")
write.csv(colData(ibd_data), "ibd_metadata.csv")

# Download CRC data
crc_data <- curatedMetagenomicData("CRC", dataType="relative_abundance", rownames="short")
write.csv(assay(crc_data), "crc_abundance.csv") 
write.csv(colData(crc_data), "crc_metadata.csv")
```

## 📋 Data Format Requirements

DAAadvisor expects data in this format:

### Count Table (`count_table.csv`)
```
             ASV_1  ASV_2  ASV_3  ...
Sample_1      10     0     45    ...
Sample_2       0    23      0    ...
Sample_3      15     5     12    ...
```

### Metadata (`metadata.csv`)
```
          condition   batch  age  disease_state
Sample_1    Control      A   25        Healthy
Sample_2  Treatment      A   30            IBD
Sample_3    Control      B   28        Healthy
```

### Ground Truth (`differential_features.txt`)
```
ASV_1
ASV_5
ASV_12
...
```

## 🎯 Using Real Data with DAAadvisor

### Standard Benchmark

```python
from daa_advisor import run_publication_benchmark

# Run with real data
results = run_publication_benchmark(
    output_dir="real_data_benchmark",
    n_bootstrap=100,
    datasets_dir="real_data/"  # Directory with your real datasets
)
```

### Custom Real Dataset

```python
import pandas as pd
from daa_advisor import DifferentialAbundanceTool

# Load your real data
count_table = pd.read_csv("your_count_table.csv", index_col=0)
metadata = pd.read_csv("your_metadata.csv", index_col=0)

# Run analysis
tool = DifferentialAbundanceTool()
results = tool.analyze(
    count_table=count_table,
    metadata=metadata,
    use_consensus=True
)

tool.summarize_results()
```

## 🏆 Publication Strategy

### Phase 1: Validation with Known Studies
1. **IBD Study**: Use Franzosa et al. 2019 data as gold standard
2. **CRC Study**: Validate with Wirbel et al. 2019 meta-analysis  
3. **Antibiotic Study**: Test with Dethlefsen & Relman 2008

### Phase 2: Cross-Study Validation
1. **Multiple IBD Studies**: Test consistency across different IBD cohorts
2. **Geographic Validation**: American vs European populations
3. **Technical Validation**: 16S vs shotgun metagenomics

### Phase 3: Novel Applications
1. **Longitudinal Studies**: Time-series antibiotic recovery
2. **Multi-omics**: Integrate metabolomics and proteomics
3. **Clinical Prediction**: Diagnostic biomarker identification

## 📚 Key Publications for Ground Truth

### Established Findings to Validate Against:

**IBD (Inflammatory Bowel Disease):**
- ↑ Enterobacteriaceae, ↓ Faecalibacterium prausnitzii
- ↑ Escherichia coli, ↓ Bifidobacterium
- Paper: [Franzosa et al. Nature Microbiology 2019](https://doi.org/10.1038/s41564-018-0306-4)

**Colorectal Cancer:**
- ↑ Fusobacterium nucleatum, ↑ Bacteroides fragilis
- ↓ Butyrate-producing bacteria
- Paper: [Wirbel et al. Nature Medicine 2019](https://doi.org/10.1038/s41591-019-0406-6)

**Type 2 Diabetes:**
- ↓ Butyrate producers, ↑ Opportunistic pathogens
- ↓ Akkermansia muciniphila
- Paper: [Qin et al. Nature 2012](https://doi.org/10.1038/nature11450)

**Antibiotic Effects:**
- Massive diversity loss (↓ Shannon index)
- ↑ Enterobacteriaceae, ↓ Bacteroidetes
- Paper: [Dethlefsen & Relman PLOS Biology 2008](https://doi.org/10.1371/journal.pbio.0060280)

## 🛠️ Troubleshooting

### Common Issues:

1. **Large File Downloads**: Use `wget -c` for resumable downloads
2. **SRA Access**: Some data requires NCBI account and dbGaP approval
3. **Format Conversion**: Use `biom convert` for BIOM format files
4. **Memory Issues**: Process large datasets in chunks

### Data Quality Checks:

```python
# Check data integrity
print(f"Count table shape: {count_table.shape}")
print(f"Metadata shape: {metadata.shape}")
print(f"Sample overlap: {len(set(count_table.index) & set(metadata.index))}")
print(f"Zero proportion: {(count_table == 0).sum().sum() / count_table.size:.2%}")
```

## 📞 Support

For help with real data download:
- **SRA Data**: [NCBI SRA Documentation](https://www.ncbi.nlm.nih.gov/sra/docs/)
- **Qiita**: [Qiita Help](https://qiita.ucsd.edu/iframe/?iframe=qiita-help)
- **curatedMetagenomicData**: [Bioconductor Vignette](https://bioconductor.org/packages/release/data/experiment/vignettes/curatedMetagenomicData/inst/doc/curatedMetagenomicData.html)

---

**Ready to run publication-quality benchmarks with real microbiome data!** 🧬✨