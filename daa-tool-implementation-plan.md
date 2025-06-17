# Differential Abundance Analysis Tool - Implementation Plan

## Executive Summary
A comprehensive tool that integrates multiple differential abundance testing methods and automatically recommends the best approach based on data characteristics, metadata properties, and data type (viral, ASVs, genes).

## Current State of the Field

### Key Findings from Literature Review

1. **Method Diversity**: At least 14+ commonly used differential abundance methods exist, including:
   - **Most Consistent**: ALDEx2, ANCOM-II
   - **High Power**: DESeq2, edgeR (especially for gene data)
   - **Compositionally Aware**: ANCOM-BC, metagenomeSeq, DACOMP
   - **Newer Methods**: ZicoSeq, LinDA, MaAsLin3
   - **Others**: LEfSe, Wilcoxon, t-test, corncob, LDM

2. **Major Challenges**:
   - Methods produce drastically different results on same data
   - No single method performs best across all scenarios
   - Performance depends heavily on data characteristics
   - No existing tool automatically selects optimal method

3. **Data-Specific Considerations**:
   - **16S/ASV data**: High sparsity, compositional effects dominant
   - **Metagenomics/Gene data**: Less sparse, more suitable for RNA-seq methods
   - **Viral data**: Extreme sparsity, presence/absence patterns important

## Tool Architecture

### Core Components

#### 1. Data Profiler Module
```python
class DataProfiler:
    """Analyzes input data characteristics"""
    
    def profile_data(self, count_table, metadata):
        return {
            'data_type': self.detect_data_type(),  # viral/ASV/gene
            'sparsity': self.calculate_sparsity(),
            'sample_size': self.get_sample_size(),
            'sequencing_depth': self.analyze_depth(),
            'zero_inflation': self.assess_zero_inflation(),
            'compositional_effects': self.detect_compositional_bias(),
            'group_balance': self.check_group_balance(),
            'effect_size_estimate': self.estimate_effect_sizes(),
            'metadata_complexity': self.analyze_metadata()
        }
```

#### 2. Method Selector Module
```python
class MethodSelector:
    """Recommends best methods based on data profile"""
    
    def recommend_methods(self, data_profile):
        # Decision tree based on benchmarking results
        recommendations = []
        
        # Primary recommendation
        primary_method = self.select_primary_method(data_profile)
        
        # Secondary methods for consensus
        secondary_methods = self.select_secondary_methods(data_profile)
        
        return {
            'primary': primary_method,
            'secondary': secondary_methods,
            'reasoning': self.generate_explanation()
        }
```

#### 3. Method Integration Layer
```python
class MethodIntegrator:
    """Wrapper for all differential abundance methods"""
    
    methods = {
        'aldex2': ALDEx2Wrapper(),
        'ancombc': ANCOMBCWrapper(),
        'deseq2': DESeq2Wrapper(),
        'edger': EdgeRWrapper(),
        'metagenomeseq': MetagenomeSeqWrapper(),
        'zicoseq': ZicoSeqWrapper(),
        'linda': LinDAWrapper(),
        'maaslin3': MaAsLin3Wrapper(),
        # ... other methods
    }
```

#### 4. Consensus Module
```python
class ConsensusAnalyzer:
    """Combines results from multiple methods"""
    
    def generate_consensus(self, results_dict):
        # Implement voting system
        # Weight by method performance for data type
        # Generate confidence scores
        pass
```

## Implementation Phases

### Phase 1: Foundation (Weeks 1-4)
- [ ] Set up project structure and dependencies
- [ ] Implement data profiler for basic characteristics
- [ ] Create wrapper classes for 5 core methods:
  - ALDEx2
  - ANCOM-BC
  - DESeq2
  - metagenomeSeq
  - ZicoSeq

### Phase 2: Method Selection Logic (Weeks 5-8)
- [ ] Implement decision tree based on benchmarking studies
- [ ] Create scoring system for method-data compatibility
- [ ] Build recommendation engine with explanations
- [ ] Add support for different data types (viral/ASV/gene)

### Phase 3: Extended Methods (Weeks 9-12)
- [ ] Add remaining method wrappers
- [ ] Implement consensus analysis module
- [ ] Create visualization components
- [ ] Build comprehensive testing framework

### Phase 4: Advanced Features (Weeks 13-16)
- [ ] Add confounder detection and adjustment
- [ ] Implement effect size estimation
- [ ] Create interactive dashboard
- [ ] Add batch processing capabilities

## Decision Logic Framework

### Primary Method Selection Rules

```yaml
data_type: ASV/16S
  if sample_size < 20:
    if sparsity > 0.8:
      primary: ANCOM-II
    else:
      primary: ALDEx2
  else:
    if compositional_effects == high:
      primary: ANCOM-BC
    else:
      primary: ALDEx2

data_type: metagenomics/genes
  if sample_size < 10:
    primary: metagenomeSeq
  else:
    if zero_inflation > 0.7:
      primary: ZicoSeq
    else:
      primary: DESeq2

data_type: viral
  if extreme_sparsity:
    primary: ZicoSeq or ANCOM-BC
  else:
    primary: DESeq2 with proper normalization
```

## Key Features

### 1. Automated Data Assessment
- Sparsity metrics
- Zero-inflation tests
- Compositional bias detection
- Sample size adequacy
- Sequencing depth variation

### 2. Smart Method Selection
- Rule-based initial selection
- Machine learning refinement (future)
- Explanation generation
- Confidence scoring

### 3. Comprehensive Reporting
- Method comparison table
- Consensus results
- Visualization suite
- Diagnostic plots
- Reproducible workflows

### 4. Data Type Specialization
- ASV-specific optimizations
- Gene-level adaptations
- Viral data handling
- Cross-domain comparisons

## Technical Stack

### Core Dependencies
```python
# Statistical Methods
- R packages: ALDEx2, ANCOMBC, DESeq2, edgeR, metagenomeSeq
- Python: scipy, statsmodels, scikit-learn

# Data Processing
- pandas, numpy
- biom-format
- rpy2 (R-Python interface)

# Visualization
- matplotlib, seaborn
- plotly (interactive plots)

# Workflow
- snakemake or nextflow
- docker/singularity
```

## Validation Strategy

### 1. Benchmark Against Published Studies
- Recreate results from Nearing et al. 2022
- Compare with Zhou et al. 2022 (ZicoSeq paper)
- Validate on Wirbel et al. 2024 datasets

### 2. Simulation Framework
- Generate synthetic data with known truth
- Test across parameter ranges
- Validate method selection accuracy

### 3. Real-World Testing
- Multiple disease cohorts
- Environmental samples
- Different sequencing platforms

## Expected Outcomes

### For Users
1. **Automated Analysis**: No need to manually select methods
2. **Increased Confidence**: Consensus approach reduces false positives
3. **Better Reproducibility**: Standardized workflows
4. **Educational Value**: Explanations for method choices

### For the Field
1. **Standardization**: Common framework for comparisons
2. **Method Development**: Platform for testing new approaches
3. **Best Practices**: Data-driven recommendations
4. **Transparency**: Open decision logic

## Future Extensions

### Version 2.0
- Machine learning for method selection
- Multi-omics integration
- Time-series support
- Batch effect correction

### Version 3.0
- Cloud deployment
- GUI interface
- Automated report generation
- Integration with major pipelines

## References

Key papers informing this design:
1. Nearing et al. 2022 - Comparison of 14 methods across 38 datasets
2. Zhou et al. 2022 - ZicoSeq development and benchmarking
3. Wirbel et al. 2024 - Realistic benchmarking framework
4. Calgaro et al. 2022 - Method performance evaluation
5. Yang et al. 2024 - Consensus approaches

## Getting Started

### Initial Development Steps
1. Fork template repository
2. Set up development environment
3. Implement data profiler
4. Create first method wrapper (suggest ALDEx2)
5. Build basic CLI interface
6. Test on example datasets

### Recommended Dataset for Testing
- American Gut Project (diverse, well-characterized)
- HMP2 (multiple body sites)
- Synthetic data from mockrobiota