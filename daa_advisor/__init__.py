"""
DAAadvisor: Differential Abundance Analysis Tool
==============================================

An intelligent tool for microbiome differential abundance analysis that 
automatically selects the best statistical method based on your data characteristics.

Main classes:
- DifferentialAbundanceTool: Main analysis orchestrator
- DataProfiler: Analyzes data characteristics 
- MethodSelector: Recommends optimal methods
"""

__version__ = "0.1.0"
__author__ = "DAAadvisor Team"
__email__ = "support@daaadvisor.org"

# Import main classes for convenient access
from .core import DifferentialAbundanceTool
from .profiler import DataProfiler
from .selector import MethodSelector
from .visualization import DAAVisualizer, create_comprehensive_report
from .benchmarking import MethodBenchmark, run_full_benchmark
from .data_generators import MicrobiomeDataGenerator, create_benchmark_datasets
from .information_theory import CompositionInformationFramework, MaximumEntropySelector, run_information_theory_analysis
from .preprocessing import InformationBasedPreprocessor, apply_information_preprocessing
from .adaptive_thresholds import AdaptiveThresholdSelector, apply_adaptive_thresholds

__all__ = [
    'DifferentialAbundanceTool',
    'DataProfiler', 
    'MethodSelector',
    'DAAVisualizer',
    'create_comprehensive_report',
    'MethodBenchmark',
    'run_full_benchmark',
    'MicrobiomeDataGenerator',
    'create_benchmark_datasets',
    'CompositionInformationFramework',
    'MaximumEntropySelector',
    'run_information_theory_analysis',
    'InformationBasedPreprocessor',
    'apply_information_preprocessing',
    'AdaptiveThresholdSelector',
    'apply_adaptive_thresholds'
]