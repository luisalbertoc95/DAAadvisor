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

__all__ = [
    'DifferentialAbundanceTool',
    'DataProfiler', 
    'MethodSelector'
]