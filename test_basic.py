#!/usr/bin/env python3
"""
Basic test of DAAadvisor without external dependencies
"""

import sys
import os

# Test the imports work
try:
    from daa_advisor.profiler import DataProfiler, DataProfile
    from daa_advisor.selector import MethodSelector, MethodRecommendation  
    from daa_advisor.methods import MethodRegistry, WilcoxonMethod
    from daa_advisor.core import DifferentialAbundanceTool
    print("✅ All imports successful!")
    
    # Test basic instantiation
    profiler = DataProfiler()
    selector = MethodSelector()
    registry = MethodRegistry()
    tool = DifferentialAbundanceTool()
    
    print("✅ All classes instantiated successfully!")
    
    # Test method registry
    available_methods = registry.list_methods()
    print(f"📋 Available methods: {available_methods}")
    
    # Test method info
    if 'wilcoxon' in available_methods:
        method_info = registry.get_method_info('wilcoxon')
        print(f"🔬 Wilcoxon method info: {method_info}")
    
    print("\n🎉 Basic tests passed! DAAadvisor is properly structured.")
    print("\n💡 To run full analysis, install dependencies:")
    print("   pip install numpy pandas scipy scikit-learn statsmodels")
    
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"❌ Error: {e}")
    sys.exit(1)