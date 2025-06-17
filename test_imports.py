#!/usr/bin/env python3
"""
Test just the import structure without numpy/pandas
"""

import sys
import os

# Test individual imports
try:
    print("Testing basic Python imports...")
    from typing import Dict, List, Optional
    from dataclasses import dataclass
    from abc import ABC, abstractmethod
    import logging
    print("✅ Basic imports successful")
    
    print("Testing DAAadvisor structure...")
    # Test if the modules exist
    import daa_advisor
    print("✅ Main package imported")
    
    # Test if __init__.py works
    print(f"DAAadvisor version: {daa_advisor.__version__}")
    print(f"Available: {daa_advisor.__all__}")
    
    print("\n🎉 Package structure is correct!")
    print("📦 Ready to install dependencies and run full analysis")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)