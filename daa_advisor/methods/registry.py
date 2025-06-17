#!/usr/bin/env python3
"""
Method registry for DAAadvisor
"""

import logging
from typing import Dict, List
from .base import DAAMethod
from .wilcoxon import WilcoxonMethod

logger = logging.getLogger(__name__)


class MethodRegistry:
    """Registry for all available differential abundance methods"""
    
    def __init__(self):
        self.methods = {}
        self._register_default_methods()
    
    def _register_default_methods(self):
        """Register built-in methods"""
        # Always available methods
        self.register_method(WilcoxonMethod())
        
        # Try to register R-based methods
        self._try_register_r_methods()
    
    def _try_register_r_methods(self):
        """Attempt to register R-based methods if available"""
        try:
            # Try importing R interface
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
            
            # Check for required R packages
            r_packages = {
                'ALDEx2': 'aldex2',
                'ANCOMBC': 'ancom-bc', 
                'DESeq2': 'deseq2',
                'edgeR': 'edger',
                'metagenomeSeq': 'metagenomeseq'
            }
            
            available_packages = []
            for r_pkg, method_name in r_packages.items():
                try:
                    robjects.r(f'library({r_pkg})')
                    available_packages.append((r_pkg, method_name))
                    logger.info(f"R package {r_pkg} available")
                except Exception:
                    logger.warning(f"R package {r_pkg} not available")
            
            # Register available R methods
            for r_pkg, method_name in available_packages:
                try:
                    if method_name == 'aldex2':
                        from .r_methods import ALDEx2Method
                        self.register_method(ALDEx2Method())
                    elif method_name == 'ancom-bc':
                        from .r_methods import ANCOMBCMethod
                        self.register_method(ANCOMBCMethod())
                    elif method_name == 'deseq2':
                        from .r_methods import DESeq2Method
                        self.register_method(DESeq2Method())
                    # Add other R methods as they are implemented
                except ImportError as e:
                    logger.warning(f"Could not import {method_name}: {e}")
                    
        except ImportError:
            logger.warning("rpy2 not available, R methods will not be registered")
        except Exception as e:
            logger.warning(f"Error registering R methods: {e}")
    
    def register_method(self, method: DAAMethod):
        """Register a new method"""
        method_name = method.name()
        self.methods[method_name] = method
        logger.info(f"Registered method: {method_name}")
    
    def get_method(self, name: str) -> DAAMethod:
        """Get a method by name"""
        if name not in self.methods:
            raise ValueError(f"Method '{name}' not available. Available methods: {list(self.methods.keys())}")
        return self.methods[name]
    
    def has_method(self, name: str) -> bool:
        """Check if method is available"""
        return name in self.methods
    
    def list_methods(self) -> List[str]:
        """List all available methods"""
        return list(self.methods.keys())
    
    def get_fallback_method(self) -> DAAMethod:
        """Get a reliable fallback method"""
        # Wilcoxon is always available as it's implemented in pure Python
        return self.methods['wilcoxon']
    
    def get_method_info(self, name: str) -> Dict:
        """Get information about a method"""
        if name not in self.methods:
            return {}
        
        method = self.methods[name]
        return {
            'name': method.name(),
            'parameters': method.get_parameters(),
            'citation': method.cite()
        }
    
    def get_available_methods_summary(self) -> Dict:
        """Get summary of all available methods"""
        summary = {}
        for name, method in self.methods.items():
            summary[name] = {
                'available': True,
                'citation': method.cite(),
                'parameters': list(method.get_parameters().keys())
            }
        return summary