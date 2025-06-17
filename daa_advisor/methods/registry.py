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
            # Import and check available R methods
            from .r_methods import check_r_package_availability
            available_methods = check_r_package_availability()
            
            # Register each available method
            for method_name, method_class in available_methods.items():
                try:
                    self.register_method(method_class())
                    logger.info(f"Successfully registered R method: {method_name}")
                except Exception as e:
                    logger.warning(f"Failed to register {method_name}: {e}")
                    
        except ImportError:
            logger.warning("R methods module not available")
        except Exception as e:
            logger.warning(f"Error checking R methods: {e}")
    
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