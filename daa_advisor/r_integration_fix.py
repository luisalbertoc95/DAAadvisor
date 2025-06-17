#!/usr/bin/env python3
"""
Enhanced R Integration with Better Error Handling and Installation Guidance

This module provides improved R integration with:
- Clear dependency checking
- Installation guidance
- Graceful degradation
- Detailed error reporting
"""

import logging
import subprocess
import sys
from typing import Dict, List, Tuple

logger = logging.getLogger(__name__)

class RIntegrationManager:
    """Manages R integration with better error handling and guidance"""
    
    def __init__(self):
        self.rpy2_available = False
        self.r_available = False
        self.available_packages = {}
        self.missing_packages = []
        self.error_messages = []
        
        self._check_dependencies()
    
    def _check_dependencies(self):
        """Comprehensive dependency checking"""
        
        # 1. Check rpy2 availability
        try:
            import rpy2.robjects as robjects
            from rpy2.robjects.packages import importr
            from rpy2.rinterface_lib.embedded import RRuntimeError
            
            self.rpy2_available = True
            logger.info("✅ rpy2 is available")
            
            # 2. Check R availability
            try:
                r_version = robjects.r('R.version.string')[0]
                self.r_available = True
                logger.info(f"✅ R is available: {r_version}")
                
                # 3. Check R packages
                self._check_r_packages(importr, RRuntimeError)
                
            except Exception as e:
                self.r_available = False
                self.error_messages.append(f"R not available: {e}")
                logger.warning(f"❌ R not available: {e}")
                
        except ImportError as e:
            self.rpy2_available = False
            self.error_messages.append(f"rpy2 not installed: {e}")
            logger.warning(f"❌ rpy2 not available: {e}")
    
    def _check_r_packages(self, importr, RRuntimeError):
        """Check availability of required R packages"""
        
        required_packages = {
            'ALDEx2': 'Differential abundance analysis with CLR transformation',
            'ANCOMBC': 'Bias correction for compositional data',
            'DESeq2': 'Negative binomial modeling for count data', 
            'edgeR': 'Empirical analysis of DGE in multifactor experiments',
            'metagenomeSeq': 'Statistical analysis for sparse high-throughput data'
        }
        
        for package, description in required_packages.items():
            try:
                pkg = importr(package)
                self.available_packages[package] = {
                    'imported': pkg,
                    'description': description
                }
                logger.info(f"✅ {package} available")
                
            except RRuntimeError as e:
                self.missing_packages.append(package)
                self.error_messages.append(f"{package} not installed: {e}")
                logger.warning(f"❌ {package} not available")
    
    def get_installation_guide(self) -> str:
        """Generate installation guide based on missing dependencies"""
        
        guide = "\n" + "="*60 + "\n"
        guide += "🔧 R METHODS INSTALLATION GUIDE\n"
        guide += "="*60 + "\n\n"
        
        if not self.rpy2_available:
            guide += "📦 STEP 1: Install rpy2\n"
            guide += "-" * 25 + "\n"
            guide += "pip install rpy2\n"
            guide += "# or with conda:\n"
            guide += "conda install -c conda-forge rpy2\n\n"
        
        if not self.r_available:
            guide += "📊 STEP 2: Install R\n" 
            guide += "-" * 18 + "\n"
            guide += "# macOS:\n"
            guide += "brew install r\n\n"
            guide += "# Or using conda:\n"
            guide += "conda install -c conda-forge r-base=4.3\n\n"
            guide += "# Or download from: https://cran.r-project.org/\n\n"
        
        if self.missing_packages:
            guide += "🧬 STEP 3: Install R Bioconductor Packages\n"
            guide += "-" * 42 + "\n"
            guide += "# Start R console and run:\n"
            guide += 'R --vanilla <<EOF\n'
            guide += 'install.packages("BiocManager")\n'
            guide += f'BiocManager::install(c({", ".join(f"\"{pkg}\"" for pkg in self.missing_packages)}))\n'
            guide += 'EOF\n\n'
            
            guide += "# Or interactively in R:\n"
            guide += "R\n"
            guide += '> install.packages("BiocManager")\n'
            for pkg in self.missing_packages:
                guide += f'> BiocManager::install("{pkg}")\n'
            guide += "\n"
        
        guide += "🧪 STEP 4: Verify Installation\n"
        guide += "-" * 30 + "\n" 
        guide += "python -c \"from daa_advisor.r_integration_fix import RIntegrationManager; RIntegrationManager().print_status()\"\n\n"
        
        guide += "💡 TROUBLESHOOTING TIPS\n"
        guide += "-" * 23 + "\n"
        guide += "• R packages require Bioconductor (not regular CRAN)\n"
        guide += "• Some packages need additional system dependencies\n"
        guide += "• On macOS, you may need Xcode command line tools: xcode-select --install\n"
        guide += "• For Ubuntu/Debian: apt-get install r-base-dev libcurl4-openssl-dev\n\n"
        
        return guide
    
    def print_status(self):
        """Print detailed status of R integration"""
        
        print("\n" + "="*60)
        print("🔍 R INTEGRATION STATUS")
        print("="*60)
        
        # rpy2 status
        status = "✅" if self.rpy2_available else "❌"
        print(f"{status} rpy2: {'Available' if self.rpy2_available else 'Not installed'}")
        
        # R status  
        status = "✅" if self.r_available else "❌"
        print(f"{status} R: {'Available' if self.r_available else 'Not available'}")
        
        # Package status
        print(f"\n📦 R PACKAGES STATUS:")
        print("-" * 25)
        
        if self.available_packages:
            print("✅ Available packages:")
            for pkg, info in self.available_packages.items():
                print(f"   • {pkg}: {info['description']}")
        
        if self.missing_packages:
            print(f"\n❌ Missing packages ({len(self.missing_packages)}):")
            for pkg in self.missing_packages:
                print(f"   • {pkg}")
        
        # Overall status
        print(f"\n🎯 OVERALL STATUS:")
        print("-" * 18)
        
        if self.rpy2_available and self.r_available and len(self.available_packages) > 0:
            print(f"✅ {len(self.available_packages)} R methods available for use")
        elif self.rpy2_available and self.r_available:
            print("⚠️  R infrastructure ready, but no R packages installed")
        else:
            print("❌ R methods not available - installation required")
        
        # Installation guide if needed
        if not self.rpy2_available or not self.r_available or self.missing_packages:
            print(self.get_installation_guide())
    
    def get_available_methods(self) -> Dict:
        """Get dictionary of available R methods for registration"""
        
        if not self.rpy2_available or not self.r_available:
            return {}
        
        # Import method classes only if R is available
        try:
            from .r_methods import (
                ALDEx2Method, ANCOMBCMethod, DESeq2Method, 
                EdgeRMethod, MetagenomeSeqMethod
            )
            
            method_mapping = {
                'ALDEx2': ALDEx2Method,
                'ANCOMBC': ANCOMBCMethod,
                'DESeq2': DESeq2Method,
                'edgeR': EdgeRMethod,
                'metagenomeSeq': MetagenomeSeqMethod
            }
            
            available_methods = {}
            for package, method_class in method_mapping.items():
                if package in self.available_packages:
                    try:
                        # Test instantiation
                        method_instance = method_class()
                        available_methods[method_instance.name()] = method_class
                        logger.info(f"✅ Method {method_instance.name()} ready for registration")
                    except Exception as e:
                        logger.warning(f"❌ Method {package} failed instantiation: {e}")
            
            return available_methods
            
        except ImportError as e:
            logger.warning(f"Failed to import R method classes: {e}")
            return {}
    
    def create_installation_script(self, filename: str = "install_r_dependencies.sh"):
        """Create a shell script for easy R dependency installation"""
        
        script_content = f"""#!/bin/bash
# DAAadvisor R Dependencies Installation Script
# Generated by R Integration Manager

set -e  # Exit on any error

echo "🚀 Installing DAAadvisor R Dependencies"
echo "======================================"

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 is required but not installed"
    exit 1
fi

# Install rpy2
echo "📦 Installing rpy2..."
python3 -m pip install rpy2

# Check if R is available
if ! command -v R &> /dev/null; then
    echo "📊 R not found. Installing R..."
    
    # macOS
    if [[ "$OSTYPE" == "darwin"* ]]; then
        if command -v brew &> /dev/null; then
            brew install r
        else
            echo "❌ Please install Homebrew or download R from https://cran.r-project.org/"
            exit 1
        fi
    
    # Linux (Ubuntu/Debian)
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        sudo apt-get update
        sudo apt-get install -y r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev
    
    else
        echo "❌ Unsupported OS. Please install R manually from https://cran.r-project.org/"
        exit 1
    fi
fi

# Install R packages
echo "🧬 Installing R Bioconductor packages..."
R --vanilla <<EOF
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org/")

BiocManager::install(c("ALDEx2", "ANCOMBC", "DESeq2", "edgeR", "metagenomeSeq"))
EOF

echo "✅ Installation complete!"
echo "🧪 Testing installation..."

python3 -c "from daa_advisor.r_integration_fix import RIntegrationManager; RIntegrationManager().print_status()"

echo "🎉 R methods are now ready for use with DAAadvisor!"
"""
        
        with open(filename, 'w') as f:
            f.write(script_content)
        
        # Make executable
        import stat
        import os
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IEXEC)
        
        print(f"✅ Installation script created: {filename}")
        print(f"📝 Run with: bash {filename}")


def quick_r_check():
    """Quick function to check R integration status"""
    manager = RIntegrationManager()
    
    if manager.rpy2_available and manager.r_available and manager.available_packages:
        print(f"✅ R integration ready: {len(manager.available_packages)} methods available")
        return True
    else:
        print("❌ R integration not ready. Run full status check:")
        print("python -c \"from daa_advisor.r_integration_fix import RIntegrationManager; RIntegrationManager().print_status()\"")
        return False


if __name__ == "__main__":
    # Run comprehensive status check
    manager = RIntegrationManager()
    manager.print_status()
    
    # Optionally create installation script
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--create-installer":
        manager.create_installation_script()