#!/usr/bin/env python3
"""
Test script for curatedMetagenomicData downloader
"""

import logging
from daa_advisor.curated_data_downloader import CuratedDataDownloader

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def test_curated_downloader():
    """Test the curated data downloader"""
    
    print("🧪 Testing curatedMetagenomicData Downloader")
    print("=" * 50)
    
    # Initialize downloader
    downloader = CuratedDataDownloader(output_dir="test_curated_data")
    
    # Test R environment check
    print("\n🔍 Testing R environment check...")
    r_available = downloader.check_r_environment()
    
    if r_available:
        print("✅ R environment ready")
        
        # Test package installation
        print("\n📦 Testing package installation...")
        packages_installed = downloader.install_required_packages()
        
        if packages_installed:
            print("✅ R packages ready")
            
            # Show available conditions
            print(f"\n📋 Available conditions:")
            for condition, info in downloader.available_conditions.items():
                print(f"   • {condition}: {info['description']}")
                print(f"     Expected features: {info['expected_features']}")
            
            print(f"\n💡 To download real data:")
            print(f"   python run_cross_validation_benchmark.py --max-conditions 1")
            
        else:
            print("❌ R packages installation failed")
            print("💡 Try manually in R:")
            print("   install.packages('BiocManager')")
            print("   BiocManager::install('curatedMetagenomicData')")
    
    else:
        print("❌ R not available")
        print("💡 Install R first:")
        print("   # macOS: brew install r")
        print("   # Linux: sudo apt-get install r-base")
        print("   # Or use conda: conda install -c conda-forge r-base")
    
    # Show synthetic data alternative
    print(f"\n🎭 Realistic synthetic data available:")
    
    from pathlib import Path
    synthetic_dir = Path("realistic_demo_data")
    
    if synthetic_dir.exists():
        files = list(synthetic_dir.glob("*.csv"))
        print(f"✅ Found {len(files)} synthetic data files:")
        for file in files:
            print(f"   📄 {file.name}")
        
        print(f"\n🚀 Ready for immediate benchmarking:")
        print(f"   python run_publication_benchmark.py --full")
    else:
        print(f"⚠️ Synthetic data not found")
        print(f"💡 Generate with: python download_real_data.py")

if __name__ == "__main__":
    test_curated_downloader()