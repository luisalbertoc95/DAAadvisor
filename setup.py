#!/usr/bin/env python3

from setuptools import setup, find_packages
import os

# Read README file
readme_path = os.path.join(os.path.dirname(__file__), 'README.md')
if os.path.exists(readme_path):
    with open(readme_path, 'r', encoding='utf-8') as f:
        long_description = f.read()
else:
    long_description = "DAAadvisor: Differential Abundance Analysis Tool"

# Read requirements
def read_requirements():
    requirements = []
    if os.path.exists('requirements.txt'):
        with open('requirements.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    # Remove version constraints for setup.py
                    package = line.split('>=')[0].split('==')[0].split('<=')[0]
                    requirements.append(package)
    return requirements

setup(
    name="daaadvisor",
    version="0.1.0",
    author="DAAadvisor Team",
    author_email="support@daaadvisor.org",
    description="Intelligent differential abundance analysis tool for microbiome data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/daaadvisor",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        'dev': [
            'pytest>=6.2.0',
            'pytest-cov>=2.12.0',
            'black>=21.0.0',
            'flake8>=3.9.0',
        ],
        'r': [
            'rpy2>=3.5.0',
        ],
        'viz': [
            'plotly>=5.0.0',
            'matplotlib>=3.4.0',
            'seaborn>=0.11.0',
        ]
    },
    entry_points={
        'console_scripts': [
            'daaadvisor=daa_advisor.cli:main',
        ],
    },
    include_package_data=True,
    keywords="microbiome differential abundance statistics bioinformatics",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/daaadvisor/issues",
        "Source": "https://github.com/yourusername/daaadvisor",
        "Documentation": "https://daaadvisor.readthedocs.io/",
    },
)