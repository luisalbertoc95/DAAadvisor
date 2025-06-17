#!/usr/bin/env python3
"""
Data generators for testing different microbiome data types
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, List, Optional
import logging
from scipy import stats

logger = logging.getLogger(__name__)


class MicrobiomeDataGenerator:
    """Generate realistic microbiome data for testing"""
    
    def __init__(self, random_seed: int = 42):
        self.random_seed = random_seed
        np.random.seed(random_seed)
    
    def generate_asv_data(self, 
                         n_samples: int = 100,
                         n_features: int = 200,
                         n_differential: int = 20,
                         sparsity: float = 0.8,
                         effect_size: float = 2.0,
                         noise_level: float = 0.1) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
        """
        Generate realistic 16S/ASV microbiome data
        
        Parameters:
        -----------
        n_samples : int
            Number of samples
        n_features : int
            Number of ASVs/features
        n_differential : int
            Number of truly differential features
        sparsity : float
            Proportion of zeros (0-1)
        effect_size : float
            Magnitude of differential abundance effect
        noise_level : float
            Amount of biological noise
        
        Returns:
        --------
        Tuple of (count_table, metadata, differential_features)
        """
        
        np.random.seed(self.random_seed)
        
        # Generate base abundance distribution (log-normal)
        base_abundances = np.random.lognormal(mean=2, sigma=2, size=n_features)
        base_abundances = base_abundances / base_abundances.sum() * 10000  # Scale to reasonable depth
        
        # Create samples with compositional variation
        count_matrix = np.zeros((n_samples, n_features))
        
        # Generate library sizes with realistic variation
        library_sizes = np.random.lognormal(mean=9, sigma=0.5, size=n_samples)  # ~8000 reads average
        
        for i in range(n_samples):
            # Add sample-specific variation
            sample_abundances = base_abundances * np.random.lognormal(mean=0, sigma=noise_level, size=n_features)
            
            # Multinomial sampling to create realistic count data
            sample_abundances = sample_abundances / sample_abundances.sum()
            counts = np.random.multinomial(int(library_sizes[i]), sample_abundances)
            count_matrix[i, :] = counts
        
        # Add sparsity by setting random elements to zero
        sparsity_mask = np.random.random((n_samples, n_features)) < sparsity
        count_matrix[sparsity_mask] = 0
        
        # Create metadata with balanced groups
        metadata = pd.DataFrame({
            'condition': ['Control'] * (n_samples // 2) + ['Treatment'] * (n_samples - n_samples // 2),
            'batch': np.random.choice(['A', 'B', 'C'], n_samples),
            'age': np.random.randint(20, 70, n_samples),
            'bmi': np.random.normal(25, 5, n_samples)
        }, index=[f"Sample_{i+1}" for i in range(n_samples)])
        
        # Add differential abundance to selected features
        differential_features = np.random.choice(n_features, n_differential, replace=False)
        treatment_mask = metadata['condition'] == 'Treatment'
        
        for feat_idx in differential_features:
            # Increase abundance in treatment group
            multiplier = np.random.uniform(effect_size, effect_size * 2)
            count_matrix[treatment_mask, feat_idx] *= multiplier
        
        # Convert to DataFrame
        count_table = pd.DataFrame(
            count_matrix.astype(int),
            index=metadata.index,
            columns=[f"ASV_{i+1}" for i in range(n_features)]
        )
        
        differential_feature_names = [f"ASV_{i+1}" for i in differential_features]
        
        logger.info(f"Generated ASV data: {n_samples} samples × {n_features} features")
        logger.info(f"Sparsity: {(count_table == 0).sum().sum() / (n_samples * n_features):.1%}")
        logger.info(f"Differential features: {len(differential_feature_names)}")
        
        return count_table, metadata, differential_feature_names
    
    def generate_gene_data(self,
                          n_samples: int = 80,
                          n_features: int = 500,
                          n_differential: int = 50,
                          sparsity: float = 0.4,
                          effect_size: float = 3.0,
                          noise_level: float = 0.2) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
        """
        Generate realistic gene/functional microbiome data (less sparse than ASV)
        """
        
        np.random.seed(self.random_seed + 1)
        
        # Gene data typically has higher counts and less sparsity
        base_abundances = np.random.gamma(shape=2, scale=1000, size=n_features)
        
        count_matrix = np.zeros((n_samples, n_features))
        library_sizes = np.random.lognormal(mean=11, sigma=0.3, size=n_samples)  # Higher depth for genes
        
        for i in range(n_samples):
            sample_abundances = base_abundances * np.random.lognormal(mean=0, sigma=noise_level, size=n_features)
            sample_abundances = sample_abundances / sample_abundances.sum()
            counts = np.random.multinomial(int(library_sizes[i]), sample_abundances)
            count_matrix[i, :] = counts
        
        # Lower sparsity for gene data
        sparsity_mask = np.random.random((n_samples, n_features)) < sparsity
        count_matrix[sparsity_mask] = 0
        
        # Create metadata
        metadata = pd.DataFrame({
            'condition': ['Healthy'] * (n_samples // 2) + ['Disease'] * (n_samples - n_samples // 2),
            'batch': np.random.choice(['Batch1', 'Batch2'], n_samples),
            'sequencing_depth': library_sizes,
            'patient_id': [f"P{i:03d}" for i in range(n_samples)]
        }, index=[f"Sample_{i+1}" for i in range(n_samples)])
        
        # Add differential features
        differential_features = np.random.choice(n_features, n_differential, replace=False)
        disease_mask = metadata['condition'] == 'Disease'
        
        for feat_idx in differential_features:
            multiplier = np.random.uniform(effect_size, effect_size * 1.5)
            if np.random.random() < 0.5:  # Half up-regulated, half down-regulated
                count_matrix[disease_mask, feat_idx] *= multiplier
            else:
                count_matrix[disease_mask, feat_idx] /= multiplier
        
        count_table = pd.DataFrame(
            count_matrix.astype(int),
            index=metadata.index,
            columns=[f"Gene_{i+1}" for i in range(n_features)]
        )
        
        differential_feature_names = [f"Gene_{i+1}" for i in differential_features]
        
        logger.info(f"Generated gene data: {n_samples} samples × {n_features} features")
        logger.info(f"Sparsity: {(count_table == 0).sum().sum() / (n_samples * n_features):.1%}")
        
        return count_table, metadata, differential_feature_names
    
    def generate_viral_data(self,
                           n_samples: int = 60,
                           n_features: int = 150,
                           n_differential: int = 10,
                           sparsity: float = 0.95,
                           effect_size: float = 5.0,
                           noise_level: float = 0.3) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
        """
        Generate realistic viral microbiome data (very sparse, presence/absence patterns)
        """
        
        np.random.seed(self.random_seed + 2)
        
        # Viral data is extremely sparse with occasional high abundance
        count_matrix = np.zeros((n_samples, n_features))
        
        # Most features have very low baseline abundance
        for i in range(n_samples):
            # Only a few features are "active" per sample
            n_active_features = np.random.poisson(5)  # Very few active viruses per sample
            active_features = np.random.choice(n_features, min(n_active_features, n_features), replace=False)
            
            for feat_idx in active_features:
                # High counts when present
                count_matrix[i, feat_idx] = np.random.negative_binomial(n=2, p=0.1)
        
        # Extreme sparsity
        sparsity_mask = np.random.random((n_samples, n_features)) < sparsity
        count_matrix[sparsity_mask] = 0
        
        # Create metadata
        metadata = pd.DataFrame({
            'condition': ['NonViral'] * (n_samples // 2) + ['Viral'] * (n_samples - n_samples // 2),
            'sample_type': np.random.choice(['Stool', 'Saliva', 'Throat'], n_samples),
            'collection_date': pd.date_range('2023-01-01', periods=n_samples, freq='D'),
            'viral_load': np.random.exponential(1000, n_samples)
        }, index=[f"Sample_{i+1}" for i in range(n_samples)])
        
        # Add differential features (viral blooms)
        differential_features = np.random.choice(n_features, n_differential, replace=False)
        viral_mask = metadata['condition'] == 'Viral'
        
        for feat_idx in differential_features:
            # Viral blooms: very high counts in viral condition
            viral_samples = np.random.choice(np.where(viral_mask)[0], 
                                           size=max(1, int(viral_mask.sum() * 0.3)), 
                                           replace=False)
            count_matrix[viral_samples, feat_idx] += np.random.negative_binomial(n=1, p=0.05, size=len(viral_samples))
        
        count_table = pd.DataFrame(
            count_matrix.astype(int),
            index=metadata.index,
            columns=[f"Virus_{i+1}" for i in range(n_features)]
        )
        
        differential_feature_names = [f"Virus_{i+1}" for i in differential_features]
        
        logger.info(f"Generated viral data: {n_samples} samples × {n_features} features")
        logger.info(f"Sparsity: {(count_table == 0).sum().sum() / (n_samples * n_features):.1%}")
        
        return count_table, metadata, differential_feature_names
    
    def generate_challenging_dataset(self,
                                   data_type: str = "asv",
                                   challenge_type: str = "unbalanced",
                                   **kwargs) -> Tuple[pd.DataFrame, pd.DataFrame, List[str]]:
        """
        Generate challenging datasets for method benchmarking
        
        Parameters:
        -----------
        data_type : str
            Type of data ('asv', 'gene', 'viral')
        challenge_type : str
            Type of challenge ('unbalanced', 'small_effect', 'confounded', 'small_sample')
        """
        
        if challenge_type == "unbalanced":
            # Highly unbalanced groups
            n_samples = kwargs.get('n_samples', 60)
            n_control = int(n_samples * 0.8)  # 80% control, 20% treatment
            kwargs['n_samples'] = n_samples
            
        elif challenge_type == "small_effect":
            # Very small effect sizes
            kwargs['effect_size'] = 1.2
            kwargs['n_differential'] = kwargs.get('n_differential', 5)
            
        elif challenge_type == "confounded":
            # Batch effects confounded with treatment
            kwargs['noise_level'] = 0.5
            
        elif challenge_type == "small_sample":
            # Very small sample size
            kwargs['n_samples'] = 20
            kwargs['n_features'] = kwargs.get('n_features', 100)
        
        # Generate base data
        if data_type == "asv":
            count_table, metadata, diff_features = self.generate_asv_data(**kwargs)
        elif data_type == "gene":
            count_table, metadata, diff_features = self.generate_gene_data(**kwargs)
        elif data_type == "viral":
            count_table, metadata, diff_features = self.generate_viral_data(**kwargs)
        else:
            raise ValueError(f"Unknown data type: {data_type}")
        
        # Apply specific challenges
        if challenge_type == "unbalanced":
            # Create unbalanced groups
            n_control = int(len(metadata) * 0.8)
            conditions = ['Control'] * n_control + ['Treatment'] * (len(metadata) - n_control)
            np.random.shuffle(conditions)
            metadata['condition'] = conditions
            
        elif challenge_type == "confounded":
            # Make batch effects confounded with treatment
            metadata['batch'] = ['A'] * (len(metadata) // 2) + ['B'] * (len(metadata) - len(metadata) // 2)
            
        logger.info(f"Generated challenging dataset: {challenge_type} {data_type} data")
        
        return count_table, metadata, diff_features


def create_benchmark_datasets() -> Dict[str, Tuple[pd.DataFrame, pd.DataFrame, List[str]]]:
    """Create a comprehensive set of benchmark datasets for testing"""
    
    generator = MicrobiomeDataGenerator()
    datasets = {}
    
    # Standard datasets
    datasets['asv_standard'] = generator.generate_asv_data(
        n_samples=100, n_features=200, sparsity=0.8, effect_size=2.0
    )
    
    datasets['gene_standard'] = generator.generate_gene_data(
        n_samples=80, n_features=500, sparsity=0.4, effect_size=3.0
    )
    
    datasets['viral_standard'] = generator.generate_viral_data(
        n_samples=60, n_features=150, sparsity=0.95, effect_size=5.0
    )
    
    # Challenging datasets
    challenge_types = ['unbalanced', 'small_effect', 'confounded', 'small_sample']
    data_types = ['asv', 'gene', 'viral']
    
    for data_type in data_types:
        for challenge in challenge_types:
            key = f"{data_type}_{challenge}"
            datasets[key] = generator.generate_challenging_dataset(
                data_type=data_type, 
                challenge_type=challenge
            )
    
    logger.info(f"Created {len(datasets)} benchmark datasets")
    return datasets