#!/usr/bin/env python3
"""
Comprehensive test suite for DAAadvisor
"""

import unittest
import numpy as np
import pandas as pd
import tempfile
import shutil
from pathlib import Path
import logging

# Import DAAadvisor components
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from daa_advisor import DifferentialAbundanceTool, DataProfiler, MethodSelector
from daa_advisor.data_generators import MicrobiomeDataGenerator, create_benchmark_datasets
from daa_advisor.methods import MethodRegistry
from daa_advisor.benchmarking import MethodBenchmark

# Suppress warnings for cleaner test output
logging.getLogger().setLevel(logging.ERROR)


class TestDataProfiler(unittest.TestCase):
    """Test data profiling functionality"""
    
    def setUp(self):
        self.generator = MicrobiomeDataGenerator(random_seed=42)
        self.profiler = DataProfiler()
    
    def test_asv_data_profiling(self):
        """Test profiling of ASV data"""
        count_table, metadata, _ = self.generator.generate_asv_data(
            n_samples=50, n_features=100, sparsity=0.8
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        
        # Test basic properties
        self.assertEqual(profile.data_type, 'asv')
        self.assertEqual(profile.n_samples, 50)
        self.assertEqual(profile.n_features, 100)
        self.assertGreater(profile.sparsity, 0.7)  # Should be high sparsity
        self.assertGreater(profile.zero_inflation, 0.5)
        
        # Test group sizes
        self.assertEqual(len(profile.group_sizes), 2)
        self.assertIn('Control', profile.group_sizes)
        self.assertIn('Treatment', profile.group_sizes)
    
    def test_gene_data_profiling(self):
        """Test profiling of gene data"""
        count_table, metadata, _ = self.generator.generate_gene_data(
            n_samples=40, n_features=200, sparsity=0.4
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        
        self.assertEqual(profile.data_type, 'gene')
        self.assertEqual(profile.n_samples, 40)
        self.assertEqual(profile.n_features, 200)
        self.assertLess(profile.sparsity, 0.6)  # Should be lower sparsity than ASV
    
    def test_viral_data_profiling(self):
        """Test profiling of viral data"""
        count_table, metadata, _ = self.generator.generate_viral_data(
            n_samples=30, n_features=150, sparsity=0.95
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        
        self.assertEqual(profile.data_type, 'viral')
        self.assertEqual(profile.n_samples, 30)
        self.assertEqual(profile.n_features, 150)
        self.assertGreater(profile.sparsity, 0.9)  # Should be very high sparsity
    
    def test_data_type_detection(self):
        """Test automatic data type detection"""
        # Test ASV detection
        asv_features = [f"ASV_{i}" for i in range(50)]
        asv_data = pd.DataFrame(
            np.random.poisson(2, (20, 50)),
            columns=asv_features,
            index=[f"Sample_{i}" for i in range(20)]
        )
        metadata = pd.DataFrame({
            'condition': ['A'] * 10 + ['B'] * 10
        }, index=asv_data.index)
        
        profile = self.profiler.profile_data(asv_data, metadata)
        self.assertEqual(profile.data_type, 'asv')
        
        # Test gene detection
        gene_features = [f"gene_{i}" for i in range(50)]
        gene_data = asv_data.copy()
        gene_data.columns = gene_features
        
        profile = self.profiler.profile_data(gene_data, metadata)
        self.assertEqual(profile.data_type, 'gene')


class TestMethodSelector(unittest.TestCase):
    """Test method selection logic"""
    
    def setUp(self):
        self.generator = MicrobiomeDataGenerator(random_seed=42)
        self.profiler = DataProfiler()
        self.selector = MethodSelector()
    
    def test_asv_method_selection(self):
        """Test method selection for ASV data"""
        count_table, metadata, _ = self.generator.generate_asv_data(
            n_samples=100, n_features=200, sparsity=0.8
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        recommendations = self.selector.recommend_methods(profile)
        
        # ALDEx2 should be recommended for compositional ASV data
        self.assertEqual(recommendations.primary_method, 'aldex2')
        self.assertGreater(recommendations.confidence, 0.8)
        self.assertIn('ancom-bc', recommendations.secondary_methods)
    
    def test_gene_method_selection(self):
        """Test method selection for gene data"""
        count_table, metadata, _ = self.generator.generate_gene_data(
            n_samples=80, n_features=300, sparsity=0.4
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        recommendations = self.selector.recommend_methods(profile)
        
        # Should recommend methods good for gene data
        self.assertIn(recommendations.primary_method, ['deseq2', 'linda', 'ancom-bc'])
        self.assertGreater(recommendations.confidence, 0.5)
    
    def test_viral_method_selection(self):
        """Test method selection for viral data"""
        count_table, metadata, _ = self.generator.generate_viral_data(
            n_samples=60, n_features=150, sparsity=0.95
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        recommendations = self.selector.recommend_methods(profile)
        
        # Should recommend methods that handle high sparsity
        self.assertIn(recommendations.primary_method, ['zicoseq', 'aldex2', 'ancom-bc'])
    
    def test_small_sample_size_handling(self):
        """Test method selection with small sample sizes"""
        count_table, metadata, _ = self.generator.generate_asv_data(
            n_samples=12, n_features=100, sparsity=0.7
        )
        
        profile = self.profiler.profile_data(count_table, metadata)
        recommendations = self.selector.recommend_methods(profile)
        
        # Should recommend methods that work with small samples
        self.assertNotEqual(recommendations.primary_method, 'maaslin3')  # Requires more samples
        self.assertGreater(len(recommendations.secondary_methods), 0)


class TestMethodRegistry(unittest.TestCase):
    """Test method registry functionality"""
    
    def setUp(self):
        self.registry = MethodRegistry()
    
    def test_method_availability(self):
        """Test that basic methods are available"""
        available_methods = self.registry.list_methods()
        
        # Wilcoxon should always be available
        self.assertIn('wilcoxon', available_methods)
        self.assertTrue(self.registry.has_method('wilcoxon'))
        
        # Should have fallback method
        fallback = self.registry.get_fallback_method()
        self.assertEqual(fallback.name(), 'wilcoxon')
    
    def test_method_info(self):
        """Test method information retrieval"""
        wilcoxon_info = self.registry.get_method_info('wilcoxon')
        
        self.assertIn('name', wilcoxon_info)
        self.assertIn('citation', wilcoxon_info)
        self.assertEqual(wilcoxon_info['name'], 'wilcoxon')
    
    def test_invalid_method(self):
        """Test handling of invalid method names"""
        self.assertFalse(self.registry.has_method('nonexistent_method'))
        
        with self.assertRaises(ValueError):
            self.registry.get_method('nonexistent_method')


class TestDifferentialAbundanceTool(unittest.TestCase):
    """Test main analysis tool"""
    
    def setUp(self):
        self.generator = MicrobiomeDataGenerator(random_seed=42)
        self.tool = DifferentialAbundanceTool()
    
    def test_basic_analysis(self):
        """Test basic differential abundance analysis"""
        count_table, metadata, true_features = self.generator.generate_asv_data(
            n_samples=50, n_features=100, n_differential=10, effect_size=3.0
        )
        
        results = self.tool.analyze(
            count_table=count_table,
            metadata=metadata,
            use_consensus=False
        )
        
        # Check results structure
        self.assertIn('profile', results)
        self.assertIn('recommendations', results)
        self.assertIn('analyses', results)
        
        # Check that analysis was performed
        self.assertGreater(len(results['analyses']), 0)
        
        # Check that significant features are found
        analysis_results = list(results['analyses'].values())[0]
        self.assertIn('pvalue', analysis_results.columns)
        self.assertIn('padj', analysis_results.columns)
        self.assertIn('feature', analysis_results.columns)
    
    def test_consensus_analysis(self):
        """Test consensus analysis with multiple methods"""
        count_table, metadata, _ = self.generator.generate_asv_data(
            n_samples=60, n_features=80, n_differential=15
        )
        
        results = self.tool.analyze(
            count_table=count_table,
            metadata=metadata,
            use_consensus=True
        )
        
        # Should have consensus results if multiple methods were run
        if len(results['analyses']) > 1:
            self.assertIn('consensus', results)
            consensus = results['consensus']
            self.assertIn('n_significant', consensus.columns)
            self.assertIn('consensus_significant', consensus.columns)
    
    def test_get_significant_features(self):
        """Test significant feature extraction"""
        count_table, metadata, true_features = self.generator.generate_asv_data(
            n_samples=40, n_features=80, n_differential=8, effect_size=4.0
        )
        
        results = self.tool.analyze(count_table, metadata)
        significant = self.tool.get_significant_features(alpha=0.05)
        
        self.assertIsInstance(significant, pd.DataFrame)
        self.assertIn('feature', significant.columns)
        self.assertIn('pvalue', significant.columns)
        
        # All features should be significant
        if len(significant) > 0:
            self.assertTrue((significant['padj'] < 0.05).all())


class TestDataGeneration(unittest.TestCase):
    """Test data generation functionality"""
    
    def setUp(self):
        self.generator = MicrobiomeDataGenerator(random_seed=42)
    
    def test_asv_data_generation(self):
        """Test ASV data generation"""
        count_table, metadata, diff_features = self.generator.generate_asv_data(
            n_samples=30, n_features=50, n_differential=5
        )
        
        # Check dimensions
        self.assertEqual(count_table.shape, (30, 50))
        self.assertEqual(metadata.shape[0], 30)
        self.assertEqual(len(diff_features), 5)
        
        # Check data types
        self.assertTrue(all(count_table.dtypes == 'int64'))
        self.assertTrue((count_table.values >= 0).all())
        
        # Check sparsity
        sparsity = (count_table == 0).sum().sum() / (30 * 50)
        self.assertGreater(sparsity, 0.5)  # Should be reasonably sparse
    
    def test_gene_data_generation(self):
        """Test gene data generation"""
        count_table, metadata, diff_features = self.generator.generate_gene_data(
            n_samples=25, n_features=100, n_differential=10
        )
        
        self.assertEqual(count_table.shape, (25, 100))
        self.assertEqual(len(diff_features), 10)
        
        # Gene data should be less sparse than ASV data
        sparsity = (count_table == 0).sum().sum() / (25 * 100)
        self.assertLess(sparsity, 0.7)
    
    def test_viral_data_generation(self):
        """Test viral data generation"""
        count_table, metadata, diff_features = self.generator.generate_viral_data(
            n_samples=20, n_features=80, n_differential=8
        )
        
        self.assertEqual(count_table.shape, (20, 80))
        self.assertEqual(len(diff_features), 8)
        
        # Viral data should be very sparse
        sparsity = (count_table == 0).sum().sum() / (20 * 80)
        self.assertGreater(sparsity, 0.8)
    
    def test_challenging_datasets(self):
        """Test challenging dataset generation"""
        challenges = ['unbalanced', 'small_effect', 'confounded', 'small_sample']
        
        for challenge in challenges:
            count_table, metadata, diff_features = self.generator.generate_challenging_dataset(
                data_type='asv', challenge_type=challenge, n_samples=30, n_features=50
            )
            
            self.assertIsInstance(count_table, pd.DataFrame)
            self.assertIsInstance(metadata, pd.DataFrame)
            self.assertIsInstance(diff_features, list)
            
            if challenge == 'unbalanced':
                # Check that groups are unbalanced
                group_sizes = metadata['condition'].value_counts()
                self.assertNotEqual(group_sizes.iloc[0], group_sizes.iloc[1])
            
            elif challenge == 'small_sample':
                # Check small sample size
                self.assertLessEqual(len(count_table), 25)


class TestBenchmarkFramework(unittest.TestCase):
    """Test benchmarking functionality"""
    
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.benchmark = MethodBenchmark(self.temp_dir)
    
    def tearDown(self):
        shutil.rmtree(self.temp_dir)
    
    def test_benchmark_dataset_creation(self):
        """Test benchmark dataset creation"""
        datasets = create_benchmark_datasets()
        
        # Should have multiple datasets
        self.assertGreater(len(datasets), 5)
        
        # Check standard datasets
        self.assertIn('asv_standard', datasets)
        self.assertIn('gene_standard', datasets)
        self.assertIn('viral_standard', datasets)
        
        # Each dataset should have the correct structure
        for name, (count_table, metadata, diff_features) in datasets.items():
            self.assertIsInstance(count_table, pd.DataFrame)
            self.assertIsInstance(metadata, pd.DataFrame)
            self.assertIsInstance(diff_features, list)
            self.assertEqual(len(count_table), len(metadata))
    
    def test_single_method_run(self):
        """Test running a single method"""
        count_table, metadata, true_features = MicrobiomeDataGenerator().generate_asv_data(
            n_samples=20, n_features=30, n_differential=3
        )
        
        result = self.benchmark._run_single_method(
            'wilcoxon', count_table, metadata, true_features
        )
        
        self.assertIn('runtime', result)
        self.assertIn('failed', result)
        self.assertFalse(result['failed'])
        self.assertIn('significant_features', result)
        self.assertIsInstance(result['significant_features'], list)
    
    def test_performance_metrics_calculation(self):
        """Test performance metrics calculation"""
        # Mock method results
        true_features = ['Feature_1', 'Feature_2', 'Feature_3']
        method_results = {
            'method1': {
                'significant_features': ['Feature_1', 'Feature_2', 'Feature_4'],
                'runtime': 1.0,
                'failed': False
            },
            'method2': {
                'significant_features': ['Feature_1'],
                'runtime': 0.5,
                'failed': False
            }
        }
        
        metrics = self.benchmark._calculate_performance_metrics(method_results, true_features)
        
        # Check method1 metrics
        self.assertIn('method1', metrics)
        m1 = metrics['method1']
        self.assertEqual(m1['true_positives'], 2)  # Feature_1, Feature_2
        self.assertEqual(m1['false_positives'], 1)  # Feature_4
        self.assertEqual(m1['false_negatives'], 1)  # Feature_3
        self.assertAlmostEqual(m1['precision'], 2/3, places=3)
        self.assertAlmostEqual(m1['recall'], 2/3, places=3)


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete pipeline"""
    
    def test_end_to_end_workflow(self):
        """Test complete end-to-end workflow"""
        # Generate test data
        generator = MicrobiomeDataGenerator(random_seed=42)
        count_table, metadata, true_features = generator.generate_asv_data(
            n_samples=40, n_features=60, n_differential=6, effect_size=3.0
        )
        
        # Run complete analysis
        tool = DifferentialAbundanceTool()
        results = tool.analyze(count_table, metadata, use_consensus=False)
        
        # Check all components worked
        self.assertIn('profile', results)
        self.assertIn('recommendations', results)
        self.assertIn('analyses', results)
        
        # Test result summary
        tool.summarize_results()  # Should not raise exception
        
        # Test significant feature extraction
        significant = tool.get_significant_features(alpha=0.1)  # Relaxed alpha for test
        self.assertIsInstance(significant, pd.DataFrame)
    
    def test_different_data_types_workflow(self):
        """Test workflow with different data types"""
        generator = MicrobiomeDataGenerator(random_seed=42)
        data_types = ['asv', 'gene', 'viral']
        
        for data_type in data_types:
            with self.subTest(data_type=data_type):
                if data_type == 'asv':
                    count_table, metadata, _ = generator.generate_asv_data(n_samples=30, n_features=40, n_differential=5)
                elif data_type == 'gene':
                    count_table, metadata, _ = generator.generate_gene_data(n_samples=30, n_features=40, n_differential=5)
                else:  # viral
                    count_table, metadata, _ = generator.generate_viral_data(n_samples=30, n_features=40, n_differential=5)
                
                tool = DifferentialAbundanceTool()
                results = tool.analyze(count_table, metadata, data_type=data_type)
                
                # Check that data type was correctly identified/set
                self.assertEqual(results['profile'].data_type, data_type)
                
                # Check that method recommendations make sense for data type
                primary_method = results['recommendations'].primary_method
                self.assertIsInstance(primary_method, str)
                self.assertGreater(len(primary_method), 0)


def run_comprehensive_tests():
    """Run all comprehensive tests"""
    
    # Create test suite
    test_classes = [
        TestDataProfiler,
        TestMethodSelector,
        TestMethodRegistry,
        TestDifferentialAbundanceTool,
        TestDataGeneration,
        TestBenchmarkFramework,
        TestIntegration
    ]
    
    suite = unittest.TestSuite()
    
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_comprehensive_tests()
    if success:
        print("\n🎉 All tests passed!")
    else:
        print("\n❌ Some tests failed!")
        exit(1)