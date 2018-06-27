"""
This file contains all testing functions for the main program execution.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'BSD'

import os
import random
import unittest

import biom
from massoc.scripts.batch import Batch
from massoc.scripts.netwrap import Nets

import massoc
from massoc.scripts.main import combine_data, run_networks, run_jobs, run_massoc, run_parallel, get_joblist

random.seed(7)

# replace with proper test file!
# looks like there is an issue with agglomerated files
testloc = list()
testloc.append(os.path.dirname(massoc.__file__)[:-6] + 'data\\coral_ID10798_42708.biom')
testloc = [x.replace('\\', '/') for x in testloc]
testbiom = {"test": biom.load_table(testloc[0])}

inputs = {'biom_file': None,
          'cluster': None,
          'otu_meta': None,
          'prefix': None,
          'sample_data': None,
          'split': None,
          'tax_table': None,
          'fp': [testloc[0][:-25]],
          'name': ['test'],
          'otu_table': None,
          'tools': ['spiec-easi'],
          'spiec': None,
          'conet': None,
          'spar_pval': None,
          'spar_boot': None,
          'levels': ['otu', 'order'],
          'prev': ['20'],
          'cores': ['4']}
batch = Batch(testbiom, inputs)
netbatch = Nets(batch)
netbatch.write_bioms()

class TestMain(unittest.TestCase):
    """"
    Tests whether the main pipeline functions
    are able to execute the appropriate commands.
    Note that the current functions only check whether network inference runs at all;
    we should have a checked dataset with known outputs to make sure inference is run correctly.
    """

    def test_combine_data(self):
        """
        Checks whether combine_data returns
        a batch object if inputs are supplied.
        """
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None}
        bioms = combine_data(inputs)
        self.assertEqual(len(bioms.otu), 1)


    def test_combine_data_error(self):
        """
        Checks whether combine_data returns an error
        if the split variable is not in the medata.
        """
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': None,
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None}
        with self.assertRaises(ValueError):
            bioms = combine_data(inputs)

    def test_run_networks(self):
        """Checks if the run_networks function
        works appropriately."""
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['spiec-easi'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['otu'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None}
        batch = Batch(testbiom, inputs)
        testnets = Nets(batch)
        testnets = run_networks(testnets)
        self.assertEqual(len(testnets.networks), 1)

    def test_run_massoc(self):
        """Checks if the run_massoc function can run."""
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['spiec-easi'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['otu'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None,
                  'min': None,
                  'rar': None}
        run_massoc(inputs)

    def test_get_joblist(self):
        """
        Checks whether the joblist function
        returns a joblist in the appropriate format:
        list of dicts with each only 1 key.
        """
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-13] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-13] + 'otu_otus.txt')],
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': ['somefile.txt'],
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family', 'class'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None}
        batch = Batch(testbiom, inputs)
        netbatch = Nets(batch)
        jobs = get_joblist(netbatch)
        self.assertEqual(len(jobs), 6)

    def test_run_jobs(self):
        """
        Checks whether run_jobs really returns only 1 network.
        """
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['spiec-easi'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None}
        batch = Batch(testbiom, inputs)
        netbatch = Nets(batch)
        jobs = get_joblist(netbatch)
        network = run_jobs(netbatch, jobs[0])
        self.assertEqual(len(network), 1)

    def test_run_parallel(self):
        """Checks if the run_parallel function works without raising an error."""
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': [(testloc[0][:-17] + 'otu_tax.txt')],
                  'fp': [testloc[0][:-25]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'min': None,
                  'name': ['test'],
                  'cores': None,
                  'rar': None}
        batch = Batch(testbiom, inputs)
        netbatch = Nets(batch)
        netbatch = run_parallel(netbatch)


if __name__ == '__main__':
    unittest.main()
