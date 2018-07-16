"""
This file contains all testing functions for the main program execution.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

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
          'levels': ['family'],
          'prev': ['20'],
          'cores': ['4']}
batch = Batch(testbiom, inputs)
netbatch = Nets(batch)
netbatch.write_bioms()
netbatch = run_networks(netbatch)

class TestMain(unittest.TestCase):
    """"
    Tests whether the main pipeline functions
    are able to execute the appropriate commands.
    Note that the current functions only check whether network inference runs at all;
    we should have a checked dataset with known outputs to make sure inference is run correctly.
    """

    def test_create_base(self):
        """
        Checks whether an object generated through netwrap
        can be converted into a neo4j database.
        """

if __name__ == '__main__':
    unittest.main()