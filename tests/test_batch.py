"""
This file contains all testing functions for the Batch object.

Note: not all errors related to agglomeration were caught.
It may be important to expand current test functions to include
agglomerated files.
"""

import unittest
from copy import deepcopy
import os
import biom
import numpy as np
import massoc

from massoc.scripts.batch import Batch

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

tabotu = '[[ 243  567  112   45   2]\n ' \
         '[ 235   56  788  232    1]\n ' \
         '[4545   22    0    1    0]\n ' \
         '[  41   20    2    4    0]]'

tabtax = "[['k__Bacteria' 'p__Firmicutes' 'c__Clostridia' 'o__Clostridiales'\n  " \
         "'f__Clostridiaceae' 'g__Anaerococcus' 's__']\n " \
         "['k__Bacteria' 'p__Bacteroidetes' 'c__Bacteroidia' 'o__Bacteroidales'\n  " \
         "'f__Prevotellaceae' 'g__Prevotella' 's__']\n " \
         "['k__Bacteria' 'p__Proteobacteria' 'c__Alphaproteobacteria'\n  " \
         "'o__Sphingomonadales' 'f__Sphingomonadaceae' 'g__Sphingomonas' 's__']\n " \
         "['k__Bacteria' 'p__Verrucomicrobia' 'c__Verrucomicrobiae'\n  " \
         "'o__Verrucomicrobiales' 'f__Verrucomicrobiaceae' 'g__Luteolibacter' 's__']]"

tabmeta = "[['Australia' 'Hot']\n " \
          "['Antarctica' 'Cold']\n " \
          "['Netherlands' 'Rainy']\n " \
          "['Belgium' 'Rainy']\n " \
          "['Iceland' 'Cold']]"

sample_ids = ['S%d' % i for i in range(1, 6)]
observ_ids = ['O%d' % i for i in range(1, 5)]

testraw = """{
     "id":null,
     "format": "Biological Observation Matrix 1.0.0-dev",
     "format_url": "http://biom-format.org",
     "type": "OTU table",
     "generated_by": "QIIME revision XYZ",
     "date": "2011-12-19T19:00:00",
     "rows":[
        {"id":"GG_OTU_1", "metadata":{"taxonomy":["k__Bacteria", "p__Proteoba\
cteria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriac\
eae", "g__Escherichia", "s__"]}},
        {"id":"GG_OTU_2", "metadata":{"taxonomy":["k__Bacteria", "p__Cyanobact\
eria", "c__Nostocophycideae", "o__Nostocales", "f__Nostocaceae", "g__Dolichosp\
ermum", "s__"]}},
        {"id":"GG_OTU_3", "metadata":{"taxonomy":["k__Archaea", "p__Euryarchae\
ota", "c__Methanomicrobia", "o__Methanosarcinales", "f__Methanosarcinaceae", "\
g__Methanosarcina", "s__"]}},
        {"id":"GG_OTU_4", "metadata":{"taxonomy":["k__Bacteria", "p__Firmicute\
s", "c__Clostridia", "o__Halanaerobiales", "f__Halanaerobiaceae", "g__Halanaer\
obium", "s__Halanaerobiumsaccharolyticum"]}},
        {"id":"GG_OTU_5", "metadata":{"taxonomy":["k__Bacteria", "p__Proteobac\
teria", "c__Gammaproteobacteria", "o__Enterobacteriales", "f__Enterobacteriace\
ae", "g__Escherichia", "s__"]}}
        ],
     "columns":[
        {"id":"Sample1", "metadata":{
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample5", "metadata":{
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample6", "metadata":{
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
        ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 6],
     "data":[[0,2,1],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,5,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [4,1,1],
             [4,2,1]
            ]
    }
"""

testbiom = {'otu': {"test": biom.parse.parse_biom_table(testraw)}}


class TestBatch(unittest.TestCase):
    """Tests Batch methods.
    """

    def test_split_biom(self):
        """Does 'split_biom' correctly split a biom file
        according to sample data properties?"""
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'BODY_SITE',
                  'tax_table': ['tax_bananas.txt'],
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        batch.split_biom()
        self.assertEqual(len(batch.otu), 3)

    def test_prev_filter(self):
        """Does the prevalence filter correctly
        reduce the number of taxa in a table?"""
        inputs = {'biom_file': None,
                  'cluster': 'Affinity',
                  'nclust': ['4'],
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'Rocket Science',
                  'tax_table': ['tax_bananas.txt'],
                  'prev': 40,
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        batch.prev_filter()
        self.assertEqual(batch.otu['test'].shape[0], 4)

    def test_collapse_tax(self):
        """Does the function for collapsing by taxonomy correctly
        add additional dictionaries?"""
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': ['BODY_SITE'],
                  'tax_table': ['tax_bananas.txt'],
                  'levels': ['otu', 'genus'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests'),
                  'name': ['test']}
        batch = Batch(testbiom, inputs)
        batch.collapse_tax()
        self.assertEqual(len(batch.genus), 1)

    def test_normalize_transform(self):
        """Is the transformed batch file different from the original one?"""
        inputs = {'biom_file': None,
                  'cluster': ['Affinity'],
                  'nclust': ['4'],
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'Rocket Science',
                  'tax_table': ['tax_bananas.txt'],
                  'name': ['test'],
                  'fp': os.path.dirname(massoc.__file__)[:-7].replace('\\', '/')}
        batch = Batch(testbiom, inputs)
        clrbatch = batch.normalize_transform(mode="clr")
        self.assertFalse(batch.otu['test'] == clrbatch.otu['test'])

    def test_cluster_bioms_kmeans_split(self):
        """Does 'cluster_bioms.py' correctly cluster
        a biom file and split the file into multiple
        subsets of the data?"""
        inputs = {'biom_file': None,
                  'cluster': 'K-means',
                  'nclust': 2,
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'TRUE',
                  'tax_table': ['tax_bananas.txt'],
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        np.random.seed(8888)
        batch = Batch(deepcopy(testbiom), inputs)
        batch.cluster_biom()
        self.assertEqual(len(batch.otu), 3)

    def test_cluster_bioms_spectral(self):
        """Does 'cluster_bioms.py' correctly cluster
        a biom file and split the file into multiple
        subsets of the data?"""
        inputs = {'biom_file': None,
                  'cluster': 'Spectral',
                  'otu_meta': None,
                  'nclust': 4,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'TRUE',
                  'tax_table': ['tax_bananas.txt'],
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        np.random.seed(8888)
        batch.cluster_biom()
        self.assertEqual(len(batch.otu), 4)

    def test_norm_machine(self):
        """While data is normalized in the machine
        learning function, it should NOT be returned
        as normalized count data."""
        inputs = {'biom_file': None,
                  'cluster': 'K-means',
                  'otu_meta': None,
                  'nclust': 4,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': ['tax_bananas.txt'],
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        batch.cluster_biom()
        self.assertEqual(batch.otu['test']._data[1, 1], testbiom['otu']['test']._data[1, 1])

    def test_rarefy(self):
        """The rarefaction function should
        remove samples below a certain read count and then perform rarefaction.
        """
        inputs = {'biom_file': None,
                  'cluster': 'K-means',
                  'otu_meta': None,
                  'nclust': 4,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': ['tax_bananas.txt'],
                  'rar': 'True',
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        rawsums = batch.otu['test'].sum(axis='sample')
        batch.rarefy()
        newsums = batch.otu['test'].sum(axis='sample')
        self.assertGreater(np.mean(rawsums), np.mean(newsums))

    def test_min(self):
        """The prevalence function should correctly filter taxa with
        mean abundances below the specified threshold."""
        inputs = {'biom_file': None,
                  'cluster': ['K-means'],
                  'otu_meta': None,
                  'nclust': ['4'],
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': ['tax_bananas.txt'],
                  'rar': ['True'],
                  'min': 3,
                  'prev': None,
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        rawsums = batch.otu['test'].sum(axis='observation')
        batch.prev_filter(mode='min')
        newsums = batch.otu['test'].sum(axis='observation')
        self.assertEqual((rawsums[0] + rawsums[4]), newsums[3])

    def test_prev_filter_qual(self):
        """Does the prevalence filter remove the correct taxon?"""
        inputs = {'biom_file': None,
                  'cluster': 'Affinity',
                  'nclust': 4,
                  'otu_meta': None,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': 'Rocket Science',
                  'tax_table': ['tax_bananas.txt'],
                  'prev': 40,
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        rawtable = batch.otu['test'].matrix_data
        batch.prev_filter(mode='prev')
        newtable = batch.otu['test'].matrix_data
        self.assertEqual((rawtable.sum(axis=1)[0] + rawtable.sum(axis=1)[4]), newtable.sum(axis=1)[3])

    def test_rarefy_qual(self):
        """The rarefaction function should
        remove samples below a certain read count and then perform rarefaction.
        Are the lowest values of the table equal to the specified rarefication
        number?
        """
        inputs = {'biom_file': None,
                  'cluster': 'K-means',
                  'otu_meta': None,
                  'nclust': 4,
                  'otu_table': ['otu_bananas.txt'],
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': ['tax_bananas.txt'],
                  'rar': 3,
                  'name': ['test'],
                  'fp': (os.path.dirname(massoc.__file__)[:-6] + 'tests')}
        batch = Batch(deepcopy(testbiom), inputs)
        batch.rarefy()
        newsums = batch.otu['test'].sum(axis='sample')
        self.assertEqual(np.mean(newsums), 3)


if __name__ == '__main__':
    unittest.main()
