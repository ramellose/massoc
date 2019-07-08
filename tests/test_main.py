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
from pathlib import Path
from biom.cli.util import write_biom_table
from subprocess import call
import massoc
from massoc.scripts.main import get_input, run_network, \
    run_neo4j, run_netstats
random.seed(7)

testloc = list()
testloc.append((os.path.dirname(massoc.__file__)[:-7]).replace('\\', '/'))

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

testbiom = {"test": biom.parse.parse_biom_table(testraw)}


class TestMain(unittest.TestCase):
    """"
    Tests whether the main pipeline functions
    are able to execute the appropriate commands.
    Note that the current functions only check whether network inference runs at all;
    we should have a checked dataset with known outputs to make sure inference is run correctly.
    """

    def test_get_input(self):
        """
        Checks whether get_input writes a settings file
        if inputs are supplied.
        """
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': testloc[0],
                  'otu_table': None,
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'conet': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': 20,
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None}
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        get_input(inputs)
        test = Path(inputs['fp'] + "/settings.json")
        self.assertTrue(test.is_file())
        call(("rm " + inputs['biom_file'][0]))
        call(("rm " + inputs['fp'] + "/settings.json"))
        call(("rm " + inputs['fp'] + "/test_family.hdf5"))
        call(("rm " + inputs['fp'] + "/test_otu.hdf5"))

    def test_get_input_error(self):
        """
        Checks whether the get_input pipe reports an error
        if no files are submitted.
        """
        inputs = {'biom_file': None,
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': [testloc[0]],
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
            get_input(inputs)

    def test_run_network(self):
        """
        Checks whether combine_data returns
        a batch object if inputs are supplied.
        """
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': testloc[0],
                  'otu_table': None,
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'spar': None,
                  'conet': (os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3'),
                  'conet_bash': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': 20,
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None}
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        get_input(inputs)
        inputs['settings'] = inputs['fp'] + '/settings.json'
        run_network(inputs)
        test = Path(inputs['fp'] + "/conet_family_test.txt")
        self.assertTrue(test.is_file())
        call(("rm " + inputs['biom_file'][0]))
        call(("rm " + inputs['fp'] + "/settings.json"))
        call(("rm " + inputs['fp'] + "/test_family.hdf5"))
        call(("rm " + inputs['fp'] + "/test_otu.hdf5"))
        call(("rm " + inputs['fp'] + "/conet_test_family.txt"))
        call(("rm " + inputs['fp'] + "/spiec-easi_test_family.txt"))

    def test_run_neo4j(self):
        """
        Checks whether combine_data returns
        a batch object if inputs are supplied.
        """
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': testloc[0],
                  'otu_table': None,
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'spar': None,
                  'conet_bash': None,
                  'conet': (os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3'),
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': 20,
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None,
                  'address': 'bolt://localhost:7687',
                  'username': 'neo4j',
                  'password': 'test',
                  'quit': False,
                  'clear': False,
                  'write': False,
                  'add': False,
                  'output': 'network',
                  'neo4j': (os.path.dirname(massoc.__file__)[:-6] + 'tests\\neo4j')}
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        get_input(inputs)
        inputs['settings'] = inputs['fp'] + '/settings.json'
        run_network(inputs)
        inputs['job'] = 'start'
        run_neo4j(inputs)
        inputs['job'] = 'upload'
        run_neo4j(inputs)
        inputs['job'] = 'write'
        run_neo4j(inputs)
        inputs['job'] = 'clear'
        run_neo4j(inputs)
        inputs['job'] = 'quit'
        run_neo4j(inputs)
        test = Path(inputs['fp'] + "/network.graphml")
        self.assertTrue(test.is_file())
        call(("rm " + inputs['biom_file'][0]))
        call(("rm " + inputs['fp'] + "/settings.json"))
        call(("rm " + inputs['fp'] + "/test_family.hdf5"))
        call(("rm " + inputs['fp'] + "/test_otu.hdf5"))
        call(("rm " + inputs['fp'] + "/conet_test_family.txt"))
        call(("rm " + inputs['fp'] + "/spiec-easi_test_family.txt"))
        call(("rm " + inputs['fp'] + "/network.graphml"))

    def test_run_netstats(self):
        """
        Checks whether combine_data returns
        a batch object if inputs are supplied.
        """
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': testloc[0],
                  'otu_table': None,
                  'tools': ['spiec-easi', 'conet'],
                  'spiec': None,
                  'spar': None,
                  'conet': (os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3'),
                  'conet_bash': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': 20,
                  'min': None,
                  'rar': None,
                  'name': ['test'],
                  'cores': None,
                  'address': 'bolt://localhost:7687',
                  'username': 'neo4j',
                  'password': 'test',
                  'quit': False,
                  'clear': False,
                  'write': False,
                  'add': False,
                  'output': 'network',
                  'logic': ['union', 'difference', 'intersection'],
                  'neo4j': (os.path.dirname(massoc.__file__)[:-6] + 'tests\\neo4j')}
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        get_input(inputs)
        inputs['settings'] = inputs['fp'] + '/settings.json'
        run_network(inputs)
        inputs['job'] = 'start'
        run_neo4j(inputs)
        inputs['job'] = 'upload'
        run_neo4j(inputs)
        run_netstats(inputs)
        inputs['job'] = 'clear'
        run_neo4j(inputs)
        inputs['job'] = 'quit'
        run_neo4j(inputs)
        test = Path(inputs['fp'] + "/difference_network.graphml")
        self.assertTrue(test.is_file())
        call(("rm " + inputs['biom_file'][0]))
        call(("rm " + inputs['fp'] + "/settings.json"))
        call(("rm " + inputs['fp'] + "/test_family.hdf5"))
        call(("rm " + inputs['fp'] + "/test_otu.hdf5"))
        call(("rm " + inputs['fp'] + "/conet_test_family.txt"))
        call(("rm " + inputs['fp'] + "/spiec-easi_test_family.txt"))
        call(("rm " + inputs['fp'] + "/network.graphml"))
        call(("rm " + inputs['fp'] + "/difference_network.graphml"))
        call(("rm " + inputs['fp'] + "/union_network.graphml"))
        call(("rm " + inputs['fp'] + "/intersection_network.graphml"))

if __name__ == '__main__':
    unittest.main()
