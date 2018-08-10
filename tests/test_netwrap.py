"""
This file contains all testing functions for netwrap.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import os
import random
import unittest
from copy import deepcopy
from subprocess import call

import biom
from massoc.scripts.batch import Batch
from massoc.scripts.netwrap import Nets

import massoc
from massoc.scripts.main import run_parallel

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
inputs = {'biom_file': None,
          'cluster': None,
          'otu_meta': None,
          'prefix': None,
          'sample_data': None,
          'split': None,
          'tax_table': None,
          'fp': testloc,
          'name': ['test'],
          'otu_table': None,
          'tools': ['spiec-easi'],
          'spiec': None,
          'conet': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3')],
          'spar': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\SparCC')],
          'spar_pval': None,
          'spar_boot': None,
          'levels': ['otu', 'order'],
          'prev': ['20'],
          'cores': ['4'],
          'neo4j': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\neo4j')]}
netbatch =Nets(Batch(testbiom, inputs))

filenames = list()
for x in inputs['name']:
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_otu.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_species.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_genus.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_family.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_order.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_class.hdf5')
    filenames.append(netbatch.inputs['fp'][0] + '/' + x + '_phylum.hdf5')

class TestNetWrap(unittest.TestCase):
    """Tests netwrap.
    More specifically, checks ability to call network inference tools.
    Netwrap requires access to these tools, and they need to be installed.
    """

    def test_spar(self):
        """Check if the SparCC function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.collapse_tax()
        testnets.write_bioms()
        testnets.run_spar(boots=10)
        for name in filenames:
            call(("rm " + name))
        self.assertEqual(len(testnets.networks), 2)

    def test_conet(self):
        """Check if the CoNet function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.collapse_tax()
        testnets.write_bioms()
        testnets.run_conet()
        for name in filenames:
            call(("rm " + name))
        self.assertEqual(len(testnets.networks), 2)

    def test_spiec(self):
        """Check if the SPIEC-EASI function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.collapse_tax()
        testnets.write_bioms()
        testnets.run_spiec()
        for name in filenames:
            call(("rm " + name))
        self.assertEqual(len(testnets.networks), 2)

if __name__ == '__main__':
    unittest.main()
