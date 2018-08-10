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
from biom.cli.util import write_biom_table
from subprocess import call
import massoc
from massoc.scripts.main import combine_data, run_jobs, run_massoc, run_parallel, get_joblist

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
          'conet': None,
          'spar': None,
          'spar_pval': None,
          'spar_boot': None,
          'levels': ['otu', 'order'],
          'prev': ['20'],
          'cores': ['4'],
          'neo4j': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\neo4j')]}


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
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
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
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        bioms = combine_data(inputs)
        call(("rm " + inputs['biom_file'][0]))
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
            bioms = combine_data(inputs)

    def test_run_massoc(self):
        """Checks if the run_massoc function can run."""
        inputs = {'biom_file': [(testloc[0]+'/data/test.biom')],
                  'cluster': None,
                  'otu_meta': None,
                  'prefix': None,
                  'sample_data': None,
                  'split': None,
                  'tax_table': None,
                  'fp': [testloc[0]],
                  'otu_table': None,
                  'tools': ['conet'],
                  'spiec': None,
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['otu'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None,
                  'min': ['10'],
                  'rar': None,
                  'conet': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3')]}
        write_biom_table(testbiom['test'], fmt='hdf5', filepath=(inputs['biom_file'][0]))
        networks = run_massoc(inputs, mode=None)
        call(("rm " + inputs['biom_file'][0]))
        self.assertEqual(len(networks.networks), 1)


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
                  'tax_table': [(testloc[0] + 'otu_tax.txt')],
                  'fp': [testloc[0]],
                  'otu_table': [(testloc[0] + 'otu_otus.txt')],
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
                  'fp': [testloc[0]+'/data'],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['conet'],
                  'spiec': None,
                  'conet': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3')],
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'name': ['test'],
                  'cores': None,
                  'min': ['10']}
        batch = Batch(testbiom, inputs)
        netbatch = Nets(batch)
        jobs = get_joblist(netbatch)
        netbatch.collapse_tax()
        netbatch.write_bioms()
        network = run_jobs(netbatch, jobs[0])
        x = inputs['name'][0]
        filename = netbatch.inputs['fp'][0] + '/' + x + '_species.hdf5'
        call("rm " + filename)
        filename = netbatch.inputs['fp'][0] + '/' + x + '_genus.hdf5'
        call("rm " + filename)
        filename = netbatch.inputs['fp'][0] + '/' + x + '_family.hdf5'
        call("rm " + filename)
        filename = netbatch.inputs['fp'][0] + '/' + x + '_order.hdf5'
        call("rm " + filename)
        filename = netbatch.inputs['fp'][0] + '/' + x + '_class.hdf5'
        call("rm " + filename)
        filename = netbatch.inputs['fp'][0] + '/' + x + '_phylum.hdf5'
        call("rm " + filename)
        call(("rm " + inputs['fp'][0] + '/' + inputs['tools'][0] + '_' + inputs['name'][0] + '_' + inputs['levels'][0] +'.hdf5'))
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
                  'fp': [testloc[0]],
                  'otu_table': [(testloc[0][:-17] + 'otu_otus.txt')],
                  'tools': ['conet'],
                  'spiec': None,
                  'conet': [(os.path.dirname(massoc.__file__)[:-6] + 'tests\\CoNet3')],  # cannot be used in general testing
                  'spar_pval': None,
                  'spar_boot': None,
                  'levels': ['family'],
                  'prev': ['20'],
                  'min': ['10'],
                  'name': ['test'],
                  'cores': None,
                  'rar': None}
        batch = Batch(testbiom, inputs)
        netbatch = Nets(batch)
        netbatch = run_parallel(netbatch)
        self.assertEqual(len(netbatch.networks), 1)


if __name__ == '__main__':
    unittest.main()
