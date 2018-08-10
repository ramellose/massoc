"""
This file contains all testing functions for the main program execution.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import unittest
import biom
from copy import deepcopy
import os
import massoc
from massoc.scripts.netwrap import Nets
from massoc.scripts.batch import Batch
import networkx as nx
from subprocess import Popen, PIPE

testloc = list()
testloc.append(os.path.dirname(massoc.__file__)[:-6])

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
networks = Nets(Batch(deepcopy(testbiom), inputs))
g = nx.Graph()
nodes = ["GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"]
g.add_nodes_from(nodes)
g.add_edges_from([("GG_OTU_1", "GG_OTU_2"),
                  ("GG_OTU_2", "GG_OTU_5"), ("GG_OTU_3", "GG_OTU_4")])
g["GG_OTU_1"]["GG_OTU_2"]['weight'] = 1.0
g["GG_OTU_2"]["GG_OTU_5"]['weight'] = 1.0
g["GG_OTU_3"]["GG_OTU_4"]['weight'] = -1.0
networks.networks['test_network'] = g

# run database
filepath= inputs['neo4j'][0] + '/bin/neo4j.bat console'
filepath = filepath.replace("\\", "/")
p = Popen(filepath, shell=True, stdout = PIPE)

from massoc.scripts.netbase import ImportDriver
from massoc.scripts.netstats import Driver

drive = ImportDriver(user='neo4j', password='test', uri='bolt://localhost:7687')
drive.clear_database()


class TestMain(unittest.TestCase):
    """"
    Tests whether the main pipeline functions
    are able to execute the appropriate commands.
    Note that the current functions only check whether network inference runs at all;
    we should have a checked dataset with known outputs to make sure inference is run correctly.
    """

    def test_create_association_relationships(self):
        """Checks whether the conversion of Association nodes
        to relationships works."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        drive.create_association_relationships()
        nodes = drive.custom_query(query="MATCH (n:Association) RETURN n")
        rels = drive.custom_query(query="MATCH ()-[r:ASSOCIATES_WITH]->() RETURN r")
        self.assertEqual(len(nodes), len(rels))

    def test_create_base(self):
        """Checks whether an object generated through netwrap
        can be converted into a neo4j database without raising an error."""
        drive.convert_nets(networks, 'weight')


    def test_pull_net(self):
        """After the database has been initialized,
        it should be possible to pull a network node."""
        drive.clear_database()
        drive.convert_nets(networks, 'weight')
        network = drive.custom_query(query=("MATCH (n:Network) RETURN n"))
        self.assertEqual(len(network), 1)

    def convert_biom(self):
        """BIOM files can be added to the database."""
        drive.clear_database()
        drive.convert_biom(testbiom['test'], 'test_biomfile')
        exp_id = drive.custom_query(query=("MATCH (n:Experiment {name: 'test_biomfile'}) RETURN n"))
        self.assertEqual(len(exp_id), 1)




if __name__ == '__main__':
    unittest.main()
