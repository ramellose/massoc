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
                                "pH":"2.0",
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample2", "metadata":{
                                "pH":"1.8",       
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample3", "metadata":{
                                "pH":"2.3",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample4", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CGCTTATCGAGA",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample5", "metadata":{
                                "pH":"2.0",        
                                "BarcodeSequence":"CATACCAGTAGC",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample6", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample7", "metadata":{
                                "pH":"1.9",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample8", "metadata":{
                                "pH":"1.9",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},    
        {"id":"Sample9", "metadata":{
                                "pH":"1.8",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},
        {"id":"Sample10", "metadata":{
                                "pH":"2.1",        
                                "BarcodeSequence":"CTCTCTACCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"gut",
                                "Description":"human gut"}},                                    
        {"id":"Sample11", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample12", "metadata":{
                                "pH":"6.9",        
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample13", "metadata":{
                                "pH":"7.1",        
                                "BarcodeSequence":"CTCTCGGCCTGT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample14", "metadata":{
                                "pH":"7.0",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample15", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample16", "metadata":{
                                "pH":"6.9",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample17", "metadata":{
                                "pH":"6.7",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},   
        {"id":"Sample18", "metadata":{
                                "pH":"7.2",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},
        {"id":"Sample19", "metadata":{
                                "pH":"6.8",        
                                "BarcodeSequence":"CTCTCTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}},                                                                                                                         
        {"id":"Sample20", "metadata":{
                                "pH":"7.0",        
                                "BarcodeSequence":"CTAACTACCAAT",
                                "LinkerPrimerSequence":"CATGCTGCCTCCCGTAGGAGT",
                                "BODY_SITE":"skin",
                                "Description":"human skin"}}
        ],
     "matrix_type": "sparse",
     "matrix_element_type": "int",
     "shape": [5, 20],
     "data":[[0,10,5],
             [0,11,5],
             [0,12,6],
             [0,13,5],
             [0,14,5],
             [0,15,5],
             [0,16,6],
             [0,17,5],
             [0,18,5],
             [0,19,6],
             [0,9,6],
             [1,0,5],
             [1,1,1],
             [1,3,2],
             [1,4,3],
             [1,8,5],
             [1,10,1],
             [1,11,2],
             [1,2,3],
             [1,14,5],
             [1,17,1],
             [1,12,2],
             [1,19,1],
             [2,2,1],
             [2,3,4],
             [2,5,2],
             [2,6,1],
             [2,8,4],
             [2,10,2],
             [2,14,4],
             [2,16,2],
             [3,0,2],
             [3,1,1],
             [3,2,1],
             [3,5,1],
             [3,7,2],
             [3,12,1],
             [3,15,2],
             [3,7,1],
             [3,10,1],
             [3,11,1],
             [4,1,1],
             [4,2,1],
             [4,4,1],
             [4,14,1],
             [4,6,1]
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
networks.networks['test_g'] = g

f = nx.Graph()
f.add_nodes_from(nodes)
f.add_edges_from([("GG_OTU_1", "GG_OTU_2"),
                  ("GG_OTU_2", "GG_OTU_3"), ("GG_OTU_3", "GG_OTU_4")])
f["GG_OTU_1"]["GG_OTU_2"]['weight'] = 1.0
f["GG_OTU_2"]["GG_OTU_3"]['weight'] = 1.0
f["GG_OTU_3"]["GG_OTU_4"]['weight'] = -1.0
networks.networks['test_f'] = f

# run database
filepath= inputs['neo4j'][0] + '/bin/neo4j.bat console'
filepath = filepath.replace("\\", "/")
p = Popen(filepath, shell=True, stdout = PIPE)

from massoc.scripts.netbase import ImportDriver
from massoc.scripts.netstats import Driver

drive = ImportDriver(user='neo4j', password='test', uri='bolt://localhost:7687')
drive.clear_database()
drive.convert_nets(networks, mode='weight')

statdrive = Driver(user='neo4j', password='test', uri='bolt://localhost:7687')


class TestMain(unittest.TestCase):
    """"
    Tests whether the main pipeline functions
    are able to execute the appropriate commands.
    Note that the current functions only check whether network inference runs at all;
    we should have a checked dataset with known outputs to make sure inference is run correctly.
    """

    def test_get_union(self):
        """There are two networks with 6 edges, but only 4 unique ones;
        the union function should collect these. """
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        union = statdrive.graph_union()
        self.assertEqual(len(union), 4)

    def test_get_intersection(self):
        """There are two networks with 6 edges, but only 4 unique ones;
        the intersection function should only collect the 3 present in both. """
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        intersection = statdrive.graph_intersection()
        self.assertEqual(len(intersection), 2)

    def test_get_difference(self):
        """When no networks are specified, only nodes that are present in 1 network are returned.
        Otherwise, nodes specific to the specified network are returned."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        difference = statdrive.graph_difference()
        self.assertEqual(len(difference), 2)

    def test_get_pairlist(self):
        """Checks if pairlists of similar associations are returned."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        pairs = statdrive.get_pairlist(level='Genus', mode='weight')
        name1 = pairs[0]['p'].nodes[0]
        name2 = pairs[0]['r'].nodes[0]
        self.assertEqual(name1, name2)

    def test_agglomerate_network(self):
        """Tests if the associations between OTU1-OTU2 and OTU2-OTU5
        are agglomerated. """
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        num_assocs = drive.custom_query("MATCH (n:Association) RETURN n")
        statdrive.agglomerate_network(level='Genus', mode='weight')
        new_assocs = drive.custom_query("MATCH (n:Association) RETURN n")
        self.assertGreater(len(num_assocs), len(new_assocs))

    def test_associate_samples_spearman(self):
        """Tests if associations between continuous sample variables and taxa are added.
        The first OTU was adjusted to only occur in skin samples."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        statdrive.associate_samples(properties=['BarcodeSequence',
                                                'LinkerPrimerSequence',
                                                'BODY_SITE',
                                                'Description', 'pH'])
        spear = drive.custom_query("MATCH (n:Taxon)-[r:SPEARMAN]-() RETURN n")
        self.assertEqual(len(spear), 1)

    def test_associate_samples_hypergeom(self):
        """Tests if associations between categorical sample variables and taxa are added.
        The first OTU was adjusted to only occur in skin samples."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        statdrive.associate_samples(properties=['BarcodeSequence',
                                                'LinkerPrimerSequence',
                                                'BODY_SITE',
                                                'Description'])
        hypergeom = drive.custom_query("MATCH p=(n:Taxon)-[r:HYPERGEOM]-() RETURN p")
        self.assertEqual(len(hypergeom), 3)

    def test_associate_pair(self):
        """Tests if associations between sample variables and taxa are added for 1 taxon.
        The first OTU was adjusted to only occur in skin samples."""
        drive.clear_database()
        drive.convert_nets(networks, mode='weight')
        taxon = drive.custom_query("MATCH (n:Taxon {name: 'GG_OTU_1'}) RETURN n")[0]['n'].get('name')
        statdrive.associate_taxon(taxon=taxon, mode='Taxon', null_input=None,
                                  properties=['BarcodeSequence',
                                              'LinkerPrimerSequence',
                                              'BODY_SITE',
                                              'Description'])
        hypergeom = drive.custom_query("MATCH p=(n:Taxon)-[r:HYPERGEOM]-() RETURN p")
        self.assertEqual(len(hypergeom), 3)

if __name__ == '__main__':
    unittest.main()
