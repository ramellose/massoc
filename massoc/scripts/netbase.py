"""
These functions convert files in the original BIOM format into a Neo4J format,
and stores this in a local database.
To reduce the size of the database, it is important to perform prevalence filtering + rarefaction first.
Helpful: find taxonomy + associations of taxa
MATCH p=(:Association)-->(n:Taxon) MATCH r=(n:Taxon)-->() RETURN p,r LIMIT 500

Functions to export to GraphML are also part of this module.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

from neo4j.v1 import GraphDatabase
from uuid import uuid4  # generates unique IDs for associations + observations
import networkx as nx
from massoc.scripts.netstats import _get_unique
import numpy as np
import logging
import sys
import os
import json
import requests
import re

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)
import logging.handlers


class ImportDriver(object):

    def __init__(self, uri, user, password, filepath):
        """
        Initializes a driver for accessing the Neo4j database.
        This driver constructs the Neo4j database and uploads extra data.

        :param uri: Adress of Neo4j database
        :param user: Username for Neo4j database
        :param password: Password for Neo4j database
        :param filepath: Filepath where logs will be written.
        """
        _create_logger(filepath)
        try:
            self._driver = GraphDatabase.driver(uri, auth=(user, password))
        except Exception:
            logger.error("Unable to start driver. \n", exc_info=True)
            sys.exit()

    def close(self):
        """
        Closes the connection to the database.
        :return:
        """
        self._driver.close()

    def clear_database(self):
        """
        Clears the entire database.
        :return:
        """
        try:
            with self._driver.session() as session:
                session.write_transaction(self._delete_all)
        except Exception:
            logger.error("Could not clear database. \n", exc_info=True)

    def custom_query(self, query):
        """
        Accepts a query and provides the results.
        :param query: String containing Cypher query
        :return: Results of transaction with Cypher query
        """
        output = None
        try:
            with self._driver.session() as session:
                output = session.read_transaction(self._query, query)
        except Exception:
            logger.error("Unable to execute query: " + query + '\n', exc_info=True)
        return output

    def convert_nets(self, nets):
        """
        Converts Nets object from netwrap.py to a Neo4J graph.
        This graph can be stored more easily in a database,
        and supports future metadata-based utilities.

        :param nets: Nets object
        :return:
        """
        try:
            for exp_id in nets.inputs['name']:
                biomfile = nets.otu[exp_id]
                self.convert_biom(biomfile, exp_id)
                with self._driver.session() as session:
                    for net in nets.networks:
                        name = net.split('.')[0]
                        self.convert_networkx(network=nets.networks[net],
                                              network_id=name, mode='weight', exp_id=exp_id)
        except Exception:
            logger.error("Could not port network object to database. \n", exc_info=True)

    def convert_networkx(self, network_id, network, exp_id=None, log=None, mode=None):
        """
        Uploads NetworkX object to Neo4j database.
        :param network_id: Name for network node.
        :param network: NetworkX object.
        :param exp_id: Name of experiment used to generate network.
        :param log: Log of steps carried out to generate network
        :param mode: if 'weight, weighted associations are uploaded
        :return:
        """
        try:
            with self._driver.session() as session:
                session.write_transaction(self._create_network, network_id, exp_id, log)
                session.write_transaction(self._create_associations, network_id, network, mode)
        except Exception:
            logger.error("Could not write networkx object to database. \n", exc_info=True)

    def convert_biom(self, biomfile, exp_id):
        """
        Stores a BIOM object in the database.

        :param biomfile: BIOM file.
        :param exp_id: Label of experiment used to generate BIOM file.
        :return:
        """
        try:
            # first check if sample metadata exists
            tax_meta = biomfile.metadata(axis='observation')
            sample_meta = biomfile.metadata(axis='sample')
            with self._driver.session() as session:
                session.write_transaction(self._create_experiment, exp_id)
                for taxon in biomfile.ids(axis='observation'):
                    session.write_transaction(self._create_taxon, taxon, biomfile)
                    tax_index = biomfile.index(axis='observation', id=taxon)
                    if tax_meta:
                        meta = biomfile.metadata(axis='observation')[tax_index]
                        for key in meta:
                            if key != 'taxonomy' and type(meta[key]) == str:
                                session.write_transaction(self._create_property,
                                                          source=taxon, sourcetype='Taxon',
                                                          target=meta[key], name=key)
                for sample in biomfile.ids(axis='sample'):
                    session.write_transaction(self._create_sample, sample, exp_id)
                    sample_index = biomfile.index(axis='sample', id=sample)
                    if sample_meta:
                        meta = biomfile.metadata(axis='sample')[sample_index]
                        # need to clean up these 'if' conditions to catch None properties
                        # there is also a problem with commas + quotation marks here
                        for key in meta:
                            # meta[key] = re.sub(r'\W+', '', str(meta[key]))
                            session.write_transaction(self._create_property,
                                                      source=sample, sourcetype='Sample',
                                                      target=meta[key], name=key, weight=None)
            obs_data = biomfile.to_dataframe()
            rows, cols = np.where(obs_data.values != 0)
            observations = list()
            for taxon, sample in list(zip(obs_data.index[rows], obs_data.columns[cols])):
                value = obs_data[sample][taxon]
                observations.append((taxon, sample, value))
            with self._driver.session() as session:
                for observation in observations:
                    session.write_transaction(self._create_observations, observation)
        except Exception:
            logger.error("Could not write BIOM file to database. \n", exc_info=True)

    def include_nodes(self, nodes, name, label, check=True):
        """
        Given a dictionary, this function tries to upload
        the file to the Neo4j database.
        The first column of the edgelist should reflect nodes
        already present in the Neo4j graph database,
        while the second column reflects node names that will be added.
        The column names are used to assign node types to the new metadata.

        The dictionary should contain another dictionary of target nodes and edge weights.

        :param nodes: Dictionary of existing nodes as values with node names as keys
        :param name: Name of variable, inserted in Neo4j graph database as type
        :param label: Label of source node (e.g. Taxon, Sample, Property, Experiment etc)
        :param check: If True, checks if all source nodes appear in the database.
        :return:
        """
        # first step:
        # check whether key values in node dictionary exist in network
        if check:
            with self._driver.session() as session:
                matches = session.read_transaction(self._find_nodes, list(nodes.keys()))
                if not matches:
                    logger.warning('No source nodes are present in the network. \n')
                    sys.exit()
        with self._driver.session() as session:
            for node in nodes:
                session.write_transaction(self._create_property,
                                          source=node, sourcetype=label,
                                          target=nodes[node]['target'], name=name, weight=nodes[node]['weight'])

    def find_nodes(self, nodes):
        """
        Returns 'True' if all names in the node list are found in the database.

        :param nodes: Dictionary of existing nodes as values with node names as keys
        :return: Boolean
        """
        with self._driver.session() as session:
            match = session.read_transaction(self._find_nodes, nodes)
        return match

    def export_network(self, path, networks=None):
        """
        Writes networks to graphML file.
        If no path is given, the network is returned as a NetworkX object.
        :param path: Filepath where network is written to.
        :param networks: Names of networks to write to disk.
        :return:
        """
        results = None
        try:
            results = self._return_networks(networks)
            if path:
                for network in results:
                    name = path + '/' + network + '.graphml'
                    nx.write_graphml(results[network], name)
        except Exception:
            logger.error("Could not write database graph to GraphML file. \n", exc_info=True)
        return results

    def export_cyto(self, networks=None):
        """
        Writes networks to Cytoscape.
        :param networks: Names of networks to write to disk.
        :return:
        """
        results = self._return_networks(networks)
        # Basic Setup
        PORT_NUMBER = 1234
        BASE = 'http://localhost:' + str(PORT_NUMBER) + '/v1/'
        HEADERS = {'Content-Type': 'application/json'}
        for network in results:
            # Define dictionary from networkx
            network_dict = {
                'data':
                    {"node_default": {}, "edge_default": {}, 'name': network},
                'elements': {
                    'nodes': [],
                    'edges': []
                }
            }
            i = 1
            id_dict = {}
            for node in results[network].nodes:
                id_dict[node] = i
                data = {'data': {'name': node, 'id': str(i),
                                 'SUID': int(i), 'selected': False,  'shared_name': node}}
                for property in results[network].nodes[node]:
                    data['data'][property] = results[network].nodes[node][property]
                network_dict['elements']['nodes'].append(data)
                i += 1
            i = 1
            for edge in results[network].edges:
                data = {'data': {'shared_name': edge[0] + '->' + edge[1],
                                 'name': edge[0] + '->' + edge[1],
                                 'source': str(id_dict[edge[0]]), 'target': str(id_dict[edge[1]]),
                                 'id': str(i), 'SUID': int(i), 'selected': False},
                        'selected': False}
                for property in results[network].edges[edge]:
                    if property == 'source':
                        # source is source node in Cytoscape
                        # but network source in Neo4j
                        data['data']['networks'] = results[network].edges[edge][property]
                    else:
                        data['data'][property] = results[network].edges[edge][property]
                network_dict['elements']['edges'].append(data)
                i += 1
            res = requests.post(BASE + 'networks?collection=Neo4jexport', data=json.dumps(network_dict),
                                headers=HEADERS)
            new_network_id = res.json()['networkSUID']
            print('Network created for ' + network + ': SUID = ' + str(new_network_id))

    def create_association_relationships(self):
        """
        For each association, a link is created from 1 taxon to another,
        with the weight as edge weight. The returned subgraph
        is similar to graphs that are returned by CoNet and SPIEC-EASI directly.
        :return:
        """
        try:
            with self._driver.session() as session:
                network = session.write_transaction(self._create_rel_assoc)
            return network
        except Exception:
            logger.error("Could not create association shortcuts. \n", exc_info=True)

    def include_sequences(self, location):
        """
        This function opens a folder of FASTA sequences with identifiers
        matching to OTU identifiers in the Neo4j database.
        The FASTA sequences are converted to a dictionary and uploaded to
        the database with the Neo4j driver function include_nodes.

        :param location: Folder containing FASTA sequences matching to OTU identifiers.
        I.e. GreenGenes FASTA files are accepted.
        :param driver: ImportDriver object
        :return: Updates database with 16S sequences.
        """
        # get list of taxa in database
        with self._driver.session() as session:
            taxa = session.read_transaction(self._get_list, 'Taxon')
        sequence_dict = dict()
        logger.info("Found " + str(len(os.listdir(location))) + " files.")
        for filename in os.listdir(location):
            with open(location + '//' + filename, 'r') as file:
                lines = file.readlines()
                logger.info("16S file " + filename + " contains " + str(int(len(lines)/2)) + " sequences.")
            for i in range(0, len(lines), 2):
                otu = lines[i].rstrip()[1:]  # remove > and \n
                sequence = lines[i + 1].rstrip()
                sequence_dict[otu] = sequence
        # with the sequence list, run include_nodes
        seqs_in_database = taxa.intersection(sequence_dict.keys())
        sequence_dict = {k: {'target': v, 'weight': None} for k, v in sequence_dict.items() if k in seqs_in_database}
        logger.info("Uploading " + str(len(sequence_dict)) + " sequences.")
        self.include_nodes(sequence_dict, name="16S", label="Taxon", check=False)

    def _return_networks(self, networks):
        """
        Returns NetworkX networks from the Neo4j database.

        :param networks: Names of networks to return.
        :return: Dictionary of networks
        """
        results = dict()
        with self._driver.session() as session:
            tax_dict = session.read_transaction(self._tax_dict)
        with self._driver.session() as session:
            tax_properties = session.read_transaction(self._tax_properties)
        for item in tax_properties:
            for taxon in tax_properties[item]:
                tax_properties[item][taxon] = str(tax_properties[item][taxon])
        if not networks:
            with self._driver.session() as session:
                networks = session.read_transaction(self._query,
                                                    "MATCH (n:Network) RETURN n")
                networks.extend(session.read_transaction(self._query,
                                                         "MATCH (n:Set) RETURN n"))
            networks = list(_get_unique(networks, key='n'))
        # create 1 network per database
        for network in networks:
            g = nx.MultiGraph()
            with self._driver.session() as session:
                edge_list = session.read_transaction(self._association_list, network)
            for edge in edge_list[0]:
                index_1 = edge[0]
                index_2 = edge[1]
                all_weights = []
                try:
                    all_weights = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",
                                             edge_list[1][edge][0])
                except TypeError:
                    if type(all_weights) == list:
                        all_weights = edge_list[1][edge][0]
                    else:
                        all_weights = []
                all_weights = [float(x) for x in all_weights]
                weight = float(np.mean(all_weights))
                g.add_edge(index_1, index_2, source=str(edge_list[0][edge]),
                           weight=weight, all_weights=str(all_weights))
            # necessary for networkx indexing
            for item in tax_dict:
                nx.set_node_attributes(g, tax_dict[item], item)
            for item in tax_properties:
                nx.set_node_attributes(g, tax_properties[item], item)
            g = g.to_undirected()
            results[network] = g
        return results

    @staticmethod
    def _delete_all(tx):
        """
        Deletes all nodes and their relationships from the database.
        :param tx: Neo4j transaction
        :return:
        """
        tx.run("MATCH (n) DETACH DELETE n")

    @staticmethod
    def _query(tx, query):
        """
        Processes custom queries.
        :param tx: Neo4j transaction
        :param query: String of Cypher query
        :return: Outcome of transaction
        """
        results = tx.run(query).data()
        return results

    @staticmethod
    def _create_rel_assoc(tx):
        """
        Requests all associations, and then generates a new relationship
        that 'shortcuts' between the taxa. This relationship is
        directed, because the database is directed, so queries
        should be treating them like undirected relationships.

        :param tx: Neo4j transaction
        :return:
        """
        assocs = tx.run("MATCH (n:Association) RETURN n").data()
        for assoc in assocs:
            nodes = tx.run(("MATCH (m)--(:Association {name: '" +
                            assoc['n'].get('name') +
                            "'})--(n) WHERE (m:Taxon OR m:Agglom_Taxon) AND"
                            "(n:Taxon OR n:Agglom_taxon) RETURN m, n LIMIT 1")).data()
            tax1 = nodes[0]['m'].get('name')
            tax2 = nodes[0]['n'].get('name')
            tx.run(("MATCH (a),(b) WHERE a.name = '" + tax1 + "' AND b.name = '" +
                    tax2 + "' CREATE (a)-[r:ASSOCIATES_WITH]->(b) RETURN r"))
        new_assocs = tx.run("MATCH p=()-[r:ASSOCIATES_WITH]-() RETURN p").data()
        return new_assocs

    @staticmethod
    def _create_experiment(tx, exp_id):
        """
        Creates a node that represents the Experiment ID.

        :param tx: Neo4j transaction
        :param exp_id: Label for experiment
        :return:
        """
        tx.run("CREATE (a:Experiment) SET a.name = $id", id=exp_id)

    @staticmethod
    def _create_taxon(tx, taxon, biomfile):
        """
        Creates a node that represents a taxon.
        Also generates taxonomy nodes + connects them, and
        includes metadata.

        :param tx: Neo4j transaction
        :param taxon: ID for taxon
        :param biomfile: BIOM file containing count data.
        :return:
        """
        # first check if OTU already exists
        hit = tx.run(("MATCH (a:Taxon {name: '" + taxon +
                      "'}) RETURN a")).data()
        if len(hit) == 0:
            tx.run("CREATE (a:Taxon) SET a.name = $id", id=taxon)
            tax_index = biomfile.index(axis='observation', id=taxon)
            if biomfile.metadata(axis='observation'):
                # it is possible that there is no metadata
                tax_dict = biomfile.metadata(axis='observation')[tax_index]['taxonomy']
                tax_levels = ['Kingdom', 'Phylum', 'Class',
                              'Order', 'Family', 'Genus', 'Species']
                # rel_list = ['IS_KINGDOM', 'IS_PHYLUM', 'IS_CLASS',
                #             'IS_ORDER', 'IS_FAMILY', 'IS_GENUS', 'IS_SPECIES']
                # tree_list = ['PART_OF_KINGDOM', 'PART_OF_PHYLUM', 'PART_OF_CLASS',
                #             'PART_OF_ORDER', 'PART_OF_FAMILY', 'PART_OF_GENUS']
                if str(taxon) is not 'Bin':
                    # define range for which taxonomy needs to be added
                    j = 0
                    if tax_dict:
                        for i in range(len(tax_dict)):
                            level = tax_dict[i]
                            if sum(c.isalpha() for c in level) > 1:
                                # only consider adding as node if there is more
                                # than 1 character in the taxonomy
                                # first request ID to see if the taxonomy node already exists
                                j += 1
                        for i in range(0, j):
                            if tax_dict[i] != 'NA':  # maybe allow user input to specify missing values
                                hit = tx.run(("MATCH (a:" + tax_levels[i] +
                                              " {name: '" + tax_dict[i] +
                                              "'}) RETURN a")).data()
                                if len(hit) == 0:
                                    tx.run(("CREATE (a:" + tax_levels[i] +
                                            " {name: '" + tax_dict[i] +
                                            "', type: 'Taxonomy'}) RETURN a")).data()
                                    if i > 0:
                                        tx.run(("MATCH (a:" + tax_levels[i] +
                                                "), (b:" + tax_levels[i-1] +
                                                ") WHERE a.name = '" + tax_dict[i] +
                                                "' AND b.name = '" + tax_dict[i-1] +
                                                "' CREATE (a)-[r: PART_OF]->(b) "
                                                "RETURN type(r)"))
                                hit = tx.run(("MATCH (a:Taxon)-[r]-(b:" + tax_levels[i] +
                                              ") WHERE a.name = '" + taxon +
                                              "' AND b.name = '" + tax_dict[i] +
                                              "' RETURN type(r)")).data()
                                if len(hit) == 0:
                                    tx.run(("MATCH (a:Taxon), (b:" + tax_levels[i] +
                                            ") WHERE a.name = '" + taxon +
                                            "' AND b.name = '" + tax_dict[i] +
                                            "' CREATE (a)-[r: BELONGS_TO]->(b) "
                                            "RETURN type(r)"))
            else:
                tx.run("CREATE (a:Taxon) SET a.name = $id", id='Bin')

    @staticmethod
    def _create_sample(tx, sample, exp_id):
        """
        Creates sample nodes and related metadata.

        :param tx: Neo4j transaction
        :param sample: Sample name
        :param exp_id: Experiment name
        :return:
        """
        tx.run("CREATE (a:Sample) SET a.name = $id", id=sample)
        tx.run(("MATCH (a:Sample), (b:Experiment) "
                "WHERE a.name = '" + sample +
                "' AND b.name = '" + exp_id +
                "' CREATE (a)-[r:IN_EXPERIMENT]->(b) "
                "RETURN type(r)"))

    @staticmethod
    def _create_property(tx, source, target, name, weight, sourcetype=''):
        """
        Creates target node if it does not exist yet
        and adds the relationship between target and source.

        :param tx: Neo4j transaction
        :param source: Source node, should exist in database
        :param target: Target node
        :param name: Type variable of target node
        :param weight: Weight of relationship
        :param sourcetype: Type variable of source node (not required)
        :return:
        """
        hit = tx.run(("MATCH (a:Property) WHERE a.name = '" +
                      target + "' AND a.type = '" +
                      name + "' RETURN a")).data()
        if len(hit) == 0:
            tx.run(("CREATE (a:Property) "
                    "SET a.name = '" + target +
                    "' SET a.type = '" + name + "' "
                                                "RETURN a")).data()
        if len(sourcetype) > 0:
            sourcetype = ':' + sourcetype
        matching_rel = tx.run(("MATCH (a" + sourcetype + ")-[r:HAS_PROPERTY]-(b:Property) "
                "WHERE a.name = '" + source +
                "' AND b.name = '" + target +
                "' AND b.type = '" + name +
                "' RETURN r")).data()
        if weight:
            rel = " {weight: [" + weight + "]}"
        else:
            rel = ""
        if len(matching_rel) == 0:
            tx.run(("MATCH (a" + sourcetype + "), (b:Property) "
                    "WHERE a.name = '" + source +
                    "' AND b.name = '" + target +
                    "' AND b.type = '" + name +
                    "' CREATE (a)-[r:HAS_PROPERTY" + rel + "]->(b) "
                    "RETURN type(r)"))


    @staticmethod
    def _create_observations(tx, observation):
        """
        Creates relationships between taxa and samples
        that represent the count number of that taxon in a sample.

        :param tx: Neo4j transaction
        :param observation: An observation (count) of a taxon in a sample.
        :return:
        """
        taxon, sample, value = observation
        tx.run(("MATCH (a:Taxon), (b:Sample) "
                "WHERE a.name = '" + taxon +
                "' AND b.name = '" + sample +
                "' CREATE (a)-[r:FOUND_IN]->(b) "
                "SET r.count = '" + str(value) +
                "' RETURN type(r)"))

    @staticmethod
    def _create_network(tx, network, exp_id=None, log=None):
        """
        Generates a network node with provenance for every network
        stored in a Nets object.

        :param tx: Neo4j transaction
        :param network: Network name
        :param exp_id: Experiment name
        :param log: Dictionary of operations carried out to generate network
        :return:
        """
        tx.run("CREATE (a:Network) "
               "SET a.name = $id "
               "RETURN a", id=network)
        if exp_id:
            tx.run(("MATCH (a:Network), (b:Experiment) "
                    "WHERE a.name = '" + network +
                    "' AND b.name = '" + exp_id +
                    "' CREATE (a)-[r:CREATED_FROM]->(b) "
                    "RETURN type(r)"))
        if log:
            for metadata in log:
                if metadata in network:
                    tx.run(("MATCH (a:Network)"
                            "WHERE a.name = '" + network +
                            "' SET a.tool = '" + metadata +
                            "' RETURN a"))
                    for network_property in log[metadata]:
                        tx.run(("MATCH (a:Network)"
                                "WHERE a.name = '" + network +
                                "' SET a." + network_property +
                                " = '" + log[metadata][network_property] +
                                "' RETURN a"))
                else:
                    if type(log[metadata]) is not dict:
                        # ensures metadata for other tools is not included
                        tx.run(("MATCH (a:Network)"
                                " WHERE a.name = '" + network +
                                "' SET a." + metadata +
                                " = '" + log[metadata] +
                                "' RETURN a"))

    @staticmethod
    def _create_associations(tx, name, network, mode=None):
        """
        Generates all the associations contained in a network and
        connects them to the related network node.
        This function uses NetworkX networks as source.

        :param tx: Neo4j transaction
        :param name: Network name
        :param network: NetworkX object
        :param mode: if 'weight', edge weights are uploaded
        :return:
        """
        # creates metadata for eventual CoNet feature associations
        for edge in network.edges:
            taxon1 = edge[0]
            taxon2 = edge[1]
            attr = network.get_edge_data(taxon1, taxon2)
            # networkx files imported from graphml will have an index
            # the 'name' property changes the name to the index
            # should probably standardize graphml imports or something
            if 'name' in network.nodes[taxon1]:
                taxon1 = network.nodes[taxon1]['name']
                taxon2 = network.nodes[taxon2]['name']
            # for CoNet imports, the taxa can actually be features
            # need to check this, because these taxa will NOT be in the dataset
            # in that case, we need to create nodes that represent
            # the features
            feature = False
            if 'isafeature' in network.nodes[edge[0]]:
                if network.nodes[edge[0]]['isafeature'] == 'yes':
                    feature = True
                    # find if TaxonProperty already exists
                    # otherwise, create it and link to Taxon
                    hit = tx.run(("MATCH (a:Property) WHERE a.type = '" +
                                  taxon1 + "' AND a.name = '" +
                                  attr['interactionType'] + "' RETURN a")).data()
                    if len(hit) == 0:
                        tx.run(("CREATE (a:Property) SET a.type = '" +
                                taxon1 + "' SET a.name = '" +
                                attr['interactionType'] + "'"))
                    tx.run(("MATCH (a:Taxon), (b:Property) WHERE a.name = '" + taxon2 +
                            "' AND b.type ='" + taxon1 +
                            "' AND b.name = '" + attr['interactionType'] +
                            "' CREATE (a)-[r:HAS_PROPERTY]->(b) RETURN type(r)"))
                elif network.nodes[edge[1]]['isafeature'] == 'yes':
                    feature = True
                    # find if TaxonProperty already exists
                    # otherwise, create it and link to Taxon
                    hit = tx.run(("MATCH (a:Property) WHERE a.type = '" +
                                  taxon2 + "' AND a.name = '" +
                                  attr['interactionType'] + "' RETURN a")).data()
                    if len(hit) == 0:
                        tx.run(("CREATE (a:Property) SET a.type = '" +
                                taxon2 + "' SET a.name = '" +
                                attr['interactionType'] + "'"))
                    tx.run(("MATCH (a:Taxon), (b:Property) WHERE a.name = '" + taxon1 +
                            "' AND b.type ='" + taxon2 +
                            "' AND b.name = '" + attr['interactionType'] +
                            "' CREATE (a)-[r:HAS_PROPERTY]->(b) RETURN type(r)"))
            if not feature:
                if mode == 'weight':
                    network_weight = str(attr['weight'])
                hit = tx.run(("MATCH p=(a)<--(:Association)-->(b) "
                              "WHERE a.name = '" + taxon1 +
                              "' AND b.name = '" + taxon2 +
                              "' RETURN p")).data()
                # first check if association is already present)
                if len(hit) > 0:
                    weights = [network_weight]
                    for association in hit:
                        uid = association['p'].nodes[1].get('name')
                        # first check if there is already a link between the association and network
                        network_hit = tx.run(("MATCH p=(a:Association)--(b:Network) "
                                              "WHERE a.name = '" +
                                              uid +
                                              "' AND b.name = '" + name +
                                              "' RETURN p")).data()
                        database_weight = association['p'].nodes[1].get('weight')
                        try:
                            database_weight = re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?",
                                                         database_weight)
                        except TypeError:
                            if type(database_weight) == list:
                                pass
                            else:
                                database_weight = []
                        weights.extend(database_weight)
                        if len(network_hit) == 0:
                            tx.run(("MATCH (a:Association), (b:Network) "
                                    "WHERE a.name = '" +
                                    uid +
                                    "' AND b.name = '" + name +
                                    "' CREATE (a)-[r:IN_NETWORK]->(b) "
                                    "RETURN type(r)"))
                        tx.run(("MATCH (a:Association) WHERE a.name = '" +
                                uid +
                                "' SET a.weight = " +
                                str(weights) +
                                " RETURN a"))
                else:
                    uid = str(uuid4())
                    # non alphanumeric chars break networkx
                    if mode == 'weight':
                        tx.run("CREATE (a:Association {name: $id}) "
                               "SET a.weight = $weight "
                               "RETURN a", id=uid, weight=str([network_weight]))
                    else:
                        tx.run("CREATE (a:Association {name: $id}) "
                               "RETURN a", id=uid)
                    tx.run(("MATCH (a:Association), (b:Taxon) "
                            "WHERE a.name = '" +
                            uid +
                            "' AND b.name = '" + taxon1 +
                            "' CREATE (a)-[r:WITH_TAXON]->(b) "
                            "RETURN type(r)"))
                    tx.run(("MATCH (a:Association), (b:Taxon) "
                            "WHERE a.name = '" +
                            uid +
                            "' AND b.name = '" + taxon2 +
                            "' CREATE (a)-[r:WITH_TAXON]->(b) "
                            "RETURN type(r)"))
                    tx.run(("MATCH (a:Association), (b:Network) "
                            "WHERE a.name = '" +
                            uid +
                            "' AND b.name = '" + name +
                            "' CREATE (a)-[r:IN_NETWORK]->(b) "
                            "RETURN type(r)"))

    @staticmethod
    def _get_list(tx, label):
        """
        Returns a list of nodes with the specified label.

        :param tx: Neo4j transaction
        :param label: Neo4j database label of nodes
        :return: List of nodes with specified label.
        """
        results = tx.run(("MATCH (n:" + label + ") RETURN n")).data()
        results = _get_unique(results, key="n")
        return results

    @staticmethod
    def _get_fasta(tx):
        """
        Generates a string of FASTA sequences.

        :param tx: Neo4j transaction
        :return: String of FASTA sequences.
        """
        results = tx.run("MATCH (n:Taxon)--(m:Property {type: '16S'}) RETURN n,m").data()
        fasta_dict = {}
        for result in results:
            fasta_dict[result['n']['name']] = result['m']['name']
        fasta_string = str()
        for key in fasta_dict:
            fasta_string += '>' + key + '\n' + fasta_dict[key] + '\n'
        return fasta_string

    @staticmethod
    def _association_list(tx, network):
        """
        Returns a list of associations, as taxon1, taxon2, and, if present, weight.

        :param tx: Neo4j transaction
        :param network: Name of network node
        :return: List of lists with source and target nodes, source networks and edge weights.
        """
        associations = tx.run(("MATCH (n:Association)--(b {name: '" + network +
                               "'}) RETURN n")).data()
        networks = dict()
        weights = dict()
        for assoc in associations:
            taxa = tx.run(("MATCH (m)--(:Association {name: '" + assoc['n'].get('name') +
                           "'})--(n) "
                           "WHERE (m:Taxon OR m:Agglom_Taxon) AND (n:Taxon OR n:Agglom_Taxon) "
                           "AND m.name <> n.name "
                           "RETURN m, n LIMIT 1")).data()
            if len(taxa) == 0:
                pass  # apparently this can happen. Need to figure out why!!
            else:
                edge = (taxa[0]['m'].get('name'), taxa[0]['n'].get('name'))
                network = tx.run(("MATCH (:Association {name: '" + assoc['n'].get('name') +
                                  "'})-->(n:Network) RETURN n"))
                network = _get_unique(network, key='n')
                network_list = list()
                for item in network:
                    network_list.append(item)
                weight = [assoc['n'].get('weight')]
                # it is possible for sets to contain associations with different weights
                if edge in networks.keys():
                    network_list.extend(networks[edge])
                    networks[edge] = set(network_list)
                    weight.extend(weights[edge])
                    weights[edge] = weight
                elif (edge[1], edge[0]) in networks.keys():
                    network_list.extend(networks[(edge[1], edge[0])])
                    networks[(edge[1], edge[0])] = set(network_list)
                    weight.extend(weights[(edge[1], edge[0])])
                    weights[(edge[1], edge[0])] = weight
                else:
                    networks[edge] = network_list
                    weights[edge] = weight
        edge_list = (networks, weights)
        return edge_list

    @staticmethod
    def _sample_list(tx):
        """
        Returns a list of sample occurrences, as taxon, sample and count.

        :param tx: Neo4j transaction
        :return: List of samples
        """
        tax_nodes = tx.run("MATCH (n)--(:Association) WHERE n:Taxon RETURN n")
        tax_nodes = _get_unique(tax_nodes, 'n')
        edge_list = list()
        for node in tax_nodes:
            samples = tx.run(("MATCH (:Taxon {name: '" + node +
                              "'})-[r:Sample]-(n:Property) RETURN n, r"))
            for sample in samples:
                sublist = list()
                sublist.append(node)
                sublist.append(sample['n'].get('name'))
                sublist.append(sample['n'].get('type'))
                sublist.append(str(type(sample['r'])))
                sublist.append(sample['r'].get('correlation'))
                edge_list.append(sublist)
        return edge_list

    @staticmethod
    def _agglom_list(tx):
        """
        Returns a list of relationships between agglomerated taxa and taxa.

        :param tx: Neo4j transaction
        :return: List of nodes labeled Agglom_Taxon
        """
        agglom_nodes = tx.run("MATCH (n:Agglom_Taxon) RETURN n")
        agglom_nodes = _get_unique(agglom_nodes, 'n')
        edge_list = list()
        for node in agglom_nodes:
            taxa = tx.run("MATCH (:Agglom_Taxon {name: '" + node +
                          "'})--(n:Taxon) RETURN n")
            for taxon in taxa:
                sublist = list()
                sublist.append(node)
                sublist.append(taxon['n'].get('name'))
                edge_list.append(sublist)
        return edge_list

    @staticmethod
    def _find_nodes(tx, names):
        """
        Returns True if all nodes in the 'names' list are found in the database.

        :param tx: Neo4j transaction
        :param names: List of names of nodes
        :return:
        """
        for name in names:
            netname = tx.run("MATCH (n {name: '" + name +
                             "'}) RETURN n").data()
            netname = _get_unique(netname, key='n')
            # only checking node name; should be unique in database!
            found = True
            if len(netname) == 0:
                found = False
            elif len(netname) > 1:
                logger.warning("Duplicated node name in database! \n")
        return found

    @staticmethod
    def _tax_dict(tx):
        """
        Returns a dictionary of taxonomic values for each node.

        :param tx: Neo4j transaction
        :return: Dictionary of taxonomy separated by taxon
        """
        taxa = tx.run("MATCH (n)--(:Association) WHERE n:Taxon OR n:Agglom_Taxon RETURN n").data()
        taxa = _get_unique(taxa, 'n')
        tax_dict = dict()
        tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        for item in tax_levels:
            tax_dict[item] = dict()
        for item in taxa:
            for level in tax_levels:
                tax = None
                level_name = tx.run("MATCH (b {name: '" + item +
                                 "'})--(n:"+ level + ") RETURN n").data()
                if len(level_name) != 0:
                    tax = level_name[0]['n'].get('name')
                if tax:
                    tax_dict[level][item] = tax
        return tax_dict

    @staticmethod
    def _tax_properties(tx):
        """
        Returns a dictionary of taxon / sample properties, to be included as taxon metadata.

        :param tx: Neo4j transaction
        :return: Dictionary of dictionary of taxon properties
        """
        nodes = tx.run("MATCH (n)--(m:Property) WHERE n:Taxon OR n:Agglom_Taxon RETURN m").data()
        nodes = _get_unique(nodes, 'm')
        properties = dict()
        for node in nodes:
            property = tx.run("MATCH (m:Property) RETURN m").data()
            property_key = property[0]['m']['type']
            properties[property_key] = dict()
            hits = tx.run("MATCH (b)--(n {name: '" + node +
                          "'}) WHERE b:Taxon OR b:Agglom_Taxon RETURN b").data()
            if hits:
                for hit in hits:
                    properties[property_key][hit['b'].get('name')] = property[0]['m']['name']
            for taxon in properties[property_key]:
                if len(properties[property_key][taxon]) == 1:
                    # tries exporting property as float instead of list
                    try:
                        properties[property_key][taxon] = np.round(float(properties[property_key][property]), 4)
                    except ValueError:
                        pass
        return properties


def _create_logger(filepath):
    """
    After a filepath has become available, loggers can be created
    when required to report on errors.

    :param filepath: Filepath where logs will be written.
    :return:
    """
    logpath = filepath + '/massoc.log'
    # filelog path is one folder above massoc
    # pyinstaller creates a temporary folder, so log would be deleted
    fh = logging.handlers.RotatingFileHandler(maxBytes=500,
                                              filename=logpath, mode='a')
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)