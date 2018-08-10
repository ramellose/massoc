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
import re
import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

class ImportDriver(object):

    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        """Closes the connection to the database."""
        self._driver.close()

    def clear_database(self):
        """Clears the entire database."""
        try:
            with self._driver.session() as session:
                session.write_transaction(self._delete_all)
        except Exception:
            logger.error("Could not clear database.", exc_info=True)

    def custom_query(self, query):
        """Accepts a query and provides the results."""
        try:
            with self._driver.session() as session:
                output = session.read_transaction(self._query, query)
            return(output)
        except Exception:
            logger.error("Unable to execute query: " + query, exc_info=True)


    def convert_nets(self, nets, mode=None):
        """
        Converts Nets object from netwrap.py to a Neo4J graph.
        This graph can be stored more easily in a database,
        and supports future metadata-based utilities.
        """
        try:
            for exp_id in nets.inputs['name']:
                biomfile = nets.otu[exp_id]
                self.convert_biom(biomfile, exp_id)
                with self._driver.session() as session:
                    for net in nets.networks:
                        self.convert_networkx(network=nets.networks[net],
                                              network_id=net, mode='weight', exp_id=exp_id, log=nets.log)
        except Exception:
            logger.error("Could not port network object to database. ", exc_info=True)

    def convert_networkx(self, network_id, network, exp_id=None, log=None, mode=None):
        try:
            with self._driver.session() as session:
                session.write_transaction(self._create_network, network_id, exp_id, log)
                session.write_transaction(self._create_associations, network_id, network, mode)
        except Exception:
            logger.error("Could not write networkx object to database. ", exc_info=True)

    def convert_biom(self, biomfile, exp_id):
        """
        Stores a BIOM object in the database.
        """
        try:
            with self._driver.session() as session:
                session.write_transaction(self._create_experiment, exp_id)
                for taxon in biomfile.ids(axis='observation'):
                    session.write_transaction(self._create_taxon, taxon, biomfile)
                for sample in biomfile.ids(axis='sample'):
                    session.write_transaction(self._create_sample, sample, exp_id, biomfile)
                session.write_transaction(self._create_observations, biomfile)
        except Exception:
            logger.error("Could not write BIOM file to database. ", exc_info=True)

    def export_network(self, path, pairlist=None, mode='basic'):
        try:
            g = nx.Graph()
            with self._driver.session() as session:
                taxon_list = session.read_transaction(self._query,
                                                      "MATCH (n) WHERE n:Taxon OR n:Agglom_Taxon RETURN n")
            taxon_dict = dict()
            for i in range(len(taxon_list)):
                taxon_dict[taxon_list[i]['n'].get('name')] = str(i)
                g.add_node(str(i), name=taxon_list[i]['n'].get('name'))
            edge_list = pairlist
            if not pairlist:
                with self._driver.session() as session:
                    edge_list = session.read_transaction(self._association_list)
            for i in range(len(edge_list)):
                if len(edge_list[i]) == 4:
                    index_1 = taxon_dict[edge_list[i][0]]
                    index_2 = taxon_dict[edge_list[i][1]]
                    g.add_edge(index_1, index_2, source=edge_list[i][2], weight=edge_list[i][3])
                else:
                    index_1 = taxon_dict[edge_list[i][0]]
                    index_2 = taxon_dict[edge_list[i][1]]
                    g.add_edge(index_1, index_2, source=edge_list[i][2])
            if mode == 'sample':
                with self._driver.session() as session:
                    sample_list = session.read_transaction(self._query,
                                                           "MATCH (n:Sample) RETURN n")
                sample_dict = dict()
                for i in range(len(sample_list)):
                    sample_dict[sample_list[i]['n'].get('name')] = 's' + str(i)
                    g.add_node(('s' + str(i)), name=sample_list[i]['n'].get('name'))
                with self._driver.session() as session:
                    assoc_list = session.read_transaction(self._sample_list)
                for i in range(len(assoc_list)):
                    taxon = taxon_dict[assoc_list[i][0]]
                    sample_dict[assoc_list[i][1]] = 's' + str(i)
                    g.add_edge(taxon, sample_dict[assoc_list[i][1]],
                               type=assoc_list[i][2], test=assoc_list[i][3], correlation=assoc_list[i][4])
                with self._driver.session() as session:
                    agglom_list = session.read_transaction(self._agglom_list)
                for edge in agglom_list:
                    g.add_edge(taxon_dict[edge[0]], taxon_dict[edge[1]])
            with self._driver.session() as session:
                tax_dict = session.read_transaction(self._tax_dict)
            # necessary for networkx indexing
            new_tax_dict = dict()
            for item in tax_dict:
                new_dict = dict()
                for subdict in tax_dict[item]:
                    if len(subdict) > 0:
                        if list(subdict.keys())[0] in taxon_dict:
                            new_dict[taxon_dict[list(subdict.keys())[0]]] = \
                                subdict[list(subdict.keys())[0]]
                new_tax_dict[item] = new_dict
                nx.set_node_attributes(g, new_tax_dict[item], item)
            with self._driver.session() as session:
                tax_properties = session.read_transaction(self._tax_properties, mode='TaxonProperty')
            attribute_dict = dict()
            attrlist = list()
            for item in tax_properties:
                taxon = taxon_dict[item]
                attributes = tax_properties[item]
                attrlist.extend(list(attributes.keys()))
            attrlist = set(attrlist)
            for type in attrlist:
                attribute_dict[type] = dict()
                for taxon in taxon_dict:
                    attribute_dict[type][taxon_dict[taxon]] = None
            for item in tax_properties:
                taxon = taxon_dict[item]
                attributes = tax_properties[item]
                for type in attributes:
                    attribute_dict[type][taxon] = attributes[type]
            for item in attribute_dict:
                nx.set_node_attributes(g, attribute_dict[item], item)
            if mode != 'sample':
                with self._driver.session() as session:
                    sample_properties = session.read_transaction(self._tax_properties, mode='SampleProperty')
                attribute_dict = dict()
                attrlist = list()
                for item in sample_properties:
                    taxon = taxon_dict[item]
                    attributes = sample_properties[item]
                    attrlist.extend(list(attributes.keys()))
                attrlist = set(attrlist)
                for type in attrlist:
                    attribute_dict[type] = dict()
                    for taxon in taxon_dict:
                        attribute_dict[type][taxon_dict[taxon]] = 'None'
                for item in sample_properties:
                    taxon = taxon_dict[item]
                    attributes = sample_properties[item]
                    for type in attributes:
                        attribute_dict[type][taxon] = attributes[type]
                for item in attribute_dict:
                    nx.set_node_attributes(g, attribute_dict[item], item)
                nx.set_node_attributes(g, sample_properties[item], 'name')
            g = g.to_undirected()
            nx.write_graphml(g, path)
        except Exception:
            logger.error("Could not write database graph to GraphML file. ", exc_info=True)

    def create_association_relationships(self):
        """For each association, a link is created from 1 taxon to another,
        with the weight as edge weight. The returned subgraph
        is similar to graphs that are returned by CoNet and SPIEC-EASI directly. """
        try:
            with self._driver.session() as session:
                network = session.write_transaction(self._create_rel_assoc)
            return network
        except Exception:
            logger.error("Could not create association shortcuts. ", exc_info=True)


    @staticmethod
    def _delete_all(tx):
        """Deletes all nodes and their relationships from the database."""
        tx.run("MATCH (n) DETACH DELETE n")

    @staticmethod
    def _query(tx, query):
        """Processes custom queries."""
        results = tx.run(query).data()
        return results

    @staticmethod
    def _create_rel_assoc(tx):
        """Requests all associations, and then generates a new relationship
        that 'shortcuts' between the taxa. This relationship is
        directed, because the database is directed, so queries
        should be treating them like undirected relationships."""
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
        """Creates a node that represents the Experiment ID."""
        tx.run("CREATE (a:Experiment) SET a.name = $id", id=exp_id)

    @staticmethod
    def _create_taxon(tx, taxon, biomfile):
        """Creates a node that represents a taxon.
        Also generates taxonomy nodes + connects them, and
        includes metadata. """
        tx.run("CREATE (a:Taxon) SET a.name = $id", id=taxon)
        tax_index = biomfile.index(axis='observation', id=taxon)
        tax_dict = biomfile.metadata(axis='observation')[tax_index]['taxonomy']
        tax_levels = ['Kingdom', 'Phylum', 'Class',
                      'Order', 'Family', 'Genus', 'Species']
        rel_list = ['IS_KINGDOM', 'IS_PHYLUM', 'IS_CLASS',
                    'IS_ORDER', 'IS_FAMILY', 'IS_GENUS', 'IS_SPECIES']
        tree_list = ['PART_OF_KINGDOM', 'PART_OF_PHYLUM', 'PART_OF_CLASS',
                    'PART_OF_ORDER', 'PART_OF_FAMILY', 'PART_OF_GENUS']
        if str(taxon) is not 'Bin':
            meta = biomfile.metadata(axis='observation')[tax_index]
            for key in meta:
                if key != 'taxonomy' and type(meta[key]) == str:
                    hit = tx.run(("MATCH (a:TaxonProperty) WHERE a.name = '" +
                                  meta[key] + "' AND a.type = '" +
                                  key + "' RETURN a")).data()
                    if len(hit) == 0:
                        if key and meta[key]:
                            tx.run("CREATE (a:TaxonProperty) "
                                   "SET a.name = $value AND a.type = $type "
                                   "RETURN a", value=meta[key], type=key).data()
                    if key and meta[key]:
                        tx.run(("MATCH (a:Taxon), (b:TaxonProperty) "
                                "WHERE a.name = '" + taxon +
                                "' AND b.name = '" + meta[key] +
                                "' CREATE (a)-[r:HAS_PROPERTY]->(b) "
                                "RETURN type(r)"))
            # define range for which taxonomy needs to be added
            j = 0
            if tax_dict:
                for i in range(0, 7):
                    level = tax_dict[i]
                    if sum(c.isalpha() for c in level) > 1:
                        # only consider adding as node if there is more
                        # than 1 character in the taxonomy
                        # first request ID to see if the taxonomy node already exists
                        j += 1
                for i in range(0, j):
                    hit = tx.run(("MATCH (a:" + tax_levels[i] +
                                  " {name: '" + tax_dict[i] +
                                  "'}) RETURN a")).data()
                    if len(hit) == 0:
                        tx.run(("CREATE (a:" + tax_levels[i] +
                                " {name: '" + tax_dict[i] +
                                "'}) RETURN a")).data()
                        if i > 0:
                            tx.run(("MATCH (a:" + tax_levels[i] +
                                    "), (b:" + tax_levels[i-1] +
                                    ") WHERE a.name = '" + tax_dict[i] +
                                    "' AND b.name = '" + tax_dict[i-1] +
                                    "' CREATE (a)-[r:" + tree_list[i-1] +
                                    "]->(b) RETURN type(r)"))
                    tx.run(("MATCH (a:Taxon), (b:" + tax_levels[i] +
                            ") WHERE a.name = '" + taxon +
                            "' AND b.name = '" + tax_dict[i] +
                            "' CREATE (a)-[r:" + rel_list[i] +
                            "]->(b) RETURN type(r)"))
        else:
            tx.run("CREATE (a:Taxon) SET a.name = $id", id='Bin')

    @staticmethod
    def _create_sample(tx, sample, exp_id, biomfile):
        """Creates sample nodes and related metadata."""
        tx.run("CREATE (a:Sample) SET a.name = $id", id=sample)
        tx.run(("MATCH (a:Sample), (b:Experiment) "
                "WHERE a.name = '" + sample +
                "' AND b.name = '" + exp_id +
                "' CREATE (a)-[r:IN_EXPERIMENT]->(b) "
                "RETURN type(r)"))
        sample_index = biomfile.index(axis='sample', id=sample)
        meta = biomfile.metadata(axis='sample')[sample_index]
        # need to clean up these 'if' conditions to catch None properties
        # there is also a problem with commas + quotation marks here
        for key in meta:
            meta[key] = re.sub(r'\W+', '', str(meta[key]))
            hit = tx.run(("MATCH (a:SampleProperty) WHERE a.name = '" +
                          meta[key] + "' AND a.type = '" +
                          key + "' RETURN a")).data()
            if len(hit) == 0:
                if key and meta[key]:
                    tx.run(("CREATE (a:SampleProperty) "
                           "SET a.name = '" + meta[key] +
                           "' SET a.type = '" + key + "' "
                           "RETURN a")).data()
            if key and meta[key]:
                tx.run(("MATCH (a:Sample), (b:SampleProperty) "
                        "WHERE a.name = '" + sample +
                        "' AND b.name = '" + meta[key] +
                        "' AND b.type = '" + key +
                        "' CREATE (a)-[r:HAS_SAMPLE_PROPERTY]->(b) "
                        "RETURN type(r)"))

    @staticmethod
    def _create_observations(tx, biomfile):
        """Creates relationships between taxa and samples
        that represent the count number of that taxon in a sample."""
        obs_data = biomfile.to_dataframe()
        for sample in biomfile.ids(axis='sample'):
            values = obs_data[sample]
            for taxon in biomfile.ids(axis='observation'):
                value = values[taxon]
                #  only add non-zero entries
                if value > 0 and taxon is not 'Bin':
                    tx.run(("MATCH (a:Taxon), (b:Sample) "
                            "WHERE a.name = '" + taxon +
                            "' AND b.name = '" + sample +
                            "' CREATE (a)-[r:FOUND_IN]->(b) "
                            "SET r.count = '" + str(value) +
                            "' RETURN type(r)"))

    @staticmethod
    def _create_network(tx, network, exp_id=None, log=None):
        """Generates a network node with provenance for every network
        stored in a Nets object."""
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
        """Generates all the associations contained in a network and
        connects them to the related network node.
        This function uses NetworkX networks as source."""
        for edge in network.edges:
            taxon1 = edge[0]
            taxon2 = edge[1]
            attr = network.get_edge_data(taxon1, taxon2)
            if mode == 'weight':
                network_weight = str(attr['weight'])
            hit = tx.run(("MATCH p=(:Taxon {name: '" +
                          taxon1 + "'})<--(:Association)-->(:Taxon {name: '" +
                          taxon2 + "'}) RETURN p")).data()
            if mode == 'weight' and len(hit)>0:
                for node in hit:
                    database_weight = node['p'].nodes[1].get('weight')
                    if database_weight == network_weight:
                        matched_hit = node
                    else:
                        hit = []
                hit[0] = matched_hit
            # first check if association is already present)
            if len(hit) > 0:
                for association in hit:
                    uid = association['p'].nodes[1].get('name')
                    tx.run(("MATCH (a:Association), (b:Network) "
                            "WHERE a.name = '" +
                            uid +
                            "' AND b.name = '" + name +
                            "' CREATE (a)-[r:IN_NETWORK]->(b) "
                            "RETURN type(r)"))
            else:
                uid = str(uuid4())
                # non alphanumeric chars break networkx
                if mode == 'weight':
                    tx.run("CREATE (a:Association {name: $id}) "
                           "SET a.weight = $weight "
                           "RETURN a", id=uid, weight=str(network_weight))
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
    def _association_list(tx):
        """Returns a list of associations, as taxon1, taxon2, and, if present, weight. """
        associations = tx.run(("MATCH (n:Association) RETURN n")).data()
        edge_list = list()
        for assoc in associations:
            sublist = list()
            taxa = tx.run(("MATCH (m)--(:Association {name: '" + assoc['n'].get('name') +
                           "'})--(n) "
                           "WHERE (m:Taxon OR m:Agglom_Taxon) AND (n:Taxon OR n:Agglom_Taxon)"
                           "RETURN m, n LIMIT 1")).data()
            if len(taxa)== 0:
                pass # apparently this can happen. Need to figure out why!!
            else:
                sublist.append(taxa[0]['m'].get('name'))
                sublist.append(taxa[0]['n'].get('name'))
                network = tx.run(("MATCH (:Association {name: '" + assoc['n'].get('name') +
                               "'})--(n:Network) RETURN n"))
                network_list = list()
                for item in network:
                    network_list.append(item['n'].get('name'))
                sublist.append(str(network_list))
                weight = assoc['n'].get('weight')
                if weight:
                    sublist.append(weight)
                edge_list.append(sublist)
        return edge_list

    @staticmethod
    def _sample_list(tx):
        """Returns a list of sample occurrences, as taxon, sample and count."""
        tax_nodes = tx.run("MATCH (n)--(:Association) WHERE n:Taxon RETURN n")
        tax_nodes = _get_unique(tax_nodes, 'n')
        edge_list = list()
        for node in tax_nodes:
            samples = tx.run(("MATCH (:Taxon {name: '" + node +
                              "'})-[r]-(n:SampleProperty) RETURN n, r"))
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
        """Returns a list of relationships between agglomerated taxa and taxa."""
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
    def _tax_dict(tx):
        """Returns a dictionary of taxonomic values for each node. """
        taxa = tx.run("MATCH (n)--(:Association) WHERE n:Taxon OR n:Agglom_Taxon RETURN n").data()
        taxa = _get_unique(taxa, 'n')
        tax_dict = dict()
        tax_levels = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        for item in tax_levels:
            tax_dict[item] = list()
        for item in taxa:
            for level in tax_levels:
                subdict = dict()
                level_name = tx.run("MATCH (b {name: '" + item +
                                 "'})--(n:"+ level + ") RETURN n").data()
                if len(level_name) != 0:
                    subdict[item] = level_name[0]['n'].get('name')
                if len(subdict) != 0:
                    tax_dict[level].append(subdict)
        return tax_dict

    @staticmethod
    def _tax_properties(tx, mode='TaxonProperty'):
        """Returns a dictionary of taxon / sample properties, to be included as taxon metadata."""
        taxa = tx.run("MATCH (n)--(:Association) WHERE n:Taxon OR n:Agglom_Taxon RETURN n").data()
        taxa = _get_unique(taxa, 'n')
        properties = dict()
        for taxon in taxa:
            properties[taxon] = dict()
            hits = tx.run("MATCH (n:" + mode + ")--(b {name: '" + taxon +
                          "'}) WHERE b:Taxon OR b:Agglom_Taxon RETURN n").data()
            if hits:
                for hit in hits:
                    properties[taxon][hit['n'].get('type')] = list()
                for hit in hits:
                    properties[taxon][hit['n'].get('type')].append(hit['n'].get('name'))
            for property in properties[taxon]:
                properties[taxon][property] = list(set(properties[taxon][property]))
        clean_properties = dict()
        for taxon in properties:
            if len(properties[taxon]) != 0:
                clean_properties[taxon] = properties[taxon]
        return clean_properties

