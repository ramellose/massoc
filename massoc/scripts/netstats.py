"""
This module describes a driver for carrying out analyses on Neo4j graph databases.
For example, the module can extract union, difference or intersecton of uploaded networks,
or it can cluster on the networks inside the Neo4j database.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'


from neo4j.v1 import GraphDatabase
import sys
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


class NetDriver(object):

    def __init__(self, uri, user, password, filepath):
        """
        Initializes a driver for accessing the Neo4j database.
        This driver carries out analyses (logic operations and statistics) on the database.

        Note that weight is automatically taken into account when carrying out logic operations;
        the functions look for associations that have specific relationships with networks.
        Because association nodes are separated if they have different weights,
        they will only be returned if they match the pattern for an assocation with a specific weight.

        I.e. the intersection of CoNet and SPIEC-EASI will only return an association
        if CoNet and SPIEC-EASI predicted the same sign.

        The difference of CoNet and SPIEC-EASI will return assocations with the same taxa
        if CoNet and SPIEC-EASI predicted a contrasting sign.

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
            exit()

    def close(self):
        """
        Closes the connection to the database.

        :return:
        """
        self._driver.close()

    def custom_query(self, query):
        """
        Accepts a query and provides the results.

        :param query: String containing Cypher query
        :return:
        """
        try:
            with self._driver.session() as session:
                output = session.read_transaction(self._query, query)
            return(output)
        except Exception:
            logger.error("Unable to execute query: " + query + '', exc_info=True)

    def graph_union(self, networks=None):
        """
        Returns a subgraph that contains all nodes present in all networks.
        If network names are specified as a list, all nodes present
        in these two networks are returned.

        :param networks: List of network names
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        union = None
        try:
            with self._driver.session() as session:
                union = session.read_transaction(self._get_union, networks)
        except Exception:
            logger.error("Could not obtain graph union. ", exc_info=True)
        return union

    def graph_intersection(self, networks=None):
        """
        Returns a subgraph that contains all nodes present in both specified networks.
        If no networks are specified, the function returns only nodes that are
        connected to all nodes in the network.

        :param networks: List of network names
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        intersection = None
        try:
            with self._driver.session() as session:
                intersection = session.read_transaction(self._get_intersection, networks)
        except Exception:
            logger.error("Could not obtain graph intersection. ", exc_info=True)
        return intersection

    def graph_difference(self, network=None):
        """
        Returns a subgraph that contains all nodes only present in one of the selected networks.
        If no networks are specified, returns all associations that are unique across multiple networks.

        :param networks: List of network names
        :param weight: Whether weight should be taken into account
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        difference = None
        try:
            with self._driver.session() as session:
                difference = session.read_transaction(self._get_difference, network)
        except Exception:
            logger.error("Could not obtain graph difference. ", exc_info=True)
        return difference

    @staticmethod
    def _query(tx, query):
        """
        Processes custom queries.

        :param tx: Neo4j transaction
        :param query: String containing Cypher query
        :return:
        """
        results = tx.run(query).data()
        return results

    @staticmethod
    def _get_weight(tx, node):
        """
        Returns the weight of an Association node.

        :param tx: Neo4j transaction
        :param node: Returns the weight of the specified node
        :return:
        """
        weight = tx.run("MATCH (n:Association {name: '" + node.get('name') +
                        "'}) RETURN n").data()[0]['n'].get('weight')
        return weight

    @staticmethod
    def _get_union(tx, networks):
        """
        Accesses database to return edge list of union of networks.

        :param tx: Neo4j transaction
        :param networks: List of network names
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        if not networks:
            assocs = tx.run("MATCH (n:Association) RETURN n").data()
        else:
            assocs = tx.run(("WITH " + str(networks) +
                             " as names MATCH (n:Association)-->(b:Network) "
                             "WHERE b.name in names RETURN n")).data()
        assocs = _get_unique(assocs, 'n')
        edge_list = list()
        for assoc in assocs:
            sublist = list()
            taxa = tx.run(("MATCH (m)--(:Association {name: '" + assoc +
                           "'})--(n) "
                           "WHERE (m:Taxon OR m:Agglom_Taxon) AND (n:Taxon OR n:Agglom_Taxon)"
                           "RETURN m, n LIMIT 1")).data()
            sublist.append(taxa[0]['m'].get('name'))
            sublist.append(taxa[0]['n'].get('name'))
            network = tx.run(("MATCH (:Association {name: '" + assoc +
                              "'})--(n:Network) RETURN n")).data()
            network_list = list()
            for item in network:
                network_list.append(item['n'].get('name'))
            sublist.append(str(network_list))
            weight = tx.run(("MATCH (n:Association {name: '" + assoc +
                             "'}) RETURN n")).data()[0]['n'].get('weight')
            if weight:
                sublist.append(weight)
            edge_list.append(sublist)
        return edge_list

    @staticmethod
    def _get_intersection(tx, networks):
        """
        Accesses database to return edge list of intersection of networks.

        :param tx: Neo4j transaction
        :param networks: List of network names
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        if not networks:
            networks = list()
            hits = tx.run("MATCH (n:Network) RETURN n").data()
            for hit in hits:
                networks.append(hit['n'].get('name'))
        queries = list()
        for node in networks:
            queries.append(("MATCH (n:Association)-->(:Network {name: '" +
                            node + "'}) "))
        query = " ".join(queries) + "RETURN n"
        assocs = tx.run(query).data()
        assocs = _get_unique(assocs, 'n')
        edge_list = list()
        for assoc in assocs:
            sublist = list()
            taxa = tx.run(("MATCH (m)--(:Association {name: '" + assoc +
                           "'})--(n) "
                           "WHERE (m:Taxon OR m:Agglom_Taxon) AND (n:Taxon OR n:Agglom_Taxon)"
                           "RETURN m, n LIMIT 1")).data()
            sublist.append(taxa[0]['m'].get('name'))
            sublist.append(taxa[0]['n'].get('name'))
            network = tx.run(("MATCH (:Association {name: '" + assoc +
                              "'})--(n:Network) RETURN n")).data()
            network_list = list()
            for item in network:
                network_list.append(item['n'].get('name'))
            sublist.append(str(network_list))
            weight = tx.run(("MATCH (n:Association {name: '" + assoc +
                             "'}) RETURN n")).data()[0]['n'].get('weight')
            if weight:
                sublist.append(weight)
            edge_list.append(sublist)
        return edge_list

    @staticmethod
    def _get_difference(tx, network):
        """
        Accesses database to return edge list of difference of networks.

        :param tx: Neo4j transaction
        :param networks: List of network names
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        if not network:
            assocs = tx.run("MATCH p=(n:Association)-[r]->(:Network) "
                            "WITH n, count(r) as num "
                            "WHERE num=1 RETURN n").data()
        else:
            assocs = tx.run(("MATCH (n:Association)-->(:Network {name: '" + network +
                             "'}) WITH n MATCH (n)-[r]->(:Network) WITH n, count(r) "
                             "as num WHERE num=1 RETURN n")).data()
        assocs = _get_unique(assocs, 'n')
        edge_list = list()
        for assoc in assocs:
            sublist = list()
            taxa = tx.run(("MATCH (m)--(:Association {name: '" + assoc +
                           "'})--(n) "
                           "WHERE (m:Taxon OR m:Agglom_Taxon) AND (n:Taxon OR n:Agglom_Taxon)"
                           "RETURN m, n LIMIT 1")).data()
            sublist.append(taxa[0]['m'].get('name'))
            sublist.append(taxa[0]['n'].get('name'))
            network = tx.run(("MATCH (:Association {name: '" + assoc +
                              "'})--(n:Network) RETURN n")).data()
            network_list = list()
            for item in network:
                network_list.append(item['n'].get('name'))
            sublist.append(str(network_list))
            weight = tx.run(("MATCH (n:Association {name: '" + assoc +
                             "'}) RETURN n")).data()[0]['n'].get('weight')
            if weight:
                sublist.append(weight)
            edge_list.append(sublist)
        return edge_list


def _get_unique(node_list, key, mode=None):
    """
    Returns number or names of unique nodes in a list.

    :param node_list: List of dictionaries returned by Neo4j transactions.
    :param key: Key accessing specific node in dictionary.
    :param mode: If 'num', the number of unique nodes is returned.
    :return: Unique nodes (list of nodes) or node number
    """
    unique_samples = list()
    for item in node_list:
        unique_samples.append(item[key].get('name'))
    unique_samples = set(unique_samples)
    if mode == 'num':
        unique_samples = len(unique_samples)
    return unique_samples


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