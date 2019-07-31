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
from itertools import combinations

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

    def graph_intersection(self, networks=None, weight=True):
        """
        Returns a subgraph that contains all nodes present in both specified networks.
        If no networks are specified, the function returns only nodes that are
        connected to all nodes in the network.

        :param networks: List of network names
        :param weight: If false, the intersection includes associations with matching partners but different weights
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        intersection = None
        try:
            with self._driver.session() as session:
                intersection = session.read_transaction(self._get_intersection, networks, weight=weight)
        except Exception:
            logger.error("Could not obtain graph intersection. ", exc_info=True)
        return intersection

    def graph_difference(self, networks=None, weight=True):
        """
        Returns a subgraph that contains all nodes only present in one of the selected networks.
        If no networks are specified, returns all associations that are unique across multiple networks.

        :param networks: List of network names
        :param weight: If false, the difference excludes associations with matching partners but different weights
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        difference = None
        try:
            with self._driver.session() as session:
                difference = session.read_transaction(self._get_difference, networks, weight=weight)
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
        assocs = tx.run(("WITH " + str(networks) +
                         " as names MATCH (n:Association)-->(b:Network) "
                         "WHERE b.name in names RETURN n")).data()
        assocs = _get_unique(assocs, 'n')
        _write_logic(tx, operation='Union', networks=networks, assocs=assocs)

    @staticmethod
    def _get_intersection(tx, networks, weight, n):
        """
        Accesses database to return edge list of intersection of networks.

        :param tx: Neo4j transaction
        :param networks: List of network names
        :param weight: If false, the intersection includes associations with matching partners but different weights
        :param n: If specified, number of networks that the intersecting node should be in
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        queries = list()
        if not n:
            for node in networks:
                queries.append(("MATCH (n:Association)-->(:Network {name: '" +
                                node + "'}) "))
            query = " ".join(queries) + "RETURN n"
            assocs = tx.run(query).data()
            assocs = list(_get_unique(assocs, 'n'))
        else:
            assocs = list()
            combos = combinations(networks, n)
            for combo in combos:
                for node in combo:
                    queries.append(("MATCH (n:Association)-->(:Network {name: '" +
                                    node + "'}) "))
                query = " ".join(queries) + "RETURN n"
                combo_assocs = tx.run(query).data()
                combo_assocs = list(_get_unique(combo_assocs, 'n'))
                assocs.append(combo_assocs)
        if weight:
            query = ("MATCH (a)-[:WITH_TAXON]-(n:Association)-[:WITH_TAXON]-(b) "
                     "MATCH (a)-[:WITH_TAXON]-(m:Association)-[:WITH_TAXON]-(b) "
                     "WHERE (n.name <> m.name) RETURN n, m")
            weighted = tx.run(query).data()
            filter_weighted = list()
            for assoc in weighted:
                # check whether associations are in all networks
                in_networks = list()
                nets = tx.run(("MATCH (a:Association {name: '"
                               + assoc['n'].get('name') +
                               "'})-->(n:Network) RETURN n")).data()
                nets = _get_unique(nets, 'n')
                in_networks.extend(nets)
                nets = tx.run(("MATCH (a:Association {name: '"
                               + assoc['m'].get('name') +
                               "'})-->(n:Network) RETURN n")).data()
                nets = _get_unique(nets, 'n')
                in_networks.extend(nets)
                if all(x in in_networks for x in networks):
                    filter_weighted.append(assoc)
            assocs.extend(_get_unique(filter_weighted, 'n'))
            assocs.extend(_get_unique(filter_weighted, 'm'))
        if weight:
            name = 'Intersection_weight'
        else:
            name = 'Intersection'
        _write_logic(tx, operation=name, networks=networks, assocs=assocs)

    @staticmethod
    def _get_difference(tx, networks, weight):
        """
        Accesses database to return edge list of difference of networks.

        :param tx: Neo4j transaction
        :param networks: List of network names
        :param weight: If false, the difference excludes associations with matching partners but different weights
        :return: Edge list of lists containing source, target, network and weight of each edge.
        """
        assocs = list()
        for network in networks:
            assocs.extend(tx.run(("MATCH (n:Association)-->(:Network {name: '" + network +
                                  "'}) WITH n MATCH (n)-[r]->(:Network) WITH n, count(r) "
                                  "as num WHERE num=1 RETURN n")).data())
        assocs = _get_unique(assocs, 'n')
        if weight:
            cleaned = list()
            for assoc in assocs:
                query = ("MATCH (a)-[:WITH_TAXON]-(n:Association {name: '" + assoc +
                         "'})-[:WITH_TAXON]-(b) "
                         "MATCH (a)-[:WITH_TAXON]-(m:Association)-[:WITH_TAXON]-(b) "
                         "WHERE (n.name <> m.name) RETURN n, m")
                check = tx.run(query).data()
                if len(check) == 0:
                    cleaned.append(assoc)
            assocs = cleaned
        if weight:
            name = 'Difference_weight'
        else:
            name = 'Difference'
        _write_logic(tx, operation=name, networks=networks, assocs=assocs)


def _write_logic(tx, operation, networks, assocs):
    """
    Accesses database to return edge list of intersection of networks.

    :param tx: Neo4j transaction
    :param operation: Type of logic operation
    :param networks: List of network names
    :param assocs: List of associations returned by logic operation
    :return:
    """
    name = operation + '_' + '_'.join(networks)
    tx.run("CREATE (n:Set {name: $id}) "
           "RETURN n", id=name)
    for assoc in assocs:
        tx.run(("MATCH (a:Association), (b:Set) WHERE a.name = '" +
                assoc +
                "' AND b.name = '" + name +
                "' CREATE (a)-[r:IN_SET]->(b) "
                "RETURN type(r)"))


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