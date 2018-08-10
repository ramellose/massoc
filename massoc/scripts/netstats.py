"""
These functions perform operations on inferred networks to
provide users with more informative networks.
For example, find_netsets calculates sets such as union, intersection and difference.
build_assoctree builds a tree of taxa based on their predicted associations.
confident_centrality calculates centrality measures for a network and assigns confidence
intervals to the values through permutation.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'


from neo4j.v1 import GraphDatabase
from uuid import uuid4
from scipy.stats import hypergeom, spearmanr
import re
import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

class Driver(object):

    def __init__(self, uri, user, password):
        self._driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        """Closes the connection to the database."""
        self._driver.close()

    def clear_database(self):
        """Clears the entire database."""
        with self._driver.session() as session:
            session.write_transaction(self._delete_all)

    def graph_union(self, networks=None):
        """Returns a subgraph that contains all nodes present in all networks.
        If network names are specified as a list, all nodes present
        in these two networks are returned."""
        try:
            with self._driver.session() as session:
                union = session.read_transaction(self._get_union, networks)
            return union
        except Exception:
            logger.error("Could not obtain graph union. ", exc_info=True)

    def graph_intersection(self, networks=None):
        """Returns a subgraph that contains all nodes present in both specified networks.
        If no networks are specified, the function returns only nodes that are
        connected to all nodes in the network."""
        try:
            with self._driver.session() as session:
                intersection = session.read_transaction(self._get_intersection, networks)
            return intersection
        except Exception:
            logger.error("Could not obtain graph intersection. ", exc_info=True)

    def graph_difference(self, network=None):
        """Returns a subgraph that contains all nodes only present in one of the selected networks.
        If no networks are specified, returns all associations that are unique across multiple networks."""
        try:
            with self._driver.session() as session:
                difference = session.read_transaction(self._get_difference, network)
            return difference
        except Exception:
            logger.error("Could not obtain graph difference. ", exc_info=True)

    def agglomerate_network(self, level=None, mode='weight'):
        """Agglomerates to specified level, or, if no level is specified,
        over all levels. Associations are agglomerated based on similarity
        at the specified taxonomic level. If the mode is set to 'weight',
        associations are only agglomerated if their weight matches.
        The stop condition is the length of the pair list;
        as soon as no pair meets the qualification, agglomeration is terminated."""
        try:
            stop_condition = False
            while not stop_condition:
                # limit is necessary to prevent excessively long queries
                pairs = self.get_pairlist(level=level, mode=mode)
                if len(pairs) > 0:
                    pair = pairs[0]
                    self.agglomerate_pair(pair, level=level, mode=mode)
                else:
                    stop_condition = True
        except Exception:
            logger.error("Could not agglomerate associations to higher taxonomic levels. ", exc_info=True)

    def get_pairlist(self, level, mode):
        """Starts a new transaction for every pair list request."""
        try:
            with self._driver.session() as session:
                pairs = session.read_transaction(self._pair_list, level, mode)
                return pairs
        except Exception:
            logger.error("Could not obtain list of matching associations. ", exc_info=True)

    def agglomerate_pair(self, pair, level, mode):
        """For one pair, as returned by get_pairlist,
        this function creates new agglomerated nodes,
        deletes old agglomerated nodes, and chains taxonomic nodes
        to the new agglomerated nodes. Morever, the two old associations
        are deleted and replaced by a new association."""
        try:
            with self._driver.session() as session:
                agglom_1 = session.write_transaction(self._create_agglom)
                agglom_2 = session.write_transaction(self._create_agglom)
                session.write_transaction(self._chainlinks, agglom_1, pair['p'].nodes[1], pair['r'].nodes[1])
                session.write_transaction(self._chainlinks, agglom_2, pair['p'].nodes[3], pair['r'].nodes[3])
                session.write_transaction(self._taxonomy, agglom_1, pair['p'].nodes[0], level)
                session.write_transaction(self._taxonomy, agglom_2, pair['r'].nodes[4], level)
                networks = session.read_transaction(self._get_network, [pair['p'].nodes[2], pair['r'].nodes[2]])
                weight = session.read_transaction(self._get_weight, [pair['p'].nodes[2], pair['r'].nodes[2]])
                session.write_transaction(self._create_association, agglom_1, agglom_2, networks, weight, mode)
                session.write_transaction(self._delete_old_associations, [pair['p'].nodes[2], pair['r'].nodes[2]])
                session.write_transaction(self._delete_old_agglomerations, (pair['p'].nodes + pair['r'].nodes))
        except Exception:
            logger.error("Could not agglomerate a pair of matching associations. ", exc_info=True)

    def associate_samples(self, properties, null_input=None):
        """Sample identities themselves are not that informative,
        but the properties associated with them are.
        To test the hypothesis that taxa are associated with specific sample properties,
        the following tests are performed:
        1. For qualitative variables, a hypergeometric test is performed;
        how many associations do we expect by chance?
        2. For quantitative variables, Spearman correlation is performed.
        Because this is a hypothesis-generating tool,
        multiple-testing correction should be applied with care."""
        try:
            with self._driver.session() as session:
                tax_nodes = session.read_transaction(self._query, "MATCH (n)--(:Association) WHERE n:Taxon RETURN n")
                tax_nodes = _get_unique(tax_nodes, 'n')
            for node in tax_nodes:
                self.associate_taxon(mode='Taxon', taxon=node, null_input=null_input, properties=properties)
            with self._driver.session() as session:
                tax_nodes = session.read_transaction(self._query, "MATCH (n)--(:Association) WHERE n:Agglom_Taxon RETURN n")
                tax_nodes = _get_unique(tax_nodes, 'n')
            for node in tax_nodes:
                self.associate_taxon(taxon=node, mode='Agglom_Taxon', null_input=null_input, properties=properties)
        except Exception:
            logger.error("Could not associate sample variables to taxa. ", exc_info=True)

    def associate_taxon(self, taxon, mode, null_input, properties):
        """Tests whether specific sample properties can be associated to a taxon."""
        try:
            conts = list()
            categs = list()
            with self._driver.session() as session:
                if mode is 'Taxon':
                    query = "WITH " + str(properties) + \
                            " as types MATCH (:Taxon {name: '" + taxon + \
                            "'})-->(:Sample)-->(n:SampleProperty) WHERE n.type in types RETURN n"
                if mode is 'Agglom_Taxon':
                    query = "WITH " + str(properties) + \
                            " as types MATCH (:Agglom_Taxon {name: '" + taxon + \
                            "'})-[:GENERATED_FROM]-(:Taxon)--" \
                            "(:Sample)-->(n:SampleProperty) WHERE n.type in types RETURN n"
                sample_properties = session.read_transaction(self._query, query)
                for item in sample_properties:
                    value = item['n'].get('name')
                    if value == null_input:
                        break
                    # try to convert value to float; if successful, adds type to continous vars
                    try:
                        value = float(value)
                    except ValueError:
                        pass
                    if type(value) == float:
                        conts.append(item['n'].get('type'))
                    else:
                        categs.append([item['n'].get('type'), item['n'].get('name')])
            conts = set(conts)
            categs = set(tuple(categ) for categ in categs)
            for categ_val in categs:
                with self._driver.session() as session:
                    hypergeom_vals = session.read_transaction(self._hypergeom_population, taxon, categ_val, mode)
                    prob = hypergeom.cdf(hypergeom_vals['success_taxon'], hypergeom_vals['total_pop'],
                                         hypergeom_vals['success_pop'], hypergeom_vals['total_taxon'])
                    if prob < 0.05:
                        session.write_transaction(self._shortcut_categorical, taxon, categ_val, mode, prob)
            for cont_val in conts:
                with self._driver.session() as session:
                    spearman_result = session.read_transaction(self._spearman_test,
                                                               taxon, cont_val, mode)
                    if spearman_result.pvalue < 0.05:
                        var_dict = {cont_val: spearman_result.correlation}
                        session.write_transaction(self._shortcut_continuous, taxon, var_dict, mode)
        except Exception:
            logger.error("Could not associate a specific taxon to sample variables. ", exc_info=True)

    @staticmethod
    def _query(tx, query):
        """Processes custom queries."""
        results = tx.run(query).data()
        return results


    @staticmethod
    def _pair_list(tx, level, mode):
        """Returns a list of association pairs, where the
        taxonomic levels at both ends match, and the name of
        the associations are different. If 'weight' is specified as mode,
        only associations with identical weight are returned."""
        if mode is 'weight':
            result = tx.run(("MATCH p=(e:" +
                             level + ")<--()<--(a:Association)-->()-->(g:" +
                             level + ") MATCH r=(h:" + level +
                             ")<--()<--(b:Association)-->()-->(f:" +
                             level +
                             ") WHERE (a.name <> b.name) AND (a.weight = b.weight) AND "
                             "(e.name = h.name) AND (g.name = f.name) RETURN p,r LIMIT 1"))
        else:
            result = tx.run(("MATCH p=(e:" + level +
                             ")<--()<--(a:Association)-->()-->(g:" + level +
                             ") MATCH r=(h:" + level +
                             ")<--()<--(b:Association)-->()-->(f:" + level +
                             ") WHERE (a.name <> b.name) AND "
                             "(e.name = h.name) AND (g.name = f.name) "
                             "RETURN p,r LIMIT 1"))
        return result.data()

    @staticmethod
    def _create_agglom(tx):
        """Creates an Agglom_Taxon node and returns its id."""
        uid = str(uuid4())
        # non alphanumeric chars break networkx
        tx.run("CREATE (a:Agglom_Taxon) SET a.name = $id", id=uid)
        return uid

    @staticmethod
    def _chainlinks(tx, node, source1, source2):
        """Each Agglom_Taxon node is linked to the Taxon node
        it originated from. If it was generated from an Agglom_Taxon node,
        that source node's relationships to Taxon nodes are copied to the new node."""
        names = [source1.get('name'), source2.get('name')]
        for name in names:
            if len(name) == 36:
                hits = tx.run(("MATCH (:Agglom_Taxon {name: '" + name +
                               "'})-[:GENERATED_FROM]->(g) RETURN g"))
                for hit in hits.data():
                    old_link = tx.run(("MATCH p=(a:Agglom_Taxon)-->(b:Taxon) WHERE a.name = '" +
                                       node + "' AND b.name ='" + hit['g'].get('name') +
                                       "' RETURN p")).data()
                    if len(old_link) == 0:
                        tx.run(("MATCH (a:Agglom_Taxon),(b:Taxon) WHERE a.name = '" +
                                node + "' AND b.name = '" + hit['g'].get('name') +
                                "' CREATE (a)-[r:GENERATED_FROM]->(b) RETURN type(r)"))
            else:
                old_link = tx.run(("MATCH p=(a:Agglom_Taxon)-->(b:Taxon) WHERE a.name = '" +
                                   node + "' AND b.name ='" + name +
                                   "' RETURN p")).data()
                if len(old_link) == 0:
                    tx.run(("MATCH (a:Agglom_Taxon),(b:Taxon) WHERE a.name = '" +
                            node + "' AND b.name = '" + name +
                            "' CREATE (a)-[r:GENERATED_FROM]->(b) RETURN type(r)"))

    @staticmethod
    def _taxonomy(tx, node, tax, level):
        """Adds appropriate taxonomic relationships to taxonomic nodes.
        Generally, if this function returns an error because the 'tree' query
        came up empty, this means the phylogenetic tree was discontinous."""
        tax_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
        rel_list = ['IS_SPECIES', 'IS_GENUS', 'IS_FAMILY', 'IS_ORDER', 'IS_CLASS', 'IS_PHYLUM', 'IS_KINGDOM']
        level_id = tax_list.index(level)
        query_list = ["MATCH p=(:Species {name: '" + tax.get('name') +
                      "'})-->(:Genus)-->(:Family)-->(:Order)-->(:Class)-->(:Phylum)-->(:Kingdom) RETURN p",
                      "MATCH p=(:Genus {name: '" + tax.get('name') +
                      "'})-->(:Family)-->(:Order)-->(:Class)-->(:Phylum)-->(:Kingdom) RETURN p",
                      "MATCH p=(:Family {name: '" + tax.get('name') +
                      "'})-->(:Order)-->(:Class)-->(:Phylum)-->(:Kingdom) RETURN p",
                      "MATCH p=(:Order {name: '" + tax.get('name') +
                      "'})-->(:Class)-->(:Phylum)-->(:Kingdom) RETURN p",
                      "MATCH p=(:Class {name: '" + tax.get('name') +
                      "'})-->(:Phylum)-->(:Kingdom) RETURN p",
                      "MATCH p=(:Phylum {name: '" + tax.get('name') +
                      "'})-->(:Kingdom) RETURN p",
                      "MATCH p=(:Kingdom {name: '" + tax.get('name') + "'}) RETURN p"]
        query = query_list[level_id]
        tree = tx.run(query).data()[0]['p']
        for i in range(7-level_id):
            rel = rel_list[i+level_id]
            tx.run(("MATCH (a:Agglom_Taxon),(b:" + tax_list[i+level_id] + ") "
                    "WHERE a.name = '" + node + "' AND b.name = '" +
                    tree.nodes[i].get('name') + "' CREATE (a)-[r:" + rel + "]->(b) RETURN type(r)"))

    @staticmethod
    def _create_association(tx, agglom_1, agglom_2, networks, weight, mode):
        """Creates new associations between agglomerated nodes, with
        the appropriate weight and Network node connections."""
        uid = str(uuid4())
        # non alphanumeric chars break networkx
        if mode is 'weight':
            tx.run("CREATE (a:Association) SET a.name = $id SET a.weight = $weight", id=uid, weight=weight)
        else:
            tx.run("CREATE (a:Association) SET a.name = $id", id=uid)
        tx.run(("MATCH (a:Association),(b:Agglom_Taxon) "
                "WHERE a.name = '" + uid + "' AND b.name = '" +
                agglom_1 + "' CREATE (a)-[r:WITH_TAXON]->(b) RETURN type(r)"))
        tx.run(("MATCH (a:Association),(b:Agglom_Taxon) "
                "WHERE a.name = '" + uid + "' AND b.name = '" +
                agglom_2 + "' CREATE (a)-[r:WITH_TAXON]->(b) RETURN type(r)"))
        for node in networks:
            tx.run(("MATCH (a:Association),(b:Network) "
                    "WHERE a.name = '" + uid + "' AND b.name = '" +
                    node.get('name') + "' CREATE (a)-[r:IN_NETWORK]->(b) RETURN type(r)"))

    @staticmethod
    def _get_network(tx, nodes):
        """When a new association is generated to replace two old ones,
        all Network nodes those were connected to are returned by this function."""
        networks = list()
        for node in nodes:
            network = tx.run("MATCH (:Association {name: '" + node.get('name') +
                             "'})-->(n:Network) RETURN n").data()[0]['n']
            networks.append(network)
        return networks

    @staticmethod
    def _get_weight(tx, nodes):
        """Returns the weight of an Association node."""
        weight = tx.run("MATCH (n:Association {name: '" + nodes[0].get('name') +
                        "'}) RETURN n").data()[0]['n'].get('weight')
        return weight

    @staticmethod
    def _delete_old_associations(tx, associations):
        """Deletes specific associations and their relationships."""
        for node in associations:
            tx.run(("MATCH (n {name: '" + node.get('name') + "'}) DETACH DELETE n"))

    @staticmethod
    def _delete_old_agglomerations(tx, nodes):
        """Deletes old Agglom_Taxon nodes and their relationships."""
        for node in nodes:
            result = tx.run(("MATCH (n:Agglom_Taxon {name: '" + node.get('name') + "'}) RETURN n")).data()
            if len(result) > 0:
                tx.run(("MATCH (n:Agglom_Taxon {name: '" + node.get('name') + "'}) DETACH DELETE n"))

    @staticmethod
    def _delete_all(tx):
        """Deletes all nodes and their relationships from the database."""
        tx.run("MATCH (n) DETACH DELETE n")

    @staticmethod
    def _get_union(tx, networks):
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

    @staticmethod
    def _hypergeom_population(tx, taxon, categ, mode):
        """
        Returns 4 numbers:
        The number of samples in the database that is linked to the specified type,
        the number of samples in the database that is linked to a success,
        and the same values for the number of samples linked to the taxon.
        Only presence / absence is tested for, not differential abundance.
        """
        type_val = categ[0]
        success = categ[1]
        hypergeom_vals = dict()
        query = "MATCH (n:Sample)-->(:SampleProperty {type: '" + type_val + \
                "'}) RETURN n"
        total_samples = tx.run(query).data()
        hypergeom_vals['total_pop'] = _get_unique(total_samples, 'n', 'num')
        query = "MATCH (n:Sample)-->(:SampleProperty {type: '" + type_val + \
                "', name: '" + success + "'}) RETURN n"
        total_samples = tx.run(query).data()
        hypergeom_vals['success_pop'] = _get_unique(total_samples, 'n', 'num')
        if mode is 'Taxon':
            query = "MATCH (:Taxon {name: '" + taxon +\
                    "'})-->(n:Sample)-->(:SampleProperty {type: '" + type_val + \
                    "'}) RETURN n"
            total_samples = tx.run(query).data()
            hypergeom_vals['total_taxon'] = _get_unique(total_samples, 'n', 'num')
        if mode is 'Agglom_Taxon':
            query = "MATCH (:Agglom_Taxon {name: '" + taxon +\
                    "'})-[:GENERATED_FROM]-(:Taxon)--(n:Sample)-->" \
                    "(:SampleProperty {type: '" + type_val + \
                    "'}) RETURN n"
            total_samples = tx.run(query).data()
            hypergeom_vals['total_taxon'] = _get_unique(total_samples, 'n', 'num')
        if mode is 'Taxon':
            query = "MATCH (:Taxon {name: '" + taxon +\
                    "'})-->(n:Sample)-->(:SampleProperty {type: '" + type_val + \
                    "', name: '" + success + "'}) RETURN n"
            total_samples = tx.run(query).data()
            hypergeom_vals['success_taxon'] = _get_unique(total_samples, 'n', 'num')
        if mode is 'Agglom_Taxon':
            query = "MATCH (:Agglom_Taxon {name: '" + taxon +\
                    "'})-[:GENERATED_FROM]-(:Taxon)-->(n:Sample)-->" \
                    "(:SampleProperty {type: '" + type_val + \
                    "', name: '" + success + "'}) RETURN n"
            total_samples = tx.run(query).data()
            hypergeom_vals['success_taxon'] = _get_unique(total_samples, 'n', 'num')
        return hypergeom_vals

    @staticmethod
    def _spearman_test(tx, taxon, type_val, mode):
        """
        Returns p-value of Spearman correlation.
        """
        # get vector of sample values
        sample_values = list()
        sample_names = list()
        taxon_values = list()
        query = "MATCH (n:Sample)-->(:SampleProperty {type: '" + type_val + \
                "'}) RETURN n"
        samples = _get_unique(tx.run(query).data(), 'n')
        for item in samples:
            query = "MATCH (:Sample {name: '" + item + \
                    "'})-->(n:SampleProperty {type: '" + type_val + \
                    "'}) RETURN n"
            sample_value = tx.run(query).data()[0]['n'].get('name')
            try:
                sample_value = float(sample_value)
            except ValueError:
                pass
            if type(sample_value) == float:
                sample_values.append(sample_value)
                sample_names.append(item)
        for sample in sample_names:
            if mode is 'Taxon':
                query = "MATCH (:Sample {name: '" + sample + \
                        "'})<-[r:FOUND_IN]-(:Taxon {name: '" + taxon + \
                        "'}) RETURN r"
                counts = tx.run(query).data()
                if len(counts) == 0:
                    count = 0
                else:
                    count = float(counts[0]['r'].get('count'))
            if mode is 'Agglom_Taxon':
                query = "MATCH (:Sample {name: '" + sample + \
                        "'})<-[r:FOUND_IN]-(:Taxon)-[:GENERATED_FROM]-" \
                        "(:Agglom_Taxon {name: '" + taxon + \
                        "'}) RETURN r"
                counts = tx.run(query).data()
                if len(counts) == 0:
                    count = 0
                else:
                    count = 0
                    for item in counts:
                        count += float(item['r'].get('count'))
            taxon_values.append(count)
        result = spearmanr(taxon_values, sample_values)
        return result

    @staticmethod
    def _shortcut_categorical(tx, taxon, categ, mode, prob):
        """Creates relationship between categorical variable and taxon."""
        tx.run(("MATCH (a:" + mode +
                "),(b:SampleProperty) "
                "WHERE a.name = '" + taxon +
                "' AND b.name = '" + categ[1] +
                "' AND b.type = '" + categ[0] +
                "' CREATE (a)-[r:HYPERGEOM]->(b) "
                "SET r.correlation = '" + str(prob) +
                "' RETURN type(r)"))

    @staticmethod
    def _shortcut_continuous(tx, taxon, cont_var, mode):
        """Creates relationship between categorical variable and taxon."""
        var_id = list(cont_var.keys())[0]
        tx.run(("CREATE (a:SampleProperty {name: '" + var_id +
                "', correlation: '" + str(cont_var[var_id]) + "'}) RETURN a"))
        tx.run(("MATCH (a:" + mode +
                "),(b:SampleProperty) "
                "WHERE a.name = '" + taxon +
                "' AND b.name = '" + var_id +
                "' CREATE (a)-[r:SPEARMAN]->(b) "
                "SET r.correlation = '" + str(cont_var[var_id]) +
                "' RETURN type(r)"))


def _get_unique(node_list, key, mode=None):
    """Returns number or names of unique nodes in a list."""
    unique_samples = list()
    for item in node_list:
        unique_samples.append(item[key].get('name'))
    unique_samples = set(unique_samples)
    if mode == 'num':
        unique_samples = len(unique_samples)
    return unique_samples