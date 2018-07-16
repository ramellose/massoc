"""
These functions convert files in the original BIOM format into a Neo4J format,
and stores this as a local database.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import py2neo

# assumes a local database is running
from py2neo import Database
db = Database("http://localhost:7474")

def convert_neo(nets):
    """
    Converts Nets object from netwrap.py to a Neo4J graph.
    This graph can be stored more easily in a database,
    and supports future metadata-based utilities.
    """
    for name in nets.otu:
        biomfile = nets.otu[name]
