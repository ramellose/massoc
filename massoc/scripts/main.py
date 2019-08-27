"""
This file contains parsers and functions that call on other functionality defined
in the rest of massoc's scripts directory.

The command line interface is intended to be called sequentially;
files are written to disk as intermediates,
while a settings file is used to transfer logs and other information
between the modules. These modules are contained in this file.
This modular design allows users to leave out parts of massoc that are not required,
and reduces the number of parameters that need to be defined in the function calls.

"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import sys
import os
from biom import load_table
from biom.parse import MetadataMap
from massoc.scripts.batch import Batch, write_settings, read_settings, read_bioms
from massoc.scripts.netwrap import Nets, run_parallel
from copy import deepcopy
from platform import system
from subprocess import Popen
from psutil import Process, pid_exists
from time import sleep
from massoc.scripts.netbase import ImportDriver
from massoc.scripts.netstats import NetDriver
from massoc.scripts.metastats import MetaDriver
import networkx as nx
from wx.lib.pubsub import pub
import logging
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


def get_input(inputs, publish=False):
    """
    Takes all input and returns a dictionary of biom files.
    If tab-delimited files are supplied, these are combined
    into a biom file. File names are used as keys.
    This is mostly a utility wrapper, as all biom-related functions
    are from biom-format.org.

    At the moment, rarefaction is performed after sample splitting.
    This means that samples with uneven sequence counts will not
    be rarefied to equal depths.

    All files are written to BIOM files, while a settings file is also written to disk
    for use by other massoc commands.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    # handler to file
    # construct logger after filepath is provided
    _create_logger(inputs['fp'])
    if inputs['biom_file'] is not None:
        logger.info('BIOM file(s) to process: ' + ", ".join(inputs['biom_file']))
    if inputs['otu_table'] is not None:
        logger.info('Tab-delimited OTU table(s) to process: ' + ", ".join(inputs['otu_table']))
    if inputs['tax_table'] is not None:
        if len(inputs['otu_table']) is not len(inputs['tax_table']):
            logger.error("Add a taxonomy table for every OTU table!", exc_info=True)
            raise ValueError("Add a taxonomy table for every OTU table!")
    if inputs['sample_data'] is not None:
        if len(inputs['otu_table']) is not len(inputs['sample_data']):
            logger.error("Add a sample data table for every OTU table!", exc_info=True)
            raise ValueError("Add a sample data table for every OTU table!")
    if inputs['otu_meta'] is not None:
        if len(inputs['otu_table']) is not len(inputs['otu_meta']):
            logger.error("Add a metadata table for every OTU table!", exc_info=True)
            raise ValueError("Add a metadata table for every OTU table!")
    filestore = {}
    if inputs['biom_file'] is None and inputs['network'] is None:
        if inputs['otu_table'] is None and inputs['network'] is None:
            logger.error("Please supply either a biom file"
                         ", a tab-delimited OTU table or a network!", exc_info=True)
            raise ValueError("Please supply either a biom file"
                             ", a tab-delimited OTU table or a network!")
    # Only process count files if present
    i = 0
    if inputs['name'] is None:
        inputs['name'] = list()
        inputs['name'].append('file_')
    if inputs['biom_file'] is not None:
        try:
            for x in inputs['biom_file']:
                biomtab = load_table(x)
                filestore[inputs['name'][i]] = biomtab
                i += 1
        except Exception:
            logger.error("Failed to import BIOM files.", exc_info=True)
    if inputs['otu_table'] is not None:
        try:
            j = 0  # j is used to match sample + tax data to OTU data
            for x in inputs['otu_table']:
                input_fp = x
                sample_metadata_fp = None
                observation_metadata_fp = None
                obs_data = None
                sample_data = None
                biomtab = load_table(input_fp)
                try:
                    sample_metadata_fp = inputs['sample_data'][j]
                    observation_metadata_fp = inputs['tax_table'][j]
                except TypeError or KeyError:
                    pass
                if sample_metadata_fp is not None:
                    sample_f = open(sample_metadata_fp, 'r')
                    sample_data = MetadataMap.from_file(sample_f)
                    sample_f.close()
                    biomtab.add_metadata(sample_data, axis='sample')
                if observation_metadata_fp is not None:
                    obs_f = open(observation_metadata_fp, 'r')
                    obs_data = MetadataMap.from_file(obs_f)
                    obs_f.close()
                    # for taxonomy collapsing,
                    # metadata variable needs to be a complete list
                    # not separate entries for each tax level
                    for b in list(obs_data):
                        tax = list()
                        for l in list(obs_data[b]):
                            tax.append(obs_data[b][l])
                            obs_data[b].pop(l, None)
                        obs_data[b]['taxonomy'] = tax
                    biomtab.add_metadata(obs_data, axis='observation')
                filestore[inputs['name'][j]] = biomtab
                j += 1
        except Exception:
            logger.warning("Failed to combine input files.", exc_info=True)
    bioms = Batch({'otu': filestore}, inputs)
    # it is possible that there are forbidden characters in the OTU identifiers
    # we can forbid people from using those, or replace those with an underscore
    if inputs['biom_file'] or inputs['otu_table']:
        for name in bioms.otu:
            biomfile = bioms.otu[name]
            taxon_ids = biomfile._observation_ids  # need to be careful with these operations
            taxon_index = biomfile._obs_index      # likely to corrupt BIOM file if done wrong
            new_ids = deepcopy(taxon_ids)
            new_indexes = deepcopy(taxon_index)
            for i in range(0, len(taxon_ids)):
                id = taxon_ids[i]
                new_id = id.replace(" ", "_")
                new_ids[i] = new_id
                new_indexes[new_id] = new_indexes.pop(id)
            biomfile._observation_ids = new_ids
            biomfile._obs_index = new_indexes
            bioms.otu[name] = biomfile
        logger.info('Collapsing taxonomy... ')
        bioms.collapse_tax()
        if inputs['cluster'] is not None:
            if publish:
                pub.sendMessage('update', msg='Clustering BIOM files...')
            logger.info('Clustering BIOM files... ')
            bioms.cluster_biom()
        if inputs['split'] is not None and inputs['split'] is not 'TRUE':
            bioms.split_biom()
        if inputs['min'] is not None:
            if publish:
                pub.sendMessage('update', msg='Setting minimum mean abundance...')
            logger.info('Removing taxa below minimum count... ')
            bioms.prev_filter(mode='min')
        if inputs['prev'] is not None:
            if publish:
                pub.sendMessage('update', msg='Setting prevalence filter...')
            logger.info('Setting prevalence filter... ')
            bioms.prev_filter(mode='prev')
        if inputs['rar'] is not None:
            if publish:
                pub.sendMessage('update', msg='Rarefying counts...')
            logger.info('Rarefying counts... ')
            bioms.rarefy()
    bioms.inputs['procbioms'] = dict()
    if inputs['biom_file'] or inputs['otu_table']:
        if 'otu' not in bioms.inputs['levels']: # add otu level always
            bioms.inputs['procbioms']['otu'] = dict()
            for name in bioms.inputs['name']:
                biomname = bioms.inputs['fp'] + '/' + name + '_' + 'otu' + '.hdf5'
                bioms.inputs['procbioms']['otu'][name] = biomname
        for level in bioms.inputs['levels']:
            bioms.inputs['procbioms'][level] = dict()
            for name in bioms.inputs['name']:
                biomname = bioms.inputs['fp'] + '/' + name + '_' + level + '.hdf5'
                bioms.inputs['procbioms'][level][name] = biomname
        all_bioms = {**bioms.otu, **bioms.genus, **bioms.family, **bioms.order,
                     **bioms.class_, **bioms.phylum}
        for biomfile in all_bioms:
            if all_bioms[biomfile].shape[0] == 1:
                logger.error("The current preprocessing steps resulted in BIOM files with only 1 row.", exc_info=True)
    if inputs['network'] is not None:
        if publish:
            pub.sendMessage('update', msg='Checking previously generated networks...')
        logger.info('Checking previously generated networks...')
        filelist = deepcopy(inputs['network'])
        for file in filelist:
            network = _read_network(file)
            nodes = len(network.nodes)
            edges = len(network.edges)
            logger.info("This network has " + str(nodes) + \
                           " nodes and " + str(edges) + " edges.")
            weight = nx.get_edge_attributes(network, 'weight')
            if len(weight) > 0:
                logger.info('This is a weighted network.')
            else:
                logger.info('This is an unweighted network.')
    try:
        if inputs['biom_file'] or inputs['otu_table']:
            bioms.write_bioms()
            logger.info('BIOM files written to disk.  ')
    except Exception:
        logger.warning('Failed to write BIOM files to disk.  ', exc_info=True)
    write_settings(bioms.inputs)
    logger.info('Settings file written to disk.  ')


def run_network(inputs, publish=False):
    """
    Pipes functions from the different massoc modules to run complete network inference.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    _create_logger(inputs['fp'])
    old_inputs = read_settings(inputs['fp'] + '/settings.json')
    old_inputs.update(inputs)
    inputs = old_inputs
    # handler to file
    filestore = read_bioms(inputs['procbioms'])
    bioms = Batch(filestore, inputs)
    bioms = Nets(bioms)
    if inputs['tools'] is not None:
        logger.info('Tools to run with default settings: ' + str(inputs['tools']) + ' ')
    bioms.inputs['network'] = list()
    network_names = list()
    for tool in bioms.inputs['tools']:
        for level in bioms.inputs['levels']:
            for name in bioms.inputs['name']:
                filename = bioms.inputs['fp'] + '/' + tool + '_' + level + '_' + name + '.txt'
                network_names.append(filename)
    bioms.inputs['network'] = network_names
    if publish:
        pub.sendMessage('update', msg='Starting network inference. This may take some time!')
    try:
        logger.info('Running network inference...  ')
        networks = run_parallel(bioms)
        networks.write_networks()
    except Exception:
        logger.warning('Failed to complete network inference.  ', exc_info=True)
    write_settings(networks.inputs)
    if publish:
        pub.sendMessage('update', msg="Finished running network inference!")
    logger.info('Finished running network inference.  ')


def run_neo4j(inputs, publish=False):
    """
    Starts and carries out operations on the Neo4j database.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    _create_logger(inputs['fp'])
    # overwritten settings should be retained
    old_inputs = read_settings(inputs['fp'] + '/settings.json')
    # handler to file
    # check if password etc is already there
    if 'username' in old_inputs:
        logins = dict((k, old_inputs[k]) for k in ('username', 'password', 'address', 'neo4j'))
    old_inputs.update(inputs)
    inputs = old_inputs
    if 'pid' in inputs:
        existing_pid = pid_exists(inputs['pid'])
    else:
        existing_pid = False
    if not inputs['neo4j']:
        inputs.update(logins)
    checks = str()
    if inputs['job'] == 'start':
        if not existing_pid:
            start_database(inputs, publish)
            existing_pid = True
        else:
            logger.info("Database is already running.  ")
    elif inputs['job'] == 'quit':
        if not existing_pid:
            logger.info("No database open.  ")
        else:
            try:
                if publish:
                    pub.sendMessage('update', msg='Getting PID...')
                # there is a lingering Java process that places a lock on the database.
                # terminating the subprocess does NOT terminate the Java process,
                # so the store lock has to be deleted manually.
                # This is different for Linux & Windows machines and may not be trivial
                # however, PID solution may be platform-independent
                # CURRENT SOLUTION:
                # get parent PID of subprocess
                # use psutil to get child PIDs
                # kill child PIDs too
                parent_pid = inputs['pid']
                parent = Process(parent_pid)
                children = parent.children(recursive=True)
                for child in children:
                    child.kill()
                # apparently killing the children also kills the parent
            except Exception:
                logger.warning("Failed to close database.  ", exc_info=True)
    elif inputs['job'] == 'clear':
        if not existing_pid:
            start_database(inputs, publish)
            existing_pid = True
        try:
            if publish:
                pub.sendMessage('update', msg='Clearing database...')
            importdriver = ImportDriver(user=inputs['username'],
                                        password=inputs['password'],
                                        uri=inputs['address'], filepath=inputs['fp'])
            importdriver.clear_database()
            importdriver.close()
        except Exception:
            logger.warning("Failed to clear database.  ", exc_info=True)
    elif inputs['job'] == 'write':
        if not existing_pid:
            start_database(inputs, publish)
            existing_pid = True
        try:
            if publish:
                pub.sendMessage('update', msg='Accessing database...')
            importdriver = ImportDriver(user=inputs['username'],
                                        password=inputs['password'],
                                        uri=inputs['address'], filepath=inputs['fp'])
            importdriver.export_network(path=inputs['fp'])
            importdriver.close()
        except Exception:
            logger.warning("Failed to write database to graphml file.  ", exc_info=True)
    else:
        if not existing_pid:
            start_database(inputs, publish)
            existing_pid = True
        if publish:
            pub.sendMessage('update', msg='Uploading files to database...')
        filestore = None
        if inputs['procbioms']:
            filestore = read_bioms(inputs['procbioms'])
        # ask users for additional input
        bioms = Batch(filestore, inputs)
        bioms = Nets(bioms)
        for file in inputs['network']:
            network = _read_network(file)
            bioms.add_networks(network, file)
        importdriver = None
        sleep(12)
        importdriver = ImportDriver(user=inputs['username'],
                                    password=inputs['password'],
                                    uri=inputs['address'], filepath=inputs['fp'])
        # importdriver.clear_database()
        try:
            # pub.sendMessage('update', msg='Uploading BIOM files...')
            logger.info("Uploading BIOM files...")
            itemlist = list()
            for level in inputs['procbioms']:
                for item in inputs['procbioms'][level]:
                    name = inputs['procbioms'][level][item]
                    biomfile = load_table(name)
                    importdriver.convert_biom(biomfile=biomfile, exp_id=name)
                    itemlist.append(name)
            checks += 'Successfully uploaded the following items and networks to the database: \n'
            for item in itemlist:
                checks += (item + '\n')
            checks += '\n'
            logger.info(checks)
        except Exception:
            logger.warning("Failed to upload BIOM files to Neo4j database.  ", exc_info=True)
        try:
            # pub.sendMessage('update', msg='Uploading network files...')
            logger.info('Uploading network files...  ')
            for item in inputs['network']:
                network = nx.read_weighted_edgelist(item)
                # try to split filename to make a nicer network id
                subnames = item.split('/')
                if len(subnames) == 1:
                    subnames = item.split('\\')
                name = subnames[-1].split('.')[0]
                importdriver.convert_networkx(network=network, network_id=name, mode='weight')
                itemlist.append(item)
        except Exception:
            logger.warning('Unable to upload network files to Neo4j database. ', exc_info=True)
            checks += 'Unable to upload network files to Neo4j database.\n'
        if publish:
            pub.sendMessage('database_log', msg=checks)
        importdriver.close()
    logger.info('Completed database operations!  ')
    write_settings(inputs)


def run_netstats(inputs, publish=False):
    """
    Runs statistical analyses on the Neo4j database, as well as logic operations.
    To do: null models.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    old_inputs = read_settings(inputs['fp'] + '/settings.json')
    old_inputs.update(inputs)
    inputs = old_inputs
    # handler to file
    _create_logger(inputs['fp'])
    checks = str()
    if 'pid' in inputs:
        existing_pid = pid_exists(inputs['pid'])
    else:
        existing_pid = False
    if not existing_pid:
        start_database(inputs, publish)
        existing_pid = True
    try:
        if publish:
            pub.sendMessage('update', msg='Starting database drivers.')
            # sys.stdout.write('Starting database drivers.')
        netdriver = NetDriver(user=inputs['username'],
                              password=inputs['password'],
                              uri=inputs['address'], filepath=inputs['fp'])
        importdriver = ImportDriver(user=inputs['username'],
                                    password=inputs['password'],
                                    uri=inputs['address'], filepath=inputs['fp'])

    except Exception:
        logger.warning("Failed to start database worker.  ", exc_info=True)
    try:
        # write operations here
        if inputs['logic']:
            if not inputs['networks']:
                networks = list()
                hits = importdriver.custom_query("MATCH (n:Network) RETURN n")
                for hit in hits:
                    networks.append(hit['n'].get('name'))
            else:
                networks = inputs['networks']
            if 'union' in inputs['logic']:
                netdriver.graph_union(networks=networks)
            if 'intersection' in inputs['logic']:
                netdriver.graph_intersection(networks=networks,
                                             weight=inputs['weight'], n=inputs['num'])
            if 'difference' in inputs['logic']:
                netdriver.graph_difference(networks=networks,
                                           weight=inputs['weight'])
            checks += 'Logic operations completed. \n'
            if publish:
                pub.sendMessage('update', msg="Exporting network...")
                if inputs['networks'] is not None:
                    names = [x.split('.')[0] for x in inputs['networks']]
                    importdriver.export_network(path=inputs['fp'] + '/' +
                                                "_".join(names) + '.graphml')
                    logger.info("Exporting networks to: " + inputs['fp'] + '/' +
                                "_".join(names) + '.graphml')
                    checks += "Exporting networks to: " + inputs['fp'] + '/' +\
                              "_".join(names) + '.graphml' "\n"
                else:
                    importdriver.export_network(path=inputs['fp'] + '/' +
                                                      '_complete.graphml')
                    logger.info("Exporting networks to: " + inputs['fp'] + '/' +
                                '_complete.graphml')
                    checks += "Exporting networks to: " + inputs['fp'] + '/' +\
                              '_complete.graphml' "\n"
        else:
            logger.warning("No logic operation specified!")
        if publish:
            pub.sendMessage('update', msg="Completed database operations!")
        # sys.stdout.write("Completed database operations!")
        checks += 'Completed database operations! \n'
    except Exception:
        logger.warning("Failed to run database worker.  ", exc_info=True)
        checks += 'Failed to run database worker. \n'
    if publish:
        pub.sendMessage('database_log', msg=checks)
    importdriver.close()
    netdriver.close()
    logger.info('Completed netstats operations!  ')
    write_settings(inputs)


def run_metastats(inputs, publish=False):
    """
    Module that carries out analysis of metadata on the database.
    This module also interfaces with external APIs to pull in additional metadata.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    old_inputs = read_settings(inputs['fp'] + '/settings.json')
    old_inputs.update(inputs)
    inputs = old_inputs
    # handler to file
    _create_logger(inputs['fp'])
    checks = str()
    try:
        if publish:
            pub.sendMessage('update', msg='Starting database drivers.')
        # sys.stdout.write('Starting database drivers.')
        metadriver = MetaDriver(user=inputs['username'],
                                password=inputs['password'],
                                uri=inputs['address'],
                                filepath=inputs['fp'])
        importdriver = ImportDriver(user=inputs['username'],
                              password=inputs['password'],
                              uri=inputs['address'],
                              filepath=inputs['fp'])
    except Exception:
        logger.warning("Failed to start database worker.  ", exc_info=True)
    if inputs['sequence']:
        try:
            logger.info('Uploading sequences to database...')
            if publish:
                pub.sendMessage('update', msg='Uploading sequences to database...')
            importdriver.include_sequences(inputs['sequence'])
        except Exception:
            logger.warning("Failed to upload sequences to database.  ", exc_info=True)
    if inputs['add']:
        try:
            logger.info('Uploading additional properties...  ')
            if publish:
                pub.sendMessage('update', msg='Uploading files to database...')
            # create dictionary from file
            # first check if this is an abundance table
            for k in range(len(inputs['add'])):
                filepath = inputs['add'][k]
                with open(filepath, 'r') as file:
                    # Second column name is type
                    # Newline is cutoff
                    colnames = file.readline().split(sep="\t")
                    lines = file.readlines()[1:]
                    if not inputs['type']:
                        label = colnames[0].rstrip()
                    else:
                        label = inputs['type']
                    # if the supplied file is a dataframe,
                    # treat first column as source and rest as target
                    logger.info('Found ' + str(len(colnames)) + ' properties.')
                    for i in range(1, len(colnames)):
                        # give a logger update every 5th property
                        node_dict = dict()
                        name = colnames[i].rstrip()
                        if i % 5 == 0:
                            logger.info('Working on the ' + str(i) + 'th property.')
                        for line in lines:
                            source = line.split(sep="\t")[0].rstrip()
                            weight = None
                            if inputs['abundance']:
                                target = colnames[i].rstrip()
                                name = inputs['abundance'][k]
                                weight = line.split(sep="\t")[i].rstrip()
                            else:
                                target = line.split(sep="\t")[i].rstrip()
                            if weight != 0:
                                node_dict[source] = {'target': target, 'weight': weight}
                        importdriver.include_nodes(nodes=node_dict, name=name, label=label)
        except Exception:
            logger.warning("Failed to upload properties to database.  ", exc_info=True)
    inputs['add'] = None
    inputs['type'] = None
    # prevents reuploading
    try:
        # write operations here
        if inputs['agglom']:
            tax_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
            level_id = tax_list.index(inputs['agglom'].capitalize())
            if inputs['weight']:
                mode = inputs['weight']
            else:
                mode = 'Ignore weight'
            for level in range(0, level_id+1):
                # pub.sendMessage('update', msg="Agglomerating edges...")
                logger.info("Agglomerating edges...")
                metadriver.agglomerate_network(level=tax_list[level], mode=mode)
            checks += 'Successfully agglomerated edges. \n'
    except Exception:
        logger.warning("Failed to carry out edge agglomeration.  ", exc_info=True)
        checks += 'Failed to carry out edge agglomeration. \n'
    try:
        if inputs['variable']:
            logger.info("Associating samples...  ")
            pub.sendMessage('update', msg="Associating samples...")
            # sys.stdout.write("Associating samples...")
            if inputs['variable'][0] == 'all':
                properties = set([x[y] for x in metadriver.custom_query("MATCH (n:Property) RETURN n.type") for y in x])
                for prop in properties:
                    metadriver.associate_samples(label=prop)
            else:
                for var in inputs['variable']:
                    metadriver.associate_samples(label=var)
            checks += 'Completed associations. \n'
    except Exception:
        logger.warning("Failed to compute metadata associations.  ", exc_info=True)
        checks += 'Failed to compute metadata associations. \n'
    if publish:
        pub.sendMessage('database_log', msg=checks)
    # functions to include:
    # include_sequences
    metadriver.close()
    importdriver.close()
    logger.info('Completed metastats operations!  ')
    write_settings(inputs)


def resource_path(relative_path):
    """
    Get absolute path to resource, works for dev and for PyInstaller.
    Source: https://stackoverflow.com/questions/7674790/bundling-data-files-with-pyinstaller-onefile

    :param relative_path: Path to MEI location.
    :return:
    """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


def start_database(inputs, publish):
    """
    Starts Neo4j database.

    :param inputs: Dictionary of inputs.
    :param publish: If True, publishes messages to be received by GUI.
    :return:
    """
    try:
        if publish:
            pub.sendMessage('update', msg='Starting database...')
        if system() == 'Windows':
            filepath = inputs['neo4j'] + '/bin/neo4j.bat console'
        else:
            filepath = inputs['neo4j'] + '/bin/neo4j console'
        filepath = filepath.replace("\\", "/")
        if system() == 'Windows' or system() == 'Darwin':
            p = Popen(filepath, shell=True)
        else:
            # note: old version used gnome-terminal, worked with Ubuntu
            # new version uses xterm to work with macOS
            # check if this conflicts!
            p = Popen(["gnome-terminal", "-e", filepath])  # x-term compatible alternative terminal
        inputs['pid'] = p.pid
        if publish:
            pub.sendMessage('pid', msg=inputs['pid'])
        sleep(12)
    except Exception:
        logger.warning("Failed to start database.  ", exc_info=True)


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


def _read_network(filepath):
    """
    Imports network file according to extension.

    :param filepath: Network filepath
    :return: NetworkX object
    """
    filename = filepath.split(sep=".")
    extension = filename[len(filename) - 1]
    network = None
    try:
        if extension == 'graphml':
            network = nx.read_graphml(filepath)
        elif extension == 'txt':
            network = nx.read_weighted_edgelist(filepath)
        elif extension == 'gml':
            network = nx.read_gml(filepath)
        else:
            logger.warning('Format not accepted. '
                           'Please specify the filename including extension (e.g. test.graphml).', exc_info=True)
            exit()
    except Exception:
        logger.error('Could not import network file!', exc_info=True)
    return network
