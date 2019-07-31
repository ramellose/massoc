"""
This function accepts a Batch object and
runs network inference on the biom files in the object.
The output of the function is a Network object,
which contains inferred networks for each of the keys
in the biom dictionary.

The settings of the network inference tools are parsed
from txt documents, and not input directly.
This is because the number of tools and their wide range
of settings would make analysis irreproducible,
as well as the command line interface obfuscated.

Default settings are included in the project.
New settings can be added by copying these files,
and providing their filenames.

Networks are output as NetworkX graph files,
and can be exported to Cytoscape-compatible graphml files.

At the moment, only one location file per tool can be accepted.
May be adjusted in the future if it is really desirable to run
with multiple settings at once.

Parameters
----------

inputs : dictionary
    dictionary of parameters
    inputs.['location']: File location for batch object and intermediate results.
    inputs.['tools']: List of tools to use for network inference.
    inputs.['settings']: List of names of the setting files.
    inputs.['save']: Saves intermediate results to the specified location.

"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import ast
import csv
import statistics
import massoc
import biom
import networkx as nx
import pandas
import sys
from copy import deepcopy
from massoc.scripts.batch import Batch
import multiprocessing as mp
from functools import partial
from subprocess import call
import os
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


class Nets(Batch):

    """Container for multiple network files.
    The container is initialized from one or more
    BIOM files, and constructs networks
    based on the provided settings.
    The Nets object logs all steps of network
    construction, and provides options for processing of
    the networks.

    Parameters
    ----------
    otu : dict
        Dictionary of biom OTU files in hdf5 format
    genus : dict, optional
        Dictionary of genus-agglomerated biom files
    family : dict, optional
        Dictionary of family-agglomerated biom files
    order : dict, optional
        Dictionary of order-agglomerated biom files
    class : dict, optional
        Dictionary of class-agglomerated biom files
    phylum : dict, optional
        Dictionary of phylum-agglomerated biom files

    networks: dict
        Dictionary of network files

    """

    def __init__(self, batch=None):
        """
        Initialization function for Nets object.
        This object can inherit values from the Batch object.

        :param batch: Batch object.
        """
        super(Nets, self).__init__()
        if batch:
            self.otu = batch.otu
            self.species = batch.species
            self.genus = batch.genus
            self.family = batch.family
            self.order = batch.order
            self.class_ = batch.class_
            self.phylum = batch.phylum
            self.inputs = batch.inputs
        else:
            otu = dict()
            for value in self.inputs['biom']:
                otutab = biom.load_table(value)
                otu[value] = otutab
            self.otu = otu
        self.networks = dict()
        if self.inputs:
            _create_logger(self.inputs['fp'])
        if type(self.otu) is not dict:
            logger.error("Please supply a dictionary of biom files. ", exc_info=True)
            raise ValueError("Please supply a dictionary of biom files.")


    def write_networks(self):
        """
        Writes all networks in a Nets file to graphml files.

        :return:
        """
        try:
            for network in self.networks:
                path = self.inputs['fp'] + '/' + network + '.txt'
                nx.write_weighted_edgelist(G=self.networks[network], path=path)
        except Exception:
            logger.error("Unable to write networks to disk. ", exc_info=True)

    def add_networks(self, filename):
        """
        In case users want to manually import a network,
        this function adds the network file to the Nets object and checks
        whether the identifiers specified in the file match those in included BIOM files.
        Currently, only edge lists are supported.

        :param filename: Filename with network object.
        :return:
        """
        network = nx.read_weighted_edgelist(filename)
        try:
            network_nodes = list(network.nodes)
        except TypeError:
            logger.error("Unable to read edge list. ", exc_info=True)
        taxon_ids = list()
        biomfiles = [self.otu, self.species, self.genus,
                     self.family, self.order, self.class_, self.phylum]
        biomlist = list()
        for subset in biomfiles:
            for file in subset:
                biomlist.append(subset[file])
        for biomfile in biomlist:
            taxon_ids.extend(list(biomfile.ids(axis='observation')))
        missing_node = any(x not in taxon_ids for x in network_nodes)
        if missing_node:
            logger.error("Imported network node not found in taxon identifiers. ", exc_info=True)
        else:
            self.networks[filename] = network

    def _prepare_conet(self):
        """
        Carries out initial work before actually running CoNet.
        The initial writing function cannot be carried out
        in a multiprocessing operation because the Biom object cannot be pickled.
        However, the bash calls can be pickled; therefore, initial data prep
        is done first, then the CoNet calls are in parallel.

        :return:
        """
        filenames = self.get_filenames()
        ids = dict()
        obs_ids = dict()
        for x in filenames:
            ids[x] = dict()
            obs_ids[x] = dict()
            for y in filenames[x]:
                tempname = filenames[x][y][:-5] + '_counts_conet.txt'
                file = biom.load_table(filenames[x][y])
                obs_ids[x][y] = deepcopy(file._observation_ids)
                # code below is necessary to fix an issue where CoNet cannot read numerical OTU ids
                orig_ids = dict()
                for i in range(len(file.ids(axis='observation'))):
                    id = file.ids(axis='observation')[i]
                    orig_ids[("otu-" + str(i))] = id
                    file.ids(axis='observation')[i] = "otu_" + str(i)
                otu = file.to_tsv()
                text_file = open(tempname, 'w')
                text_file.write(otu[34:])
                text_file.close()
                ids[x][y] = orig_ids
        return ids, obs_ids

    def _prepare_spar(self):
        """
        Carries out initial work before actually running SparCC.
        The initial writing function cannot be carried out
        in a multiprocessing operation because the Biom object cannot be pickled.
        However, the bash calls can be pickled; therefore, initial data prep
        is done first, then the SparCC calls are in parallel.

        :return:
        """
        filenames = self.get_filenames()
        for x in filenames:
            for y in filenames[x]:
                file = biom.load_table(filenames[x][y])
                otu = file.to_tsv()
                tempname = filenames[x][y][:-5] + '_otus_sparcc.txt'
                text_file = open(tempname, 'w')
                text_file.write(otu[29:])
                text_file.close()

def _add_tax(network, file):
    """
    Adds taxon names from filename.

    :param network: NetworkX object
    :param file: File with taxonomy
    :return: Taxonomically annotated network
    """
    file = biom.load_table(file)
    tax = file._observation_metadata
    try:
        if tax is not None:
            for species in tax:
                species.pop('Genus (Aggregated)', None)
                species.pop('collapsed_ids', None)
            tax = file.metadata_to_dataframe('observation')
            taxnames = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
            taxdict = {}
            i = 0
            for name in tax.columns:
                taxdict[name] = taxnames[i]
                i = i + 1
            tax = tax.rename(index=str, columns=taxdict)
            for column in tax.columns:
                taxdict = tax[column].to_dict()
                nx.set_node_attributes(network, values=taxdict, name=column)
    except Exception:
        logger.error("Unable to collect taxonomy for agglomerated files. ", exc_info=True)
    return network


def run_conet(filenames, conet, orig_ids, obs_ids, settings=None):
    """
    Runs a Bash script containing the CoNet Bash commands.
    Unfortunately, errors produced by CoNet cannot be caught,
    because the exit status of the script is 0 regardless
    of CoNet producing a network or not.
    conet = nets.inputs['conet']

    :param filenames: Location of BIOM files written to disk.
    :param conet: Location of CoNet folder.
    :param orig_ids: OTU ids before annotation
    :param obs_ids: OTU ids with forbidden characters removed
    :param settings: Dictionary containing settings for CoNet
    :return: CoNet networks as NetworkX objects
    """
    if settings:
        path = settings
        if path[-3:] != '.sh':
            logger.error("Please supply a .sh executable to run CoNet. ", exc_info=True)
            raise ValueError("Please supply a .sh executable to run CoNet.")

    else:
        path = resource_path('CoNet.sh')
    path = path.replace('\\', '/')
    libpath = conet + '\\lib\\CoNet.jar'
    libpath = libpath.replace('\\', '/')
    results = dict()
    for x in filenames:
        for y in filenames[x]:
            graphname = filenames[x][y][:-5] + '_conet.tsv'
            tempname = filenames[x][y][:-5] + '_counts_conet.txt'
            # Code below removed because taxonomic information is not necessary
            # tax = file._observation_metadata
            # for species in tax:
            #    species.pop('Genus (Aggregated)', None)
            #    species.pop('collapsed_ids', None)
            # tax = file.metadata_to_dataframe('observation')
            # num = tax.shape[1]
            # for i in range(num, 7):
            #    level = 'taxonomy_' + str(i)
            #    tax[level] = 'Merged'
            #tax.to_csv(taxname, sep="\t")
            #f = open(taxname, 'r')
            #lines = f.read()
            #f.close()
            #lines = 'id' + lines
            #f = open(taxname, 'w')
            #f.write(lines)
            #f.close()
            # solving issue where guessingparam is higher than maximum edge number
            n_otus = len(orig_ids[x][y])
            guessingparam = str(n_otus * n_otus -1)
            if int(guessingparam) > 1000:
                guessingparam = str(1000)
            cmd = path + ' ' + tempname + ' ' + ' ' + graphname + ' ' + libpath + \
                  ' ' + resource_path("") + str(x) + '_' + str(y) + ' ' + guessingparam
            call(cmd, shell=True)
            call("rm " + tempname, shell=True)
            call("rm " + resource_path("") + str(x) + '_' + str(y) + "_threshold", shell=True)
            call("rm " + resource_path("") + str(x) + '_' + str(y) + "_permnet", shell=True)
            try:
                with open(graphname, 'r') as fin:
                    data = fin.read().splitlines(True)
                    fin.close()
                with open(graphname, 'w') as fout:
                    fout.writelines(data[2:])
                    fout.close()
            except FileNotFoundError:
                logger.error("Warning: CoNet did not complete network inference on: " + str(x) + "_" + str(y) + ' ', exc_info=True)
            signs = [b[0] for b in csv.reader(open(graphname, 'r'), delimiter='\t')]
            signs = [word.replace('mutualExclusion', '-1') for word in signs]
            signs = [word.replace('copresence', '1') for word in signs]
            signs = [word.replace('unknown', 'None') for word in signs]
            signs = [ast.literal_eval(b) for b in signs]
            clean_signs = list()  # None values need to be removed to make sure median is not 0.5.
            for sublist in signs:
                cleaned = [elem for elem in sublist if elem is not None]
                clean_signs.append(cleaned)
            signs = [statistics.median(x) for x in clean_signs]
            # methods = [x[2] for x in csv.reader(open(graphname, 'r'), delimiter='\t')]
            names = [b[15] for b in csv.reader(open(graphname, 'r'), delimiter='\t')]
            names = [b.split('->') for b in names]
            new_names = list()
            for item in names:
                new_item = [y.replace(y, orig_ids[x][y][b]) for b in item]
                new_names.append(new_item)
            i = 0
            adj = pandas.DataFrame(index=obs_ids[x][y], columns=obs_ids[x][y])
            adj = adj.fillna(0)
            for name in new_names:
                id1 = adj.columns.get_loc(name[0])
                id2 = adj.columns.get_loc(name[1])
                sign = signs[i]
                i = i+1
                adj.iloc[id1, id2] = sign
                adj.iloc[id2, id1] = sign
            net = nx.from_pandas_adjacency(adj)
            net = _add_tax(net, filenames[x][y])
            results[("conet_" + x + "_" + y)] = net
            call("rm " + graphname, shell=True)
    return results


def run_spiec(filenames, settings=None):
    """
    Runs a R executable containing settings for SPIEC-EASI network inference.

    :param filenames: Location of BIOM files written to disk.
    :param settings: Dictionary containing settings for SPIEC-EASI
    :return: SPIEC-EASI networks as NetworkX objects
    """
    results = dict()
    if settings:
        path = settings
        if path[-2:] != '.R':
            logger.error("Please supply an R executable to run SPIEC-EASI. ", exc_info=True)
            raise ValueError("Please supply an R executable to run SPIEC-EASI.")
    else:
        path = resource_path('spieceasi.r')
    path = path.replace('\\', '/')
    for x in filenames:
        for y in filenames[x]:
            graphname = filenames[x][y][:-5] + '_spiec'
            cmd = "Rscript " + path + " -i " + filenames[x][y] + " -o " + graphname
            call(cmd, shell=True)
            try:
                corrtab = pandas.read_csv(graphname, sep='\t', index_col=0)
            except FileNotFoundError:
                logger.error("Warning: SPIEC-EASI did not complete network inference. " + str(x) + "_" + str(y) + ' ', exc_info=True)
                exit(1)
            corrtab.columns = corrtab.index
            corrtab[corrtab > 0] = 1
            corrtab[corrtab < 0] = -1
            net = nx.from_pandas_adjacency(corrtab)
            net = _add_tax(net, filenames[x][y])
            results[("spiec-easi_" + x + "_" + y)] = net
            call("rm " + graphname, shell=True)
    return results


def run_spar(filenames, spar, boots=100, pval_threshold=0.001):
    """
    Runs python 2.7 SparCC code.
    spar = nets.inputs['spar'][0]

    :param filenames: Location of BIOM files written to disk.
    :param spar: Location of SparCC Python code
    :param boots: Number of bootstraps
    :param pval_threshold: p-value threshold for SparCC
    :return: SparCC networks as NetworkX objects
    """
    path = list()
    path.append(spar + '\\SparCC.py')
    path.append(spar + '\\MakeBootstraps.py')
    path.append(spar + '\\PseudoPvals.py')
    path = [x.replace('\\', '/') for x in path]
    results = dict()
    for x in filenames:
        for y in filenames[x]:
            tempname = filenames[x][y][:-5] + '_otus_sparcc.txt'
            corrs = filenames[x][y][:-5] + '_spar_corrs.tsv'
            cov = filenames[x][y][:-5] + '_spar_cov.tsv'
            pvals = filenames[x][y][:-5] + '_spar_pvals.tsv'
            bootstraps = filenames[x][y][:-(5 + len(x))] + 'bootstraps'
            cmd = "python2 " + path[0] + " " + tempname + " -i 5 " +\
                  " --cor_file " + corrs + " --cov_file " + cov
            call(cmd, shell=True)
            call("mkdir " + bootstraps, shell=True)
            n_bootstraps = str(boots)
            cmd = "python2 " + path[1] + " " + tempname + " -n " + n_bootstraps + \
                  " -t /permutation_#.txt -p " + bootstraps
            call(cmd, shell=True)
            for i in range(0, int(n_bootstraps)):
                permpath = bootstraps + '/permutation_' + str(i) + '.txt'
                pvalpath = bootstraps + '/perm_cor_' + str(i) + '.txt'
                cmd = "python2 " + path[0] + " " + permpath + " -i 5 " + \
                      " --cor_file " + pvalpath + " --cov_file " + cov
                call(cmd, shell=True)
            cmd = "python2 " + path[2] + ' ' + corrs + ' ' + bootstraps + \
                  '/perm_cor_#.txt 5 -o ' + pvals + ' -t two_sided'
            call(cmd, shell=True)
            call("rm -rf " + bootstraps, shell=True)
            call("rm " + tempname, shell=True)
            call("rm " + cov, shell=True)
            try:
                corrtab = pandas.read_csv(corrs, sep='\t', index_col=0)
            except FileNotFoundError:
                logger.error("Warning: SparCC did not complete network inference. " + str(x) + "_" + str(y) + ' ', exc_info=True)
                exit(1)
            corrtab.columns = corrtab.index
            pvaltab = pandas.read_csv(pvals, sep='\t', index_col=0)
            pvaltab = pvaltab < pval_threshold  # p value threshold for SparCC pseudo p-values
            corrtab = corrtab.where(pvaltab)
            corrtab = corrtab.fillna(0)
            corrtab[corrtab > 0] = 1
            corrtab[corrtab < 0] = -1
            net = nx.from_pandas_adjacency(corrtab)
            net = _add_tax(net, filenames[x][y])
            results[("sparcc_" + x + "_" + y)] = net
            call("rm " + corrs + " " + pvals +
                 " " + os.path.dirname(massoc.__file__)[:-6] +
                 "\cov_mat_SparCC.out", shell=True)
    return results


def resource_path(relative_path):
    """
     Get absolute path to resource, works for dev and for PyInstaller.
     Source: https://stackoverflow.com/questions/7674790/bundling-data-files-with-pyinstaller-onefile

    :param relative_path: Path to MEI location
    :return:
    """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


def run_jobs(job, spar, conet, orig_ids, obs_ids, filenames,
             spiec_settings=None, conet_settings=None):
    """
    Accepts a job from a joblist to run network inference in parallel.

    :param job: dictionary of dictionary with tools as keys and taxonomic levels as second-layer keys
    :param spar: Location of SparCC folder
    :param conet: Location of CoNet folder
    :param orig_ids: Original OTU IDs
    :param obs_ids: OTU IDs with forbidden characters removed
    :param filenames: Locations of BIOM files
    :param spiec_settings: Location of alternative Rscript for SPIEC-EASI
    :param conet_settings: Location of alternative Bash script for CoNet
    :return: NetworkX networks
    """
    select_filenames = {job[0]: {job[2]: filenames[job[0]][job[2]]}}
    # only filenames with the same taxonomic level are included
    if 'spiec-easi' in job:
        logger.info('Running SPIEC-EASI... ')
        networks = run_spiec(select_filenames, settings=spiec_settings)
    if 'sparcc' in job:
        logger.info('Running SparCC... ')
        if 'spar_setting' in job['sparcc']:
            if len(job['spar_setting'][1]) == 2:
                networks = run_spar(spar=spar, filenames=select_filenames,
                                    boots=job['spar_setting'][1]['spar_boot'],
                                    pval_threshold=job['spar_setting'][1]['spar_pval'])
            else:
                if 'spar_boot' in job['spar_setting'][1]:
                    networks = run_spar(spar=spar, filenames=select_filenames,
                                        boots=job['spar_setting'][1]['spar_boot'])
                if 'spar_pval' in job['spar_setting'][1]:
                    networks = run_spar(spar=spar, filenames=select_filenames,
                                        pval_threshold=job['spar_setting'][1]['spar_pval'])
        else:
            networks = run_spar(spar=spar, filenames=select_filenames)
    if 'conet' in job:
        logger.info('Running CoNet... ')
        networks = run_conet(conet=conet, filenames=select_filenames,
                             orig_ids=orig_ids, obs_ids=obs_ids, settings=conet_settings)
    return networks


def get_joblist(nets):
    """
    Creates a list of jobs that can be distributed over multiple processes.
    Note: should be appended to handle multiple taxonomic levels + files!
    Each job is a tuple of the taxonomic level, tool and name.

    :param nets: Nets object
    :return: Dictionary of dictionary of jobs
    """
    joblist = list()
    for name in nets.inputs['name']:
        for level in nets.inputs['levels']:
            sublist = dict()
            if nets.inputs['tools']:
                for i in nets.inputs['tools']:
                    sublist[i] = level
            if nets.inputs['spiec'] is not None:
                sublist['spiec_setting'] = [level, nets.inputs['spiec']]
            if nets.inputs['spar_boot'] or nets.inputs['spar_pval'] is not None:
                sublist['spar_setting'] = [level, {'spar_boot': nets.inputs['spar_boot'],
                'spar_pval': nets.inputs['spar_pval']}]
            for value in sublist:
                joblist.append((sublist[value], value, name))
    return joblist


def run_parallel(nets):
    """
    Creates partial function to run as pool.

    :param nets: Nets object
    :return:
    """
    cores = nets.inputs['cores']
    jobs = get_joblist(nets)
    filenames = nets.inputs['procbioms']
    logger.info('Collecting jobs... ')
    pool = mp.Pool(cores)
    # multiprocess supports passing objects
    # multiprocessing does not
    # however, multiprocess cannot be frozen
    # need to rewrite netwrap as pickle-able objects!
    orig_ids = None
    obs_ids = None
    if 'conet' in nets.inputs['tools']:
        orig_ids, obs_ids = nets._prepare_conet()
    if 'sparcc' in nets.inputs['tools']:
        nets._prepare_spar()
    func = partial(run_jobs, filenames=filenames, orig_ids=orig_ids,
                   obs_ids=obs_ids, spar=nets.inputs['spar'], conet=nets.inputs['conet'],
                   spiec_settings=nets.inputs['spiec'], conet_settings=nets.inputs['conet_bash'])
    try:
        logger.info('Distributing jobs... ')
        # network_list = list()
        # for job in jobs:
            # result = run_jobs(nets, job)
            # network_list.append(result)
        results = pool.map(func, iter(jobs))
    except Exception:
        logger.error('Failed to generate workers. ', exc_info=True)
    for item in results:
        for network in item:
            nets.networks[network] = item[network]
    # for i in range(1, len(jobs)):
    #    nets.networks = {**nets.networks, **results[i]}
    # clean up old written BIOM files
    logger.info('Completed tasks! ')
    return nets


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