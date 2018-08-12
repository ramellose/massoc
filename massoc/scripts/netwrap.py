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
import os
import statistics
from datetime import datetime
from subprocess import call
import massoc
import biom
import networkx
import pandas
from copy import deepcopy
import massoc.execs
from massoc.scripts.batch import Batch
from massoc.scripts.main import resource_path
import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

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
    log : str
        Logs operations used on batch object
    networks: dict
        Dictionary of network files

    """

    def __init__(self, batch=None):
        super(Nets, self).__init__()
        if batch is not None:
            self.otu = batch.otu
            self.species = batch.species
            self.genus = batch.genus
            self.family = batch.family
            self.order = batch.order
            self.class_ = batch.class_
            self.phylum = batch.phylum
            self.inputs = batch.inputs
            self.log = batch.log
        else:
            otu = dict()
            self.log["Creation"] = "Creation: " + datetime.now().strftime('%B %d %Y %H:%M:%S') + "\n"
            for value in self.inputs['biom']:
                otutab = biom.load_table(value)
                otu[value] = otutab
            self.otu = otu
        self.networks = dict()
        self.names = list(self.otu)
        if type(self.otu) is not dict:
            raise ValueError("Please supply a dictionary of biom files.")


    def get_filenames(self):
        """
        If the standard processing functions in the Batch object are used,
        filenames follow a standard format. This function retrieves such filenames
        so BIOM files can be read without requiring them to be stated explicitly
        by the user.
        """
        try:
            filenames = {}
            for x in self.names:
                for level in self.inputs['levels']:
                    if level is 'class':
                        filenames[(x + '_') + level + '_'] = self.inputs['fp'][0] + '/' + x + '_' + level + '.hdf5'
                    else:
                        filenames[(x + '_') + level] = self.inputs['fp'][0] + '/' + x + '_' + level + '.hdf5'
        except TypeError:
            logger.error("Unable to generate filenames", exc_info=True)
        return filenames

    def run_spiec(self, settings=None):
        """
        Runs a R executable containing settings for SPIEC-EASI network inference.
        """
        if settings is None:
            path = resource_path('spieceasi.r')
            path = path.replace('\\', '/')
        else:
            path = settings[0]
            if path[-2:] != '.R':
                raise ValueError("Please supply an R executable to run SPIEC-EASI.")
        filenames = self.get_filenames()
        for x in filenames:
            try:
                graphname = filenames[x][:-5] + '_spiec'
                cmd = "Rscript " + path + " -i " + filenames[x] + " -o " + graphname
                call(cmd, shell=True)
                corrtab = pandas.read_csv(graphname, sep='\t', index_col=0)
                corrtab.columns = corrtab.index
                corrtab[corrtab > 0] = 1
                corrtab[corrtab < 0] = -1
                net = networkx.from_pandas_adjacency(corrtab)
                net = _add_tax(net, filenames[x])
                self.networks[("spiec-easi_" + x)] = net
                call("rm " + graphname, shell=True)
            except Exception:
                logger.error("Unable to finish SPIEC-EASI on " + str(x), exc_info=True)

    def run_spar(self, boots=100, pval_threshold=0.001):
        """
        Runs python 2.7 SparCC code.
        """
        path = list()
        path.append(self.inputs['spar'][0] + '\\SparCC.py')
        path.append(self.inputs['spar'][0] + '\\MakeBootstraps.py')
        path.append(self.inputs['spar'][0] + '\\PseudoPvals.py')
        path = [x.replace('\\', '/') for x in path]
        filenames = self.get_filenames()
        for x in filenames:
            try:
                file = biom.load_table(filenames[x])
                otu = file.to_tsv()
                tempname = filenames[x][:-5] + '_otus_sparcc.txt'
                text_file = open(tempname, 'w')
                text_file.write(otu[29:])
                text_file.close()
                corrs = filenames[x][:-5] + '_spar_corrs.tsv'
                cov = filenames[x][:-5] + '_spar_cov.tsv'
                pvals = filenames[x][:-5] + '_spar_pvals.tsv'
                bootstraps = filenames[x][:-(5 + len(x))] + 'bootstraps'
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
                corrtab = pandas.read_csv(corrs, sep='\t', index_col=0)
                corrtab.columns = corrtab.index
                pvaltab = pandas.read_csv(pvals, sep='\t', index_col=0)
                pvaltab = pvaltab < pval_threshold  # p value threshold for SparCC pseudo p-values
                corrtab = corrtab.where(pvaltab)
                corrtab = corrtab.fillna(0)
                corrtab[corrtab > 0] = 1
                corrtab[corrtab < 0] = -1
                net = networkx.from_pandas_adjacency(corrtab)
                net = _add_tax(net, filenames[x])
                self.networks[("sparcc_" + x)] = net
                call("rm " + corrs + " " + pvals +
                     " " + os.path.dirname(massoc.__file__)[:-6] +
                     "\cov_mat_SparCC.out", shell=True)
            except Exception:
                logger.error("Unable to finish SparCC on " + str(x), exc_info=True)

    def run_conet(self):
        """
        Runs a Bash script containing the CoNet Bash commands.
        Unfortunately, errors produced by CoNet cannot be caught,
        because the exit status of the script is 0 regardless
        of CoNet producing a network or not.
        """
        path = resource_path('CoNet.sh')
        path = path.replace('\\', '/')
        libpath = self.inputs['conet'][0] + '\\lib\\CoNet.jar'
        libpath = libpath.replace('\\', '/')
        filenames = self.get_filenames()
        for x in filenames:
            try:
                graphname = filenames[x][:-5] + '_conet.tsv'
                tempname = filenames[x][:-5] + '_counts_conet.txt'
                file = biom.load_table(filenames[x])
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
                n_otus = file.shape[0]
                guessingparam = str(n_otus * n_otus -1)
                if int(guessingparam) > 1000:
                    guessingparam = str(1000)
                cmd = path + ' ' + tempname + ' ' + ' ' + graphname + ' ' + libpath + \
                      ' ' + filenames[x][:-5] + ' ' + guessingparam
                call(cmd, shell=True)
                call("rm " + tempname + " " + " " + filenames[x][:-5] + "_threshold" + " " +
                    filenames[x][:-5] + "_permnet", shell=True)
                with open(graphname, 'r') as fin:
                    data = fin.read().splitlines(True)
                    fin.close()
                with open(graphname, 'w') as fout:
                    fout.writelines(data[2:])
                    fout.close()
                signs = [x[0] for x in csv.reader(open(graphname, 'r'), delimiter='\t')]
                signs = [word.replace('mutualExclusion', '-1') for word in signs]
                signs = [word.replace('copresence', '1') for word in signs]
                signs = [word.replace('unknown', 'None') for word in signs]
                signs = [ast.literal_eval(x) for x in signs]
                clean_signs = list()  # None values need to be removed to make sure median is not 0.5.
                for sublist in signs:
                    cleaned = [elem for elem in sublist if elem is not None]
                    clean_signs.append(cleaned)
                signs = [statistics.median(x) for x in clean_signs]
                # methods = [x[2] for x in csv.reader(open(graphname, 'r'), delimiter='\t')]
                names = [x[15] for x in csv.reader(open(graphname, 'r'), delimiter='\t')]
                names = [x.split('->') for x in names]
                new_names = list()
                for item in names:
                    new_item = [x.replace(x, orig_ids[x]) for x in item]
                    new_names.append(new_item)
                file = biom.load_table(filenames[x])
                i = 0
                allspecies = file._observation_ids
                adj = pandas.DataFrame(index=allspecies, columns=allspecies)
                adj = adj.fillna(0)
                for name in new_names:
                    id1 = adj.columns.get_loc(name[0])
                    id2 = adj.columns.get_loc(name[1])
                    sign = signs[i]
                    i = i+1
                    adj.iloc[id1, id2] = sign
                    adj.iloc[id2, id1] = sign
                net = networkx.from_pandas_adjacency(adj)
                net = _add_tax(net, filenames[x])
                self.networks[("conet_" + x)] = net
                call ("rm " + graphname, shell=True)
            except Exception:
                logger.error("Unable to finish CoNet on " + str(x), exc_info=True)

    def write_networks(self):
        """
        Writes all networks in a Nets file to graphml files.
        """
        try:
            for network in self.networks:
                path = self.inputs['fp'][0] + '/' + network + '.txt'
                networkx.write_weighted_edgelist(G=self.networks[network], path=path)
        except Exception:
            logger.error("Unable to write networks to XML files", exc_info=True)

    def add_networks(self):
        """
        In case users want to manually import a network,
        this function adds the network file to the Nets object and checks
        whether the identifiers specified in the file match those in included BIOM files.
        Currently, only edge lists are supported.
        """
        filename = self.inputs['network']['filename']
        filetype = self.inputs['network']['filetype']
        if filetype == 'weighted':
            network = networkx.read_weighted_edgelist(filename)
        if filetype == 'unweighted':
            network = networkx.read_edgelist(filename)
        try:
            network_nodes = list(network.nodes)
        except TypeError:
            logger.error("Unable to read edge list.")
        taxon_ids = list()
        for biomfile in self.otu:
            taxon_ids.extend(list(self.otu[biomfile].ids(axis='observation')))
        missing_node = any(x not in taxon_ids for x in network_nodes)
        if missing_node:
            logger.error("Imported network node not found in taxon identifiers.")
        else:
            self.network[filename] = network


def _add_tax(network, file):
    """
    Adds taxon names from filename.
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
                networkx.set_node_attributes(network, values=taxdict, name=column)
    except Exception:
        logger.error("Unable to collect taxonomy for agglomerated files", exc_info=True)
    return network
