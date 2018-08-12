"""
These scripts allow massoc to run operations on Batch and Nets objects.
massoc accepts biom files or tab delimited files and constructs clusters.
Run_networks runs network inference sequentially, while
run_parallel runs them in parallel using the Python Pool function.
Run_networks is deprecated once run_parallel works appropriately:
running sequentially is also possible by only specifying 1 pool.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import argparse
import sys
from functools import partial
import os
from biom import load_table
from biom.parse import MetadataMap
from multiprocess import Pool
import multiprocessing
from massoc.scripts.batch import Batch
import massoc
from subprocess import call
from copy import deepcopy

import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

general_settings = {'biom_file': None, 'otu_table': None, 'tax_table': None, 'sample_data': None,
                    'otu_meta': None, 'cluster': None, 'split': None, 'prev': None, 'fp': None,
                    'levels': None, 'tools': None, 'spiec': None, 'conet': None, 'spar': None, 'spar_pval': None,
                    'spar_boot': None, 'nclust': None, 'name': None, 'cores': None, 'rar': None, 'min': None,
                    'network': None, 'assoc': None, 'agglom': None, 'logic': None,
                    'agglom_weight': None, 'export': None, 'neo4j': None, 'gml_name': None,
                    'procbioms': None, 'address': None, 'username': None, 'password': None}


def get_input(argv):
    """This parser gets inputs for BIOM and tab-delimited file processing.
    It combines these files in a Batch object.
    Moreover, it accepts arguments to perform machine learning algorithms
    on the BIOM files and can split datasets according to sample variables.
    Moreover, it accepts settings for network inference.
    It requires input biom files and constructs a Nets object.
    Multiple tools can be called at once by giving the tool names.
    However, users can also specify custom settings by providing
    alternative scripts for calling the tools."""
    parser = argparse.ArgumentParser(
        description='Import biom files'
                    'and cluster samples')

    parser.add_argument('-biom', '--biom_file',
                        dest='biom_file',
                        ninputs='+',
                        help='Input BIOM file. '
                             'Note: do not supply separate '
                             'OTU tables or taxonomy tables '
                             'if you give a BIOM file.',
                        default=None)
    parser.add_argument('-otu', '--otu_table',
                        dest='otu_table',
                        ninputs='+',
                        help='Input tab-delimited '
                             'OTU table. '
                             'If you supply multiple files, '
                             'make sure to name your files '
                             'OTU_xxx, tax_xxx, sample_xxx '
                             'and meta_xxx with xxx '
                             'being identical.',
                        default=False)
    parser.add_argument('-tax', '--tax_table',
                        dest='tax_table',
                        ninputs='+',
                        help='Input taxonomy table.',
                        default=None)
    parser.add_argument('-s', '--sample_data',
                        dest='sample_data',
                        ninputs='+',
                        help='Input sample data, '
                             'such as measurements '
                             'or geographical location.',
                        default=None)
    parser.add_argument('-od', '--otu_meta',
                        dest='otu_meta',
                        ninputs='+',
                        help='Input OTU metadata, '
                             'such as KO terms in genome '
                             'or phenotype.',
                        default=None)
    parser.add_argument('-cl', '--cluster_samples',
                        dest='cluster',
                        choices=['K-means', 'DBSCAN',
                                 'Gaussian', 'Spectral',
                                 'Affinity', None],
                        required=False,
                        help='Split data with clustering',
                        default=None)
    parser.add_argument('-nclust', '--number_of_clusters',
                        dest='nclust',
                        required=False,
                        help='Number of clusters to evaluate',
                        default=['4'])
    parser.add_argument('-split', '--split_samples',
                        dest='split',
                        required=False,
                        help='For splitting data with sample variable, '
                             'specify variable as the column name.'
                             'For splitting data by cluster: '
                             'specify TRUE.',
                        default=None)
    parser.add_argument('-prev', '--prevalence_filter',
                        dest='prev',
                        required=False,
                        help='Filter for OTUs that are present in '
                             'multiple samples. Please include prevalence'
                             'as a fraction between 0 and 1.',
                        default=None)
    parser.add_argument('-rar', '--rarefaction',
                        dest='rar',
                        required=False,
                        help='Set to TRUE to rarefy to even depth, '
                             'or specify a number to rarefy to.',
                        default=None)
    parser.add_argument('-min', '--mininum_abundance',
                        dest='min',
                        required=False,
                        help='Minimum mean abundance required'
                             'before filtering.',
                        default=None)
    parser.add_argument('-fp', '--output_filepath',
                        dest='fp',
                        help='Filepath for saving output files',
                        default=None)
    parser.add_argument('-name', '--file_name',
                        dest='name',
                        help='Prefix for saving output files. '
                             'Specify a unique name for every '
                             'BIOM file or OTU table.',
                        default='Out')
    parser.add_argument('-levels', '--tax_levels',
                        dest='levels',
                        help='Taxonomic levels used for network inference.',
                        default=['otu', 'genus', 'class'],
                        choices=['otu', 'genus', 'class',
                                 'family', 'order', 'phylum'])
    parser.add_argument('-tools', '--tool_names',
                        dest='tools',
                        required=False,
                        choices=['spiec-easi', 'sparcc', 'conet'],
                        help='Runs all listed tools with default settings.',
                        default=None)
    parser.add_argument('-spiec_settings', '--SPIEC-EASI_settings',
                        dest='spiec_settings',
                        required=False,
                        help='Location of SPIEC-EASI settings file. ',
                        default=None)
    parser.add_argument('-conet_settings', '--CoNet_settings',
                        dest='conet_settings',
                        required=False,
                        help='Location of CoNet settings file. ',
                        default=None)
    parser.add_argument('-spar', '--SparCC_executable',
                        dest='spar',
                        required=False,
                        help='Location of SparCC folder (not the .py file)',
                        default=None)
    parser.add_argument('-conet', '--CoNet_executable',
                        dest='conet_settings',
                        required=False,
                        help='Location of CoNet3 folder',
                        default=None)
    parser.add_argument('-spar_pval', '--SparCC_pval',
                        dest='spar_pval',
                        required=False,
                        help='Threshold for SparCC pseudo-pvalues. ',
                        default=None)
    parser.add_argument('-spar_boot', '--SparCC_boot',
                        dest='spar_boot',
                        required=False,
                        help='Number of bootstraps for SparCC. ',
                        default=None)
    parser.add_argument('-cores', '--number_of_processes',
                        dest='cores',
                        required=False,
                        help='Number of processes to distribute'
                             'jobs across.',
                        default=None)
    inputs = parser.parse_args(argv)
    if inputs.biom_file is not None:
        print('BIOM file(s) to process: ', inputs.biom_file)
    if inputs.otu_table is not None:
        print('Tab-delimited OTU table(s) to process: ', inputs.otu_table)
    if inputs.otu_table and inputs.tax_table is not None:
        if len(inputs.otu_table) is not len(inputs.tax_table):
            raise ValueError("Add a taxonomy table for every OTU table!")
    if inputs.otu_table and inputs.sample_data is not None:
        if len(inputs.otu_table) is not len(inputs.sample_data):
            raise ValueError("Add a sample data table for every OTU table!")
    if inputs.otu_table and inputs.otu_meta is not None:
        if len(inputs.otu_table) is not len(inputs.otu_meta):
            raise ValueError("Add a metadata table for every OTU table!")
    if inputs.tools is not None:
        print('Tools to run with default settings: ', inputs.tools)
    if inputs.spiec:
        print('Running SPIEC-EASI with these settings: ', inputs.spiec)
    if inputs.conet:
        print('Running CoNet with these settings: ', inputs.conet)
    if inputs.spar_boot or inputs.spar_pval:
        print('Running SparCC with these settings: ', inputs.spar_boot, inputs.spar_pval)
    if inputs.gcoda:
        print('Running gCoda with these settings: ', inputs.gcoda)
    print(inputs)
    return inputs

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller.
     Source: https://stackoverflow.com/questions/7674790/bundling-data-files-with-pyinstaller-onefile"""
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)


def combine_data(inputs):
    """Takes all input and returns a dictionary of biom files.
    If tab-delimited files are supplied, these are combined
    into a biom file. File names are used as keys.
    This is mostly a utility wrapper, as all biom-related functions
    are from biom-format.org.

    At the moment, rarefaction is performed after sample splitting.
    This means that samples with uneven sequence counts will not
    be rarefied to equal depths."""
    filestore = {}
    if inputs['biom_file'] is None:
        if inputs['otu_table'] is None:
            raise ValueError("Please supply either a biom file "
                             "or a tab-delimited OTU table!")
    i = 0  # i is used to assign correct 'name' vars
    if inputs['biom_file'] is not None:
        for x in inputs['biom_file']:
            biomtab = load_table(x)
            filestore[inputs['name'][i]] = biomtab
            i += 1
    if inputs['otu_table'] is not None:
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
            except TypeError:
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
                for i in list(obs_data):
                    tax = list()
                    for j in list(obs_data[i]):
                        tax.append(obs_data[i][j])
                        obs_data[i].pop(j, None)
                    obs_data[i]['taxonomy'] = tax
                biomtab.add_metadata(obs_data, axis='observation')
            filestore[inputs['name'][i]] = biomtab
            i += 1
            j += 1
    bioms = Batch(filestore, inputs)
    # it is possible that there are forbidden characters in the OTU identifiers
    # we can forbid people from using those, or replace those with an underscore
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
    if inputs['cluster'] is not None:
        sys.stdout.write('Clustering BIOM files...')
        sys.stdout.flush()
        bioms.cluster_biom()
    if inputs['split'] is not None and inputs['split'] is not 'TRUE':
        bioms.split_biom()
    if inputs['min'] is not None:
        sys.stdout.write('Removing taxa below minimum count...')
        sys.stdout.flush()
        bioms.prev_filter(mode='min')
    if inputs['rar'] is not None:
        sys.stdout.write('Rarefying samples...')
        sys.stdout.flush()
        bioms.rarefy()
    if inputs['prev'] is not None:
        sys.stdout.write('Setting prevalence filter...')
        sys.stdout.flush()
        bioms.prev_filter(mode='prev')
    bioms = massoc.scripts.netwrap.Nets(bioms)
    return bioms

def run_jobs(nets, job):
    """
    Accepts a job from a joblist to run network inference in parallel.
    """
    nets.inputs['levels'] = [list(job.values())[0]]
    if 'spiec-easi' in job:
        sys.stdout.write('Running SPIEC-EASI...')
        sys.stdout.flush()
        nets.run_spiec()
    if 'sparcc' in job:
        sys.stdout.write('Running SparCC...')
        sys.stdout.flush()
        if 'spar_setting' in job['sparcc']:
            if len(job['spar_setting'][1]) == 2:
                nets.run_spar(boots=job['spar_setting'][1]['spar_boot'], pval_threshold=job['spar_setting'][1]['spar_pval'])
            else:
                if 'spar_boot' in job['spar_setting'][1]:
                    nets.run_spar(boots=job['spar_setting'][1]['spar_boot'])
                if 'spar_pval' in job['spar_setting'][1]:
                    nets.run_spar(pval_threshold=job['spar_setting'][1]['spar_pval'])
        else:
            nets.run_spar()
    if 'conet' in job:
        sys.stdout.write('Running CoNet...')
        sys.stdout.flush()
        nets.run_conet()
    if 'spiec_setting' in job:
        sys.stdout.write('Running SPIEC-EASI with custom settings...')
        sys.stdout.flush()
        nets.run_spiec(list(job.values())[0][1])
    return nets.networks


def get_joblist(nets):
    """
    Creates a list of jobs that can be distributed over multiple processes.
    Note: should be appended to handle multiple taxonomic levels + files!
    Each job is a dictionary: the key represents a tool, while the
    value represents taxonomic level and file location.
    """
    joblist = list()
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
            joblist.append({value: sublist[value]})
    return joblist


def run_parallel(nets):
    """Creates partial function to run as pool."""
    if nets.inputs['levels'] is not None:
        if len(nets.inputs['levels']) > 1 or nets.inputs['levels'][0] is not 'otu':
            sys.stdout.write('Collapsing taxonomy...')
            sys.stdout.flush()
            nets.collapse_tax()
    else:
        nets.inputs['levels'] = list()
        nets.inputs['levels'].append('otu')
    cores = 4
    if nets.inputs['cores'] is not None:
        cores = int(nets.inputs['cores'][0])
    try:
        sys.stdout.write('Writing files to disk...')
        sys.stdout.flush()
        nets.write_bioms()
    except Exception:
        logger.error('Could not write ' + str(nets.inputs['name'][0]) + ' to disk', exc_info=True)
    jobs = get_joblist(nets)
    sys.stdout.write('Collecting jobs...')
    sys.stdout.flush()
    for item in jobs:
        for key in item.keys():
            if key == 'conet':
                # copies CoNet script path
                path = resource_path('CoNet.sh')
                path = path.replace('\\', '/')
                nets.log['conet'] = {'path': path}
                nets.log['conet']['level'] = item[key]
            if key == 'spar':
                nets.log['sparcc'] = {'bootstraps': 100, 'pvalue': 0.001}
                nets.log['sparcc']['level'] = item[key]
            if nets.inputs['spar_boot'] is not None:
                nets.log['sparcc']['bootstraps'] = nets.inputs['spar_boot']
            if nets.inputs['spar_pval'] is not None:
                nets.log['sparcc']['pvalue'] = nets.inputs['spar_pval']
            if key == 'spiec-easi':
                path = resource_path('spieceasi.r')
                path = path.replace('\\', '/')
                file = open(path, 'r')
                txt = file.read().splitlines()
                nets.log['spiec-easi'] = {'method': txt[31], 'stars': txt[33][-36:]}
                nets.log['spiec-easi']['level'] = item[key]
    # func = partial(run_jobs, nets)
    # pool = Pool(cores)
    # multiprocess supports passing objects
    # multiprocessing does not
    # however, multiprocess cannot be frozen
    # need to rewrite netwrap as pickle-able objects!
    try:
        sys.stdout.write('Distributing jobs...')
        sys.stdout.flush()
        network_list = list()
        for job in jobs:
            result = run_jobs(nets, job)
            network_list.append(result)
        # results = pool.map(func, iter(jobs))
    except Exception:
        logger.error('Failed to generate workers', exc_info=True)
    for item in network_list:
        nets.networks[list(item.keys())[0]] = item[list(item.keys())[0]]
    logfile = open(resource_path("massoc.log"), 'r')
    logtext = logfile.read()
    logfile.close()
    # for i in range(1, len(jobs)):
    #    nets.networks = {**nets.networks, **results[i]}
    # clean up old written BIOM files
    sys.stdout.write('Cleaning up old files...')
    sys.stdout.flush()
    for x in nets.inputs['name']:
        filename = nets.inputs['fp'][0] + '/' + x + '_otu.hdf5'
        call("rm " + filename, shell=True)
        if len(nets.genus) != 0:
            filename = nets.inputs['fp'][0] + '/' + x + '_species.hdf5'
            call("rm " + filename, shell=True)
            filename = nets.inputs['fp'][0] + '/' + x + '_genus.hdf5'
            call("rm " + filename, shell=True)
            filename = nets.inputs['fp'][0] + '/' + x + '_family.hdf5'
            call("rm " + filename, shell=True)
            filename = nets.inputs['fp'][0] + '/' + x + '_order.hdf5'
            call("rm " + filename, shell=True)
            filename = nets.inputs['fp'][0] + '/' + x + '_class.hdf5'
            call("rm " + filename, shell=True)
            filename = nets.inputs['fp'][0] + '/' + x + '_phylum.hdf5'
            call("rm " + filename, shell=True)
    sys.stdout.write('Completed tasks!')
    sys.stdout.flush()
    return(nets)


def run_massoc(settings, mode='write'):
    """
    Pipes functions from the different massoc modules to run complete network inference.
    """
    if type(settings) is tuple:
        settings = settings[0]
    allbioms = combine_data(settings) # combines data, runs the preprocessing and writes to disk
    networks = run_parallel(allbioms)
    networks.summary()
    if mode == 'write':
        networks.write_networks()
    else:
        return networks


if __name__ == '__main__':
    multiprocessing.freeze_support()
    options = get_input(sys.argv[1:])
    run_massoc(vars(options))
