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
from massoc.scripts.netwrap import Nets
from multiprocess import Pool
from massoc.scripts.batch import Batch
import massoc

import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

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
    if inputs['cluster'] is not None:
        sys.stdout.write('Clustering BIOM files...')
        sys.stdout.flush()
        bioms.cluster_biom()
    if inputs['split'] is not None and inputs['split'] is not 'TRUE':
        bioms.split_biom()
    if inputs['min'] is not None:
        sys.stdout.write('Removing samples below minimum count...')
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
    bioms.write_bioms()
    return bioms


def run_networks(nets):
    """
    If network names are specified in the input files,
    this function calls the appropriate networks.
    This function is run sequentially when parallel is set to False.
    Otherwise, it pickles the network object
    and passes this to a function
    that can process the pickled object in parallel.
    Currently only works with default setting files.
    """
    if nets.inputs['levels'] is not None:
        if len(nets.inputs['levels']) > 1 or nets.inputs['levels'][0] is not 'otu':
            sys.stdout.write('Collapsing taxonomy...')
            sys.stdout.flush()
            nets.collapse_tax()
    else:
        nets.inputs['levels'] = list()
        nets.inputs['levels'].append('otu')
    if nets.inputs['tools'] is not None:
        if 'spiec-easi' in nets.inputs['tools']:
            sys.stdout.write('Running SPIEC-EASI...')
            sys.stdout.flush()
            nets.run_spiec()
        if 'sparcc' in nets.inputs['tools']:
            sys.stdout.write('Running SparCC...')
            sys.stdout.flush()
            nets.run_spar()
        if 'conet' in nets.inputs['tools']:
            sys.stdout.write('Running CoNet...')
            sys.stdout.flush()
            nets.run_conet()
    if nets.inputs['spiec'] is not None:
        sys.stdout.write('Running SPIEC-EASI with custom settings...')
        # sys.stdout.flush()
        nets.run_spiec(settings=nets.inputs['spiec'])
    if nets.inputs['conet'] is not None:
        sys.stdout.write('Running CoNet with custom settings...')
        sys.stdout.flush()
    if nets.inputs['spar_boot'] is not None and nets.inputs['spar_pval'] is None:
        sys.stdout.write('Running SparCC with custom settings...')
        sys.stdout.flush()
        nets.run_spar(boots=nets.inputs['spar_boot'])
    if nets.inputs['spar_pval'] is not None and nets.inputs['spar_boot'] is None:
        sys.stdout.write('Running SparCC with custom settings...')
        sys.stdout.flush()
        nets.run_spar(pval_threshold=nets.inputs['spar_pval'])
    if nets.inputs['spar_boot'] and nets.inputs['spar_pval'] is not None:
        sys.stdout.write('Running SparCC with custom settings...')
        sys.stdout.flush()
        nets.run_spar(boots=nets.inputs['spar_boot'], pval_threshold=nets.inputs['spar_pval'])
    sys.stdout.write('Finished running network inference!')
    sys.stdout.flush()
    return(nets)


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
        nets.run_spar()
    if 'conet' in job:
        sys.stdout.write('Running CoNet...')
        sys.stdout.flush()
        nets.run_conet()
    if 'spiec_setting' in job:
        sys.stdout.write('Running SPIEC-EASI with custom settings...')
        sys.stdout.flush()
        nets.run_spiec(list(job.values())[0][1])
    if 'conet_setting' in job:
        sys.stdout.write('Running CoNet with custom settings...')
        sys.stdout.flush()
        nets.run_conet(list(job.values())[0][1])
    if 'spar_setting' in job:
        sys.stdout.write('Running SparCC with custom settings...')
        sys.stdout.flush()
        if len(job['spar_setting'][1]) == 2:
            nets.run_spar(boots=job['spar_setting'][1]['spar_boot'], pval_threshold=job['spar_setting'][1]['spar_pval'])
        else:
            if 'spar_boot' in job['spar_setting'][1]:
                nets.run_spar(boots=job['spar_setting'][1]['spar_boot'])
            if 'spar_pval' in job['spar_setting'][1]:
                nets.run_spar(pval_threshold=job['spar_setting'][1]['spar_pval'])
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
        for i in nets.inputs['tools']:
            sublist[i] = level
        if nets.inputs['spiec'] is not None:
            sublist['spiec_setting'] = [level, nets.inputs['spiec']]
        if nets.inputs['conet'] is not None:
            sublist['conet_setting'] = [level, nets.inputs['conet']]
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
        nets.write_bioms()
    except Exception:
        logger.error('Could not write ' + str(nets.inputs['name'][0]) + ' to disk', exc_info=True)
    pool = Pool(cores)
    jobs = get_joblist(nets)
    func = partial(run_jobs, nets)
    try:
        results = pool.map(func, iter(jobs))
        logfile = open(resource_path("massoc.log"), 'r')
        logtext = logfile.read()
        logfile.close()
        dump = open(nets.inputs['fp'], 'w')
        dump.write(logtext)
        dump.close()
    except Exception:
        logger.error('Failed to generate workers', exc_info=True)
    nets.networks = results[0]
    for i in range(1, len(jobs)):
        nets.networks = {**nets.networks, **results[i]}
    return(nets)


def run_massoc(settings, mode="parallel"):
    """
    Pipes functions from the different massoc modules to run complete network inference.
    """
    if type(settings) is tuple:
        settings = settings[0]
    allbioms = combine_data(settings) # combines data, runs the preprocessing and writes to disk
    networks = Nets(allbioms)
    if mode != "parallel":
        networks = run_networks(networks)
    if mode == "parallel":
        networks = run_parallel(networks)
    networks.write_networks()
    networks.summary()

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller.
     Source: https://stackoverflow.com/questions/7674790/bundling-data-files-with-pyinstaller-onefile"""
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

if __name__ == '__main__':
    options = get_input(sys.argv[1:])
    run_massoc(vars(options))
