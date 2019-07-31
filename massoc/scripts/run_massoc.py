"""
This file contains the CLI for massoc.

The command line interface is intended to be called sequentially;
files are written to disk as intermediates,
while a settings file is used to transfer logs and other information
between the modules.
This modular design allows users to leave out parts of massoc that are not required,
and reduces the number of parameters that need to be defined in the function calls.

"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import sys
import multiprocessing as mp
from massoc.scripts.main import get_input, run_network, \
    run_neo4j, run_netstats, run_metastats
import os
import argparse
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


def massoc(massoc_args):
    """
    Main function for running massoc. 
    :param massoc_args: Parameters generated through massoc_parser.
    :return: 
    """
    if 'input' in massoc_args:
        logger.info('Running input module. ')
        get_input(massoc_args)
    if 'network' in massoc_args:
        logger.info('Running network inference module. ')
        run_network(massoc_args)
    if 'database' in massoc_args:
        logger.info('Working on Neo4j database. ')
        run_neo4j(massoc_args)
    if 'netstats' in massoc_args:
        logger.info('Performing network analysis on Neo4j database. ')
        run_netstats(massoc_args)
    if 'metastats' in massoc_args:
        logger.info('Performing metadata analysis on Neo4j database. ')
        run_metastats(massoc_args)
    logger.info('Completed tasks! ')


massoc_parser = argparse.ArgumentParser(description='massoc pipeline')
massoc_parser.add_argument('-s', '--set',
                           dest='settings',
                           help='Settings txt file containing '
                                'filepaths to processed data, '
                                'as well as other settings. ',
                           default=os.getcwd() + '\\settings.json')
subparsers = massoc_parser.add_subparsers(title="massoc modules",
                                          description="Each module carries out a part of massoc. "
                                                      "Modules can be used independently "
                                                      "as long as the correct settings file "
                                                      "is provided. ")
inputparser = subparsers.add_parser('input', description='Import files and preprocesses them.',
                                    help='The input module accepts BIOM files or tab-delimited files'
                                         'and preprocesses these. This can include sample separation,'
                                         'clustering, rarefaction, and prevalence filtering. '
                                         'Additionally, in case users want to run network inference on '
                                         'different taxonomic levels, this module agglomerates the data.')
inputparser.add_argument('-fp', '--output_filepath',
                         dest='fp',
                         help='Filepath for saving output files',
                         default=os.getcwd())
inputparser.add_argument('-biom', '--biom_file',
                         dest='biom_file',
                         nargs='+',
                         help='Input BIOM file. '
                              'Note: do not supply separate '
                              'OTU tables or taxonomy tables '
                              'if you give a BIOM file.',
                         default=None)
inputparser.add_argument('-otu', '--otu_table',
                         dest='otu_table',
                         nargs='+',
                         help='Input tab-delimited '
                              'OTU table. '
                              'If you supply multiple files, '
                              'make sure to name your files '
                              'OTU_xxx, tax_xxx, sample_xxx '
                              'and meta_xxx with xxx '
                              'being identical.',
                         default=None)
inputparser.add_argument('-tax', '--tax_table',
                         dest='tax_table',
                         nargs='+',
                         help='Input taxonomy table.',
                         default=None)
inputparser.add_argument('-sample', '--sample_data',
                         dest='sample_data',
                         nargs='+',
                         help='Input sample data, '
                              'such as measurements '
                              'or geographical location.',
                         default=None)
inputparser.add_argument('-od', '--otu_meta',
                         dest='otu_meta',
                         nargs='+',
                         help='Input OTU metadata, '
                              'such as KO terms in genome '
                              'or phenotype.',
                         default=None)
inputparser.add_argument('-cl', '--cluster_samples',
                         dest='cluster',
                         choices=['K-means', 'DBSCAN',
                                  'Gaussian', 'Spectral',
                                  'Affinity', None],
                         required=False,
                         type=str,
                         help='Split data with clustering',
                         default=None)
inputparser.add_argument('-nclust', '--number_of_clusters',
                         dest='nclust',
                         required=False,
                         help='Number of clusters to evaluate',
                         default=4,
                         type=int)
inputparser.add_argument('-split', '--split_samples',
                         dest='split',
                         type=str,
                         required=False,
                         help='For splitting data with sample variable, '
                              'specify variable as the column name.'
                              'For splitting data by cluster: '
                              'specify TRUE.',
                         default=None)
inputparser.add_argument('-prev', '--prevalence_filter',
                         dest='prev',
                         required=False,
                         help='Filter for OTUs that are present in '
                              'multiple samples. Please include prevalence'
                              'as a fraction between 0 and 1.',
                         type=float,
                         default=None)
inputparser.add_argument('-rar', '--rarefaction',
                         dest='rar',
                         required=False,
                         help='Set to True to rarefy to even depth, '
                              'or specify a number to rarefy to.',
                         default=None)
inputparser.add_argument('-min', '--mininum_abundance',
                         dest='min',
                         required=False,
                         help='Minimum mean abundance required'
                         'before filtering.',
                         type=int,
                         default=None)
inputparser.add_argument('-name', '--file_name',
                         dest='name',
                         help='Prefix for saving output files. '
                              'Specify a unique name for every '
                              'BIOM file or OTU table.',
                         nargs='+',
                         default=['Out'])
inputparser.add_argument('-levels', '--tax_levels',
                         dest='levels',
                         nargs='+',
                         help='Taxonomic levels used for network inference.',
                         default=['otu'],
                         choices=['otu', 'species', 'genus', 'family', 'order', 'class', 'phylum'])
inputparser.add_argument('-net', '--networks',
                         dest='network',
                         nargs='+',
                         help='Weighted edge lists of previously generated networks.',
                         default=None)
inputparser.set_defaults(input=True)

networkparser = subparsers.add_parser('network', description='Runs network inference.',
                                      help='Given a settings file with preprocessed biom files,'
                                           'this module carries out network construction. '
                                           'Currently, SPIEC-EASI, CoNet and SparCC are supported. '
                                           'If you have difficulties running the tools through massoc, '
                                           'consider importing completed networks through the neo4j module. ')
networkparser.add_argument('-tools', '--tool_names',
                           dest='tools',
                           required=False,
                           choices=['spiec-easi', 'sparcc', 'conet'],
                           nargs='+',
                           help='Runs all listed tools with default settings.',
                           default=None)
networkparser.add_argument('-spiec_settings', '--SPIEC-EASI_settings',
                           dest='spiec',
                           required=False,
                           help='Location of SPIEC-EASI Rscript file. ',
                           default=None)
networkparser.add_argument('-conet_bash', '--CoNet_bashscript',
                           dest='conet_bash',
                           required=False,
                           help='Location of CoNet bash file. ',
                           default=None)
networkparser.add_argument('-spar', '--SparCC_executable',
                           dest='spar',
                           required=False,
                           help='Location of SparCC folder (not the .py file)',
                           default=None)
networkparser.add_argument('-conet', '--CoNet_executable',
                           dest='conet',
                           required=False,
                           help='Location of CoNet3 folder',
                           default=None)
networkparser.add_argument('-spar_pval', '--SparCC_pval',
                           dest='spar_pval',
                           required=False,
                           help='Threshold for SparCC pseudo-pvalues. ',
                           type=float,
                           default=None)
networkparser.add_argument('-spar_boot', '--SparCC_boot',
                           dest='spar_boot',
                           required=False,
                           help='Number of bootstraps for SparCC. ',
                           type=int,
                           default=None)
networkparser.add_argument('-cores', '--number_of_processes',
                           dest='cores',
                           required=False,
                           help='Number of processes to distribute'
                                'jobs across.',
                           type=int,
                           default=4)
networkparser.add_argument('-fp', '--output_filepath',
                           dest='fp',
                           help='Filepath for saving output files and reading settings.',
                           default=os.getcwd())
networkparser.set_defaults(network=True)


neo4jparser = subparsers.add_parser('neo4j', description='Sets up a Neo4j graph database.',
                                    help='If the user provides a settings file with '
                                         'filenames that need to be imported '
                                         '(those include processed BIOM files and network files),'
                                         'this module uses the supplied username, password and '
                                         'filepath to start up a Neo4j database. The database model is '
                                         'described in more detail in the manual. '
                                         'Other files can also be provided here, '
                                         'as long as they are given as edge lists or other network files.')
neo4jparser.add_argument('-n', '--neo4j',
                         dest='neo4j',
                         help='Filepath to neo4j folder. ',
                         required=False,
                         type=str,
                         default=None)
neo4jparser.add_argument('-u', '--username',
                         dest='username',
                         required=False,
                         help='Username for neo4j database access. ',
                         type=str,
                         default='neo4j')
neo4jparser.add_argument('-p', '--password',
                         dest='password',
                         required=False,
                         type=str,
                         help='Password for neo4j database access. ')
neo4jparser.add_argument('-a', '--address',
                         dest='address',
                         required=False,
                         help='Address for neo4j database. ',
                         type=str,
                         default='bolt://localhost:7687')
neo4jparser.add_argument('-j', '--job',
                         dest='job',
                         required=False,
                         choices=['clear', 'quit', 'start', 'upload', 'write'],
                         default='upload',
                         help='Operation to carry out on Neo4j database. \n'
                              'clear: Remove all edges and nodes from the database. \n'
                              'quit: Safely shut down the database. Retains data. \n'
                              'start: Start local database. \n'
                              'upload: Upload BIOM file(s) in settings to the database. \n'
                              'write: Export Cytoscape-compatible file from the database.')
neo4jparser.add_argument('-o', '-output',
                         dest='output',
                         required=False,
                         default='graph',
                         type=str,
                         help='Filename for exporting Cytoscape-compatible graphs from the Neo4j database.')
neo4jparser.add_argument('-fp', '--output_filepath',
                         dest='fp',
                         help='Filepath for saving output files and reading settings.',
                         default=os.getcwd())
neo4jparser.set_defaults(database=True)


netstatsparser = subparsers.add_parser('netstats', description='Network analysis module.',
                                       help='If the user has previously set up a Neo4j graph database,'
                                            'this module carries out cluster analysis and '
                                            'identifies central nodes. '
                                            'In the future, the module will also generates null models'
                                            ' for statistical analyses. ')
netstatsparser.add_argument('-l', '--logic',
                            dest='logic',
                            required=False,
                            nargs='+',
                            help='Logic operations to carry out on Neo4j database with multiple networks. ',
                            choices=['union', 'intersection', 'difference'],
                            default=None)
netstatsparser.add_argument('-w', '--weight',
                            dest='weight',
                            action='store_true',
                            required=False,
                            help='If flagged, intersections include associations with\n '
                                 'matching partners but different weights, and differences exclude these. ',
                            default=None)
netstatsparser.add_argument('-net', '--networks',
                            dest='networks',
                            required=False,
                            action='store_true',
                            default=None,
                            help='Intersection-only option that returns the intersection for \n'
                                 'any n-combination of the specified networks.')
netstatsparser.add_argument('-n', '--num',
                            dest='num',
                            required=False,
                            default=None,
                            help='Carry out set operations on a subset of networks in the database.'
                                 '\n Make sure supplied names match the network node labels.')
netstatsparser.add_argument('-fp', '--output_filepath',
                            dest='fp',
                            help='Filepath for saving output files and reading settings.',
                            default=os.getcwd())
netstatsparser.set_defaults(netstats=True)


metastatsparser = subparsers.add_parser('metastats', description='Network meta-analysis module.',
                                        help='If the user has previously set up a Neo4j graph database,'
                                             'this module carries out network analyses that incorporate'
                                             ' provided metadata. Such metadata can include taxonomic data,'
                                             ' or it can include nodes with specific labels. ')
metastatsparser.add_argument('-add', '--additional_data',
                             dest='add',
                             required=False,
                             default=None,
                             nargs='+',
                             help='Filepath to edge list or table that will be uploaded to the Neo4j database. \n'
                             'If the file is a table, the first column should be of nodes inside the database. \n'
                             'If the column header does not match the network label, (e.g. #SampleID instead of Sample) \n'
                             'specify this with the type argument.')
metastatsparser.add_argument('-abn', '--abundance_data',
                             dest='abundance',
                             required=False,
                             default=None,
                             nargs='+',
                             help='If you are adding metadata that is an abundance table,'
                             ' specify the type of data (e.g. KO terms) here. ')
metastatsparser.add_argument('-type', '--annotation_type',
                             dest='type',
                             help='Node label used in place of 1st column header. \n'
                             'This label (e.g. Sample) should be specified when using QIITA files,\n'
                             'or other tables that have a # in front of the column names. ',
                             default=None,
                             type=str, )
metastatsparser.add_argument('-tax', '--tax_agglomeration',
                             dest='agglom',
                             help='Taxonomic level used for network agglomeration.',
                             default=None,
                             choices=['species', 'genus', 'family', 'order', 'class', 'phylum'])
metastatsparser.add_argument('-w', '--weight',
                             dest='weight',
                             required=False,
                             action='store_true',
                             help='Take weight into account for taxonomic agglomeration.')
metastatsparser.add_argument('-v', '--variable',
                             dest='variable',
                             required=False,
                             type=str,
                             nargs='+',
                             help='Sample metadata variables to associate taxa to. Can be more than one variable. \n'
                                  'Specify "all" if you want to test all variables.')
metastatsparser.add_argument('-fp', '--output_filepath',
                             dest='fp',
                             help='Filepath for saving output files and reading settings.',
                             default=os.getcwd())
metastatsparser.add_argument('-s', '--sequence',
                             dest='sequence',
                             required=False,
                             type=str,
                             help='Location of 16S sequences (e.g. GreenGenes folder containing FASTA files).')
metastatsparser.set_defaults(metastats=True)


def main():
    mp.freeze_support()
    options = massoc_parser.parse_args()
    massoc(vars(options))


if __name__ == '__main__':
    main()
