"""
BIOM batch

The biom-format provides a hdf5 format for storage
and annotation of microbiome data.
The Batch class expands on this by storing multiple
biom files and performing processes on these files in batch.
Effectively, the Batch class ensures that all biom files
analysed by a researcher undergo the same processing steps.
Moreover, the Batch class contains a log file that
documents these processing steps.

Note that the split + cluster code creates a deep copy of the OTU dict
and changes the batch.otu dictionary to this deep copy.
This is to prevent strange behaviour, e.g. iterating over the dict
and transforming it needs to be possible.

Kaufman, L., & Rousseeuw, P. J. (1990).
Finding groups in data: an introduction to cluster analysis, 68-125.

"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

from biom.cli.util import write_biom_table
from skbio.stats.composition import multiplicative_replacement, clr, ilr
from scipy.sparse import csr_matrix, csc_matrix
import copy
import numpy as np
from copy import deepcopy
from sklearn.cluster import KMeans, DBSCAN, SpectralClustering, AffinityPropagation
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
import json
import sys
from biom import load_table
import logging.handlers

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)

# handler to file
# only handler with 'w' mode, rest is 'a'
# once this handler is started, the file writing is cleared
# other handlers append to the file


class Batch(object):

    """Container for multiple BIOM files.
    The container is initialized from an absolute
    abundance count table, and then agglomerates this
    to higher taxonomic levels.
    These are all accessible from the Batch object,
    so no need to agglomerate biom files one by one.
    The Batch object can write files simultaneously,
    stores command line arguments (and appends these)
    and prints out the arguments if required.
    The Batch object also contains methods for
    normalization, log transforms and prevalence filtering.

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

    Raises
    ------
    TableException
        batch only operates on BIOM OTU tables
    """

    def __init__(self, counts=None, inputs=None):
        """
        Initializes the Batch object.

        :param counts: Dictionary of dictionary of BIOM files.
        The first-level dictionary has taxonomic levels as keys, the second names.
        :param inputs: Dictionary of inputs.
        :return:
        """
        self.otu = {}
        self.species = {}
        self.genus = {}
        self.family = {}
        self.order = {}
        self.class_ = {}
        self.phylum = {}
        self.inputs = inputs
        if inputs:
            _create_logger(self.inputs['fp'])
        if counts is not None:
            if 'otu' in counts:
                self.otu = counts['otu']
            if 'species' in counts:
                self.species = counts['species']
            if 'genus' in counts:
                self.genus = counts['genus']
            if 'family' in counts:
                self.family = counts['family']
            if 'order' in counts:
                self.order = counts['order']
            if 'class' in counts:
                self.class_ = counts['class']
            if 'phylum' in counts:
                self.phylum = counts['phylum']
        self.levels = {'otu': self.otu,
                       'species': self.species,
                       'genus': self.genus,
                       'family': self.family,
                       'order': self.order,
                       'class': self.class_,
                       'phylum': self.phylum}
        # necessary for collapse func
        self.n = {'species': 7,
                  'genus': 6,
                  'family': 5,
                  'order': 4,
                  'class': 3,
                  'phylum': 2}

    def collapse_tax(self):
        """
        The collapse_tax function allows users to generate BIOM files
        containing agglomerated data. This means network analysis can be
        performed simultaneously on genus, order and other taxonomic levels.

        :return:
        """
        try:
            for level in self.inputs['levels']:
                if level != 'otu':
                    for x in list(self.levels['otu']):
                        self.levels[level][x] = _data_bin(self.otu[x], self.n[level], level + '_' + x)
            self.write_bioms()
        except TypeError:
            logger.error("Could not collapse taxonomy", exc_info=True)

    def get_filenames(self):
        """
        If the standard processing functions in the Batch object are used,
        filenames follow a standard format. This function retrieves such filenames
        so BIOM files can be read without requiring them to be stated explicitly
        by the user.

        :return:
        """
        try:
            filenames = dict()
            for level in self.inputs['levels']:
                filenames[level] = dict()
                for x in self.inputs['name']:
                    if level is 'class':
                        filenames[level][x] = self.inputs['fp'] + '/' + x + '_' + level + '.hdf5'
                    else:
                        filenames[level][x] = self.inputs['fp'] + '/' + x + '_' + level + '.hdf5'
        except TypeError:
            logger.error('Unable to generate filenames. \n', exc_info=True)
        return filenames

    def write_bioms(self, fmt='hdf5'):
        """
        Utility function that writes BIOM files
        in a Batch object to HDF5 files.
        OTU files are always written to disk,
        the rest only if required.

        :param fmt: Format for writing; 'hdf5' or 'json'.
        :return:
        """
        for x in self.inputs['name']:
            for level in self.inputs['levels']:
                filename = self.inputs['fp'] + '/' + x + '_' + level + '.hdf5'
                try:
                    write_biom_table(self.levels[level][x], fmt, filename)
                except Exception:
                    logger.error("Cannot write " + str(x) + " to disk", exc_info=True)

    def normalize_transform(self, mode='clr'):
        """
        Some operations may require transformed data.
        This function performs normalization and
        a clr transform on all OTU tables in a Batch object.
        It returns a deep copy of the original Batch object,
        so the original file is not modified.

        :param mode: transformation mode; clr (centered log-ratio) or ilr (isometric log-ratio)
        :return: Transformed copy of Batch object.
        """
        batchcopy = copy.deepcopy(self)
        try:
            for x in list(self.otu):
                # normalizes the data by samples
                normbiom = batchcopy.otu[x].norm(axis='sample', inplace=False)
                mat = csr_matrix.toarray(normbiom.matrix_data)
                # replaces all zeros with a small value
                # multiplicative replacement preserves ratios between values
                mat = multiplicative_replacement(mat)
                if mode is 'clr':
                    mat = clr(mat)
                elif mode is 'ilr':
                    mat = ilr(mat)
                else:
                    raise ValueError("Only CLR and ILR transformations are currently supported.")
                normbiom._data = csc_matrix(mat)
                batchcopy.otu[x] = normbiom
        except Exception:
            logger.error("Failed to normalize data", exc_info=True)
        return batchcopy

    def prev_filter(self, mode='prev'):
        """
        Some operations may require transformed data.
        This function performs normalization and
        a clr transform on all OTU tables in a Batch object.
        It returns a deep copy of the original Batch object,
        so the original file is not modified.

        :param mode: prev or min, specifies whether taxa should be filtered
        based on prevalence or minimum abundance. The values are stored in the batch.inputs dictionary.
        :return:
        """
        for level in self.levels:
            for name in self.levels[level]:
                data = self.levels[level][name].matrix_data
                data = csr_matrix.todense(data)
                keep_otus = list()
                binotu = None
                try:
                    if mode == 'prev':  # calculates prevalence
                        fracs = np.count_nonzero(data, axis=1)
                        nsamples = data.shape[1]
                        fracs = fracs / nsamples
                        for y in range(0, len(fracs)):
                            if fracs[y] >= (float(self.inputs['prev'])/100):
                                keep_otus.append(self.levels[level][name]._observation_ids[y])
                            else:
                                binotu = self.levels[level][name]._observation_ids[y]
                        if binotu is not None and 'Bin' not in keep_otus:
                            keep_otus.append(binotu)
                except Exception:
                    logger.error("Could not set prevalence filter", exc_info=True)
                try:
                    if mode == 'min':
                        mincount = np.sum(data, axis=1)
                        for y in range(0, len(mincount)):
                            if mincount[y] >= (int(self.inputs['min'])):
                                keep_otus.append(self.levels[level][name]._observation_ids[y])
                            else:
                                binotu = self.levels[level][name]._observation_ids[y]
                        if binotu is not None:
                            keep_otus.append(binotu)
                except Exception:
                    logger.error("Could not set a minimum count filter", exc_info=True)
                keep = self.levels[level][name].filter(keep_otus, axis="observation", inplace=False)
                try:
                    if binotu is not None:
                        bin = self.levels[level][name].filter(keep_otus[:-1], axis="observation", inplace=False, invert=True)
                        binsums = np.sum(bin.matrix_data, axis=0) # sums all binned OTUs
                        # need to recreate keep._data as lil matrix, is more efficient
                        orig = keep._data.tolil(copy=True)
                        if 'Bin' not in keep_otus:
                            bin_id = keep._obs_index[binotu]
                            orig[bin_id] = binsums
                            keep._observation_ids[bin_id] = "Bin"
                            keep._obs_index["Bin"] = keep._obs_index.pop(binotu)
                        if 'Bin' in keep_otus:  # necessary to prevent duplicate Bin ID
                            old_bin_id = keep._obs_index["Bin"]
                            old_bin_sums = keep._data[old_bin_id]
                            new_bin_sums = binsums + old_bin_sums
                            orig[old_bin_id] = new_bin_sums
                        # update keep._data with orig
                        keep._data = orig.tocsr()
                except Exception:
                    logger.error("Could not preserve binned taxa", exc_info=True)
                self.levels[level][name] = keep

    def rarefy(self):
        """
        For each BIOM file, a rarefaction filter is applied.
        A mininum read depth can be specified;
        samples with reads lower than this read depth are removed,
        and then samples are rarefied to equal depth.

        :return:
        """
        all_bioms = {'otu': self.otu, 'genus': self.genus,
                     'family': self.family, 'order': self.order,
                     'class': self.class_, 'phylum': self.phylum}
        batchcopy = deepcopy(all_bioms)
        for level in all_bioms:
            for name in all_bioms[level]:
                try:
                    if self.inputs['rar'] == 'True':
                        lowest_count = int(min(all_bioms[level][name].sum(axis='sample')))
                    else:
                        lowest_count = int(self.inputs['rar'])
                    data = all_bioms[level][name].matrix_data
                    data = csr_matrix.todense(data)
                    keep_samples = list()
                    mincount = np.sum(data, axis=0)
                    for y in range(mincount.shape[1]):
                        if mincount.item(y) >= lowest_count:
                            keep_samples.append(all_bioms[level][name]._sample_ids[y])
                    keep = all_bioms[level][name].filter(keep_samples, axis="sample", inplace=False)
                    batchcopy[level][name] = keep.subsample(n=lowest_count, axis='sample')
                except Exception:
                    logger.error("Unable to rarefy file", exc_info=True)
                for name in list(all_bioms[level]):
                    all_bioms[level][name] = batchcopy[level][name]
        self.otu = all_bioms['otu']
        self.genus = all_bioms['genus']
        self.family = all_bioms['family']
        self.order = all_bioms['order']
        self.class_ = all_bioms['class']
        self.phylum = all_bioms['phylum']

    def split_biom(self):
        """
        Splits bioms into several subfiles according to
        sample metadata variable. Source: biom-format.org.
        The original file is preserved, so returned files
        include the split- and non-split files.

        :return:
        """
        inputs = self.inputs
        part_f = lambda id_, md: md[inputs['split']]
        for level in self.inputs['levels']:
            new_dict = deepcopy(self.levels[level])
            if type(self.levels[level]) is not dict:
                logger.warning('Split_biom requires a dictionary of biom files to be supplied. \n', exc_info=True)
                raise ValueError("Split_biom requires a dictionary of biom files to be supplied.")
            for x in list(self.levels[level]):
                try:
                    biomtab = new_dict[x]
                    if inputs['split'] not in biomtab._sample_metadata[1]:
                        if inputs['split'] is not 'TRUE':
                            raise Warning("Sample metadata of " + x + "does not contain this header!")
                    new_tables = biomtab.partition(part_f, axis='sample')
                    for new in new_tables:
                        key = x + '_' + new[0]
                        new_dict[key] = new[1]
                        self.inputs['name'].append(key)
                except Exception:
                    logger.error("Failed to split files", exc_info=True)
            self.levels[level] = new_dict

    def cluster_biom(self):
        """
        First normalizes bioms so clustering is not affected,
        performs transformation and then applies clustering.
        Note that the returned biom files are not normalized,
        this is just for the clustering process.
        Many network inference tools require absolute counts.
        Silhouette score is used to determine the optimal
        number of clusters.
        Clustering adds metadata info to the samples.
        Splitting according to cluster ID is done
        by wrapping the split_biom function.

        :return:
        """
        inputs = self.inputs
        if inputs['nclust'] is not None:
            nums = list(range(2, (int(inputs['nclust']) + 1)))
        else:
            nums = list(range(2,5))
        new_dict = {}
        if type(self.otu) is not dict:
            logger.warning('Cluster_biom requires a dictionary of biom files to be supplied. \n', exc_info=True)
            raise ValueError("Cluster_biom requires a dictionary of biom files to be supplied.")
        normbatch = self.normalize_transform(mode='clr')
        # CLR transform places data in Euclidean space
        for x in list(self.otu):
            try:
                # define topscore and bestcluster for no cluster
                norm_table = normbatch.otu[x]
                topscore = 0
                bestcluster = [1] * len(norm_table.ids())
                data = csr_matrix.todense(norm_table.matrix_data)
                data = np.matrix.transpose(data)
                data = PCA(n_components=2).fit_transform(data)
                randomclust = np.random.randint(2, size=len(data))
                sh_score = [silhouette_score(data, randomclust)]
                # K-means clustering, tests 2-4 clusters
                if inputs['cluster'] == 'K-means':
                    for i in nums:
                        clusters = KMeans(i).fit_predict(data)
                        silhouette_avg = silhouette_score(data, clusters)
                        sh_score.append(silhouette_avg)
                    topscore = int(np.argmax(sh_score) + 1)
                    bestcluster = KMeans(topscore).fit_predict(data)
                # DBSCAN clustering, automatically finds optimal cluster size
                if inputs['cluster'] == 'DBSCAN':
                    bestcluster = DBSCAN().fit_predict(data)
                    topscore = len(set(bestcluster)) - (1 if -1 in bestcluster else 0)
                # Gaussian Mixture Model (gmm) probability distribution
                if inputs['cluster'] == 'Gaussian':
                    for i in nums:
                        fit = GaussianMixture(i).fit(data)
                        clusters = fit.predict(data)
                        silhouette_avg = silhouette_score(data, clusters)
                        sh_score.append(silhouette_avg)
                    topscore = int(np.argmax(sh_score) + 1)
                    bestfit = GaussianMixture(topscore).fit(data)
                    bestcluster = bestfit.predict(data)
                # Spectral Clustering
                if inputs['cluster'] == 'Spectral':
                    for i in nums:
                        clusters = SpectralClustering(i).fit_predict(data)
                        silhouette_avg = silhouette_score(data, clusters)
                        sh_score.append(silhouette_avg)
                    topscore = int(np.argmax(sh_score) + 1)
                    bestcluster = SpectralClustering(topscore).fit_predict(data)
                # Affinity Propagation clustering
                if inputs['cluster'] == 'Affinity':
                    bestcluster = AffinityPropagation().fit_predict(data)
                    topscore = len(set(bestcluster)) - (1 if -1 in bestcluster else 0)
                if max(sh_score) < 0.25:
                    raise ValueError("Silhouette score too low: please try a different algorithm. "
                                     "Your data may not be suitable for clustering.")
                new_dict[x] = deepcopy(self.otu[x])
                for i in range(topscore):
                    mask, = np.where(bestcluster == i)
                    for j in mask:
                        new_dict[x]._sample_metadata[j]['cluster'] = inputs['cluster'] + '_' + str(i)
                self.otu = new_dict
                if inputs['split'] is not None:
                    if inputs['split'] == 'TRUE':
                        inputs['split'] = 'cluster'
                        self.split_biom()
            except Exception:
                logger.error("Error occurred when clustering samples", exc_info=True)


def _data_bin(biomfile, taxnum, key):
    """
    While the BIOM file collapse function collapses counts, it stores taxonomy in a manner
    that is not compatible with the OTU taxonomy; taxonomy is stored as observation ID and not as
    observation metadata.
    Here, new observation IDs are generated for agglomerated data that are a combination
    of the old IDs. Moreover, the taxonomy of the agglomerated data is concatenated to
    the agglomeration level and added to the observation metadata.

    :param biomfile: BIOM file according to the biom-format standards.
    :param taxnum: Number indicating taxonomic level (e.g. 6 for Genus).
    :param key: Value used for naming agglomerated taxa.
    :return: Taxonomically collapsed BIOM file.
    """
    collapse_f = lambda id_, md: '; '.join(md['taxonomy'][:taxnum])
    collapsed = biomfile.collapse(collapse_f, norm = False, min_group_size = 1,
                                  axis='observation', include_collapsed_metadata=True)
    agglom_ids = collapsed._observation_metadata
    obs_ids = collapsed._observation_ids
    for id in agglom_ids:
        keep_id = id['collapsed_ids'][0]  # gets OTU ID from original biom file to construct taxonomy
        keep_obs = biomfile._obs_index[keep_id]  # translates OTU ID to index
        keep_tax = biomfile._observation_metadata[keep_obs]['taxonomy'][:taxnum]  # removes tax up to agglom
        id['taxonomy'] = keep_tax
    new_obs_ids = deepcopy(obs_ids)
    new_obs_index = deepcopy(collapsed._obs_index)
    for number in range(len(obs_ids)):
        id = obs_ids[number]
        num_id = collapsed._obs_index[id]
        old_ids = agglom_ids[num_id]
        new_id = key + "-agglom-" + str(number)
        new_obs_ids[number] = new_id
        # need to update index as well
        new_obs_index[new_id] = new_obs_index.pop(id)
    collapsed._observation_ids = new_obs_ids
    collapsed._obs_index = new_obs_index
    return collapsed


def write_settings(settings, path=None):
    """
    Writes a dictionary of settings to a json file.

    :param settings: Dictionary of settings (i.e. 'inputs' in above functions).
    :param settings: Filepath to settings file.
    :return:
    """
    if not path:
        path = settings['fp'] + '/settings.json'
    # as requested in comment
    with open(path, 'w') as file:
        file.write(json.dumps(settings))
    file.close()
    logger.info("Wrote settings to: " + path)


def read_settings(path):
    """
    Writes a dictionary of settings to a json file.

    :param path: Filepath to settings file.
    :return: Dictionary of settings
    """
    # as requested in comment
    with open(path, 'r') as file:
        settings = json.load(file)
    file.close()
    return settings


def read_bioms(counts):
    """
    Given a dictionary of filenames per taxonomic level,
    this function generates a dictionary of biom files.

    :param counts: Dictionary of filenames leading to BIOM files
    :return: Dictionary of BIOM files
    """
    bioms = deepcopy(counts)
    for level in counts:
        for filename in counts[level]:
            bioms[level][filename] = load_table(counts[level][filename])
    return bioms


def _create_logger(filepath):
    """
    After a filepath has become available, loggers can be created
    when required to report on errors.

    :param filepath: Path where logger files are stored.
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
