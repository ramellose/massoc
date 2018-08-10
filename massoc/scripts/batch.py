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

from datetime import datetime
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
import os
import massoc
import logging
import logging.handlers as handlers
logger = logging.getLogger()
logger.setLevel(logging.WARNING)

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

    def __init__(self, otu=None, inputs=None):
        self.log = dict()
        self.log["Creation"] = datetime.now().strftime('%B %d %Y %H:%M:%S')
        self.otu = otu
        self.species = {}
        self.genus = {}
        self.family = {}
        self.order = {}
        self.class_ = {}
        self.phylum = {}
        self.inputs = inputs

    def collapse_tax(self):
        """
        The collapse_tax function allows users to generate BIOM files
        containing agglomerated data. This means network analysis can be
        performed simultaneously on genus, order and other taxonomic levels.
        """
        try:
            for x in list(self.otu):
                self.species[x + '_species'] = _data_bin(self.otu[x], 7)
                self.genus[x + '_genus'] = _data_bin(self.otu[x], 6)
                self.family[x + '_family'] = _data_bin(self.otu[x], 5)
                self.order[x + '_order'] = _data_bin(self.otu[x], 4)
                self.class_[x + '_class'] = _data_bin(self.otu[x], 3)
                self.phylum[x + '_phylum'] = _data_bin(self.otu[x], 2)
            self.log['collaped_tax'] = "True"
            self.write_bioms()
        except TypeError:
            logger.error("Could not collapse taxonomy", exc_info=True)

    def summary(self):
        """
        This function prints a summary of the settings,
        as well as performed operations.
        """
        for key in self.log:
            if type(self.log[key]) is str:
                print(key + ': ' + self.log[key])
            else:
                for item in self.log[key]:
                    print(item + ': ' + self.log[key][item])


    def write_bioms(self, fmt='hdf5'):
        """
        Utility function that writes all BIOM files
        in a Batch object to HDF5 files.
        """
        for x in self.inputs['name']:
            filename = self.inputs['fp'][0] + '/' + x + '_otu.hdf5'
            try:
                write_biom_table(self.otu[x], fmt, filename)
            except Exception:
                logger.error("Cannot write " + str(x) + " to disk", exc_info=True)
            if len(self.genus) != 0:
                filename = self.inputs['fp'][0] + '/' + x + '_species.hdf5'
                write_biom_table(self.species[x + '_species'], fmt, filename)
                filename = self.inputs['fp'][0] + '/' + x + '_genus.hdf5'
                write_biom_table(self.genus[x + '_genus'], fmt, filename)
                filename = self.inputs['fp'][0] + '/' + x + '_family.hdf5'
                write_biom_table(self.family[x + '_family'], fmt, filename)
                filename = self.inputs['fp'][0] + '/' + x + '_order.hdf5'
                write_biom_table(self.order[x + '_order'], fmt, filename)
                filename = self.inputs['fp'][0] + '/' + x + '_class.hdf5'
                write_biom_table(self.class_[x + '_class'], fmt, filename)
                filename = self.inputs['fp'][0] + '/' + x + '_phylum.hdf5'
                write_biom_table(self.phylum[x + '_phylum'], fmt, filename)

    def normalize_transform(self, mode='clr'):
        """
        Some operations may require transformed data.
        This function performs normalization and
        a clr transform on all OTU tables in a Batch object.
        It returns a deep copy of the original Batch object,
        so the original file is not modified.
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
        Filters all OTUs in a Batch.otu dictionary for prevalence.
        prev should be a value between 0 and 1, and is the
        minimum fraction of non-zero values required per OTU.
        After OTUs are filtered for prevalence,
        higher taxonomic levels are updated to reflect this.
        """
        for x in list(self.otu):
            data = self.otu[x].matrix_data
            data = csr_matrix.todense(data)
            keep_otus = list()
            binotu = None
            try:
                if mode is 'prev':  # calculates prevalence
                    fracs = np.count_nonzero(data, axis=1)
                    nsamples = data.shape[1]
                    fracs = fracs / nsamples
                    for y in range(0, len(fracs)):
                        if fracs[y] >= (float(self.inputs['prev'][0])/100):
                            keep_otus.append(self.otu[x]._observation_ids[y])
                        else:
                            binotu = self.otu[x]._observation_ids[y]
                    if binotu is not None and 'Bin' not in keep_otus:
                        keep_otus.append(binotu)
                    self.log['prevalence_filter'] = self.inputs['prev'][0] + "%"
            except Exception:
                logger.error("Could not set prevalence filter", exc_info=True)
            try:
                if mode is 'min':
                    mincount = np.sum(data, axis=1)
                    for y in range(0, len(mincount)):
                        if mincount[y] >= (int(self.inputs['min'][0])):
                            keep_otus.append(self.otu[x]._observation_ids[y])
                        else:
                            binotu = self.otu[x]._observation_ids[y]
                    if binotu is not None:
                        keep_otus.append(binotu)
                    self.log['abundance_filter'] = self.inputs['min'][0] + " counts"
            except Exception:
                logger.error("Could not set a minimum count filter", exc_info=True)
            keep = self.otu[x].filter(keep_otus, axis="observation", inplace=False)
            try:
                if binotu is not None:
                    bin = self.otu[x].filter(keep_otus[:-1], axis="observation", inplace=False, invert=True)
                    data = bin.matrix_data
                    data = csr_matrix.todense(data)
                    binsums = np.sum(data, axis=0) # sums all binned OTUs
                    if 'Bin' not in keep_otus:
                        bin_id = keep._obs_index[binotu]
                        keep._data[bin_id] = binsums
                        keep._observation_ids[bin_id] = "Bin"
                        keep._obs_index["Bin"] = keep._obs_index.pop(binotu)
                    if 'Bin' in keep_otus:  # necessary to prevent duplicate Bin ID
                        old_bin_id = keep._obs_index["Bin"]
                        old_bin_sums = keep._data[old_bin_id]
                        new_bin_sums = binsums + old_bin_sums
                        keep._data[old_bin_id] = new_bin_sums
            except Exception:
                logger.error("Could not preserve binned taxa", exc_info=True)
            self.otu[x] = keep
        if len(self.genus) > 0:
            self.collapse_tax()


    def rarefy(self):
        """
        For each BIOM file, a rarefaction filter is applied.
        A mininum read depth can be specified;
        samples with reads lower than this read depth are removed,
        and then samples are rarefied to equal depth.
        """
        batchcopy = deepcopy(self)
        for x in list(self.otu):
            try:
                if self.inputs['rar'][0] == 'True':
                    lowest_count = int(min(self.otu[x].sum(axis='sample')))
                else:
                    lowest_count = int(self.inputs['rar'][0])
                data = self.otu[x].matrix_data
                data = csr_matrix.todense(data)
                keep_samples = list()
                mincount = np.sum(data, axis=0)
                for y in range(mincount.shape[1]):
                    if mincount.item(y) >= lowest_count:
                        keep_samples.append(self.otu[x]._sample_ids[y])
                keep = self.otu[x].filter(keep_samples, axis="sample", inplace=False)
                batchcopy.otu[x] = keep.subsample(n=lowest_count, axis='sample')
                self.log['rarefaction'] = str(lowest_count) + " counts"
            except Exception:
                logger.error("Unable to rarefy file", exc_info=True)
            for x in list(self.otu):
                self.otu[x] = batchcopy.otu[x]

    def split_biom(self, mode="sample", *args):
        """
        Splits bioms into several subfiles according to
        sample metadata variable. Source: biom-format.org.
        The original file is preserved, so returned files
        include the split- and non-split files.
        """
        inputs = self.inputs
        part_f = lambda id_, md: md[inputs['split'][0]]
        new_dict = deepcopy(self.otu)
        if type(self.otu) is not dict:
            raise ValueError("Split_biom requires a dictionary of biom files to be supplied.")
        for x in list(self.otu):
            try:
                biomtab = new_dict[x]
                if inputs['split'][0] not in biomtab._sample_metadata[1]:
                    if inputs['split'][0] is not 'TRUE':
                        raise ValueError("Sample metadata does not contain this header!")
                new_tables = biomtab.partition(part_f, axis='sample')
                for new in new_tables:
                    key = x + '_' + new[0]
                    new_dict[key] = new[1]
            except Exception:
                logger.error("Failed to split files", exc_info=True)
        self.otu = new_dict
        self.log['split_variables'] = inputs['split'][0]

    def cluster_biom(self):
        """First normalizes bioms so clustering is not affected,
        performs transformation and then applies clustering.
        Note that the returned biom files are not normalized,
        this is just for the clustering process.
        Many network inference tools require absolute counts.
        Silhouette score is used to determine the optimal
        number of clusters.
        Clustering adds metadata info to the samples.
        Splitting according to cluster ID is done
        by wrapping the split_biom function.
        """
        inputs = self.inputs
        if inputs['nclust'] is not None:
            nums = list(range(2, (int(inputs['nclust'][0]) + 1)))
        else:
            nums = list(range(2,5))
        new_dict = {}
        if type(self.otu) is not dict:
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
                if inputs['cluster'][0] == 'K-means':
                    for i in nums:
                        clusters = KMeans(i).fit_predict(data)
                        silhouette_avg = silhouette_score(data, clusters)
                        sh_score.append(silhouette_avg)
                    topscore = int(np.argmax(sh_score) + 1)
                    bestcluster = KMeans(topscore).fit_predict(data)
                # DBSCAN clustering, automatically finds optimal cluster size
                if inputs['cluster'][0] == 'DBSCAN':
                    bestcluster = DBSCAN().fit_predict(data)
                    topscore = len(set(bestcluster)) - (1 if -1 in bestcluster else 0)
                # Gaussian Mixture Model (gmm) probability distribution
                if inputs['cluster'][0] == 'Gaussian':
                    for i in nums:
                        fit = GaussianMixture(i).fit(data)
                        clusters = fit.predict(data)
                        silhouette_avg = silhouette_score(data, clusters)
                        sh_score.append(silhouette_avg)
                    topscore = int(np.argmax(sh_score) + 1)
                    bestfit = GaussianMixture(topscore).fit(data)
                    bestcluster = bestfit.predict(data)
                # Spectral Clustering
                if inputs['cluster'][0] == 'Spectral':
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
                        new_dict[x]._sample_metadata[j]['cluster'] = inputs['cluster'][0] + '_' + str(i)
                self.otu = new_dict
                self.log['cluster'] = str(topscore) + " clusters with " + inputs['cluster'][0]
                if inputs['split'] is not None:
                    if inputs['split'][0] == 'TRUE':
                        inputs['split'][0] = 'cluster'
                        self.split_biom()
            except Exception:
                logger.error("Error occurred when clustering samples", exc_info=True)


def _data_bin(biomfile, taxnum):
    """While the BIOM file collapse function collapses counts, it stores taxonomy in a manner
    that is not compatible with the OTU taxonomy; taxonomy is stored as observation ID and not as
    observation metadata.
    Here, new observation IDs are generated for agglomerated data that are a combination
    of the old IDs. Moreover, the taxonomy of the agglomerated data is concatenated to
    the agglomeration level and added to the observation metadata. """
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
    for number in range(len(obs_ids)):
        id = obs_ids[number]
        num_id = collapsed._obs_index[id]
        old_ids = agglom_ids[num_id]
        new_obs_ids[number] = "agglom-" + str(number)
    collapsed._observation_ids = new_obs_ids
    return collapsed

