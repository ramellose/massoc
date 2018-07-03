"""
This file contains all testing functions for netwrap.
"""

__author__ = 'Lisa Rottjers'
__maintainer__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import os
import random
import unittest
from copy import deepcopy
from subprocess import call

import biom
from massoc.scripts.batch import Batch
from massoc.scripts.netwrap import Nets

import massoc
from massoc.scripts.main import run_networks

random.seed(7)

# replace with proper test file!
# looks like there is an issue with agglomerated files
testloc = list()
testloc.append(os.path.dirname(massoc.__file__)[:-6] + 'data\\arctic_soils.biom')
testloc = [x.replace('\\', '/') for x in testloc]
testbiom = {"test": biom.load_table(testloc[0])}

inputs = {'biom_file': None,
          'cluster': None,
          'otu_meta': None,
          'otu_table': ['otu_bananas.txt'],
          'prefix': None,
          'sample_data': None,
          'split': 'BODY_SITE',
          'tax_table': ['tax_bananas.txt'],
          'fp': [testloc[0][:-18]],
          'tools': ['spiec-easi'],
          'spiec': None,
          'conet': None,
          'spar_pval': None,
          'spar_boot': None,
          'levels': ['family'],
          'name': ['test'],
          'min': None,
          'rar': None}
batch = Batch(testbiom, inputs)
netbatch = Nets(batch)
netbatch.write_bioms()

class TestNetWrap(unittest.TestCase):
    """Tests netwrap.
    More specifically, checks ability to call network inference tools.
    Netwrap requires access to these tools, and they need to be installed.
    """

    def test_spar(self):
        """Check if the SparCC function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.run_spar(boots=10)
        self.assertEqual(len(testnets.networks), 1)

    def test_conet(self):
        """Check if the CoNet function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.run_conet()
        self.assertEqual(len(testnets.networks), 1)

    def test_spiec(self):
        """Check if the SPIEC-EASI function call works
        by testing length of Nets.networks."""
        testnets = deepcopy(netbatch)
        testnets.run_spiec()
        self.assertEqual(len(testnets.networks), 1)

    def test_run_spiec_settings(self):
        """Checks if a supplied .R file can be used to run SPIEC-EASI.
        """
        location = os.path.dirname(massoc.__file__)[:-6] + 'tests\\spieceasi_test.R'
        location = location.replace('\\', '/')
        execspiec = '#!/usr/bin/Rscript\n\n#\' @title Calls SPIEC-EASI from command line\n#\' ' \
               '@description Calls SPIEC-EASI and saves output as a graphml file.\n#\' ' \
               '@details Authors: Zachary D. Kurtz et al. SPIEC-EASI is available at: ' \
               'https://github.com/zdk123/SpiecEasi.' \
               '\n#\' Note that the settings for SPIEC-EASI can be adjusted from within this file.\n#\' ' \
               'The most important settings are: Meinshausen-Buhlmann vs graphical lasso algorithm\n#\' ' \
               'icov.select.params ' \
               '(increasing this setting will improve the sparsity of the output network)\n#\'\n#\' ' \
               '@param i input biom file\n#\' @param o output filename + path\n#\' @example spieceasi.R -fp ' \
               '<filepath> -o <filepath>\n#\' @export\n\nlibrary(biom)\nlibrary(docopt)\nlibrary(SpiecEasi)\n\n' \
               'doc = \'Usage:\n     spieceasi.r [-i input] [-o network]\n\nOptions:\n      ' \
               '-i Filepath for input biom file\n      ' \
               '-o Filepath for output network\n\'\n\nopts = docopt(doc)\n\n' \
               'file = read_biom(opts$i)\ncounttab = t(as.matrix(biom_data(file)))\n' \
               '# change SPIEC-EASI  method: Meinshausen-Buhlmann (mb) or graphical lasso (gl)\n' \
               'method = "glasso"\n' \
               '# number of STARS iterations is set with icov.select.params\n' \
               'spiec.out = spiec.easi(counttab, method, icov.select.params=list(rep.num=30))\n' \
               'if (method == "mb"){\n  adj = as.matrix(getOptBeta(spiec.out))\n}\n' \
               'if (method == "glasso"){\n  adj = as.matrix(getOptMerge(spiec.out))\n}\n' \
               'colnames(adj) = colnames(counttab)\nrownames(adj) = colnames(counttab)\n' \
               'write.table(adj, opts$o, sep="\\t")\n'
        file = open(location, "w")
        file.write(execspiec)
        file.close()
        newsets = {'biom_file': [testloc[0]],
                   'cluster': None,
                   'otu_meta': None,
                   'otu_table': ['otu_bananas.txt'],
                   'prefix': None,
                   'sample_data': None,
                   'split': ['BODY_SITE'],
                   'tax_table': ['tax_bananas.txt'],
                   'fp': [testloc[0][:-18]],
                   'tools': None,
                   'spiec': [location],
                   'conet': None,
                   'spar_pval': None,
                   'spar_boot': None,
                   'levels': ['family'],
                   'name': ['test']}
        spiecbatch = Batch(testbiom, newsets)
        spiecbatch = Nets(spiecbatch)
        run_networks(spiecbatch)
        call('rm ' + location[0])
        self.assertEqual(len(spiecbatch.networks), 1)


if __name__ == '__main__':
    unittest.main()
