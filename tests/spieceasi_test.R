#!/usr/bin/Rscript

#' @title Calls SPIEC-EASI from command line
#' @description Calls SPIEC-EASI and saves output as a graphml file.
#' @details Authors: Zachary D. Kurtz et al. SPIEC-EASI is available at: https://github.com/zdk123/SpiecEasi.
#' Note that the settings for SPIEC-EASI can be adjusted from within this file.
#' The most important settings are: Meinshausen-Buhlmann vs graphical lasso algorithm
#' icov.select.params (increasing this setting will improve the sparsity of the output network)
#'
#' @param i input biom file
#' @param o output filename + path
#' @example spieceasi.R -fp <filepath> -o <filepath>
#' @export

require(biom)
require(docopt)
require(SpiecEasi)

doc = 'Usage:
     spieceasi.r [-i input] [-o network]

Options:
      -i Filepath for input biom file
      -o Filepath for output network
'

opts = docopt(doc)

file = read_biom(opts$i)
counttab = t(as.matrix(biom_data(file)))
# change SPIEC-EASI  method: Meinshausen-Buhlmann (mb) or graphical lasso (gl)
method = "glasso"
# number of STARS iterations is set with icov.select.params
spiec.out = spiec.easi(counttab, method, icov.select.params=list(rep.num=30))
if (method == "mb"){
  adj = as.matrix(getOptBeta(spiec.out))
}
if (method == "glasso"){
  adj = as.matrix(getOptMerge(spiec.out))
}
colnames(adj) = colnames(counttab)
rownames(adj) = colnames(counttab)
write.table(adj, opts$o, sep="\t")
