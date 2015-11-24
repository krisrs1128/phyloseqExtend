
################################################################################
# Tools for normalization of genomic data.
################################################################################

#' @title Compute RPKGs from RPKMs
#' @param rpkm A matrix of RPKM's, with genes on rows and samples on columns.
#' @param lib A vector of library sizes for each sample.
#' @param ags A vector of estimated average genome sizes for each sample.
#' @return rpkg The matrix of reads normalized according to RPKG. The definition
#' is [{# mapped reads} / {gene length}] / [library size / AGS]. The quantity
#' library size / AGS is called the number of genome equivalents for that
#' sample.
#' @references "Average genome size estimation improves comparative metagenomics
#' and sheds light on the functional ecology of the human microbiome" by
#' Nayfach and Pollard.
#' @export
rpkg_from_rpkm <- function(rpkm, lib, ags) {
  pk <- rpkm %*% diag(lib) / 1e6
  rpkg <- pk %*% diag(lib / ags)
  colnames(rpkg) <- colnames(rpkm)
  rpkg
}
