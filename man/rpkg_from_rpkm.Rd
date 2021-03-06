% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/normalize.R
\name{rpkg_from_rpkm}
\alias{rpkg_from_rpkm}
\title{Compute RPKGs from RPKMs}
\usage{
rpkg_from_rpkm(rpkm, lib, ags)
}
\arguments{
\item{rpkm}{A matrix of RPKM's, with genes on rows and samples on columns.}

\item{lib}{A vector of library sizes for each sample.}

\item{ags}{A vector of estimated average genome sizes for each sample.}
}
\value{
rpkg The matrix of reads normalized according to RPKG. The definition
is [{# mapped reads} / {gene length}] / [library size / AGS]. The quantity
library size / AGS is called the number of genome equivalents for that
sample.
}
\description{
Compute RPKGs from RPKMs
}
\references{
"Average genome size estimation improves comparative metagenomics
and sheds light on the functional ecology of the human microbiome" by
Nayfach and Pollard.
}

