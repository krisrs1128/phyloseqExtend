% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{plot_spectra}
\alias{plot_spectra}
\title{Plot Spectra Intensities}
\usage{
plot_spectra(physeq, method = "speaq", log_scale = FALSE, x_min = NULL,
  x_max = NULL, plot_title = NULL, ...)
}
\arguments{
\item{physeq}{A phyloseq object containing the spectrum to plot}

\item{method}{The R package to use for the plot. Currently only supports
"speaq" or "ggplot2".}

\item{log_scale}{Should the intensities be plotted on a log scale?}

\item{x_min}{What is the minimum index to display?}

\item{x_max}{What is the maximum index to display?}
}
\value{
If "ggplot2" is used, then return the ggplot2 object containing the
 plot. Otherwise returns null.
}
\description{
If spectra(physeq)[["peaks"]] is FALSE, then this will plot the
raw intensities of each spectrum in the sample. Otherwise, we only plot
the positions and heights of peaks for each sample.
}
