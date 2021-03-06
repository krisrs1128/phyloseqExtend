% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/spectra-utils.R
\name{get_peaks_list}
\alias{get_peaks_list}
\title{Wrapper for MALDIquant detectPeaks}
\usage{
get_peaks_list(specmat, peak_opts)
}
\arguments{
\item{spectra_object}{A spectra object with a nonempty specmat slot. If this
slot has column names, they will be used in the peak detection step.
Otherwise, placeholder intensities will be used.}

\item{...}{Additional options to detecPeaks in MALDIquant.}
}
\value{
peaks_list A list whose i^th element contains the indices
of peaks in the i^th sample
}
\description{
Wrapper for MALDIquant detectPeaks
}

