
################################################################################
# Preprocess different components before they are input to the spectra object.
################################################################################

# contingency-tables -----------------------------------------------------------
#' @title Merge default contingency table processing options
#' @param opts A partially filled list specifying options to contingency table
#' processing. Those that are not specified will be merged according to their
#' defaults below. These are,
#'  $id_col: Is the first column of the data a nonnumeric identifier column?
#'   Defaults to TRUE.
#'  $filter_ix: A function that can be applied to the full matrix (excluding a
#'   potential id column), and which returns the indices of the row filtered
#'   version. Defaults to filter_nonzero().
#'  $transform: A function that can be applied to the full matrix, transforming
#'   its values. Defaults to log().
#' @return A modified version of opts, with defaults filled in.
#' @export
merge_ctable_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$id_col <- TRUE
  default_opts$filter_ix <- filter_nonzero
  default_opts$transform <- function(x) { x }
  modifyList(default_opts, opts)
}

#' @title Wrapper for contingency table processing
#' @param X A matrix / data.frame / data.table with the samples as rows, and
#' with ids as either row names or as the first column [opts$id_col in
#' merge_ctable_opts()].
#' @param opts A list of processing options. Choices and their defaults are
#' detailed in merge_ctable_opts().
#' @return X A processed version of the input data X.
#' @importFrom data.table data.table
#' @export
process_ctable <- function(X, opts = list()) {
  init_class <- class(X)
  X <- data.table(X, keep.rownames = TRUE)

  opts <- merge_ctable_opts(opts)

  # extract just the numeric part of the data
  if(opts$id_col) {
    ids <- X[, 1, with = F]
    X <- X[, -1, with = F]
  }

  # perform filtering / transformations
  ix <- opts$filter_ix(X)
  X <- opts$transform(X[ix, ])

  # put back into original form
  if(opts$id_col) {
    X <- cbind(ids[ix], X)
  }
  class(X) <- init_class
  X
}

# spectra ----------------------------------------------------------------------
#' @title Merge default spectra processing options
#' @param opts A partially filled list specifying options to spectra processing.
#' Those that are not specified will be merged according to their defaults
#' below. These are,
#'  $align_opts Arguments passed to alignSpectra() in the MALDIquant package.
#'   Defaults to the default options for that function. \cr
#'  $baseline_opts Arguments passed to the removeBaseline() function in the
#'   MALDIquant package. Defaults to the default options for that function. \cr
#'  $binary Return a sample x potential-peak indicator matrix? Or just the peaks
#'   intensities [with zeros, if it was not a peak]? Defaults to FALSE. \cr
#'  $peak_opts Arguments passed to detectPeak() in the MALDIquant package.
#'   Defaults to the default options for that function. \cr
#'  $thresh_max The maximum value for any sample, above which we discard it as
#'   an outlier. Defaults to Inf (don't discard anything) \cr
#' @return opts The version of opts with unspecified options filled in.
#' @export
merge_spectra_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$align_opts <- list()
  default_opts$baseline_opts <- list()
  default_opts$binary <- FALSE
  default_opts$peak_opts <- list()
  default_opts$thresh_max <- Inf
  modifyList(default_opts, opts)
}

#' @title Wrapper for processing spectra matrices
#' @description This function 1) removes outliers, 2) removes spectra baseline,
#' 3) aligns spectra, and 3) detects peaks, across samples in a spectrum.
#' @param X A spectra matrix with samples on rows.
#' @param opts A potentially partially specified list of processing options. See
#' merge_spectra_opts() for options and defaults.
#' @return A list with the followign objects
#'   $X The spectra matrix with peaks aligned. \cr
#'   $peaks_ix A list whose i^th element give the indices of peaks in the i^th
#'    sample. \cr
#'   $X_peaks A matrix whose columns are the positions where a peak was detected
#'    in at least one sample. \cr
#'   $X_peaks_zeros The same as X_peaks, but with zeros for samples where that
#'    given position was not identified as a peak.
#' @export
process_spectra <- function(X, opts = list()) {
  opts <- merge_spectra_opts(opts)

  # remove any outliers
  row_maxes <- apply(X, 1, max)
  X <- X[row_maxes < opts$thresh_max, ]
  if(any(row_maxes > opts$thresh_max)) {
    message(sprintf("Discarding rows %s as outliers \n",
                    paste0(which(row_maxes > opts$thresh_max), collapse = ", ")))
  }

  # do all alignment / peak detection
  X <- remove_baseline(X, opts$baseline_opts)
  X <- align_spectra(X, opts$align_opts)
  peaks_ix <- get_peaks_list(X, opts$peak_opts)

  # get zeros matrices
  X_peaks <- get_spectra_at_peaks(X, peaks_ix)
  X_peaks_zeros <- get_spectra_at_peaks_zeros(X, peaks_ix, opts$binary)
  list(X = X, peaks_ix = peaks_ix, X_peaks = X_peaks,
       X_peaks_zeros = X_peaks_zeros)
}
