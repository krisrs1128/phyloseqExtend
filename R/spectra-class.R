# build-spectra-object ----------------------------------------------------
#' @title Build or access an object of class spectra
#'
#' @param spectra_matrix An object of class matrix, containing the spectral
#'  samples as rows.
#' @param peaks A logical indicating whether the matrix columns corresponds to
#'  wavelengths (FALSE, for raw spectra) or just peaks (TRUE, for spectra
#'  that have been processed). Defaults to FALSE.
#'
#' @return spectra An object of class spectra, containing the input spectra
#' matrix as the main data component.
#'
#' @docType methods
#' @export

setGeneric("spectra", function(object, peaks = FALSE) {
  standardGeneric("spectra")
})

# If a phyloseqExtend object with a spectra slot is input, return that slot
setMethod("spectra", "phyloseqExtend", function(object) {
  spectra_object <- slot(object, "spectra")
  if(is.null(spectra_object)) {
    warning("No spectra available.")
  }
  return (spectra_object)
})

# If a matrix is input, create a new spectra object from that matrix
setMethod("spectra", "matrix", function(object, peaks = FALSE) {
  new("spectra", object, peaks = peaks)
})

# If a data.frame or data.table is input, convert to a matrix
setMethod("spectra", "data.frame", function(object, peaks = FALSE) {
  spectra(as(object, "matrix"), peaks = peaks)
})

# call-peaks --------------------------------------------------------------
#' @title Wrapper for speaq::detectSpecPeaks
#'
#' @param spectra A spectra object containing raw spectra reads
#' @return peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @importFrom speaq detectSpecPeaks
get_peaks_list <- function(spectra, ...) {
  peaks_list <- speaq::detectSpecPeaks(spectra@.Data, ...)
  return (peaks_list)
}

# get-spectra-at-peaks ----------------------------------------------------
#' @title Extract Spectra at Peak Positions
#'
#' @param peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#'
#' @return peaks A matrix containing spectra values filtered to samples
#'  and measurement indices containining peaks.
#'
#' @importFrom reshape2 dcast melt
#' @importFrom magrittr %>%
get_spectra_at_peaks <- function(peaks_list) {
  peaks_ix  <- reshape2::melt(peaks_list)
  colnames(peaks_ix) <- c("index", "sample")
  peaks <- peaks_ix %>%
    cbind(value = t(spectra)[as.matrix(peaks_ix)]) %>%
    reshape2::dcast(sample ~ index, fill = 0)

  # Give appropriate sample and spectra names
  rownames(peaks) <- peaks$sample
  peaks$sample <- NULL
  colnames(peaks) <- colnames(spectra)[as.numeric(colnames(peaks))]
  peaks <- as.matrix(peaks)
  return (peaks)
}

#' @title Wrapper for dohCluster in speaq
#'
#' @param spectra A spectra object containing raw spectra reads
#' @param peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#'
#' @return aligned_peaks A matrix whose ij^th element is the value of the j^th
#' unique peak in the i^th unique sample.
align_peaks <- function(spectra, peaks_list, ...) {
  ref_peak <- speaq::findRef(peaks_list)
  aligned_peaks <- speaq::dohCluster(spectra, peakList = peaks_list,
                              refInd = ref_peak$refInd, ...)
  return (aligned_peaks)
}

#' @title Remove Outliers
#'
#' @param physeq A phyloseqExtend object with a nonempty spectra slot
#' @param thresh Any spectra with maximum value above thresh will be discarded.
#'
#' @return None. Removes the outlier samples from the spectra component of the
#' input physeq object.
remove_outlier_spectra <- function(physeq, thresh) {
  spectra <- spectra(physeq)
  stopifnot(!is.null(spectra))
  max_vals <- apply(spectra, 1, max)
  physeq@spectra@.Data <- spectra[which(max_vals <= thresh), ]
  return (physeq)
}

#' @title Subsample spectra and filter spectra to range based on column names
#'
#' @description For long spectra we may want to either 1) subsample to lower
#' resolution, so that plots can be made more quickly, or 2) filter down to
#' a specified range of indices, to understand smaller scale phenomena. The
#' first task here can be accomplished by setting subsample_frac; for example,
#' setting subsample_frac = 1/2 will return the spectrum matrix filtered to
#' every other column. For the second task, we require the columns of
#' \code{spectra} to be numeric, so we can filter cols on whether they are
#' between x_min and x_max.
#'
#' @param spectra An object of class matrix, containing the spectral samples as
#'  rows.
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#'  the spectrum is very long, but can lead to missed peaks.
#'
#' @return spectra A matrix with the same number of rows as the input, but with
#'  filtered columns.
subsample_spectra_cols <- function(spectra, subsample_frac = 1, x_min = NULL,
                                   x_max = NULL) {
  spectra <- spectra[, seq(1, ncol(spectra), by = 1 / subsample_frac)]
  if(!is.null(x_min)) {
    cols_val <- as.numeric(colnames(spectra))
    spectra <- spectra[, which(cols_val >= x_min)]
  }
  if(!is.null(x_max)) {
    # have to recompute cols_val, in case cols were dropped above
    cols_val <- as.numeric(colnames(spectra))
    spectra <- spectra[, which(cols_val <= x_max)]
  }
  return (spectra)
}
