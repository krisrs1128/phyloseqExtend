
################################################################################
# Utilities for processing spectra
################################################################################

# align-and-detect -------------------------------------------------------------

#' @title Wrapper for MALDIquant detectPeaks
#' @param spectra_object A spectra object with a nonempty specmat slot. If this
#' slot has column names, they will be used in the peak detection step.
#' Otherwise, placeholder intensities will be used.
#' @param ... Additional options to detecPeaks in MALDIquant.
#' @return peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @importFrom MALDIquant createMassSpectrum detectPeaks
#' @export
get_peaks_list <- function(specmat, peak_opts) {
  ppm <- get_ppm(specmat)
  apply(specmat, 1, function(x) {
    ppm_vals <- detectPeaks(createMassSpectrum(ppm, x))@mass
    which(ppm %in% ppm_vals)
  })
}

#' @title Get ppm values associated with a spectra matrix
#' @param specmat A matrix whose column names will be used to derive ppms
#' @export
get_ppm <- function(specmat) {
  ppm <- as.numeric(colnames(specmat))
  if(any(is.na(ppm)) | length(ppm) == 0) {
    ppm <- seq_len(ncol(specmat))
    warning("Using placeholder ppm values")
  }
  ppm
}

#' @title Wrapper for removeBaseline and alignSpectra in MALDIquant
#' @param spectra_object A spectra object with a nonempty specmat slot.
#' @param baseline_opts A list of arguments to pass to removeBaseline() in
#' MALDIquant.
#' @importFrom MALDIquant createMassSpectrum removeBaseline
#' @return baselined_peaks A matrix whose ij^th element is the value of the j^th
#' unique peak in the i^th unique sample.
#' @export
remove_baseline <- function(specmat, baseline_opts = list()) {
  ppms <- get_ppm(specmat)
  nmr_list <- apply(specmat, 1, function(x) createMassSpectrum(ppms, x))
  nmr_list <- do.call(removeBaseline, c(list(object = nmr_list), baseline_opts))
  baselined_peaks <- do.call(rbind, lapply(nmr_list, function(x) x@intensity))
  rownames(baselined_peaks) <- rownames(specmat)
  baselined_peaks
}

#' @title Wrapper for removeBaseline and alignSpectra in MALDIquant
#' @param spectra_object A spectra object with a nonempty specmat slot.
#' @param baseline_opts A list of arguments to pass to removeBaseline() in
#' MALDIquant.
#' @param align_opts A list of arguments to pass to alignSpectra in MALDIquant.
#' @importFrom MALDIquant createMassSpectrum alignSpectra
#' @return aligned_peaks A matrix whose ij^th element is the value of the j^th
#' unique peak in the i^th unique sample.
#' @export
align_spectra <- function(specmat, align_opts = list()) {
  ppms <- get_ppm(specmat)
  nmr_list <- apply(specmat, 1, function(x) createMassSpectrum(ppms, x))
  nmr_list <- do.call(alignSpectra, c(list(spectra = nmr_list), align_opts))
  aligned_peaks <- do.call(rbind, lapply(nmr_list, function(x) x@intensity))
  rownames(aligned_peaks) <- rownames(specmat)
  aligned_peaks
}

# outliers ---------------------------------------------------------------------
#' @title Remove Outliers
#' @param physeq A phyloseqExtend object with a nonempty spectra slot
#' @param thresh Any spectra with maximum value above thresh will be discarded.
#' @return physeq  Version of physeq object with outlier samples removed (across
#' all components in the input physeq object).
#' @importFrom phyloseq sample_data otu_table
#' @export
remove_outlier_spectra <- function(physeq, thresh) {
  spectra_matrix <- spectra(physeq)@specmat
  stopifnot(!is.null(spectra_matrix))
  max_vals <- apply(spectra_matrix, 1, max)

  keep_samples <- rownames(spectra_matrix)[which(max_vals <= thresh)]

  if(length(keep_samples) == 0) {
    stop(sprintf("Smallest spectrum max is %s. This threshold choice would remove all samples.", min(max_vals)))
  }
  physeq@spectra <- subset_samples_spectra(keep_samples, physeq@spectra)
  physeq@sam_data <- sample_data(physeq)[keep_samples,, drop = F ]
  if(physeq@otu_table@taxa_are_rows) {
    physeq@otu_table@.Data <- physeq@otu_table@.Data[, keep_samples, drop = F]
  } else {
    physeq@otu_table@.Data <- physeq@otu_table@.Data[keep_samples, , drop = F]
  }
  physeq
}

# subsampling ------------------------------------------------------------------
#' @title Extract Spectra at Peak Positions, with zeros for non-peaks
#' @param specmat A spectra matrix with samples on rows.
#' @param peaks_list A list whose i^th element contains the indices of peaks in
#' the i^th sample
#' @param binary Should the actual spectrum intensity be returned, or just an
#' indicator of whether the position was a peak?
#' @return peaks A matrix containing spectra values filtered to samples and
#' measurement indices containining peaks, and placing zeros if a
#' measurement was not called a peak.
#' @importFrom data.table dcast melt
#' @importFrom magrittr %>%
#' @export
get_spectra_at_peaks_zeros <- function(specmat, peaks_list, binary = FALSE) {
  names(peaks_list) <- NULL # so melt doesn't try to guess id variables
  peaks_ix  <- melt(peaks_list)
  colnames(peaks_ix) <- c("index", "sample")
  peaks <- peaks_ix %>%
    cbind(value = t(specmat)[as.matrix(peaks_ix)]) %>%
    dcast(sample ~ index, fill = 0)

  # Give appropriate sample and spectra names
  rownames(peaks) <- peaks$sample
  peaks$sample <- NULL
  colnames(peaks) <- colnames(specmat)[as.numeric(colnames(peaks))]
  peaks <- as.matrix(peaks)

  # convert to binary, if desired
  if(binary) {
    peaks <- peaks > 0
    class(peaks) <- "numeric"
  }
  rownames(peaks) <- rownames(specmat)
  peaks
}

#' @title Extract Spectra at Peak Positions
#' @param peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @return peaks A matrix containing spectra values filtered to samples
#' and measurement indices containining peaks.
#' @export
get_spectra_at_peaks <- function(specmat, peaks_list) {
  specmat[, unique(unlist(peaks_list))]
}

#' @title Subsample spectra and filter spectra to range based on column names
#' @description For long spectra we may want to either 1) subsample to lower
#' resolution, so that plots can be made more quickly, or 2) filter down to
#' a specified range of indices, to understand smaller scale phenomena. The
#' first task here can be accomplished by setting subsample_frac; for example,
#' setting subsample_frac = 1/2 will return the spectrum matrix filtered to
#' every other column. For the second task, we require the columns of
#' \code{spectra} to be numeric, so we can filter cols on whether they are
#' between x_min and x_max.
#' @param spectra_matrix An object of class matrix, containing the spectral
#'  samples as rows.
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#'  the spectrum is very long, but can lead to missed peaks.
#' @return spectra_matrix A matrix with the same number of rows as the input,
#'  but with filtered columns.
#' @export
subsample_spectra_cols <- function(spectra_matrix, subsample_frac = 1,
                                   x_min = NULL, x_max = NULL) {
  keep_ix <- seq(1, ncol(spectra_matrix), by = 1 / subsample_frac)
  spectra_names <- colnames(spectra_matrix)
  spectra_matrix <- spectra_matrix[, keep_ix, drop = F]
  colnames(spectra_matrix) <- spectra_names[keep_ix]
  if(!is.null(x_min)) {
    cols_val <- as.numeric(colnames(spectra_matrix))
    spectra_matrix <- spectra_matrix[, which(cols_val >= x_min), drop = F]
  }
  if(!is.null(x_max)) {
    # have to recompute cols_val, in case cols were dropped above
    cols_val <- as.numeric(colnames(spectra_matrix))
    spectra_matrix <- spectra_matrix[, which(cols_val <= x_max), drop = F]
  }
  spectra_matrix
}

# misc -------------------------------------------------------------------------
#' @title Interpolate Duplicated Values
#' @description When working with spectra, it is often necessarily to
#' interpolate duplicated ppm values. For example, c(1, 1, 2) should be turned
#' into c(1, 1.5, 2).
#' @param x The vector to interpolate.
#' @return The interpolated vector.
#' @examples
#' interpolate_duplicates(c(1, 1, 2))
#' @export
interpolate_duplicates <- function(x) {
  interpolation <- vector(length = length(x))
  cur_ix <- 1
  unique_x <- unique(x)
  for(i in seq_along(unique_x)) {
    cur_x <- unique_x[i]
    if(i < length(unique_x)) {
      next_x <- unique_x[i + 1]
    } else {
      next_x <- unique_x[i] + (unique_x[i] - unique_x[i - 1])
    }
    n_prev_ix <- sum(x == unique_x[i])
    interpolation[cur_ix:(cur_ix + n_prev_ix - 1)] <- seq(cur_x, next_x, length.out = n_prev_ix + 1)[-(n_prev_ix + 1)] # don't want a term equal to "to"
    cur_ix <- cur_ix + n_prev_ix
  }
  interpolation
}
