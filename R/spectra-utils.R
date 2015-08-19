# call-peaks --------------------------------------------------------------
#' @title Wrapper for speaq::detectSpecPeaks
#'
#' @param spectra_object A spectra object with a nonempty specmat slot
#' @return peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @importFrom speaq detectSpecPeaks
#' @export
get_peaks_list <- function(spectra_object, ...) {
  peaks_list <- speaq::detectSpecPeaks(spectra_object@spectra@specmat, ...)
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
#'
#' @export
get_spectra_at_peaks <- function(specmat, peaks_list) {
  peaks_ix  <- reshape2::melt(peaks_list)
  colnames(peaks_ix) <- c("index", "sample")
  peaks <- peaks_ix %>%
    cbind(value = t(specmat)[as.matrix(peaks_ix)]) %>%
    reshape2::dcast(sample ~ index, fill = 0)

  # Give appropriate sample and spectra names
  rownames(peaks) <- peaks$sample
  peaks$sample <- NULL
  colnames(peaks) <- colnames(specmat)[as.numeric(colnames(peaks))]
  peaks <- as.matrix(peaks)
  return (peaks)
}

# wrap-peak-alignment -----------------------------------------------------
#' @title Wrapper for dohCluster in speaq
#'
#' @param spectra_object A spectra object with a nonempty specmat slot.
#' @param peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#'
#' @return aligned_peaks A matrix whose ij^th element is the value of the j^th
#' unique peak in the i^th unique sample.
align_peaks <- function(spectra_object, peaks_list, ...) {
  ref_peak <- speaq::findRef(peaks_list)
  aligned_peaks <- speaq::dohCluster(spectra_object@specmat,
                                     peakList = peaks_list,
                                     refInd = ref_peak$refInd, ...)
  return (aligned_peaks)
}

# remove-outliers ---------------------------------------------------------
#' @title Remove Outliers
#'
#' @param physeq A phyloseqExtend object with a nonempty spectra slot
#' @param thresh Any spectra with maximum value above thresh will be discarded.
#'
#' @return None. Removes the outlier samples from the spectra, sample_data and
#'  otu_table components of the input physeq object.
#'
#' @importFrom phyloseq sample_data otu_table
#'
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
  return (physeq)
}

# subsample-cols ----------------------------------------------------------
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
#' @param spectra_matrix An object of class matrix, containing the spectral
#'  samples as rows.
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#'  the spectrum is very long, but can lead to missed peaks.
#'
#' @return spectra_matrix A matrix with the same number of rows as the input,
#'  but with filtered columns.
#' @export
subsample_spectra_cols <- function(spectra_matrix, subsample_frac = 1,
                                   x_min = NULL, x_max = NULL) {
  spectra_matrix <- spectra_matrix[, seq(1, ncol(spectra_matrix), by = 1 / subsample_frac), drop = F]
  if(!is.null(x_min)) {
    cols_val <- as.numeric(colnames(spectra_matrix))
    spectra_matrix <- spectra_matrix[, which(cols_val >= x_min), drop = F]
  }
  if(!is.null(x_max)) {
    # have to recompute cols_val, in case cols were dropped above
    cols_val <- as.numeric(colnames(spectra_matrix))
    spectra_matrix <- spectra_matrix[, which(cols_val <= x_max), drop = F]
  }
  return (spectra_matrix)
}

# detect-and-align --------------------------------------------------------
#' @title Wrapper for peak detection and alignment
#'
#' @description A common preprocessing step in work with spectra is peak
#'  alignment. Here, we wrap the peak detection, reference finding, and
#'  hierarchical clustering peak alignment steps from the speaq package.
#'
#' @param physeq A phyloseqExtend object with a nonempty spectra slot
#' @param detectSpecPeaksOpts A list containing options to pass into the
#'  detectSpecPeaks function in package speaq.
#' @param dohClusterOpts A list containing options to pass into the
#'  dohCluster function in package speaq.
#' @param binary If TRUE, returns 0/1 for whether there is a peak in the sample
#'  at a particular position.
#'
#' @return physeq The input physeq object, with the specmat element of the
#'  spectra slot modified so that all peaks are aligned, and with the peakmat
#'  element included.
#' @export
detect_and_align_peaks <- function(physeq, detectSpecPeaksOpts, dohClusterOpts,
                                   binary = F) {

  # warning in case there are peaks
  if(!is.null(spectra(physeq)@peakmat)) {
    warning("Overwriting existing peakmat slot.")
  }

  # detect peaks
  detectSpecPeaksOpts <- modifyList(
    detectSpecPeaksOpts,
    list(spectra_object = spectra(physeq))
  )
  message("Detecting peaks...")
  peaks <- do.call(get_peaks_list, detectSpecPeaksOpts)

  # align peaks
  message("Aligning peaks to reference...")
  dohClusterInput <- modifyList(
    dohClusterOpts,
    list(spectra_object = spectra(physeq), peaks_list = peaks)
  )

  # combine results
  specmat <- do.call(align_peaks, dohClusterInput)
  peakmat = get_spectra_at_peaks(spectra(physeq)@specmat, peaks)

  if(binary) {
    peakmat[peakmat > 0] <- 1
  }

  rownames(peakmat) <- rownames(specmat)
  result_list <- list(specmat = specmat,
                      peakmat = peakmat,
                      adjmat = spectra(physeq)@adjmat)
  physeq@spectra <- spectra(result_list)
  return (physeq)
}

# input_peak_slot ---------------------------------------------------------
#' @title Shrink a Full Spectrum to Just Peaks
#'
#' @description Rather than working with full spectrum reads, it is often
#'  convenient to work with a set of noteworthy peaks. In particular, we may
#'  want to first identify peaks across samples according to some peak calling
#'  algorithm, and then simplify the problem by setting
#'
#' @param physeq A phyloseqExtend object with a nonempty spectra slot.
#' @param peakmat A matrix whose rows correspond to samples in the specmat and
#'  whose columns index unique samples in the data.
#' @export
input_peaks <- function(physeq, binary = FALSE, peakmat = NULL, ...) {

  # If we don't have any peaks yet, extract them
  if(is.null(peakmat)) {
    peaks_list <- get_peaks_list(spectra(physeq), ...)
    if(length(unique(unlist(peaks_list))) == 0) {
      stop("No peaks remain with current peak detection paramters.
           Try setting baselineThresh = smaller number in argument
           (see speaq::detectSpecPeaks for more details).")
    }
    peakmat <- get_spectra_at_peaks(spectra(physeq)@specmat, peaks_list)
  }

  # extract ppms associated with peaks
  ppms <- colnames(spectra(physeq)@specmat)
  rownames(peakmat) <- rownames(spectra(physeq)@specmat)

  if(binary) {
    peakmat[peakmat > 0] <- 1
  }

  result_list <- list(specmat = specmat,
                      peakmat = peakmat,
                      adjmat = spectra(physeq)@adjmat)
  physeq@spectra <- spectra(result_list)
  return (physeq)
}

# interpolate-duplicates -------------------------------------------------------
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
  return (interpolation)
}