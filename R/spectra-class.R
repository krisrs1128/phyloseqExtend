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
#' @param spectra_object A spectra object containing raw spectra reads
#' @return peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @importFrom speaq detectSpecPeaks
#' @export
get_peaks_list <- function(spectra_object, ...) {
  peaks_list <- speaq::detectSpecPeaks(spectra_object@.Data, ...)
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


# wrap-peak-alignment -----------------------------------------------------

#' @title Wrapper for dohCluster in speaq
#'
#' @param spectra_object A spectra object containing raw spectra reads
#' @param peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#'
#' @return aligned_peaks A matrix whose ij^th element is the value of the j^th
#' unique peak in the i^th unique sample.
align_peaks <- function(spectra_object, peaks_list, ...) {
  ref_peak <- speaq::findRef(peaks_list)
  aligned_peaks <- speaq::dohCluster(spectra_object@.Data,
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
  spectra_matrix <- spectra(physeq)@.Data
  stopifnot(!is.null(spectra_matrix))
  max_vals <- apply(spectra_matrix, 1, max)

  keep_samples <- rownames(spectra_matrix)[which(max_vals <= thresh)]

  if(length(keep_samples) == 0) {
    stop(sprintf("Smallest spectrum max is %s. This threshold choice would remove all samples.", min(max_vals)))
  }
  physeq@spectra@.Data <- spectra_matrix[keep_samples,, drop = F]
  physeq@sam_data <- sample_data(physeq)[keep_samples,, drop = F ]
  physeq@otu_table@.Data <- physeq@otu_table@.Data[, keep_samples, drop = F]
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
#'
#' @return physeq_and_peaks A list containing the following elements
#'    $physeq: The input physeq object, with the spectra slot modified
#'      so that all peaks are aligned.
#'    $peaks: A list of lists whose i^th element is a list of peaks for the
#'      i^th sample.
#'
#' @export
detect_and_align_peaks <- function(physeq, detectSpecPeaksOpts, dohClusterOpts) {

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
  physeq@spectra <- spectra(do.call(align_peaks, dohClusterInput), peaks = FALSE)
  return (list(physeq = physeq, peaks = peaks))
}

# convert-to-peaks --------------------------------------------------------

#' @title Shrink a Full Spectrum to Just Peaks
#'
#' @description Rather than working with full spectrum reads, it is often
#'  convenient to work with a set of noteworthy peaks. In particular, we may
#'  want to first identify peaks across samples according to some peak calling
#'  algorithm, and then simplify the problem by setting
#'
#' @param physeq A phyloseqExtend object with a nonempty spectra slot.
#' @export
convert_to_peaks <- function(physeq, binary = FALSE, peaks_list = NULL, ...) {
  if(spectra(physeq)@peaks) {
    warning("spectra slot already in peaks form, doing nothing.")
    return(physeq)
  }

  # If we don't have any peaks yet, extract them
  if(is.null(peaks_list)) {
    peaks_list <- get_peaks_list(spectra(physeq), ...)
    if(length(unique(unlist(peaks_list))) == 0) {
      stop("No peaks remain with current peak detection paramters.
            Try setting baselineThresh = smaller number in argument
            (see speaq::detectSpecPeaks for more details).")
    }
  }

  # extract ppms associated with peaks
  ppms <- colnames(spectra(physeq))
  unique_peaks <- unique(unlist(peaks_list))
  unique_peaks <- ppms[unique_peaks]

  # Create a binary matrix where 1 indicates a peak
  peaks_matrix <- ldply(peaks_list, function(x) {
    vec  <- setNames(rep(0, length(unique_peaks)), unique_peaks)
    vec[ppms[x]] <- 1
    return (vec)
  })
  rownames(peaks_matrix) <- rownames(spectra(physeq))

  # If the user only wants binary peaks, return at this step
  if(binary) {
    physeq@spectra <- spectra(peaks_matrix, peaks = TRUE)
    return (physeq)
  }

  # If we want the actual peaks heights, search for the actual peak heights
  # in the input matrix, and substitute these
  ones_ix <- which(peaks_matrix == 1, arr.ind = T)
  named_ones_ix <- cbind(rownames(peaks_matrix)[ones_ix[, 1]],
                         colnames(peaks_matrix)[ones_ix[, 2]])

  peaks_matrix[ones_ix] <- spectra(physeq)@.Data[named_ones_ix]
  physeq@spectra <- spectra(peaks_matrix, peaks = TRUE)
  return (physeq)
}
