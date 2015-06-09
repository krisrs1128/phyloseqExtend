# build-spectra-object ----------------------------------------------------
#' @title Build or access an object of class spectra
#'
#' @param spectra_matrix An object of class matrix, containing the spectral
#'  samples as rows.
#'
#' @return spectra An object of class spectra, containing the input spectra
#' matrix as the main data component.
#'
#' @docType methods
#' @export
setGeneric("spectra", function(object) {
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

#' @title Add default row and column names for a spectra object
add_dim_names <- function(object) {
  if(is.null(rownames(object@specmat)) & nrow(object@specmat) > 0) {
    rownames(object@specmat) <- paste0("sa", 1:nrow(object@specmat))
  }
  if(is.null(colnames(object@specmat)) & ncol(object@specmat) > 0) {
    colnames(object@specmat)  <- paste0("ind", 1:ncol(object@specmat))
  }
  if(is.null(rownames(object@peakmat)) & nrow(object@peakmat) > 0) {
    rownames(object@peakmat) <- paste0("sa", 1:nrow(object@peakmat))
  }
  if(is.null(colnames(object@peakmat)) & ncol(object@peakmat) > 0) {
    rownames(object@peakmat) <- paste0("pk", 1:ncol(object@peakmat))
  }
  if(is.null(rownames(object@adjmat)) & nrow(object@adjmat) > 0) {
   rownames(object@adjmat) <- paste0("sa", 1:nrow(object@adjmat))
  }
  if(is.null(rownames(object@adjmat)) & ncol(object@adjmat) > 0) {
   rownames(object@adjmat) <- paste0("sa", 1:nrow(object@adjmat))
  }
  return (object)
}

# If a single matrix is input, assume that matrix is the specmat
setMethod("spectra", "matrix", function(object) {
  spectra_obj <- new("spectra", specmat = object)
  spectra_obj <- add_dim_names(spectra_obj)
  return (spectra_obj)
})

# If a data.frame or data.table is input, convert to a matrix and
# assume it is a specmat
setMethod("spectra", "data.frame", function(object) {
  spectra(as(object, "matrix"))
})

# If a list is input, check the names to determine the appropriate slots
setMethod("spectra", "list", function(object) {
  spectra_obj <- new("spectra", specmat = as.matrix(object[["specmat"]]),
      peakmat = as.matrix(object[["peakmat"]]),
      adjmat = as.matrix(object[["adjmat"]]))
  spectra_obj <- add_dim_names(spectra_obj)
  return (spectra_obj)
})

# validity ----------------------------------------------------------------

setValidity("spectra",
            function(object) {
              nspec <- nrow(object@specmat)
              npeak <- nrow(object@peakmat)
              if(min(nspec, npeak) > 0) {
                if(nspec != npeak) {
                  "Number of spectra samples not equal to number of peak samples."
                }
                if(!identical(rownames(object@specmat), rownames(object@peakmat))) {
                  "Spectra row names are inconsistent"
                }
              } else {
                TRUE
              }
            })

# sample-names ------------------------------------------------------------

setMethod("sample_names", "spectra", function(physeq) {
  # Have to use argument physeq for consistency for phyloseq, but it
  # actually expects a spectra object as its argument.
  return (rownames(physeq@specmat))
})

# prune-samples -----------------------------------------------------------

#' @title Subset the samples of an object of class spectra
#'
#' @param samples_ix A vector of integers giving the rows of the spectra
#'  matrices to include.
#' @param spectra An object of class spectra, whose matrices we want to filter.
subset_samples_spectra <- function(samples_ix, spectra) {
  if(nrow(spectra@specmat) > 0) spectra@specmat <- spectra@specmat[samples_ix, ]
  if(nrow(spectra@peakmat) > 0) spectra@peakmat <- spectra@peakmat[samples_ix, ]
  if(nrow(spectra@adjmat) > 0) spectra@adjmat <- spectra@adjmat[samples_ix, ]
  return (spectra)
}

#' @title Define a prune_samples methods for spectra class
#'
setMethod("prune_samples", signature("character", "spectra"), function(samples, x) {
  samples_ix <- which(sample_names(x) %in% samples)
  x <- subset_samples_spectra(samples_ix, x)
  return (x)
})

# call-peaks --------------------------------------------------------------
#' @title Wrapper for speaq::detectSpecPeaks
#'
#' @param spectra_object A spectra object with a nonempty @specmat slot
#' @return peaks_list A list whose i^th element contains the indices
#' of peaks in the i^th sample
#' @importFrom speaq detectSpecPeaks
#' @export
get_peaks_list <- function(spectra_object, ...) {
  peaks_list <- speaq::detectSpecPeaks(spectra_object@specmat, ...)
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
