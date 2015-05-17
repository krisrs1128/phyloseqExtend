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

# If a matrix is input, create a new spectra object from that matrix
setMethod("spectra", "matrix", function(object) {
  new("spectra", object)
})

# If a data.frame or data.table is input, convert to a matrix
setMethod("spectra", "data.frame", function(object) {
  spectra(as(object, "matrix"))
})

# call-peaks --------------------------------------------------------------
#' @title Extract aligned peaks from raw spectra
#'
#' @param spectra A spectra object containing raw spectra reads
#' @return spectra The original spectra object, updated so that the data
#'  object includes a binary indicator for peaks
#' @importFrom speaq detectSpecPeaks
extract_peaks <- function(spectra, ...) {
  if(spectra@peaks) {
    warning("spectra peaks are already extracted")
  } else {
    n <- nrow(spectra_object@.Data)
    p <- ncol(spectra_object@.Data)
    peaks_list <- detectSpecPeaks(spectra_object@.Data, ...)
    peaks_matrix <- matrix(0, n, p)
    for(i in 1:length(peaks_list)) {
      peaks_matrix[i, peaks_list[[i]]] <- 1
    }
    spectra@peaks <- TRUE
  }
  spectra@.Data <- peaks_matrix
  return (spectra)
}
