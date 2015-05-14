# build-spectra-object ----------------------------------------------------
#'  @title Build or access an object of class spectra
#'  @param spectra_matrix An object of class matrix, containing the spectral
#'    samples as rows.
#'  @param peaks A logical indicating whether the matrix columns corresponds to
#'    wavelengths (FALSE, for raw spectra) or just peaks (TRUE, for spectra
#'    that have been processed). Defaults to FALSE.
#'  @return spectra An object of class spectra, containing the input spectra
#'    matrix as the main data component.
#'  @docType methods
#'  @export

setGeneric("spectra", function(object) {
  standardGeneric("spectra")
})

# If a phyloseqExtend object with a spectra slot is input, return that slot
setMethod("spectra", "phyloseqExtend", function(object) {
  spectra_object <- slot(object, "spectra")
  if(is.null(spectra_object)) {
    warning("spectra slot is null.")
  }
  return (spectra_object)
})

# If a matrix is input, create a new spectra object from that matrix
setMethod("spectra", "matrix", function(object) {
  return (new("spectra", object))
})
