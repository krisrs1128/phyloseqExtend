# show-phyloseqExtend -----------------------------------------------------
#' @title Print a phyloseqExtend object
#'
#' @export
#' @docType methods
setMethod("show", "phyloseqExtend", function(object) {
   # Print usual phyloseq tables
   signature(object = "phyloseq")
   callNextMethod()

    # print spectra summary if there
   spec_obj <- spectra(object)
   if(!is.null(spec_obj)) {
    specmat_string <- paste0("\t\t\t specmat [ ", nrow(spec_obj@specmat), " samples by ", ncol(spec_obj@specmat), " indices]\n")
    peakmat_string <- paste0("\t\t\t peakmat [ ", nrow(spec_obj@peakmat), " samples by ", ncol (spec_obj@peakmat), " peaks]\n")
    adjmat_string <- paste0("\t\t\t adjmat  [ ", nrow(spec_obj@adjmat), " samples by ", ncol(spec_obj@adjmat), " samples]\n")
    cat(paste0("spectra()     Spectra: \n", specmat_string, peakmat_string, adjmat_string), fill = TRUE)
  }
})
