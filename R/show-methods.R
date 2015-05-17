# show-phyloseqExtend -----------------------------------------------------
#' @title Print a phyloseqExtend object
#'
#' @export
#' @docType methods
setMethod("show", "phyloseqExtend", function(object){
  # Print usual phyloseq tables
  signature(object = "phyloseq")
  callNextMethod()

  # print spectra summary if there
  spec_obj <- spectra(object)
  if(!is.null(spec_obj)) {
    cat(paste("spectra()     ", class(spec_obj)[1], ":           [ ", nrow(spec_obj), " samples by ", ncol(spec_obj), " intensities ]", sep = ""), fill = TRUE)
  }

})
