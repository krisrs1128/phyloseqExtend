# spectra -----------------------------------------------------------------
#' The S4 class for storing metabolomic information
#'
#'  @description A matrix containing either raw spectra measurements or
#'    processed peaks indicators or heights.
#'
#'  @name spectra-class
#'  @exportClass spectra
setClass("spectra", representation(peaks = "logical"),
         contains = "matrix")

# phyloseq-extend ---------------------------------------------------------
#' The high-level extended phyloseq object, including new slot types
#'
#'  @slot spectra An object of class spectra, used for storing metabolomic
#'    information
#'  @importClassesFrom phyloseq phyloseq
#'  @exportClass phyloseqExtend
setClass("phyloseqExtend",
         representation = representation(
           spectra = "spectra"
          ),
         contains = "phyloseq")
