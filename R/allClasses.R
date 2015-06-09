# spectra -----------------------------------------------------------------
#' @title The S4 class for storing metabolomic information
#' @description A class containing 1) a matrix with the n spectra samples
#'  as rows 2) a matrix with n samples as rows and with columns giving the
#'  union of all detected peaks, and 3) the adjacency matrix for samples
#'  using the coocurrence of peaks.
#'
#' @name spectra-class
#' @exportClass spectra
setClass("spectra", slots = list(specmat = "matrix", peakmat = "matrix",
                                 adjmat = "matrix"))

# null-unions -------------------------------------------------------------
#' @title spectra S4 class with placeholder
#' @description If a component of the phyloseq extend object is not present, we keep a NULL
#'   placeholder.
#' @keywords internal
#' @name spectraOrNull-class
#' @exportClass spectraOrNull
setClassUnion("spectraOrNull", c("spectra", "NULL"))

# phyloseq-extend ---------------------------------------------------------
#' @title The high-level extended phyloseq object, including new slot types
#' @slot spectra An object of class spectra, used for storing metabolomic
#'   information
#' @name phyloseqExtend-class
#' @importClassesFrom phyloseq phyloseq
#' @exportClass phyloseqExtend
setClass("phyloseqExtend",
         representation = representation(
           spectra = "spectraOrNull"
          ),
         contains = "phyloseq")

