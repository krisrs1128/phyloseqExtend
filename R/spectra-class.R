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
