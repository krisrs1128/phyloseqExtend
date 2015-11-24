
################################################################################
# Preprocess different components before they are input to the spectra object.
################################################################################

#' @title Merge default contingency table processing options
#' @param opts A partially filled list specifying options to contingency table
#' processing. Those that are not specified will be merged according to their
#' defaults below. These are,
#'  $id_col: Is the first column of the data a nonnumeric identifier column?
#'   Defaults to TRUE.
#'  $filter_ix: A function that can be applied to the full matrix (excluding a
#'   potential id column), and which returns the indices of the row filtered
#'   version. Defaults to filter_nonzero().
#'  $transform: A function that can be applied to the full matrix, transforming
#'   its values. Defaults to log().
#' @return A modified version of opts, with defaults filled in.
#' @export
merge_ctable_opts <- function(opts = list()) {
  default_opts <- list()
  default_opts$id_col <- TRUE
  default_opts$filter_ix <- filter_nonzero
  default_opts$transform <- function(x) { x }
  modifyList(default_opts, opts)
}

#' @title Wrapper for contingency table processing
#' @param X A matrix / data.frame / data.table with the samples as rows, and
#' with ids as either row names or as the first column [opts$id_col in
#' merge_ctable_opts()].
#' @param opts A list of processing options. Choices and their defaults are
#' detailed in merge_ctable_opts().
#' @return X A processed version of the input data X.
#' @importFrom data.table data.table
#' @export
process_ctable <- function(X, opts = list()) {
  init_class <- class(X)
  X <- data.table(X, keep.rownames = TRUE)

  opts <- merge_ctable_opts(opts)

  # extract just the numeric part of the data
  if(opts$id_col) {
    ids <- X[, 1, with = F]
    X <- X[, -1, with = F]
  }

  # perform filtering / transformations
  ix <- opts$filter_fun(X)
  X <- opts$transform(X[ix, ])

  # put back into original form
  if(opts$id_col) {
    X <- cbind(ids[ix], X)
  }
  class(X) <- init_class
  X
}
