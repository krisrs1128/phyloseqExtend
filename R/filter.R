
################################################################################
# More general / faster alternatives to phyloseq's filter_samples and
# filter_taxa.
################################################################################

# utils ------------------------------------------------------------------------
#' @title Fast row variance
#' @description Faster than apply. See http://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
row_var <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# vectorized-filters -----------------------------------------------------------
#' @title Vectorized filter function in the spirit of k over a
#' @description We want to keep only those rows / cols that have many values
#' above some threshold.
#' @param X A matrix whose rows want to filter.
#' @param a The threshold expression value.
#' @param q A quantile for the number of samples above the threshold,
#' determining which samples we keep. So, only samples with many samples above
#' a will be returned.
#' @return X A version of X filtered according to a and q.
#' @export
filter_nonzero <- function(X, a = 20, q = 0.9) {
  stats <- rowSums(X > a)
  x_q <- quantile(stats, q)
  X[stats > x_q, ]
}

#' @title Vectorized variance filtering
#' @description We want to keep only those rows / cols that have many values
#' above some threshold.
#' @param X A matrix whose rows want to filter.
#' @param q A quantile for the minimum variance of the rows we keep.
#' @return X A row-variance filtered version of X.
#' @export
filter_variance <- function(X, q = 0.9) {
  stats <- row_var(X)
  x_q <- quantile(stats, q)
  X[stats > x_q, ]
}
