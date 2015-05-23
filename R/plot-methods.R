# description -------------------------------------------------------------
# phyloseq comes with a suite of nice, automatic graphics for various
# slots. Here, we consider plots based on the data occupying the new slots
# in phyloseqExtend.

# plot_spectra ------------------------------------------------------------
#' @title Plot Spectra Intensities
#'
#' @description If spectra(physeq)[["peaks"]] is FALSE, then this will plot the
#' raw intensities of each spectrum in the sample. Otherwise, we only plot
#' the positions and heights of peaks for each sample.
#'
#' @param physeq A phyloseq object containing the spectrum to plot
#' @param method The R package to use for the plot. Currently only supports
#'  "speaq" or "ggplot2".
#' @param log_scale Should the intensities be plotted on a log scale?
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#'
#' @return If "ggplot2" is used, then return the ggplot2 object containing the
#'  plot. Otherwise returns null.
#'
#' @export
plot_spectra  <- function(physeq, method = "speaq", log_scale = FALSE,
                          x_min = NULL, x_max = NULL, plot_title = NULL, ...) {
  method <- match.arg(method, choices = c("speaq", "ggplot2"))
  spectra_mat <- spectra(physeq)
  stopifnot(!is.null(spectra_mat))
  if(spectra_mat@peaks) {
    p <- plot_spectra_peaks(spectra_mat, log_scale, x_min, x_max,
                            plot_title, ...)
  } else {
    p <- plot_raw_spectra(spectra_mat, method, log_scale, x_min,
                          x_max, plot_title, ...)
  }
  return (p)
}

#' @title Plot Raw Spectra Intensities Across All Indices
#'
#' @description Before preprocessing to a collection of peaks, it may be
#' useful to plot the full raw spectra. This function provides an interface
#' to \code{speaq} and \code{ggplot2} that can be used to plot these raw
#' spectra.
#'
#' @param spectra_mat An object of class matrix, containing the spectral
#'  samples as rows.
#' @param method The R package to use for the plot. Currently only supports
#'  "speaq" or "ggplot2".
#' @param log_scale Should the intensities be plotted on a log scale?
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#'  the spectrum is very long, but can lead to missed peaks.
#' @param plot_title A title to include for the figure.
#'
#' @importFrom speaq drawSpec
#' @importFrom ggplot2 ggplot geom_line aes_string
#' @importFrom reshape2 melt
#' @importFrom data.table data.table
#' @importFrom dplyr filter
plot_raw_spectra <- function(spectra_mat, method = "speaq", log_scale = FALSE,
                             subsample_frac = 1, x_min = NULL, x_max = NULL,
                             plot_title = NULL, ...) {
  spectra_mat <- subsample_spectra_cols(spectra_mat, subsample_frac, x_min, x_max)
  if(method == "speaq") {
    speaq::drawSpec(spectra_mat, main = plot_title)
  } else if(method == "ggplot2") {
    # Convert from "wide" to "long" format, for ggplot2
    spectra_dat <- data.table::data.table(spectra_mat)
    spectra_dat <- cbind(id = rownames(spectra_mat), spectra_dat)
    spectra_dat <- reshape2::melt(spectra_dat, id.vars = "id", variable.name = "index", value.name = "intensity")
    spectra_dat$index <- as.numeric(as.character(spectra_dat$index))

    # Construct the desired plot
    p <- ggplot(spectra_dat) +
      geom_line(aes_string(x = "index", y = "intensity", group = "id")) +
      ggtitle(plot_title)
    if(log_scale) {
      p <- p + scale_y_log10()
    }
    return (p)
  }
}
