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
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#'
#' @return If "ggplot2" is used, then return the ggplot2 object containing the
#'  plot. Otherwise returns null.
#'
#' @export
plot_spectra  <- function(physeq, method = "speaq", log_scale = FALSE,
                          subsample_frac = 1, x_min = NULL, x_max = NULL,
                          col = NULL, linetype = NULL, plot_title = NULL, ...) {
  method <- match.arg(method, choices = c("speaq", "ggplot2"))
  stopifnot(!is.null(spectra(physeq)))
  if(spectra_object@peaks) {
    stop("peak plotting is not yet supported \n
         set @peaks element in spectra slot to FALSE.")
    p <- plot_spectra_peaks(physeq, log_scale, x_min, x_max, plot_title, ...)
  } else {
    p <- plot_raw_spectra(physeq, method, log_scale, subsample_frac,
                          x_min, x_max, col, linetype, plot_title, ...)
  }
  if(!is.null(p)) {
    return (p)
  }
}

#' @title Plot Raw Spectra Intensities Across All Indices
#'
#' @description Before preprocessing to a collection of peaks, it may be
#' useful to plot the full raw spectra. This function provides an interface
#' to \code{speaq} and \code{ggplot2} that can be used to plot these raw
#' spectra.
#'
#' @param physeq A phyloseq object containing the spectrum to plot
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
#' @importFrom phyloseq sam_data
plot_raw_spectra <- function(physeq, method = "speaq", log_scale = FALSE,
                             subsample_frac = 1, x_min = NULL, x_max = NULL,
                             col = NULL, linetype = NULL, facet_cols = NULL,
                             plot_title = NULL, ...) {

  # Extract spectra
  spectra_matrix <- spectra(physeq)@.Data
  spectra_matrix <- subsample_spectra_cols(spectra_matrix, subsample_frac, x_min, x_max)

  # Extract features to annotate with
  annotation_names <- unique(c(col, linetype, facet_cols))
  if(!is.null(annotation_names)) {
    sample_data_table <- data.table(data.frame(sample_data(physeq)))
    annotation <- sample_data_table[, annotation_names, with=F]
    interaction_label <- droplevels(interaction(annotation))
  } else {
    annotation <- NULL
    interaction_label <- NULL
  }

  if(method == "speaq") {
    drawSpec(
      spectra_matrix,
      main = plot_title,
      useLog = ifelse(log_scale, 1, -1),
      groupLabel = interaction_label
    )
    return()
  } else if(method == "ggplot2") {
    # Convert from "wide" to "long" format, for ggplot2
    spectra_dat <- data.table(id = rownames(spectra_matrix), spectra_matrix)
    if(!is.null(annotation)) {
      spectra_dat <- data.table(spectra_dat, annotation)
    }

    spectra_dat <- reshape2::melt(spectra_dat, id.vars = c("id", annotation_names),
                                  variable.name = "index",
                                  value.name = "intensity")
    spectra_dat$index <- as.numeric(as.character(spectra_dat$index))

    # Construct the desired plot
    p <- ggplot(spectra_dat) +
      geom_line(aes_string(x = "index", y = "intensity", group = "id",
                           col = col, linetype = linetype)) +
      ggtitle(plot_title)
    if(!is.null(facet_cols)) {
      facet_formula <- formula(paste0(facet_cols, collapse = "~"))
      p <- p + facet_grid(facet_formula)
    }
    if(log_scale) {
      p <- p + scale_y_log10()
    }
    return (p)
  }
}
