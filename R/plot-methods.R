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
#' @param col A string giving a column name in sample_data(physeq) according
#'  to which we will color the spectra.
#' @param linetype A string giving a column name in sample_data(physeq)
#'  according to which we will modify the spectra linetypes.
#' @param facet_cols A character vector giving the column names in
#'  sample_data(physeq) according to which to facet the plot by.
#' @param plot_title A title to include for the figure.
#' @param alpha The transparency parameter for the plots, when using ggplot2.
#' @param line_thickness The thickness of the spectra lines, when using ggplot2.
#'
#' @return If "ggplot2" is used, then return the ggplot2 object containing the
#'  plot. Otherwise does not return anything.
#'
#' @export
plot_spectra  <- function(physeq, method = "speaq", plot_type = "specmat",
                          log_scale = FALSE, subsample_frac = 1, x_min = NULL,
                          x_max = NULL, col = NULL, linetype = NULL,
                          facet_cols = NULL, plot_title = NULL, alpha = 1,
                          line_thickness = 1, ...) {
  method <- match.arg(method, choices = c("speaq", "ggplot2"))
  plot_type <- match.arg(plot_type, choices = c("specmat", "peakmat"))
  stopifnot(!is.null(spectra(physeq)))
  if(plot_type == "specmat") {
    p <- plot_spectra_peaks(physeq, log_scale, subsample_frac, x_min,
                            x_max, col, linetype, facet_cols, plot_title,
                            alpha, line_thickness, ...)
  } else if(plot_type == "peakmat") {
    p <- plot_raw_spectra(physeq, method, log_scale, subsample_frac,
                          x_min, x_max, col, linetype, facet_cols,
                          plot_title, alpha, line_thickness, ...)
  }
  if(!is.null(p)) {
    return (p)
  }
}

# plot-raw-spectra --------------------------------------------------------

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
#' @param col A string giving a column name in sample_data(physeq) according
#'  to which we will color the spectra.
#' @param linetype A string giving a column name in sample_data(physeq)
#'  according to which we will modify the spectra linetypes.
#' @param facet_cols A character vector giving the column names in
#'  sample_data(physeq) according to which to facet the plot by.
#' @param plot_title A title to include for the figure.
#' @param alpha The transparency parameter for the plots, when using ggplot2.
#' @param line_thickness The thickness of the spectra lines, when using ggplot2.
#'
#' @return p The ggplot object containing plots of the raw spectra, with
#'  color / linetype / faceting annotation as desired.
#'
#' @importFrom speaq drawSpec
#' @importFrom ggplot2 ggplot geom_line aes_string facet_grid scale_y_log10
#' @importFrom reshape2 melt
#' @importFrom data.table data.table
#' @importFrom dplyr filter
#' @importFrom phyloseq sample_data
plot_raw_spectra <- function(physeq, method = "speaq", log_scale = FALSE,
                             subsample_frac = 1, x_min = NULL, x_max = NULL,
                             col = NULL, linetype = NULL, facet_cols = NULL,
                             plot_title = NULL, alpha = 1, line_thickness = 1) {

  # Extract spectra
  spectra_matrix <- spectra(physeq)@specmat@.Data
  spectra_matrix <- subsample_spectra_cols(spectra_matrix, subsample_frac, x_min, x_max)

  # Extract features to annotate with
  annotation_names <- unique(c(col, linetype, facet_cols))
  annotation_list <- get_annotation(physeq, annotation_names)
  annotation <- annotation_list$annotation
  interaction_label <- annotation_list$interaction_label

  if(method == "speaq") {
  # Plot using speaq
    drawSpec(
      spectra_matrix,
      main = plot_title,
      useLog = ifelse(log_scale, 1, -1),
      groupLabel = interaction_label
    )
    return()
  } else if(method == "ggplot2") {

    spectra_dat <- data.table(id = rownames(spectra_matrix), spectra_matrix)
    if(!is.null(annotation)) {
      spectra_dat <- data.table(spectra_dat, annotation)
    }

    # Convert from "wide" to "long" format, for ggplot2
    spectra_dat <- reshape2::melt(spectra_dat,
                                  id.vars = c("id", annotation_names),
                                  variable.name = "index",
                                  value.name = "intensity")
    spectra_dat$index <- as.numeric(as.character(spectra_dat$index))

    # Construct the desired plot
    p <- ggplot(spectra_dat) +
      geom_line(aes_string(x = "index", y = "intensity", group = "id",
                           col = col, linetype = linetype),
                alpha = alpha, size = line_thickness) +
      ggtitle(plot_title)
    if(!is.null(facet_cols)) {
      if(length(facet_cols) == 1) facet_cols <- c(facet_cols, ".")
      facet_formula <- paste0(facet_cols, collapse = "~")
      p <- p + facet_grid(facet_formula)
    }
    if(log_scale) {
      p <- p + scale_y_log10()
    }
    return (p)
  }
}


# utils -------------------------------------------------------------------

#' @title Helper to get annotation list from a physeq object
#'
#' @param physeq A phyloseq object with a nonempty sample data slot
#' @param annotation_names A character vector giving the names of sample data
#'  columns to extract
#'
#' @return annotation_interaction A list containing the following elements
#'  $annotation: The columns of the sample_data data.frame requested from
#'    annotation_names
#'  $interaction_label: The factor vector containing all combinations of
#'    levels in the annotation variables.
#' @importFrom phyloseq sample_data
get_annotation <- function(physeq, annotation_names = NULL) {
  if(!is.null(annotation_names)) {
    sample_data_table <- data.table(data.frame(sample_data(physeq)))
    annotation <- sample_data_table[, annotation_names, with=F]
    interaction_label <- droplevels(interaction(annotation))
  } else {
    annotation <- NULL
    interaction_label <- NULL
  }
  return(list(annotation = annotation, interaction_label = interaction_label))
}

#' @title Helper to add layers to ggplot plot
add_plot_layers <- function(p, plot_title = NULL, facet_cols = NULL,
                            log_scale = FALSE) {
  # add faceting
  if(!is.null(facet_cols)) {
    facet_formula <- formula(paste0(facet_cols, collapse = "~"))
    p <- p + facet_grid(facet_formula)
  }

  # add log scale
  if(log_scale) {
    p <- p + scale_y_log10()
  }

  # add title
  p <- p + ggtitle(plot_title)
  return (p)
}

# plot-spectra-peaks ------------------------------------------------------

#' @title Plot only the Peaks in a spectra matrix
#'
#' @param physeq A phyloseq object containing the spectrum to plot
#' @param log_scale Should the intensities be plotted on a log scale?
#' @param x_min What is the minimum index to display?
#' @param x_max What is the maximum index to display?
#' @param subsample_frac We will only plot ever the value at every
#'  1 / subsample_frac indices. This can accelerate plotting in the case that
#'  the spectrum is very long, but can lead to missed peaks.
#' @param col A string giving a column name in sample_data(physeq) according
#'  to which we will color the spectra.
#' @param linetype A string giving a column name in sample_data(physeq)
#'  according to which we will modify the spectra linetypes.
#' @param facet_cols A character vector giving the column names in
#'  sample_data(physeq) according to which to facet the plot by.
#' @param plot_title A title to include for the figure.
#' @param alpha The transparency parameter for the plots, when using ggplot2.
#' @param line_thickness The thickness of the spectra lines, when using ggplot2.
#'
#' @return p The ggplot object containing plots of the raw spectra, with
#'  color / linetype / faceting annotation as desired.
plot_spectra_peaks <- function(physeq, log_scale = FALSE, subsample_frac = 1,
                               x_min = NULL, x_max = NULL, col = NULL,
                               linetype = NULL, facet_cols = NULL,
                               plot_title = NULL, alpha = 1, line_thickness = 1, ...) {

  # Extract spectra
  spectra_matrix <- spectra(physeq)@peakmat@.Data
  spectra_matrix <- subsample_spectra_cols(spectra_matrix, subsample_frac, x_min, x_max)

  # Extract features to annotate with
  annotation_names <- unique(c(col, linetype, facet_cols))
  annotation <- get_annotation(physeq, annotation_names)$annotation

  # append the annotation to the spectra matrix
  spectra_dat <- cbind(spectra_matrix, annotation)
  spectra_dat <- melt(spectra_dat, id.vars = annotation_names,
                      value.name = "intensity")
  if("variable" %in% colnames(spectra_dat)) {
    spectra_dat$index <- as.numeric(as.character(spectra_dat$variable))
    spectra_dat$variable <- NULL
  } else {
    spectra_dat$index <- as.numeric(as.character(spectra_dat$Var2))
    spectra_dat$Var1 <- NULL
    spectra_dat$Var2 <- NULL
  }

  # construct the plot
  p <- ggplot(spectra_dat) +
    geom_segment(aes_string(x = "index", y = "0", xend = "index", yend = "intensity",
                     col = col, linetype = linetype), alpha = alpha, size = line_thickness) +
    scale_y_continuous("intensity")
  p <- add_plot_layers(p, plot_title, facet_cols, log_scale)
  return (p)
}
