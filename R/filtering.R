

# Filter samples based on number of reads in peaks -----------------------------

#' filterSamples
#' 
#' function to get indices of samples that pass filtters
#' @param object SummarizedExperiment with matrix of fragment counts per peak
#'  per sample, as computed by \code{\link{getCounts}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param shiny make shiny gadget?
#' @param ix_return return indices of sample to keep instead of subsetted counts
#'  object
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated
#'  based on data. min_in_peaks is set to 0.5 times the median proportion of
#'  fragments in peaks. min_depth is set to the maximum of 500 or 10% of the
#'  median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{getCounts}},  \code{\link{getPeaks}}, 
#' \code{\link{filterPeaks}}
#' @export 
#' @examples         
#' data(example_counts, package = "chromVAR")
#' 
#' counts_filtered <- filterSamples(example_counts, min_depth = 1500,
#'                                   min_in_peaks = 0.15, shiny = FALSE)
filterSamples <- function(object, min_in_peaks = NULL, min_depth = NULL, 
                          shiny = interactive(), 
  ix_return = FALSE) {
  object <- counts_check(object)
  if ("depth" %ni% colnames(colData(object))) 
    stop("colData for object must have column named depth with total reads ",
         "per sample")
  if (shiny && !interactive()) 
    shiny <- FALSE
  depths <- colData(object)$depth
  fragments_per_sample <- getFragmentsPerSample(object)
  if (is.null(min_in_peaks)) {
    min_in_peaks <- round(median(fragments_per_sample/depths) * 0.5, digits = 3)
    if (!shiny) 
      message("min_in_peaks set to ", min_in_peaks)
  }
  if (is.null(min_depth)) {
    min_depth <- max(500, median(depths) * 0.1)
    if (!shiny) 
      message("min_depth set to ", min_depth)
  }
  if (shiny) {
    shiny_returned <- filter_samples_gadget(depths, fragments_per_sample, 
                                            min_in_peaks, 
                                            min_depth, colnames(object))
    message("min_in_peaks set to ", shiny_returned$min_in_peaks)
    message("min_depth set to ", shiny_returned$min_depth)
    keep_samples <- shiny_returned$keep
  } else {
    keep_samples <- 
      intersect(which(depths >= min_depth), 
                which(fragments_per_sample/depths >= min_in_peaks))
  }
  if (ix_return) 
    return(keep_samples) else return(object[, keep_samples])
}

#' filterSamplesPlot
#' 
#' plot filtering of samples
#' @param object SummarizedExperiment with matrix of fragment counts per peak 
#' per sample, as computed by \code{\link{getCounts}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param use_plotly make interactive plot?
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated
#' based on data. min_in_peaks is set to 0.5 times the median proportion of 
#' fragments in peaks.  min_depth is set to the maximum of 500 or 10% of the 
#' median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{getCounts}},  \code{\link{getPeaks}}, 
#' \code{\link{filterPeaks}}
#' @export
#' @examples 
#' data(example_counts, package = "chromVAR")
#' 
#' counts_filtered <- filterSamples(example_counts, min_depth = 1500,
#'                                   min_in_peaks = 0.15, shiny = FALSE)
#' counts_filtered_plot <- filterSamplesPlot(counts_filtered, 
#'                                           min_in_peaks = 0.15,
#'                                           min_depth = 1500, 
#'                                           use_plotly = FALSE)
#' 
filterSamplesPlot <- function(object, min_in_peaks = NULL, min_depth = NULL, 
                                use_plotly = interactive()) {
  object <- counts_check(object)
  if ("depth" %ni% colnames(colData(object))) 
    stop("colData for object must have column named depth with total reads ",
         "per sample")
  if (use_plotly && !interactive()) 
    use_plotly <- FALSE
  depths <- colData(object)$depth
  fragments_per_sample <- getFragmentsPerSample(object)
  if (is.null(min_in_peaks)) {
    min_in_peaks <- round(median(fragments_per_sample/depths) * 0.5, digits = 3)
    message("min_in_peaks set to ", min_in_peaks)
  }
  if (is.null(min_depth)) {
    min_depth <- max(500, median(depths) * 0.1)
    message("min_depth set to ", min_depth)
  }
  keep_samples <- 
    intersect(which(depths >= min_depth), 
              which(fragments_per_sample/depths >=min_in_peaks))
  tmp_df <- data.frame(x = depths,
                       y = fragments_per_sample/depths, 
                       z = ((seq_along(fragments_per_sample)) %in% 
    keep_samples), name = colnames(object))
  p <- ggplot(tmp_df, aes_string(x = "x", y = "y", col = "z", text = "name")) + 
    geom_point() + 
    xlab("Number of fragments") + 
    ylab("Proportion of fragments in peaks") + 
    scale_x_log10() +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, min(1, 
                                         max(tmp_df$y) * 1.2))) + 
    scale_color_manual(name = "Pass?", values = c("gray", 
                                                  "black"), 
                       breaks = c(TRUE, FALSE), labels = c("Yes", "No")) + 
    chromVAR_theme() +
    geom_hline(yintercept = min_in_peaks, col = "red", lty = 2) + 
    geom_vline(xintercept = min_depth, 
               col = "red", lty = 2)
  if (use_plotly) {
    return(ggplotly(p))
  } else {
    return(p + annotation_logticks(sides = "b"))
  }
}


filter_samples_gadget <- function(depths, 
                                  fragments_per_sample, 
                                  min_in_peaks, 
                                  min_depth, 
                                  sample_names) {
  
  ui <- miniPage(
    gadgetTitleBar("Adjust parameters to change filtering"), 
    fillCol(
      flex = c(1,3), 
      fillRow(flex = c(1, 1), 
              sliderInput("min_in_peaks", 
                          "Minimum fraction of fragments in peaks:", 
                          min = 0,
                          max = 1, 
                          value = min_in_peaks), 
              numericInput("min_depth", 
                           "Minimum total reads:", 
                           min = 0, 
                           max = max(depths), 
                           value = min_depth)), 
      plotlyOutput("plot", height = "100%")
      )
    )
  
  server <- function(input, output, session) {
    
    keep_samples <- 
      reactive(
        intersect(which(depths >= input$min_depth),
                  which(fragments_per_sample/depths >= input$min_in_peaks)))
    
    # Render the plot
    output$plot <- renderPlotly({
      tmp_df <- 
        data.frame(x = depths, 
                   y = fragments_per_sample/depths,
                   pass_filter = (seq_along(fragments_per_sample) %in% 
                                    keep_samples()), name = sample_names)
      p <- ggplot(tmp_df, aes_string(x = "x", y = "y", col = "pass_filter", 
                                     text = "name")) + 
        geom_point() + 
        xlab("Number of fragments") + 
        ylab("Proportion of fragments in peaks") + 
        scale_x_log10() +
        scale_y_continuous(expand = c(0, 0), 
                           limits = c(0, min(1, max(tmp_df$y) * 1.2))) +
        scale_color_manual(name = "Pass?", 
                           values = c("gray", "black"), breaks = c(TRUE, FALSE),
                           labels = c("Yes","No")) + chromVAR_theme()
      p <- p + 
        geom_hline(yintercept = input$min_in_peaks, col = "red", lty = 2) + 
        geom_vline(xintercept = input$min_depth, col = "red", lty = 2)
      ggplotly(p)
      p
    })
    
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      stopApp(list(keep = keep_samples(), min_in_peaks = input$min_in_peaks, 
        min_depth = input$min_depth))
    })
  }
  
  runGadget(ui, server)
}




# Filter samples based on biases ----------------------------------------------


bias_skew <- function(object, nbins = 10, expectation = NULL, 
                      what = c("bias", "counts")) {
  
  object <- counts_check(object)
  
  what <- match.arg(what)
  
  fragments_per_sample <- getFragmentsPerSample(object)
  
  if (what == "bias") {
    bias <- rowRanges(object)$bias
  } else {
    bias <- getFragmentsPerPeak(object)
  }
  
  if (min(getFragmentsPerPeak(object)) <= 0) 
    stop("All peaks must have at least one fragment in one sample")
  
  bias <- bias + runif(length(bias), 0, 0.001)
  bias_quantiles <- quantile(bias, seq(0, 1, 1/nbins))
  bias_cut <- cut(bias, breaks = bias_quantiles)
  
  bias_bins <- split(seq_len(nrow(object)), bias_cut)
  
  if (is.null(expectation)) {
    expectation <- computeExpectations(object)
  } else {
    stopifnot(length(expectation) == nrow(object))
  }
  
  
  bias_mat <- sparseMatrix(j = unlist(bias_bins), 
                           i = unlist(lapply(seq_len(nbins),
                                             function(x) 
                                               rep(x, length(bias_bins[[x]])))),
                           x = 1, 
                           dims = c(nbins, nrow(object)))
  
  observed <- as.matrix(bias_mat %*% assays(object)$counts)
  
  expected <- as.vector(bias_mat %*% expectation %*% fragments_per_sample)
  
  
  out <- (observed - expected)/expected
  return(out)
}

upper_bias_limit_helper <- function(x, k) {
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  return(q[2] + k * (q[2] - q[1]))
}


# Filter peaks based on counts -------------------------------------------------

#' filterPeaks
#' 
#' function to get indices of peaks that pass filters
#' @param object SummarizedExperiment with matrix of fragment counts per peak
#'  per sample, as computed by \code{\link{getCounts}}
#' @param min_fragments_per_peak minimum number of fragmints in peaks across all
#'  samples 
#' @param non_overlapping reduce peak set to non-overlapping peaks, see details
#' @param ix_return return indices of peaks to keep instead of subsetted counts
#'  object
#' @details if non_overlapping is set to true, when peaks overlap the 
#' overlapping peak with lower counts is removed
#' @return vector of indices, representing peaks that should be kept
#' @seealso \code{\link{getPeaks}},  \code{\link{filterSamples}},
#' \code{\link{getCounts}}
#' @export   
#' @author Alicia Schep       
#' @examples 
#' data(example_counts, package = "chromVAR")
#' 
#' counts_filtered <- filterSamples(example_counts, min_depth = 1500,
#'                                   min_in_peaks = 0.15, shiny = FALSE)
#' counts_filtered <- filterPeaks(example_counts)
filterPeaks <- function(object, 
                         min_fragments_per_peak = 1, 
                         non_overlapping = TRUE, 
                         ix_return = FALSE) {
  object <- counts_check(object)
  fragments_per_peak <- getFragmentsPerPeak(object)
  peaks <- rowRanges(object)
  keep_peaks <- which(fragments_per_peak >= min_fragments_per_peak)
  if (non_overlapping) {
    strand(peaks) <- "*"
    if (!isDisjoint(peaks) && !isSorted(peaks)) {
      stop("Peaks must be sorted to be able to filter overlapping peaks!\n",
           "Please use 'sort' function to filter peaks")
    }
    while (!(isDisjoint(peaks[keep_peaks]))) {
      chr_names <- as.character(seqnames(peaks[keep_peaks]))
      starts <- start(peaks[keep_peaks])
      ends <- end(peaks[keep_peaks])
      overlap_next <- 
        intersect(which(chr_names[seq_len(length(keep_peaks) - 1)] == 
                          chr_names[seq_len(length(keep_peaks) - 1) + 1]), 
                  which(ends[seq_len(length(keep_peaks) - 1)] >= 
                          starts[seq_len(length(keep_peaks) - 1) + 1]))
      overlap_previous <- overlap_next + 1
      overlap_comparison <- fragments_per_peak[keep_peaks[overlap_previous]] > 
        fragments_per_peak[keep_peaks[overlap_next]]
      discard <- keep_peaks[c(overlap_previous[!overlap_comparison], 
                              overlap_next[overlap_comparison])]
      keep_peaks <- keep_peaks[keep_peaks %ni% discard]
    }
  }
  if (ix_return) 
    return(keep_peaks) else return(object[keep_peaks, ])
}
