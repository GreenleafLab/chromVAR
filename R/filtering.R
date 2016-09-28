

# Filter samples based on number of reads in peaks -----------------------------

#' filter_samples
#' 
#' function to get indices of samples that pass filtters
#' @param object SummarizedExperiment with matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param shiny make shiny gadget?
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated based on data.
#' min_in_peaks is set to 0.5 times the median proportion of fragments in peaks.  min_depth is 
#' set to the maximum of 500 or 10% of the median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_peaks}}
#' @import miniUI ggplot2 plotly shiny
#' @export          
filter_samples <- function(object,  min_in_peaks = NULL, min_depth = NULL, shiny = TRUE, ix_return = FALSE){
  depths <- colData(object)$depth
  fragments_per_sample = get_fragments_per_sample(object)
  if (is.null(min_in_peaks)){
    min_in_peaks = round(median(fragments_per_sample/depths)*0.5, digits = 3)
    if (!shiny) message(paste("min_in_peaks set to ",min_in_peaks,sep="",collapse=""))
  } 
  if (is.null(min_depth)){
    min_depth = max(500,median(depths)*0.1)
    if (!shiny) message(paste("min_depth set to ",min_depth,sep="",collapse=""))
  }
  if (shiny){
    shiny_returned <- filter_samples_gadget(depths, fragments_per_sample, min_in_peaks, min_depth, colnames(object))
    message(paste("min_in_peaks set to ", shiny_returned$min_in_peaks,sep="",collapse=""))
    message(paste("min_depth set to ", shiny_returned$min_depth,sep="",collapse=""))
    keep_samples <- shiny_returned$keep
  } else{
    keep_samples <- intersect(which(depths >= min_depth), 
                            which(fragments_per_sample/depths >= min_in_peaks))   
  }
  if (ix_return) return(keep_samples) else return(object[,keep_samples])    
}

#' filter_samples_plot
#' 
#' plot filtering of samples
#' @param object SummarizedExperiment with matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' @param min_in_peaks minimum fraction of samples within peaks 
#' @param min_depth minimum library size
#' @param interactive make interactive plot?
#' @details If unspecified, min_in_peaks and min_depth cutoffs will be estimated based on data.
#' min_in_peaks is set to 0.5 times the median proportion of fragments in peaks.  min_depth is 
#' set to the maximum of 500 or 10% of the median library size.
#' @return indices of samples to keep
#' @seealso \code{\link{get_counts}},  \code{\link{get_inputs}}, \code{\link{filter_peaks}}
#' @import miniUI ggplot2 plotly shiny
#' @export          
filter_samples_plot <- function(object,  min_in_peaks = NULL, min_depth = NULL, interactive = TRUE){
  depths <- colData(object)$depth
  fragments_per_sample = get_fragments_per_sample(object)
  if (is.null(min_in_peaks)){
    min_in_peaks = round(median(fragments_per_sample/depths)*0.5, digits = 3)
    message(paste("min_in_peaks set to ",min_in_peaks,sep="",collapse=""))
  } 
  if (is.null(min_depth)){
    min_depth = max(500,median(depths)*0.1)
    message(paste("min_depth set to ",min_depth,sep="",collapse=""))
  }
  keep_samples <- intersect(which(depths >= min_depth), 
                            which(fragments_per_sample/depths >= min_in_peaks))   
  tmp_df = data.frame(x = depths, y= fragments_per_sample/depths, z = ((1:length(fragments_per_sample)) %in% keep_samples),
                      name = colnames(object))
  p = ggplot(tmp_df, aes_string(x="x", y="y",col="z", text = "name")) + geom_point() +
    xlab("Number of fragments") + ylab("Proportion of fragments in peaks") + scale_x_log10() + 
    scale_y_continuous(expand =c(0,0), limits =c(0, min(1,max(tmp_df$y)*1.2)))+
    scale_color_manual(name = "Pass?",values = c("gray","black"),breaks = c(TRUE,FALSE), labels = c("Yes","No")) +
    chromVAR_theme()
  p = p + geom_hline(yintercept = min_in_peaks, col="red", lty = 2) + geom_vline(xintercept = min_depth, col="red", lty=2)
  if (interactive){
    return(ggplotly(p))
  } else{
    return(p+ annotation_logticks(sides="b"))
  }
}


filter_samples_gadget <- function(depths, fragments_per_sample,  min_in_peaks, min_depth, sample_names) {
  
  ui <- miniPage(
    gadgetTitleBar("Adjust parameters to change filtering"),
    fillCol(flex = c(1,3),
           fillRow(flex = c(1,1), sliderInput("min_in_peaks", "Minimum fraction of fragments in peaks:", 
                                min=0, max=1, value=min_in_peaks),
             numericInput("min_depth", "Minimum total reads:", 
                                 min=0, max=max(depths), value=min_depth)),
            plotlyOutput("plot", height = "100%")
    )
  )
  
  server <- function(input, output, session) {
    
    keep_samples <- reactive(intersect(which(depths >= input$min_depth), 
                                       which(fragments_per_sample/depths >= input$min_in_peaks)))
    
    # Render the plot
    output$plot <- renderPlotly({
      tmp_df = data.frame(x = depths, y= fragments_per_sample/depths, pass_filter = ((1:length(fragments_per_sample)) %in% keep_samples()),
                          name = sample_names)
      p = ggplot(tmp_df, aes_string(x="x", y="y",col="pass_filter", text = "name")) + geom_point() +
        xlab("Number of fragments") + ylab("Proportion of fragments in peaks") + scale_x_log10() + 
        scale_y_continuous(expand =c(0,0), limits =c(0, min(1,max(tmp_df$y)*1.2)))+
        scale_color_manual(name = "Pass?",values = c("gray","black"),breaks = c(TRUE,FALSE), labels = c("Yes","No")) +
        chromVAR_theme()
      p = p + geom_hline(yintercept = input$min_in_peaks, col="red", lty = 2) + geom_vline(xintercept = input$min_depth, col="red", lty=2)
      ggplotly(p)
    })
    
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      stopApp(list(keep = keep_samples(), min_in_peaks = input$min_in_peaks, min_depth = input$min_depth))
    })
  }
  
  runGadget(ui, server)
}




# Filter samples based on biases  ----------------------------------------------


bias_skew <- function(counts_mat,
                      bias,
                      nbins = 10,
                      expectation = NULL){
  
  if (inherits(counts_mat,"matrix")){
    counts_mat <- Matrix(counts_mat)
  }
  stopifnot(inherits(counts_mat,"Matrix"))
  
  if (inherits(bias, "GenomicRanges")){
    bias <- mcols(bias)$gc
  }
  
  keep  <- which(rowSums(counts_mat) > 0)
  counts_mat <- counts_mat[keep,]
  bias <- bias[keep]
  
  counts_info <- counts_summary(counts_mat)
  
  bias = bias + runif(length(bias),0,0.001)
  bias_quantiles = quantile(bias, seq(0,1,1/nbins))
  bias_cut = cut(bias, breaks = bias_quantiles)

  bias_bins = split(1:counts_info$npeak, bias_cut)
  
  if (is.null(expectation)){
    expectation <- compute_expectations(counts_mat)
  } else{
    stopifnot(length(expectation) == nrow(counts_mat))
  }
  
  sample_names <- colnames(counts_mat)
  
  
  bias_mat <- sparseMatrix(j = unlist(bias_bins), 
                           i = unlist(lapply(1:nbins,function(x) rep(x, length(bias_bins[[x]])))),
                           x = 1,
                           dims = c(nbins, counts_info$npeak))
  
  observed <- as.matrix(bias_mat %*% counts_mat)
  
  expected <- as.vector(bias_mat %*% expectation %*% counts_info$fragments_per_sample)
  
  
  out <- (observed - expected) / expected
  return(out)
}

upper_bias_limit_helper <- function(x, k){
  q <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  return(q[2] + k * (q[2] - q[1]))
}


# Filter peaks based on counts -------------------------------------------------

#' filter_peaks
#' 
#' function to get indices of peaks that pass filters
#' @param counts_mat SummarizedExperiment with matrix of fragment counts per peak per sample, as computed 
#' by \code{\link{getFragmentCounts}}
#' \code{\link{read_peaks}}
#' @param min_fragments_per_peak minimum number of fragmints in peaks across all samples 
#' @param non_overlapping reduce peak set to non-overlapping peaks, see details
#' @details if non_overlapping is set to true, when peaks overlap the overlapping peak with lower counts is removed
#' @return vector of indices, representing peaks that should be kept
#' @seealso \code{\link{get_peaks}},  \code{\link{get_inputs}}, \code{\link{filter_samples}},
#' \code{\link{get_counts}}
#' @export          
filter_peaks <- function(counts_mat, min_fragments_per_peak = 1, non_overlapping = TRUE, ix_return = FALSE){
  fragments_per_peak = get_fragments_per_peak(counts_mat)
  peaks <- rowRanges(counts_mat)
  keep_peaks <- which(fragments_per_peak >= min_fragments_per_peak)
  if (non_overlapping){
    strand(peaks) <- "*"    
    if (!isTRUE(all.equal(peaks, sort(peaks)))){
      stop("peaks must be sorted to be able to filter non-overlapping peaks!")
    }
    while (!(isDisjoint(peaks[keep_peaks]))){
      chr_names = seqnames(peaks[keep_peaks])
      starts = start(peaks[keep_peaks])
      ends = end(peaks[keep_peaks])      
      overlap_next = intersect(which(chr_names[1:(length(keep_peaks) -1)] == chr_names[2:(length(keep_peaks))]),
                               which(ends[1:(length(keep_peaks) -1)] >= starts[2:(length(keep_peaks))]))
      overlap_previous = overlap_next + 1
      overlap_comparison = fragments_per_peak[keep_peaks[overlap_previous]] > fragments_per_peak[keep_peaks[overlap_next]]
      discard = keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
      keep_peaks = keep_peaks[keep_peaks %ni% discard]
    }   
  }
  if (ix_return) return(keep_peaks) else return(counts_mat[keep_peaks,])
}
