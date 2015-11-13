#' An S4 class that holds fragment counts across peak regions
#'
#' @slot counts Matrix of fragment counts, with each row representing a peak, each column a sample
#' @slot peaks GRanges object with peaks for each row of counts
#' @slot total_fragments total fragments within peaks
#' @slot fragments_per_sample total fragments within peaks per sample
#' @slot fragments_per_peak total fragments within each peak across samples
#' @slot nsample number of samples
#' @slot npeak number of peaks
#' @slot sample_meta optional data.frame with meta data for samples, with rows corresponding to samples
#' @slot depth total read depth for each sample
#' @importClassesFrom Matrix Matrix
#' @importClassesFrom GenomicRanges GenomicRanges
#'
#' @export
fragmentCounts <- setClass("fragmentCounts",
                           slots = c(counts = 'Matrix',
                                     peaks = 'GenomicRanges',
                                     total_fragments = 'numeric',
                                     fragments_per_sample = 'numeric',
                                     fragments_per_peak = 'numeric',
                                     nsample = 'numeric',
                                     npeak = 'numeric',
                                     sample_meta = 'data.frame',
                                     depth = 'numeric'
                           ))


#' Fragment Counts Constructor
#'
#' makes fragmentCounts object using counts Matrix and peaks.  
#' generally not for use by end-user.
#' @param counts A Matrix with fragment counts-- rows should correspond to peaks, columns to samples
#' @param peaks GRanges with peaks
#' @rdname fragmentCounts
#' @export
fragmentCounts <- function(counts, peaks, depth){
  total_fragments = sum(counts)
  fragments_per_sample = colSums(counts)
  fragments_per_peak = rowSums(counts)
  nsample = ncol(counts)
  npeak = nrow(counts)
  return(new("fragmentCounts", counts = counts, peaks = peaks,
             fragments_per_sample = fragments_per_sample, 
             fragments_per_peak = fragments_per_peak,
             nsample = nsample, npeak = npeak,
             depth = depth))
}


setMethod("show",
          signature="fragmentCounts",
          definition = function(object){
            cat("An object of class ", class(object), "\n", sep = "")
            cat("Counts for ", ncol(object@counts), " cells/samples and ",
                nrow(object@counts), " peaks. \n", sep = "")
            invisible(NULL)
          })

setValidity("fragmentCounts", function(object) {
  msg <- NULL
  valid <- TRUE
  #Check dimensions
  if (nrow(object@counts) != length(object@peaks)){
    valid <- FALSE
    msg <- c(msg, "Number of rows of counts matrix must be same as length of peaks")
  }
  if (valid) TRUE else msg})


setMethod("[", signature = signature(x = "fragmentCounts", j = "missing"),
          definition = function(x, i, j) {
            x@counts = x@counts[i, ]
            x@fragments_per_sample = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@peaks = x@peaks[i]
            x@total_fragments = sum(x@counts)
            x@npeak = lenth(x@peaks)
            return(x)
          })

setMethod("[", signature = signature(x = "fragmentCounts"),
          definition = function(x, i ,j ) {
            x@counts = x@counts[i,j]
            x@fragments_per_sample = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@peaks = x@peaks[i]
            x@total_fragments = sum(x@counts)
            x@npeak = lenfth(x@peaks)
            x@nsample = ncol(x@counts)
            x@sample_meta = x@sample_meta[j,]
            x@depth = x@depth[j]
            return(x)
          })

setMethod("[", signature = signature(x = "fragmentCounts", i = "missing"),
          definition = function(x, i ,j ) {
            x@counts = x@counts[,j]
            x@fragments_per_sample = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@total_fragments = sum(x@counts)
            x@nsample = ncol(x@counts)
            x@sample_meta = x@sample_meta[j,]
            x@depth = x@depth[j]
            return(x)
})


#' make fragmentCounts object using bam file with RG tags 
#'
#' @param bam filename for bam file with aligned reads
#' @param peaks GRanges object with peaks 
#' @export
getFragmentCountsByRG <- function(bam, peaks, BPPARAM = BiocParallel::bpparam()){

  rg_fragments <- bamToFragmentsByRG(bam, BPPARAM)
  depth <- lapply(rg_fragments, length)
  counts_mat <- Matrix::Matrix(simplify2array(BiocParallel::bplapply(rg_fragments,
                                                         function(fragments) GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=T),
                                                         BPPARAM = BPPARAM)))

  counts <- fragmentCounts(counts = counts_mat,
                           peaks = peaks,
                           depth = depth)

  return(counts)
}



#' make fragmentCounts object using multiple bam files
#'
#' @param bams filenames for bam file with aligned reads
#' @param peaks GRanges object with peaks 
#' @export
getFragmentCounts <- function(bams, peaks){

  depth = rep(0, length(bams))
  mat = matrix(nrow = length(peaks), ncol = length(bams))
  
  for (i in 1:length(bams)){
    fragments = bamToFragments(bams[i])
    depth[i] = length(fragments)
    mat[,i] = GenomicRanges::countOverlaps(peaks, fragments, type="any", ignore.strand=T)
  }
  
  counts_mat = Matrix::Matrix(mat)
  
  counts <- fragmentCounts(counts = counts_mat,
                           peaks = peaks,
                           depth = depth)

  return(counts)
}

#not exported
left_right_to_grglist <- function(left, right){
  stopifnot(length(left) == length(right))
  if (length(left) == 0){
    return(GenomicRanges::GenomicRangesList())
  }
  x = c(left,right)[as.vector(matrix(seq_len(2L * length(left)), nrow=2L, byrow=TRUE))]
  p = IRanges::PartitioningByEnd(cumsum(rep(2,length(x)/2)))
  out = BiocGenerics::relist(x,p)
  return(out)
}


#' @export
bamToFragmentsByRG <- function(bamfile, BPPARAM = BiocParallel::bpparam()){

  scanned <- Rsamtools::scanBam(bamfile,
                                param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE),
                                                                         what = c("rname","pos","isize"),
                                                                         tag = "RG"))[[1]]
  RG_tags <- gtools::mixedsort(unique(scanned$tag$RG), decreasing = TRUE)

  out <- BiocParallel::bplapply(RG_tags, function(RG){
    match_RG = which(scanned$tag$RG == RG)
    scanned_left <- as(GenomicRanges::GRanges(seq = scanned$rname[match_RG],
                                              IRanges::IRanges(start = scanned$pos[match_RG], width = 1), strand = "+"),
                       "GAlignments")
    scanned_right <- as(GenomicRanges::GRanges(seq = scanned$rname[match_RG],
                                               IRanges::IRanges(start = scanned$pos[match_RG] + abs(scanned$isize[match_RG]) - 1, width = 1), strand = "-"),
                        "GAlignments")
    ##Convert to GRangesList
    left_right_to_grglist(scanned_left, scanned_right)
    #grl <- GenomicAlignments::grglist(c(scanned_left, scanned_right)[as.vector(matrix(seq_len(2L * length(scanned_left)), nrow=2L, byrow=TRUE))],
    #                                  order.as.in.query=TRUE,
    #                                  drop.D.ranges=FALSE)
    #GenomicAlignments:::shrinkByHalf(grl)
  }, BPPARAM = BPPARAM)

  names(out) = RG_tags

  return(out)
}

#' @export
bamToFragments <- function(bamfile){

  scanned <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE), what = c("rname","pos","isize")))[[1]]

  scanned_left <- as(GenomicRanges::GRanges(seq = scanned$rname, IRanges::IRanges(start = scanned$pos, width = 1), strand = "+"), "GAlignments")
  scanned_right <- as(GenomicRanges::GRanges(seq = scanned$rname, IRanges::IRanges(start = scanned$pos + abs(scanned$isize) - 1, width = 1), strand = "-"), "GAlignments")
  ##Convert to GRangesList
  left_right_to_grglist(scanned_left, scanned_right)
  #grl <- GenomicAlignments::grglist(c(scanned_left, scanned_right)[as.vector(matrix(seq_len(2L * length(scanned_left)), nrow=2L, byrow=TRUE))],
  #                                  order.as.in.query=TRUE,
  #                                  drop.D.ranges=FALSE)
  #out <- GenomicAlignments:::shrinkByHalf(grl)

  return(out)

}


#' Filter Fragment Counts by quality
#'
#' @param counts_mat fragmentCounts object
#' @param min_in_peaks minimum fraction of fragments within peaks (default = 0.25)
#' @param min_fragments minimum fragments in peaks (default = 5000)
#' @export          
filterFragmentCounts <- function(counts_mat, min_in_peaks = 0.25, min_fragments = 5000){
  stopifnot(inherits(counts_mat, "fragmentCounts"))
  stopifnot((sum(counts_mat@depth) > sum(counts_mat@fragments_per_cell)) == counts_mat@nsample)
  keep <- intersect(which(counts_mat@fragments_per_cell >= min_fragments), 
                    which(counts_mat@fragments_per_sample/counts_mat@depth >= min_in_peaks))
            
  counts_mat <- counts_mat[,keep]            
  return(counts_mat)
}






