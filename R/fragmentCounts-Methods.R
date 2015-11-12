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
            x@fragments_per_cell = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@peaks = x@peaks[i]
            x@total_fragments = sum(x@counts)
            return(x)
          })

setMethod("[", signature = signature(x = "fragmentCounts"),
          definition = function(x, i ,j ) {
            x@counts = x@counts[i,j]
            x@fragments_per_cell = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@peaks = x@peaks[i]
            x@total_fragments = sum(x@counts)
            return(x)
          })

setMethod("[", signature = signature(x = "fragmentCounts", i = "missing"),
          definition = function(x, i ,j ) {
            x@counts = x@counts[,j]
            x@fragments_per_cell = colSums(x@counts)
            x@fragments_per_peak = rowSums(x@counts)
            x@total_fragments = sum(x@counts)
            return(x)
          })



getFragmentCountsByRG <- function(bam, bed, BPPARAM = BiocParallel::bpparam()){

  rg_fragments <- bamToFragmentsByRG(bam, BPPARAM)
  counts_mat <- Matrix(simplify2array(BiocParallel::bplapply(rg_fragments,
                                                         function(fragments) GenomicRanges::countOverlaps(bed, fragments, type="any", ignore.strand=T),
                                                         BPPARAM = BPPARAM)))

  counts <- fragmentCounts(counts = counts_mat,
                           peaks = bed,
                           total_fragments = sum(counts_mat),
                           fragments_per_cell = colSums(counts_mat),
                           fragments_per_peak = rowSums(counts_mat))

  return(counts)
}




getFragmentCounts <- function(bams, bed, BPPARAM = BiocParallel::bpparam()){

  counts_mat <- Matrix(simplify2array(BiocParallel::bplapply(bams,
                                                         function(bam){
                                                           fragments = bamToFragments(bam)
                                                           GenomicRanges::countOverlaps(bed, fragments, type="any", ignore.strand=T)
                                                         },
                                                         BPPARAM = BPPARAM)))
  counts <- fragmentCounts(counts = counts_mat,
                           peaks = bed,
                           total_fragments = sum(counts_mat),
                           fragments_per_cell = colSums(counts_mat),
                           fragments_per_peak = rowSums(counts_mat))

  return(counts)
}

getTotalCountsByRG <- function(bamfile){
  scanned <- Rsamtools::scanBam(bamfile,
                                param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE),
                                                                tag = "RG"))[[1]]
  RG_tags <- gtools::mixedsort(unique(scanned$tag$RG), decreasing = TRUE)
  rg_factor <- factor(scanned$tag$RG, levels = RG_tags, ordered = TRUE)
  RG_counts <- tabulate(rg_factor)
  names(RG_counts) = RG_tags
  return(RG_counts)
}


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
    grl <- GenomicAlignments::grglist(c(scanned_left, scanned_right)[S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(scanned_left))],
                                      order.as.in.query=TRUE,
                                      drop.D.ranges=FALSE)
    GenomicAlignments:::shrinkByHalf(grl)
  }, BPPARAM = BPPARAM)

  names(out) = RG_tags

  return(out)
}

bamToFragments <- function(bam){

  scanned <- Rsamtools::scanBam(bamfile, param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isMinusStrand=FALSE, isProperPair = TRUE), what = c("rname","pos","isize")))[[1]]

  scanned_left <- as(GenomicRanges::GRanges(seq = scanned$rname, IRanges::IRanges(start = scanned$pos, width = 1), strand = "+"), "GAlignments")
  scanned_right <- as(GenomicRanges::GRanges(seq = scanned$rname, IRanges::IRanges(start = scanned$pos + abs(scanned$isize) - 1, width = 1), strand = "-"), "GAlignments")
  ##Convert to GRangesList
  grl <- GenomicAlignments::grglist(c(scanned_left, scanned_right)[S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(scanned_left))],
                                    order.as.in.query=TRUE,
                                    drop.D.ranges=FALSE)
  out <- GenomicAlignments:::shrinkByHalf(grl)

  return(out)

}

setGeneric("filterFragmentCounts", function(counts_mat, total_fragments, ...) standardGeneric("filterFragmentCounts"))

setMethod("filterFragmentCounts", "fragmentCounts",
          function(counts_mat, total_fragments, min_in_peaks = 0.25, min_fragments = 5000){
            keep <- intersect(which(total_fragments >= min_fragments), which(counts_mat@fragments_per_cell/total_fragments >= min_in_peaks))
            counts_mat <- counts_mat[,keep]
            return(counts_mat)
          })






