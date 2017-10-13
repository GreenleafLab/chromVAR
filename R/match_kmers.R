remove_rc <- function(seqs) {
  temp <- cbind(as.character(seqs), as.character(reverseComplement(seqs)))
  temp[temp[, 1] > temp[, 2], 1] <- temp[temp[, 1] > temp[, 2], 2]
  DNAStringSet(unique(temp[, 1]))
}


match_kmers_inner_helper <- function(kmers, seqs, var = FALSE) {
  if (is.character(kmers))
    kmers <- DNAStringSet(kmers)
  stopifnot(inherits(kmers, "DNAStringSet"))
  if (!all_true(width(kmers) == width(kmers[1])) || var) {
    indices <- vwhichPDict(kmers, seqs, fixed = FALSE)
    indices_rc <- vwhichPDict(kmers,
                              reverseComplement(seqs),
                              fixed = FALSE)
  } else {
    pd <- PDict(kmers)
    indices <- vwhichPDict(pd, seqs)
    indices_rc <- vwhichPDict(pd, reverseComplement(seqs))
  }
  indices <- merge_lists(indices, indices_rc, by = "order")
  indices <- lapply(indices, unique)

  out <- sparseMatrix(i = unlist(lapply(seq_along(indices),
                                        function(x)
                                          rep(x, length(indices[[x]]))),
                                 use.names = FALSE),
                      j = unlist(indices, use.names = FALSE),
                      x = TRUE,
                      dims = c(length(seqs), length(kmers)),
                      dimnames = list(NULL, as.character(kmers)))

  return(out)
}

get_kmer_positions <- function(kmer, peaks, seqs) {
  matches <- vmatchPattern(DNAString(kmer), seqs, fixed = FALSE)
  rc_matches <- vmatchPattern(reverseComplement(DNAString(kmer)),
                              seqs, fixed = FALSE)

  tmp1 <- elementNROWS(matches)
  tmp2 <- unlist(lapply(seq_along(peaks), function(x) rep(x, tmp1[x])),
                 use.names = FALSE)
  f_pos <- resize(shift(peaks[tmp2], shift = start(unlist(matches))) + 1,
                  width = 1)
  BiocGenerics::strand(f_pos) <- "+"

  tmp1 <- elementNROWS(rc_matches)
  tmp2 <- unlist(lapply(seq_along(peaks), function(x) rep(x, tmp1[x])),
                 use.names = FALSE)
  r_pos <- resize(shift(peaks[tmp2], shift = start(unlist(rc_matches)) - 1),
                  width = 1)
  BiocGenerics::strand(r_pos) <- "-"

  return(BiocGenerics::sort(c(f_pos, r_pos)))
}

match_kmers_helper <- function(seqs, kmers, out, ranges) {

  if (out == "matches") {
    kmer_ix <- match_kmers_inner_helper(kmers, seqs)
    if (is.null(ranges)) {
      out <- SummarizedExperiment(assays = list(matches = kmer_ix),
                                  colData = DataFrame(name = kmers))
    } else {
      out <- SummarizedExperiment(assays = list(matches = kmer_ix),
                                  rowRanges = ranges,
                                  colData = DataFrame(name = kmers))
    }
  } else if (out == "positions") {
    stopifnot(!is.null(ranges))
    out <- lapply(kmers, get_kmer_positions, ranges, seqs)
  }
  return(out)
}


#' matchKmers
#'
#' Find kmer matches in the DNA string-based subject
#'
#'
#' @param subject either \code{\link[GenomicRanges]{GenomicRanges}},
#' \code{\link[Biostrings]{DNAStringSet}}, \code{\link[Biostrings]{DNAString}},
#' or character vector
#' @param k k
#' @param genome BSgenome object, only used if subect is
#' \code{\link[GenomicRanges]{GenomicRanges}}
#' @param out what to return? see details
#' @param ranges if subject is not GenomicRanges, ranges to use when out is
#' positions
#' @param ... additional arguments
#' @seealso \code{\link{getAnnotations}}, \code{\link{computeDeviations}}
#' @details  Can either return a SummarizedExperiment with just sparse matrix
#' with values set to 1 for a match (if return == 'matches'),  or a
#' GenomicRanges  object with all the positions of matches
#' @return SummarizedExperiment with matches assay storing which peaks contain
#' which kmers
#' @export
#' @examples 
#' 
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' 
#' # Get peak-kmer annotation matrix for 6mers
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' kmer_ix <- matchKmers(6, mini_counts, 
#'                        genome = BSgenome.Hsapiens.UCSC.hg19)
setGeneric("matchKmers",
           function(k, subject, ...) standardGeneric("matchKmers"))

#' @describeIn matchKmers For DNAStringSet Objects
#' @export
setMethod("matchKmers", signature(k = "character", subject = "DNAStringSet"),
          function(k,
                   subject,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

            match_kmers_helper(subject, k, out, ranges)
          })

#' @describeIn matchKmers For character strings
#' @export
setMethod("matchKmers", signature(k = "character", subject = "character"),
          function(k,
                   subject,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

  match_kmers_helper(subject, k, out, ranges)
})

#' @describeIn matchKmers For DNA String objects
#' @export
setMethod("matchKmers", signature(k = "character", subject = "DNAString"),
          function(k,
                   subject,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

            seqs <- as.character(subject)

  match_kmers_helper(seqs, k, out, ranges)
})

#' @describeIn matchKmers For GenomicRanges
#' @export
setMethod("matchKmers", signature(k = "character", subject = "GenomicRanges"),
          function(k,
                   subject, 
                   genome = GenomeInfoDb::genome(subject),
                   out = c("matches", "positions")) {
            out <- match.arg(out)
            GenomicRanges::strand(subject) <- "+"
            genome <- validate_genome_input(genome)
            seqs <- getSeq(genome, subject)
            
            match_kmers_helper(seqs, k, out, subject)
          })

#' @describeIn matchKmers For RangedSummarizedExperiment (containing GRanges in
#'  rowRanges)
#' @export
setMethod("matchKmers", signature(k = "character",
                                   subject = "RangedSummarizedExperiment"),
          function(k, subject,
                   ...) {
            matchKmers(k, rowRanges(subject), ...)
          })


#' @describeIn matchKmers Catch-all for other un-documented types
#' @export
setMethod("matchKmers", signature(k = "numeric", subject = "ANY"),
          function(k, subject, ...) {
            
            kmers <- DNAStringSet(mkAllStrings(c("A", "C", "G", "T"),
                                               width = k))
            kmers <- remove_rc(kmers)
            matchKmers(kmers, subject,...)
          })


#' @describeIn matchKmers Catch-all for other un-documented types with 
#' DNAStringSet
#' @export
setMethod("matchKmers", signature(k = "DNAStringSet", subject = "ANY"),
          function(k,
                   subject,
                   ...) {
            kmers <- as.character(k)
            matchKmers(kmers, subject, ...)
})


#' @describeIn matchKmers Catch-all for other un-documented types with DNAString
#' @export
setMethod("matchKmers", signature(k = "DNAString", subject = "ANY"),
          function(k,
                   subject,
                   ...) {
            
            kmers <- as.character(k)
            matchKmers(kmers, subject, ...)
          })

