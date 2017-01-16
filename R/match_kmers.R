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

  tmp1 <- elementLengths(matches)
  tmp2 <- unlist(sapply(1:length(peaks), function(x) rep(x, tmp1[x])),
                 use.names = FALSE)
  f_pos <- resize(shift(peaks[tmp2], shift = start(unlist(matches))) + 1,
                  width = 1)
  BiocGenerics::strand(f_pos) <- "+"

  tmp1 <- elementLengths(rc_matches)
  tmp2 <- unlist(sapply(1:length(peaks), function(x) rep(x, tmp1[x])),
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


#' match_kmers
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
#' @details  Can either return a SummarizedExperiment with just sparse matrix
#' with values set to 1 for a match (if return == 'matches'),  or a
#' GenomicRanges  object with all the positions of matches
#'@export
setGeneric("match_kmers",
           function(k, subject, ...) standardGeneric("match_kmers"))

#' @describeIn match_kmers For DNAStringSet Objects
#' @export
setMethod("match_kmers", signature(k = "character", subject = "DNAStringSet"),
          function(k,
                   subject, genome = NULL,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

            match_kmers_helper(subject, k, out, ranges)
          })

#' @describeIn match_kmers For character strings
#' @export
setMethod("match_kmers", signature(k = "character", subject = "character"),
          function(k,
                   subject,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

  match_kmers_helper(subject, k, out, ranges)
})

#' @describeIn match_kmers For DNA String objects
#' @export
setMethod("match_kmers", signature(k = "character", subject = "DNAString"),
          function(k,
                   subject,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)

            seqs <- as.character(subject)

  match_kmers_helper(seqs, k, out, ranges)
})

#' @describeIn match_kmers For GenomicRanges
#' @export
setMethod("match_kmers", signature(k = "character", subject = "GenomicRanges"),
          function(k,
                   subject, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                   out = c("matches", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)
            GenomicRanges::strand(subject) <- "+"
  seqs <- getSeq(genome, subject)

  match_kmers_helper(seqs, k, out, subject)
})

#' @describeIn match_kmers For RangedSummarizedExperiment (containing GRanges in rowRanges)
#' @export
setMethod("match_kmers", signature(k = "character",
                                   subject = "RangedSummarizedExperiment"),
          function(k, subject,
                   genome = BSgenome.Hsapiens.UCSC.hg19,
                   out = c("matches", "counts", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)
            match_kmers(k, rowRanges(subject), genome, out)
          })


#' @describeIn match_kmers Catch-all for other un-documented types
#' @export
setMethod("match_kmers", signature(k = "numeric", subject = "ANY"),
          function(k, subject,
                   genome = BSgenome.Hsapiens.UCSC.hg19,
                   out = c("matches", "counts", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)
            kmers <- DNAStringSet(mkAllStrings(c("A", "C", "G", "T"),
                                     width = k))
  kmers <- remove_rc(kmers)
  match_kmers(kmers, subject, genome, out, ranges)
})


#' @describeIn match_kmers Catch-all for other un-documented types with DNAStringSet
#' @export
setMethod("match_kmers", signature(k = "DNAStringSet", subject = "ANY"),
          function(k,
                   subject,
                   genome = BSgenome.Hsapiens.UCSC.hg19,
                   out = c("matches", "counts", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)
            kmers <- as.character(k)
            match_kmers(kmers, subject, genome, out, ranges)
})


#' @describeIn match_kmers Catch-all for other un-documented types with DNAString
#' @export
setMethod("match_kmers", signature(k = "DNAString", subject = "ANY"),
          function(k,
                   subject,
                   genome = BSgenome.Hsapiens.UCSC.hg19,
                   out = c("matches", "counts", "positions"),
                   ranges = NULL) {
            out <- match.arg(out)
            kmers <- as.character(k)
  match_kmers(kmers, subject, genome, out, ranges)
})

