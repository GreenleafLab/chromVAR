context("matchKmers")

library(Biostrings)
library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)

rc_string <- function(x) {
  vapply(strsplit(chartr("ACGT", "TGCA", x), ""),
         function(y) paste(rev(y), collapse = ""), "")
}

make_seq <- function(n, a, b) {
  paste(c("A", "C", "G", "T")[floor(abs(sin(seq_len(n) * a + b)) * 1e4) %% 4 + 1],
        collapse = "")
}

seq_strings <- c(make_seq(60, 1.1, 0.3), make_seq(45, 2.3, 1.7),
                 "AAAAAAAAAA", "ACGTACGTAA", "TTTTGGGGCC", "ACG")
subject_set <- DNAStringSet(seq_strings)

count_all <- function(pattern, string) {
  hits <- gregexpr(paste0("(?=", pattern, ")"), string, perl = TRUE)[[1]]
  if (hits[1] == -1) integer(0) else as.integer(hits)
}

brute_matches <- function(kmers, seqs) {
  out <- matrix(FALSE, nrow = length(seqs), ncol = length(kmers))
  for (i in seq_along(seqs)) {
    for (j in seq_along(kmers)) {
      out[i, j] <- grepl(kmers[j], seqs[i], fixed = TRUE) ||
        grepl(rc_string(kmers[j]), seqs[i], fixed = TRUE)
    }
  }
  out
}

all_kmers <- function(k) {
  apply(expand.grid(rep(list(c("A", "C", "G", "T")), k),
                    stringsAsFactors = FALSE)[, rev(seq_len(k)), drop = FALSE],
        1, paste, collapse = "")
}

canonical_kmers <- function(k) {
  km <- all_kmers(k)
  sort(unique(pmin(km, rc_string(km))))
}

brute_positions <- function(kmer, seqs, chrs, starts) {
  k <- nchar(kmer)
  out <- character(0)
  for (i in seq_along(seqs)) {
    for (p in count_all(kmer, seqs[i])) {
      out <- c(out, paste(chrs[i], starts[i] + p - 1, "+"))
    }
    for (p in count_all(rc_string(kmer), seqs[i])) {
      out <- c(out, paste(chrs[i], starts[i] + p - 1, "-"))
    }
  }
  sort(out)
}

position_keys <- function(gr) {
  sort(paste(as.character(seqnames(gr)), start(gr), as.character(strand(gr))))
}

test_that("matches for equal-width kmers agree with brute-force search", {
  kmers <- c("ACGT", "AAAA", "GGGT", "CCCA", "TTGA", "ACGA")
  res <- matchKmers(kmers, subject_set)
  expect_equal(unname(as.matrix(assay(res, "matches"))),
               unname(brute_matches(kmers, seq_strings)))
  expect_true(any(brute_matches(kmers, seq_strings)))
  expect_false(all(brute_matches(kmers, seq_strings)))
})

test_that("matches for kmers of mixed widths agree with brute-force search", {
  kmers <- c("ACG", "AAAAA", "GGTT", "CCCAAAC", "TA")
  res <- matchKmers(kmers, subject_set)
  expect_equal(unname(as.matrix(assay(res))),
               unname(brute_matches(kmers, seq_strings)))
})

test_that("matches are found on the reverse complement strand", {
  res <- matchKmers("AAAA", DNAStringSet(c("TTTTCC", "GGGG", "CCCC")))
  expect_equal(as.vector(as.matrix(assay(res))), c(TRUE, FALSE, FALSE))
  res <- matchKmers("CCCC", DNAStringSet(c("TTTTCC", "GGGG", "CCCC")))
  expect_equal(as.vector(as.matrix(assay(res))), c(FALSE, TRUE, TRUE))
})

test_that("kmers of mixed widths are also found on the reverse complement strand", {
	res <- matchKmers(c("AAAA", "CC", "GATC"),
	                  DNAStringSet(c("TTTTGG", "ACGT", "GGGG")))
	expect_equal(unname(as.matrix(assay(res))),
	             matrix(c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE,
	                      FALSE), nrow = 3))
})

test_that("overlapping occurrences and short sequences are handled", {
  res <- matchKmers("AAA", DNAStringSet(c("AAAAA", "AA", "")))
  expect_equal(as.vector(as.matrix(assay(res))), c(TRUE, FALSE, FALSE))
})

test_that("the matches result is a SummarizedExperiment with a sparse assay", {
  kmers <- c("ACGT", "AAAA", "GGGT")
  res <- matchKmers(kmers, subject_set)
  expect_true(is(res, "SummarizedExperiment"))
  expect_false(is(res, "RangedSummarizedExperiment"))
  expect_equal(assayNames(res), "matches")
  mat <- assay(res, "matches")
  expect_true(is(mat, "lgCMatrix"))
  expect_equal(dim(mat), c(length(seq_strings), 3))
  expect_equal(dim(res), c(length(seq_strings), 3))
  expect_null(rownames(mat))
  expect_equal(colnames(mat), kmers)
  expect_equal(colnames(res), kmers)
  expect_equal(as.character(colData(res)$name), kmers)
})

test_that("a sequence without matches gives an all-FALSE row", {
  res <- matchKmers(c("CGCG", "GGGG"), DNAStringSet(c("AAAAAA", "ACGTAC")))
  expect_equal(unname(as.matrix(assay(res))[1, ]), c(FALSE, FALSE))
  expect_equal(dim(assay(res)), c(2, 2))
})

test_that("duplicated kmers give identical columns", {
  res <- matchKmers(c("ACGT", "ACGT"), subject_set)
  mat <- as.matrix(assay(res))
  expect_equal(unname(mat[, 1]), unname(mat[, 2]))
})

test_that("ranges are stored as rowRanges when given", {
  pk <- GRanges("chr1", IRanges(start = seq(1, by = 100, length.out = 6),
                                width = 10))
  res <- matchKmers(c("ACGT", "AAAA"), subject_set, ranges = pk)
  expect_true(is(res, "RangedSummarizedExperiment"))
  expect_equal(start(rowRanges(res)), start(pk))
  expect_equal(dim(res), c(6, 2))
  expect_equal(unname(as.matrix(assay(res))),
               unname(brute_matches(c("ACGT", "AAAA"), seq_strings)))
})

test_that("numeric k uses the canonical kmers of that length", {
  for (k in 1:4) {
    res <- matchKmers(k, subject_set)
    expect_equal(sort(colnames(res)), canonical_kmers(k))
    expect_equal(unname(as.matrix(assay(res))),
                 unname(brute_matches(colnames(res), seq_strings)))
  }
  expect_equal(ncol(matchKmers(3, subject_set)), 32)
  expect_equal(ncol(matchKmers(4, subject_set)), 136)
})

test_that("every sequence of length k matches exactly its canonical kmer", {
  km <- all_kmers(3)
  res <- matchKmers(3, DNAStringSet(km))
  mat <- as.matrix(assay(res))
  canon <- pmin(km, rc_string(km))
  expect_equal(unname(colSums(mat)[canon]), as.vector(table(canon)[canon]))
  expect_equal(unname(rowSums(mat)), rep(1, length(km)))
})

test_that("kmers can be given as DNAStringSet or DNAString", {
  kmers <- c("ACGT", "AAAA", "GGGT")
  expected <- brute_matches(kmers, seq_strings)
  res <- matchKmers(DNAStringSet(kmers), subject_set)
  expect_equal(unname(as.matrix(assay(res))), unname(expected))
  expect_equal(colnames(res), kmers)
  res <- matchKmers(DNAString("GGGT"), subject_set)
  expect_equal(as.vector(as.matrix(assay(res))), expected[, 3])
  expect_equal(colnames(res), "GGGT")
})

test_that("out = 'positions' returns the match positions per kmer", {
  starts <- c(1000, 5000, 20, 300, 7000, 90)
  chrs <- c("chr1", "chr2", "chr1", "chr1", "chr3", "chr2")
  pk <- GRanges(chrs, IRanges(start = starts, width = nchar(seq_strings)))
  kmers <- c("AAA", "ACGT", "GGTT", "TTCA")
  res <- matchKmers(kmers, subject_set, out = "positions", ranges = pk)
  expect_equal(length(res), length(kmers))
  for (i in seq_along(kmers)) {
    expect_true(is(res[[i]], "GRanges"))
    expect_equal(position_keys(res[[i]]),
                 brute_positions(kmers[i], seq_strings, chrs, starts))
    expect_true(all(width(res[[i]]) == 1))
    expect_false(is.unsorted(res[[i]]))
  }
  expect_true(length(res[[1]]) > 0)
})

test_that("positions of a palindromic kmer are reported on both strands", {
  pk <- GRanges("chr1", IRanges(start = 101, width = 10))
  res <- matchKmers("ACGT", DNAStringSet("ACGTACGTAA"), out = "positions",
                    ranges = pk)
  expect_equal(position_keys(res[[1]]),
               sort(c("chr1 101 +", "chr1 105 +", "chr1 101 -", "chr1 105 -")))
})

test_that("positions are empty for a kmer that does not occur", {
  pk <- GRanges("chr1", IRanges(start = 101, width = 10))
  res <- matchKmers("CCCC", DNAStringSet("ACGTACGTAA"), out = "positions",
                    ranges = pk)
  expect_equal(length(res[[1]]), 0)
})

test_that("out = 'positions' requires ranges", {
  expect_error(matchKmers("ACGT", subject_set, out = "positions"),
               "ranges")
})

test_that("out must be matches or positions", {
  expect_error(matchKmers("ACGT", subject_set, out = "counts"),
               "should be one of")
})

genome_set <- DNAStringSet(c(chr1 = make_seq(80, 1.9, 0.4),
                             chr2 = make_seq(70, 0.8, 2.2)))
genome_strings <- as.character(genome_set)
peak_ranges <- GRanges(c("chr1", "chr2", "chr1", "chr2"),
                       IRanges(start = c(1, 5, 41, 31),
                               width = c(30, 30, 30, 30)),
                       strand = c("-", "*", "+", "-"))
peak_seqs <- mapply(substr, genome_strings[as.character(seqnames(peak_ranges))],
                    start(peak_ranges), end(peak_ranges), USE.NAMES = FALSE)

test_that("GenomicRanges subjects are scanned in a DNAStringSet genome", {
  kmers <- c("ACG", "TTT", "GCA", "CAA")
  res <- matchKmers(kmers, peak_ranges, genome = genome_set)
  expect_true(is(res, "RangedSummarizedExperiment"))
  expect_equal(unname(as.matrix(assay(res))),
               unname(brute_matches(kmers, peak_seqs)))
  expect_equal(as.character(seqnames(rowRanges(res))),
               as.character(seqnames(peak_ranges)))
  expect_equal(start(rowRanges(res)), start(peak_ranges))
  expect_equal(end(rowRanges(res)), end(peak_ranges))
})

test_that("the strand of GenomicRanges subjects does not change the matches", {
  kmers <- c("ACG", "TTT", "GCA", "CAA")
  stranded <- matchKmers(kmers, peak_ranges, genome = genome_set)
  unstranded <- peak_ranges
  strand(unstranded) <- "*"
  expect_equal(as.matrix(assay(stranded)),
               as.matrix(assay(matchKmers(kmers, unstranded,
                                          genome = genome_set))))
})

test_that("GenomicRanges subjects work with numeric k and positions", {
  res <- matchKmers(3, peak_ranges, genome = genome_set)
  expect_equal(unname(as.matrix(assay(res))),
               unname(brute_matches(colnames(res), peak_seqs)))
  kmers <- c("ACG", "TT")
  pos <- matchKmers(kmers, peak_ranges, genome = genome_set,
                    out = "positions")
  for (i in seq_along(kmers)) {
    expect_equal(position_keys(pos[[i]]),
                 brute_positions(kmers[i], peak_seqs,
                                 as.character(seqnames(peak_ranges)),
                                 start(peak_ranges)))
  }
})

test_that("RangedSummarizedExperiment subjects use their rowRanges", {
  rse <- SummarizedExperiment(assays = list(counts = matrix(0, 4, 2)),
                              rowRanges = peak_ranges)
  kmers <- c("ACG", "TTT")
  res <- matchKmers(kmers, rse, genome = genome_set)
  ref <- matchKmers(kmers, peak_ranges, genome = genome_set)
  expect_equal(as.matrix(assay(res)), as.matrix(assay(ref)))
  expect_equal(dim(res), c(4, 2))
})

test_that("a GenomicRanges subject needs a usable genome", {
  expect_error(matchKmers("ACG", peak_ranges), "No genome provided")
  expect_error(matchKmers("ACG", peak_ranges, genome = 1),
               "genome input must be")
})

test_that("a BSgenome genome gives the same matches as its sequences", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
  genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  pk <- GRanges("chrM", IRanges(start = c(1, 1001, 5001), width = 200))
  seqs <- as.character(getSeq(genome, pk))
  kmers <- c("ACG", "TTTA", "GGC")
  res <- matchKmers(kmers, pk, genome = genome)
  expect_equal(unname(as.matrix(assay(res))),
               unname(brute_matches(kmers, seqs)))
})
