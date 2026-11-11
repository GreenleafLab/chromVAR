context("get_inputs: peaks, counts, depths")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(SummarizedExperiment)
})
local_serial_param()

# Helpers -----------------------------------------------------------------------

write_tmp <- function(lines, ext = ".txt") {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}

quiet <- function(x) suppressWarnings(suppressMessages(x))

sam_line <- function(name, flag, chr, pos, mate_pos, isize, len = 10,
                     rg = NULL) {
  paste(c(name, flag, chr, pos, 60, paste0(len, "M"),
          if (mate_pos == 0) "*" else "=", mate_pos, isize,
          strrep("A", len), strrep("I", len),
          if (!is.null(rg)) paste0("RG:Z:", rg)), collapse = "\t")
}

# Both reads of a proper pair whose fragment covers left..(left + size - 1).
pair_lines <- function(name, chr, left, size, rg = NULL, proper = TRUE,
                       len = 10) {
  flags <- if (proper) c(99, 147) else c(97, 145)
  right_read <- left + size - len
  c(sam_line(name, flags[1], chr, left, right_read, size, len, rg),
    sam_line(name, flags[2], chr, right_read, left, -size, len, rg))
}

make_bam <- function(lines, seqs = c(chr1 = 1000, chr2 = 1000)) {
  sam <- tempfile(fileext = ".sam")
  writeLines(c("@HD\tVN:1.0\tSO:unsorted",
               paste0("@SQ\tSN:", names(seqs), "\tLN:", seqs), lines), sam)
  Rsamtools::asBam(sam, tempfile(), indexDestination = FALSE, overwrite = TRUE)
}

dense_counts <- function(x) as.matrix(assays(x)$counts)

# Reference implementation: a fragment is counted in a peak when either of its
# two end positions lies in the peak.
ref_fragment_counts <- function(chr, left, right, peaks) {
  pc <- as.character(seqnames(peaks))
  vapply(seq_along(peaks), function(i) {
    in_peak <- function(p) p >= start(peaks)[i] & p <= end(peaks)[i]
    sum(chr == pc[i] & (in_peak(left) | in_peak(right)))
  }, 0)
}

proper_fragments <- function(bam, rg = NULL) {
  s <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(
    what = c("rname", "pos", "isize", "flag"), tag = "RG"))[[1]]
  keep <- bitwAnd(s$flag, 2L) > 0 & bitwAnd(s$flag, 16L) == 0
  if (!is.null(rg)) keep <- keep & s$tag$RG %in% rg
  data.frame(chr = as.character(s$rname)[keep], left = s$pos[keep],
             right = s$pos[keep] + abs(s$isize[keep]) - 1,
             rg = s$tag$RG[keep], stringsAsFactors = FALSE)
}

peaks_around <- function(fr, every = 5, half = 25) {
  ends <- data.frame(chr = c(fr$chr, fr$chr), pos = c(fr$left, fr$right),
                     stringsAsFactors = FALSE)
  ends <- ends[seq(1, nrow(ends), by = every), ]
  GRanges(ends$chr, IRanges(pmax(1, ends$pos - half), ends$pos + half))
}

# Five proper pairs on chr1/chr2 with fragments
#   f1 chr1 11-40, f2 chr1 101-210, f3 chr2 5-24, f4 chr1 301-330,
#   f5 chr1 11-40 (duplicate of f1)
frag_reads <- function(rgs = c(NA, NA, NA, NA, NA)) {
  r <- function(x) if (is.na(x)) NULL else x
  c(pair_lines("f1", "chr1", 11, 30, r(rgs[1])),
    pair_lines("f2", "chr1", 101, 110, r(rgs[2])),
    pair_lines("f3", "chr2", 5, 20, r(rgs[3])),
    pair_lines("f4", "chr1", 301, 30, r(rgs[4])),
    pair_lines("f5", "chr1", 11, 30, r(rgs[5])))
}

hand_peaks <- GRanges(
  c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr3"),
  IRanges(c(1, 35, 150, 300, 1, 20, 1), c(12, 45, 160, 340, 3, 30, 10)))

# getPeaks ----------------------------------------------------------------------

test_that("getPeaks converts 0-based bed starts to 1-based starts", {
  f <- write_tmp(c("chr1\t0\t10", "chr1\t20\t30"))
  peaks <- getPeaks(f)
  expect_equal(start(peaks), c(1, 21))
  expect_equal(end(peaks), c(10, 30))
  expect_equal(width(peaks), c(10, 10))
  expect_equal(as.character(seqnames(peaks)), c("chr1", "chr1"))
  expect_equal(as.character(strand(peaks)), c("*", "*"))
})

test_that("getPeaks keeps the values of extra columns in the given order", {
  f <- write_tmp(c("chr1\t0\t10\tpeakA\t3.5\t7", "chr1\t20\t30\tpeakB\t1\t8"))
  peaks <- getPeaks(f, extra_cols = c(score = 5, name = 4))
  expect_equal(colnames(mcols(peaks)), c("score", "name"))
  expect_equal(peaks$score, c(3.5, 1))
  expect_equal(peaks$name, c("peakA", "peakB"))
})

test_that("getPeaks uses an extra column named strand as the peak strand", {
  f <- write_tmp(c("chr1\t0\t10\t+", "chr1\t20\t30\t-"))
  peaks <- getPeaks(f, extra_cols = c(strand = 4))
  expect_equal(as.character(strand(peaks)), c("+", "-"))
  expect_equal(ncol(mcols(peaks)), 0)
})

test_that("getPeaks with sort_peaks orders chromosomes naturally and then by start", {
  f <- write_tmp(c("chr10\t100\t110", "chr2\t50\t60", "chr1\t300\t310",
                   "chr1\t5\t15", "chrX\t0\t10"))
  peaks <- quiet(getPeaks(f, sort_peaks = TRUE))
  expect_equal(as.character(seqnames(peaks)),
               c("chr1", "chr1", "chr2", "chr10", "chrX"))
  expect_equal(start(peaks), c(6, 301, 51, 101, 1))
})

test_that("getPeaks keeps non-standard chromosomes", {
  f <- write_tmp(c("chr1\t0\t10", "chrUn_gl000220\t0\t10", "chr1_random\t0\t10"))
  peaks <- quiet(getPeaks(f, sort_peaks = TRUE))
  expect_equal(length(peaks), 3)
  expect_true(all(c("chrUn_gl000220", "chr1_random") %in%
                    as.character(seqnames(peaks))))
  expect_equal(as.character(seqnames(peaks))[1], "chr1")
})

test_that("getPeaks without sort_peaks keeps the file order", {
  f <- write_tmp(c("chr1\t100\t110", "chr1\t0\t10"))
  peaks <- quiet(getPeaks(f))
  expect_equal(start(peaks), c(101, 1))
})

test_that("getPeaks announces when sorting changed unsorted peaks", {
  f <- write_tmp(c("chr2\t0\t10", "chr1\t0\t10"))
  expect_message(suppressWarnings(getPeaks(f, sort_peaks = TRUE)),
                 "Peaks sorted")
})

test_that("getPeaks does not announce sorting for already sorted peaks", {
  f <- write_tmp(c("chr1\t0\t10", "chr2\t0\t10"))
  expect_message(getPeaks(f, sort_peaks = TRUE), NA)
})

test_that("getPeaks warns about unsorted peaks only when not sorting", {
  f <- write_tmp(c("chr2\t0\t10", "chr1\t0\t10"))
  expect_warning(getPeaks(f), "Peaks not sorted")
  expect_warning(suppressMessages(getPeaks(f, sort_peaks = TRUE)), NA)
})

test_that("getPeaks warns about overlapping peaks", {
  f <- write_tmp(c("chr1\t0\t10", "chr1\t5\t15"))
  expect_warning(getPeaks(f), "Peaks are overlapping")
})

test_that("getPeaks warns about peaks of unequal width", {
  f <- write_tmp(c("chr1\t0\t10", "chr1\t20\t35"))
  expect_warning(getPeaks(f), "Peaks are not equal width")
})

test_that("getPeaks gives no warning for sorted, disjoint peaks of equal width", {
  f <- write_tmp(c("chr1\t0\t10", "chr1\t20\t30", "chr2\t0\t10"))
  expect_warning(getPeaks(f), NA)
})

test_that("getPeaks errors on a missing file", {
  expect_error(getPeaks(tempfile()), "does not exist|cannot open")
})

test_that("getPeaks errors when extra_cols point past the last column", {
  f <- write_tmp(c("chr1\t0\t10\tx"))
  expect_error(getPeaks(f, extra_cols = c(a = 5)),
               "subset columns|undefined columns")
})

test_that("getPeaks gives the same result without readr", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  f <- write_tmp(c("chr1\t0\t10\tpeakA\t3", "chr1\t20\t30\tpeakB\t4"))
  expected <- getPeaks(f, extra_cols = c(name = 4, score = 5))
  testthat::local_mocked_bindings(is.installed = function(pkg) FALSE,
                                  .package = "chromVAR")
  without_readr <- getPeaks(f, extra_cols = c(name = 4, score = 5))
  expect_equal(start(without_readr), start(expected))
  expect_equal(end(without_readr), end(expected))
  expect_equal(without_readr$name, expected$name)
  expect_equal(without_readr$score, expected$score)
})

# readNarrowpeaks -----------------------------------------------------------------

np_line <- function(chr, start, name, qval, summit, end = start + 200) {
  paste(chr, start, end, name, 1, ".", 5, 4, qval, summit, sep = "\t")
}

test_that("readNarrowpeaks returns peaks of the requested width", {
  f <- write_tmp(c(np_line("chr1", 100, "p1", 3, 100),
                   np_line("chr1", 1000, "p2", 4, 50),
                   np_line("chr2", 100, "p3", 5, 120)))
  for (w in c(50, 51, 300)) {
    peaks <- quiet(readNarrowpeaks(f, width = w))
    expect_equal(width(peaks), rep(w, 3))
  }
  expect_equal(width(quiet(readNarrowpeaks(f))), rep(500, 3))
})

test_that("readNarrowpeaks centers each peak within one base of its summit", {
  f <- write_tmp(c(np_line("chr1", 100, "p1", 3, 100),
                   np_line("chr1", 1000, "p2", 4, 50),
                   np_line("chr2", 100, "p3", 5, 120)))
  peaks <- quiet(readNarrowpeaks(f, width = 51))
  centers <- (start(peaks) + end(peaks)) / 2
  expect_true(all(abs(centers - c(200, 1050, 220)) <= 1))
})

test_that("readNarrowpeaks keeps name, score and qval", {
  f <- write_tmp(c(np_line("chr1", 100, "p1", 3, 100),
                   np_line("chr1", 1000, "p2", 4, 50)))
  peaks <- quiet(readNarrowpeaks(f, width = 50))
  expect_equal(sort(colnames(mcols(peaks))), c("name", "qval", "score"))
  expect_equal(peaks$name, c("p1", "p2"))
  expect_equal(peaks$qval, c(3, 4))
  expect_equal(peaks$score, c(1, 1))
})

test_that("readNarrowpeaks with non_overlapping = FALSE keeps every peak", {
  f <- write_tmp(c(np_line("chr1", 100, "p1", 3, 100),
                   np_line("chr1", 100, "p2", 4, 105),
                   np_line("chr1", 100, "p3", 5, 110)))
  peaks <- quiet(readNarrowpeaks(f, width = 100, non_overlapping = FALSE))
  expect_equal(length(peaks), 3)
  expect_equal(peaks$name, c("p1", "p2", "p3"))
})

test_that("readNarrowpeaks keeps the overlapping peak with the higher qval", {
  f <- write_tmp(c(np_line("chr1", 100, "low", 3, 100),
                   np_line("chr1", 100, "high", 6, 120)))
  peaks <- quiet(readNarrowpeaks(f, width = 100))
  expect_equal(peaks$name, "high")
  f <- write_tmp(c(np_line("chr1", 100, "high", 6, 100),
                   np_line("chr1", 100, "low", 3, 120)))
  expect_equal(quiet(readNarrowpeaks(f, width = 100))$name, "high")
})

test_that("readNarrowpeaks resolves chains of overlapping peaks", {
  f <- write_tmp(c(np_line("chr1", 0, "A", 5, 100),
                   np_line("chr1", 0, "B", 1, 140),
                   np_line("chr1", 0, "C", 3, 180),
                   np_line("chr1", 0, "D", 9, 400)))
  peaks <- quiet(readNarrowpeaks(f, width = 100))
  expect_equal(peaks$name, c("A", "D"))
  expect_true(isDisjoint(peaks))
})

test_that("readNarrowpeaks result is disjoint and sorted by chromosome", {
  f <- write_tmp(c(np_line("chr10", 0, "a", 5, 100),
                   np_line("chr10", 0, "b", 6, 130),
                   np_line("chr2", 0, "c", 5, 100),
                   np_line("chr1", 0, "d", 5, 100),
                   np_line("chr1", 0, "e", 5, 160),
                   np_line("chr1", 0, "f", 5, 900)))
  peaks <- quiet(readNarrowpeaks(f, width = 100))
  expect_true(isDisjoint(peaks))
  expect_equal(as.character(seqnames(peaks))[1], "chr1")
  expect_equal(match(c("chr1", "chr2", "chr10"),
                     as.character(seqnames(peaks))),
               sort(match(c("chr1", "chr2", "chr10"),
                          as.character(seqnames(peaks)))))
  expect_true("b" %in% peaks$name)
  expect_false("a" %in% peaks$name)
})

test_that("readNarrowpeaks keeps peaks with identical coordinates on different chromosomes", {
  f <- write_tmp(c(np_line("chr1", 0, "a", 5, 100),
                   np_line("chr2", 0, "b", 5, 100)))
  expect_equal(length(quiet(readNarrowpeaks(f, width = 100))), 2)
})

test_that("readNarrowpeaks keeps peaks that are far apart", {
  f <- write_tmp(c(np_line("chr1", 0, "a", 5, 100),
                   np_line("chr1", 5000, "b", 1, 100)))
  peaks <- quiet(readNarrowpeaks(f, width = 100))
  expect_equal(sort(peaks$name), c("a", "b"))
})

test_that("readNarrowpeaks errors on a missing file", {
  expect_error(readNarrowpeaks(tempfile()), "does not exist|cannot open")
})

test_that("readNarrowpeaks gives the same result without readr", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  f <- write_tmp(c(np_line("chr1", 0, "A", 5, 100),
                   np_line("chr1", 0, "B", 1, 140),
                   np_line("chr2", 0, "C", 3, 180)))
  expected <- quiet(readNarrowpeaks(f, width = 100))
  testthat::local_mocked_bindings(is.installed = function(pkg) FALSE,
                                  .package = "chromVAR")
  without_readr <- readNarrowpeaks(f, width = 100)
  expect_equal(start(without_readr), start(expected))
  expect_equal(end(without_readr), end(expected))
  expect_equal(without_readr$name, expected$name)
})

# counts and counts<- --------------------------------------------------------------

test_that("counts returns the counts assay", {
  m <- matrix(1:6, 3, 2)
  se <- SummarizedExperiment(assays = list(counts = m))
  expect_equal(counts(se), m)
})

test_that("counts errors when there is no counts assay", {
  se <- SummarizedExperiment(assays = list(other = matrix(1:6, 3, 2)))
  expect_error(counts(se), "counts")
})

test_that("counts<- replaces the counts assay with a matrix or a Matrix", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3, 2)))
  new_m <- matrix(6:1, 3, 2)
  counts(se) <- new_m
  expect_equal(counts(se), new_m)
  counts(se) <- Matrix::Matrix(new_m * 2, sparse = TRUE)
  expect_equal(as.matrix(counts(se)), new_m * 2)
  expect_true(is(counts(se), "Matrix"))
})

test_that("counts<- adds a counts assay when none exists and keeps other assays", {
  se <- SummarizedExperiment(assays = list(other = matrix(1:6, 3, 2)))
  counts(se) <- matrix(0, 3, 2)
  expect_equal(sort(assayNames(se)), c("counts", "other"))
  expect_equal(assays(se)$other, matrix(1:6, 3, 2))
})

test_that("counts<- rejects a value with the wrong number of rows", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3, 2)))
  expect_error(counts(se) <- matrix(1:4, 2, 2), "rows")
})

test_that("counts<- rejects a value that is neither matrix nor Matrix", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3, 2)))
  expect_error(counts(se) <- data.frame(a = 1:3, b = 1:3))
})

# left_right_to_grglist ----------------------------------------------------------------

test_that("left_right_to_grglist pairs each left end with its right end", {
  left <- GRanges(c("chr1", "chr2", "chr1"), IRanges(c(1, 5, 9), width = 1))
  right <- GRanges(c("chr1", "chr2", "chr1"), IRanges(c(10, 20, 30), width = 1))
  out <- chromVAR:::left_right_to_grglist(left, right)
  expect_equal(length(out), 3)
  expect_equal(unname(lengths(out)), c(2L, 2L, 2L))
  expect_equal(start(out[[1]]), c(1, 10))
  expect_equal(start(out[[2]]), c(5, 20))
  expect_equal(start(out[[3]]), c(9, 30))
  expect_equal(as.character(seqnames(out[[2]])), c("chr2", "chr2"))
})

test_that("left_right_to_grglist errors when the ends differ in number", {
  left <- GRanges("chr1", IRanges(1, width = 1))
  right <- GRanges(c("chr1", "chr1"), IRanges(c(5, 6), width = 1))
  expect_error(chromVAR:::left_right_to_grglist(left, right),
               "length\\(left\\) == length\\(right\\)")
})

# getCounts from bam: paired and single end --------------------------------------------

test_that("getCounts counts a paired-end fragment in each peak holding one of its ends", {
  bam <- make_bam(frag_reads())
  x <- quiet(getCounts(bam, hand_peaks, paired = TRUE))
  expect_is(x, "RangedSummarizedExperiment")
  expect_equal(unname(dense_counts(x)[, 1]), c(2, 2, 0, 1, 0, 1, 0))
  expect_equal(colData(x)$depth, 5)
  expect_equal(colnames(x), basename(bam))
  expect_equal(nrow(x), length(hand_peaks))
  expect_equal(start(rowRanges(x)), start(hand_peaks))
  expect_equal(as.character(seqnames(rowRanges(x))),
               as.character(seqnames(hand_peaks)))
  expect_true(is(assays(x)$counts, "Matrix"))
})

test_that("getCounts counts a fragment once when both ends fall in one peak", {
  bam <- make_bam(pair_lines("a", "chr1", 301, 30))
  peaks <- GRanges("chr1", IRanges(300, 340))
  x <- quiet(getCounts(bam, peaks, paired = TRUE))
  expect_equal(unname(dense_counts(x)[1, 1]), 1)
})

test_that("getCounts does not count a fragment that only spans a peak", {
  bam <- make_bam(pair_lines("a", "chr1", 101, 110))
  peaks <- GRanges("chr1", IRanges(150, 160))
  x <- quiet(getCounts(bam, peaks, paired = TRUE))
  expect_equal(unname(dense_counts(x)[1, 1]), 0)
})

test_that("getCounts counts duplicate fragments separately", {
  bam <- make_bam(c(pair_lines("a", "chr1", 11, 30),
                    pair_lines("b", "chr1", 11, 30),
                    pair_lines("c", "chr1", 11, 30)))
  x <- quiet(getCounts(bam, GRanges("chr1", IRanges(1, 12)), paired = TRUE))
  expect_equal(unname(dense_counts(x)[1, 1]), 3)
  expect_equal(colData(x)$depth, 3)
})

test_that("getCounts single-end mode counts the 5' end of every read", {
  bam <- make_bam(frag_reads())
  x <- quiet(getCounts(bam, hand_peaks, paired = FALSE))
  expect_equal(unname(dense_counts(x)[, 1]), c(2, 2, 0, 2, 0, 1, 0))
  expect_equal(colData(x)$depth, 10)
})

test_that("getCounts single-end mode uses the right end of reverse-strand reads", {
  bam <- make_bam(c(sam_line("fwd", 0, "chr1", 100, 0, 0),
                    sam_line("rev", 16, "chr1", 200, 0, 0)))
  peaks <- GRanges("chr1", IRanges(c(95, 105, 195, 205), c(101, 110, 201, 210)))
  x <- quiet(getCounts(bam, peaks, paired = FALSE))
  expect_equal(unname(dense_counts(x)[, 1]), c(1, 0, 0, 1))
})

test_that("getCounts single-end mode puts the 5' end of a reverse read on its last base", {
	bam <- make_bam(c(sam_line("fwd", 0, "chr1", 100, 0, 0),
	                  sam_line("rev", 16, "chr1", 200, 0, 0)))
	peaks <- GRanges("chr1", IRanges(start = c(99, 100, 101, 200, 209, 210),
	                                 end = c(99, 100, 101, 208, 209, 210)))
	x <- quiet(getCounts(bam, peaks, paired = FALSE))
	expect_equal(unname(dense_counts(x)[, 1]), c(0, 1, 0, 0, 1, 0))
})

test_that("getCounts paired mode counts the first and last base of a fragment", {
	bam <- make_bam(pair_lines("a", "chr1", 11, 30))
	peaks <- GRanges("chr1", IRanges(start = c(10, 11, 12, 39, 40, 41),
	                                 end = c(10, 11, 12, 39, 40, 41)))
	x <- quiet(getCounts(bam, peaks, paired = TRUE))
	expect_equal(unname(dense_counts(x)[, 1]), c(0, 1, 0, 0, 1, 0))
	bam_rg <- make_bam(pair_lines("a", "chr1", 11, 30, rg = "g1"))
	x <- quiet(getCounts(bam_rg, peaks, paired = TRUE, by_rg = TRUE))
	expect_equal(unname(dense_counts(x)[, 1]), c(0, 1, 0, 0, 1, 0))
})

test_that("getCounts paired mode skips reads that are not proper pairs", {
  bam <- make_bam(c(pair_lines("good", "chr1", 11, 30),
                    pair_lines("bad", "chr1", 101, 30, proper = FALSE)))
  peaks <- GRanges("chr1", IRanges(c(1, 100), c(12, 110)))
  paired <- quiet(getCounts(bam, peaks, paired = TRUE))
  expect_equal(unname(dense_counts(paired)[, 1]), c(1, 0))
  expect_equal(colData(paired)$depth, 1)
  single <- quiet(getCounts(bam, peaks, paired = FALSE))
  expect_equal(unname(dense_counts(single)[, 1]), c(1, 1))
  expect_equal(colData(single)$depth, 4)
})

test_that("getCounts gives zero for peaks on chromosomes missing from the bam", {
  bam <- make_bam(pair_lines("a", "chr1", 11, 30))
  peaks <- GRanges(c("chr1", "chrNew"), IRanges(c(1, 1), c(12, 12)))
  x <- quiet(getCounts(bam, peaks, paired = TRUE))
  expect_equal(unname(dense_counts(x)[, 1]), c(1, 0))
})

test_that("getCounts handles non-standard chromosome names in bam and peaks", {
  bam <- make_bam(pair_lines("a", "chrUn_gl000220", 11, 30),
                  seqs = c(chr1 = 1000, chrUn_gl000220 = 1000))
  peaks <- GRanges(c("chr1", "chrUn_gl000220"), IRanges(c(1, 1), c(12, 12)))
  x <- quiet(getCounts(bam, peaks, paired = TRUE))
  expect_equal(unname(dense_counts(x)[, 1]), c(0, 1))
})

test_that("getCounts reports one column per bam file in a multi-bam call", {
  bam1 <- make_bam(frag_reads())
  bam2 <- make_bam(pair_lines("z", "chr1", 301, 30))
  x <- quiet(getCounts(c(bam1, bam2), hand_peaks, paired = TRUE,
                       colData = S4Vectors::DataFrame(celltype = c("A", "B"))))
  expect_equal(dim(x), c(7, 2))
  expect_equal(colnames(x), basename(c(bam1, bam2)))
  expect_equal(unname(dense_counts(x)[, 1]), c(2, 2, 0, 1, 0, 1, 0))
  expect_equal(unname(dense_counts(x)[, 2]), c(0, 0, 0, 1, 0, 0, 0))
  expect_equal(colData(x)$depth, c(5, 1))
  expect_equal(colData(x)$celltype, c("A", "B"))
})

test_that("getCounts without colData has a colData holding only depth", {
  bam <- make_bam(frag_reads())
  x <- quiet(getCounts(bam, hand_peaks, paired = TRUE))
  expect_equal(colnames(colData(x)), "depth")
})

test_that("getCounts agrees with a reference count on the example paired bams", {
  bams <- vapply(c("test_single1.bam", "test_single2.bam", "test_single3.bam"),
                 function(f) system.file("extdata", f, package = "chromVAR"), "")
  fr <- do.call(rbind, lapply(bams, proper_fragments))
  fr <- fr[order(fr$chr, fr$left), ]
  peaks <- peaks_around(fr, every = 3, half = 40)
  x <- quiet(getCounts(unname(bams), peaks, paired = TRUE))
  for (i in seq_along(bams)) {
    f <- proper_fragments(bams[i])
    expect_equal(unname(dense_counts(x)[, i]),
                 ref_fragment_counts(f$chr, f$left, f$right, peaks))
    expect_equal(colData(x)$depth[i], nrow(f))
  }
  expect_gt(sum(dense_counts(x)), 0)
})

# getCounts from bam with RG tags --------------------------------------------------------

rg_assign <- c("g2", "g10", "g10", "g2", "g1")

test_that("getCounts by_rg gives one naturally ordered column per read group", {
  bam <- make_bam(frag_reads(rg_assign))
  x <- quiet(getCounts(bam, hand_peaks, paired = TRUE, by_rg = TRUE))
  expect_equal(colnames(x), c("g1", "g2", "g10"))
  expect_equal(unname(dense_counts(x)[, "g1"]), c(1, 1, 0, 0, 0, 0, 0))
  expect_equal(unname(dense_counts(x)[, "g2"]), c(1, 1, 0, 1, 0, 0, 0))
  expect_equal(unname(dense_counts(x)[, "g10"]), c(0, 0, 0, 0, 0, 1, 0))
  expect_equal(unname(colData(x)$depth), c(1, 2, 2))
  expect_true(is(assays(x)$counts, "Matrix"))
})

test_that("getCounts by_rg single-end counts the 5' end of every read per group", {
  bam <- make_bam(frag_reads(rg_assign))
  x <- quiet(getCounts(bam, hand_peaks, paired = FALSE, by_rg = TRUE))
  expect_equal(colnames(x), c("g1", "g2", "g10"))
  expect_equal(unname(dense_counts(x)[, "g1"]), c(1, 1, 0, 0, 0, 0, 0))
  expect_equal(unname(dense_counts(x)[, "g2"]), c(1, 1, 0, 2, 0, 0, 0))
  expect_equal(unname(dense_counts(x)[, "g10"]), c(0, 0, 0, 0, 0, 1, 0))
  expect_equal(unname(colData(x)$depth), c(2, 4, 4))
})

test_that("getCounts by_rg column sums add up to the counts without read groups", {
  rg_bam <- make_bam(frag_reads(rg_assign))
  plain_bam <- make_bam(frag_reads())
  by_rg <- quiet(getCounts(rg_bam, hand_peaks, paired = TRUE, by_rg = TRUE))
  plain <- quiet(getCounts(plain_bam, hand_peaks, paired = TRUE))
  expect_equal(unname(rowSums(dense_counts(by_rg))),
               unname(dense_counts(plain)[, 1]))
  expect_equal(sum(colData(by_rg)$depth), colData(plain)$depth)
})

test_that("getCounts by_rg stacks the groups of several bams and repeats colData", {
  bam1 <- make_bam(frag_reads(rg_assign))
  bam2 <- make_bam(pair_lines("z", "chr1", 11, 30, rg = "x1"))
  x <- quiet(getCounts(c(bam1, bam2), hand_peaks, paired = TRUE, by_rg = TRUE,
                       colData = S4Vectors::DataFrame(batch = c("A", "B"))))
  expect_equal(colnames(x), c("g1", "g2", "g10", "x1"))
  expect_equal(colData(x)$batch, c("A", "A", "A", "B"))
  expect_equal(unname(colData(x)$depth), c(1, 2, 2, 1))
  expect_equal(unname(dense_counts(x)[, "x1"]), c(1, 1, 0, 0, 0, 0, 0))
})

test_that("getCounts by_rg agrees with a reference count on the example RG bam", {
  bam <- system.file("extdata", "test_RG.bam", package = "chromVAR")
  fr <- proper_fragments(bam)
  tags <- unique(fr$rg)
  tags <- tags[order(as.numeric(sub(".*-", "", tags)))]
  peaks <- peaks_around(fr[order(fr$chr, fr$left), ], every = 2, half = 40)
  x <- quiet(getCounts(bam, peaks, paired = TRUE, by_rg = TRUE))
  expect_equal(colnames(x), tags)
  for (tag in tags) {
    f <- fr[fr$rg == tag, ]
    expect_equal(unname(dense_counts(x)[, tag]),
                 ref_fragment_counts(f$chr, f$left, f$right, peaks))
    expect_equal(unname(colData(x)$depth[colnames(x) == tag]), nrow(f))
  }
  expect_gt(sum(dense_counts(x)), 0)
})

# getCounts from bed ------------------------------------------------------------------------

bed_lines <- c("chr1\t10\t40\t+", "chr1\t10\t40\t-", "chr1\t100\t150\t+",
               "chr1\t100\t150\t-", "chrX\t0\t5\t-", "chrUn_gl000\t0\t20\t+")
bed_peaks <- GRanges(
  c("chr1", "chr1", "chr1", "chr1", "chr2", "chrX", "chrUn_gl000"),
  IRanges(c(1, 35, 101, 145, 1, 1, 1), c(12, 45, 110, 155, 3, 3, 5)))

test_that("getCounts bed single-end counts the 5' end according to strand", {
  f <- write_tmp(bed_lines, ".bed")
  x <- quiet(getCounts(f, bed_peaks, paired = FALSE, format = "bed"))
  expect_is(x, "RangedSummarizedExperiment")
  expect_equal(unname(dense_counts(x)[, 1]), c(1, 1, 1, 1, 0, 0, 1))
  expect_equal(colData(x)$depth, 6)
  expect_equal(colnames(x), basename(f))
})

test_that("getCounts bed paired mode counts both fragment ends", {
  f <- write_tmp(bed_lines, ".bed")
  x <- quiet(getCounts(f, bed_peaks, paired = TRUE, format = "bed"))
  expect_equal(unname(dense_counts(x)[, 1]), c(2, 2, 2, 2, 0, 1, 1))
  expect_equal(colData(x)$depth, 6)
})

test_that("getCounts bed finds the strand column in a six-column bed", {
  four <- write_tmp(bed_lines[1:4], ".bed")
  six <- write_tmp(c("chr1\t10\t40\tn1\t0\t+", "chr1\t10\t40\tn2\t0\t-",
                     "chr1\t100\t150\tn3\t0\t+", "chr1\t100\t150\tn4\t0\t-"),
                   ".bed")
  x4 <- quiet(getCounts(four, bed_peaks, paired = FALSE, format = "bed"))
  x6 <- quiet(getCounts(six, bed_peaks, paired = FALSE, format = "bed"))
  expect_equal(unname(dense_counts(x6)[, 1]), c(1, 1, 1, 1, 0, 0, 0))
  expect_equal(unname(dense_counts(x6)[, 1]), unname(dense_counts(x4)[, 1]))
})

test_that("getCounts bed handles several files and keeps colData", {
  f1 <- write_tmp(bed_lines[1], ".bed")
  f2 <- write_tmp(bed_lines[2:3], ".bed")
  x <- quiet(getCounts(c(f1, f2), bed_peaks, paired = FALSE, format = "bed",
                       colData = S4Vectors::DataFrame(g = c("a", "b"))))
  expect_equal(dim(x), c(7, 2))
  expect_equal(colnames(x), basename(c(f1, f2)))
  expect_equal(unname(dense_counts(x)[, 1]), c(1, 0, 0, 0, 0, 0, 0))
  expect_equal(unname(dense_counts(x)[, 2]), c(0, 1, 1, 0, 0, 0, 0))
  expect_equal(colData(x)$depth, c(1, 2))
  expect_equal(colData(x)$g, c("a", "b"))
  expect_true(is(assays(x)$counts, "Matrix"))
})

test_that("getCounts bed gives the same counts for fragments given in a bam", {
  bam <- make_bam(pair_lines("a", "chr1", 11, 30))
  bed <- write_tmp("chr1\t10\t40\t+", ".bed")
  peaks <- GRanges("chr1", IRanges(c(1, 35, 100), c(12, 45, 110)))
  from_bam <- quiet(getCounts(bam, peaks, paired = TRUE))
  from_bed <- quiet(getCounts(bed, peaks, paired = TRUE, format = "bed"))
  expect_equal(unname(dense_counts(from_bam)[, 1]),
               unname(dense_counts(from_bed)[, 1]))
  expect_equal(colData(from_bam)$depth, colData(from_bed)$depth)
})

# getCounts: invalid input -----------------------------------------------------------------------

test_that("getCounts rejects an unknown format", {
  f <- write_tmp(bed_lines, ".bed")
  expect_error(getCounts(f, bed_peaks, paired = FALSE, format = "sam"),
               "should be one of")
})

test_that("getCounts requires the paired argument", {
  f <- write_tmp(bed_lines, ".bed")
  expect_error(getCounts(f, bed_peaks, format = "bed"), "paired")
  expect_error(getCounts(make_bam(frag_reads()), hand_peaks), "paired")
})

test_that("getCounts errors for a missing alignment file", {
  expect_error(quiet(getCounts(tempfile(fileext = ".bam"), hand_peaks,
                               paired = TRUE)))
  expect_error(quiet(getCounts(tempfile(fileext = ".bed"), bed_peaks,
                               paired = FALSE, format = "bed")))
})

test_that("getCounts errors when colData has the wrong number of rows", {
  bam1 <- make_bam(frag_reads())
  bam2 <- make_bam(pair_lines("z", "chr1", 301, 30))
  expect_error(quiet(getCounts(c(bam1, bam2), hand_peaks, paired = TRUE,
                               colData = S4Vectors::DataFrame(a = 1))))
})

# getSampleDepths -----------------------------------------------------------------------------

test_that("getSampleDepths counts proper pairs once when paired", {
  bam <- make_bam(frag_reads())
  depth <- getSampleDepths(bam, paired = TRUE)
  expect_equal(unname(depth), 5)
  expect_equal(names(depth), basename(bam))
})

test_that("getSampleDepths counts every read when not paired", {
  bam <- make_bam(frag_reads())
  expect_equal(unname(getSampleDepths(bam, paired = FALSE)), 10)
})

test_that("getSampleDepths excludes reads that are not proper pairs when paired", {
  bam <- make_bam(c(pair_lines("good", "chr1", 11, 30),
                    pair_lines("bad", "chr1", 101, 30, proper = FALSE)))
  expect_equal(unname(getSampleDepths(bam, paired = TRUE)), 1)
  expect_equal(unname(getSampleDepths(bam, paired = FALSE)), 4)
})

test_that("getSampleDepths gives one named value per bam file", {
  bam1 <- make_bam(frag_reads())
  bam2 <- make_bam(pair_lines("z", "chr1", 301, 30))
  depth <- getSampleDepths(c(bam1, bam2), paired = TRUE)
  expect_equal(unname(depth), c(5, 1))
  expect_equal(names(depth), basename(c(bam1, bam2)))
})

test_that("getSampleDepths defaults to paired, without read groups", {
  bam <- make_bam(frag_reads())
  expect_equal(getSampleDepths(bam), getSampleDepths(bam, paired = TRUE))
})

test_that("getSampleDepths by_rg tabulates read groups in natural order", {
  bam <- make_bam(frag_reads(rg_assign))
  by_rg <- getSampleDepths(bam, by_rg = TRUE, paired = TRUE)
  expect_equal(names(by_rg), c("g1", "g2", "g10"))
  expect_equal(as.numeric(by_rg), c(1, 2, 2))
  single <- getSampleDepths(bam, by_rg = TRUE, paired = FALSE)
  expect_equal(names(single), c("g1", "g2", "g10"))
  expect_equal(as.numeric(single), c(2, 4, 4))
})

test_that("getSampleDepths by_rg concatenates the groups of several bams", {
  bam1 <- make_bam(frag_reads(rg_assign))
  bam2 <- make_bam(pair_lines("z", "chr1", 11, 30, rg = "x1"))
  depth <- getSampleDepths(c(bam1, bam2), by_rg = TRUE, paired = TRUE)
  expect_equal(names(depth), c("g1", "g2", "g10", "x1"))
  expect_equal(as.numeric(depth), c(1, 2, 2, 1))
})

test_that("getSampleDepths by_rg agrees with the example RG bam", {
  bam <- system.file("extdata", "test_RG.bam", package = "chromVAR")
  fr <- proper_fragments(bam)
  tags <- unique(fr$rg)
  tags <- tags[order(as.numeric(sub(".*-", "", tags)))]
  depth <- getSampleDepths(bam, by_rg = TRUE, paired = TRUE)
  expect_equal(names(depth), tags)
  expect_equal(as.numeric(depth), as.numeric(table(factor(fr$rg, tags))))
  expect_equal(sum(depth), nrow(fr))
})

test_that("getSampleDepths counts the lines of bed files", {
  f1 <- write_tmp(bed_lines, ".bed")
  f2 <- write_tmp(bed_lines[1:2], ".bed")
  depth <- getSampleDepths(c(f1, f2), format = "bed")
  expect_equal(as.numeric(depth), c(6, 2))
  expect_equal(names(depth), basename(c(f1, f2)))
})

test_that("getSampleDepths of bed files works without readr", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  testthat::local_mocked_bindings(is.installed = function(pkg) FALSE,
                                  .package = "chromVAR")
  f <- write_tmp(bed_lines, ".bed")
  expect_equal(as.numeric(getSampleDepths(f, format = "bed")), 6)
})

test_that("getSampleDepths rejects an unknown format", {
  expect_error(getSampleDepths(write_tmp(bed_lines), format = "cram"),
               "should be one of")
})

test_that("getSampleDepths matches the depth column of getCounts", {
  bam <- make_bam(frag_reads())
  x <- quiet(getCounts(bam, hand_peaks, paired = TRUE))
  expect_equal(unname(getSampleDepths(bam, paired = TRUE)), colData(x)$depth)
})
