context("filtering")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(SummarizedExperiment)
})

# Helpers ---------------------------------------------------------------------

# Counts matrix whose column sums are `fps`; all counts sit in peak 1 and 2.
sample_se <- function(depth, fps, nrow = 4) {
  m <- matrix(0, nrow, length(depth))
  m[1, ] <- fps - floor(fps / 2)
  m[2, ] <- floor(fps / 2)
  colnames(m) <- paste0("s", seq_along(depth))
  SummarizedExperiment(
    assays = list(counts = m),
    rowRanges = GRanges("chr1", IRanges(seq(1, by = 100, length.out = nrow),
                                        width = 50)),
    colData = DataFrame(depth = depth))
}

# fragments in peaks per sample: 500, 1000, 800, 200, 800
# proportion in peaks:           0.5, 0.5,  0.2, 0.5,  0.1
filter_se <- function() {
  sample_se(depth = c(1000, 2000, 4000, 400, 8000),
            fps = c(500, 1000, 800, 200, 800))
}

peak_se <- function(chr, start, end, row_counts, strand = "*") {
  m <- cbind(row_counts, 0)
  SummarizedExperiment(
    assays = list(counts = m),
    rowRanges = GRanges(chr, IRanges(start, end), strand = strand))
}

# filterSamples --------------------------------------------------------------------

test_that("filterSamples estimates both cutoffs from the data and reports them", {
  se <- filter_se()
  msgs <- capture_messages(
    kept <- filterSamples(se, shiny = FALSE, ix_return = TRUE))
  expect_match(msgs, "min_in_peaks set to 0\\.25", all = FALSE)
  expect_match(msgs, "min_depth set to 500", all = FALSE)
  expect_equal(kept, c(1L, 2L))
})

test_that("filterSamples sets min_depth to 10% of the median depth when that exceeds 500", {
  se <- sample_se(depth = c(10000, 20000, 40000, 4000, 80000),
                  fps = c(5000, 10000, 8000, 2000, 8000))
  kept <- suppressMessages(filterSamples(se, shiny = FALSE, ix_return = TRUE))
  expect_message(filterSamples(se, shiny = FALSE, ix_return = TRUE),
                 "min_depth set to 2000")
  expect_equal(kept, c(1L, 2L, 4L))
})

test_that("filterSamples rounds the estimated min_in_peaks to three digits", {
  se <- sample_se(depth = c(10000, 10000, 10000), fps = c(2468, 2468, 2468))
  expect_message(filterSamples(se, shiny = FALSE, ix_return = TRUE),
                 "min_in_peaks set to 0\\.123")
})

test_that("filterSamples gives no message when both cutoffs are given", {
  expect_message(
    filterSamples(filter_se(), min_in_peaks = 0.3, min_depth = 100,
                  shiny = FALSE, ix_return = TRUE), NA)
})

test_that("filterSamples reports only the cutoff that was estimated", {
  se <- filter_se()
  msgs <- character()
  withCallingHandlers(
    filterSamples(se, min_in_peaks = 0.3, shiny = FALSE, ix_return = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_equal(length(msgs), 1)
  expect_match(msgs, "min_depth set to 500")
})

test_that("filterSamples keeps samples at or above both cutoffs", {
  se <- filter_se()
  kept <- filterSamples(se, min_in_peaks = 0.3, min_depth = 100, shiny = FALSE,
                        ix_return = TRUE)
  expect_equal(kept, c(1L, 2L, 4L))
  kept <- filterSamples(se, min_in_peaks = 0.05, min_depth = 3000,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(kept, c(3L, 5L))
  kept <- filterSamples(se, min_in_peaks = 0.15, min_depth = 3000,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(kept, 3L)
  kept <- filterSamples(se, min_in_peaks = 0.3, min_depth = 1500,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(kept, 2L)
})

test_that("filterSamples cutoffs are inclusive", {
  se <- filter_se()
  kept <- filterSamples(se, min_in_peaks = 0.5, min_depth = 1000,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(kept, c(1L, 2L))
  kept <- filterSamples(se, min_in_peaks = 0.5, min_depth = 1001,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(kept, 2L)
  kept <- filterSamples(se, min_in_peaks = 0.51, min_depth = 0,
                        shiny = FALSE, ix_return = TRUE)
  expect_equal(length(kept), 0)
})

test_that("filterSamples keeps everything at zero cutoffs and nothing at huge ones", {
  se <- filter_se()
  expect_equal(filterSamples(se, min_in_peaks = 0, min_depth = 0,
                             shiny = FALSE, ix_return = TRUE), 1:5)
  none <- filterSamples(se, min_in_peaks = 0, min_depth = 1e9, shiny = FALSE)
  expect_equal(ncol(none), 0)
  expect_equal(nrow(none), nrow(se))
})

test_that("filterSamples returns the subsetted object unless ix_return is set", {
  se <- filter_se()
  out <- filterSamples(se, min_in_peaks = 0.3, min_depth = 100, shiny = FALSE)
  expect_is(out, "RangedSummarizedExperiment")
  expect_equal(colnames(out), c("s1", "s2", "s4"))
  expect_equal(colData(out)$depth, c(1000, 2000, 400))
  expect_equal(as.matrix(assays(out)$counts),
               as.matrix(assays(se)$counts)[, c(1, 2, 4)])
  expect_equal(start(rowRanges(out)), start(rowRanges(se)))
})

test_that("filterSamples ix_return gives increasing sample indices", {
  se <- filter_se()
  kept <- filterSamples(se, min_in_peaks = 0.05, min_depth = 0,
                        shiny = FALSE, ix_return = TRUE)
  expect_false(is.unsorted(kept))
  expect_true(all(kept %in% seq_len(ncol(se))))
})

test_that("filterSamples accepts a plain matrix and sparse counts", {
  se <- filter_se()
  sparse_se <- se
  assays(sparse_se)$counts <- Matrix::Matrix(assays(se)$counts, sparse = TRUE)
  expect_equal(filterSamples(sparse_se, min_in_peaks = 0.3, min_depth = 100,
                             shiny = FALSE, ix_return = TRUE), c(1L, 2L, 4L))
  out <- filterSamples(se, min_in_peaks = 0.3, min_depth = 100, shiny = FALSE)
  expect_true(is(assays(out)$counts, "Matrix"))
})

test_that("filterSamples with shiny = TRUE falls back to cutoffs outside interactive sessions", {
  skip_if(interactive())
  se <- filter_se()
  expect_equal(
    suppressMessages(filterSamples(se, shiny = TRUE, ix_return = TRUE)),
    c(1L, 2L))
  expect_message(filterSamples(se, shiny = TRUE, ix_return = TRUE),
                 "min_depth set to 500")
})

test_that("filterSamples requires a depth column in colData", {
  se <- filter_se()
  colData(se)$depth <- NULL
  expect_error(filterSamples(se, shiny = FALSE),
               "colData for object must have column named depth")
})

test_that("filterSamples validates the counts object", {
  expect_error(filterSamples(1:3, shiny = FALSE), "SummarizedExperiment")
  se <- filter_se()
  assayNames(se) <- "other"
  expect_error(filterSamples(se, shiny = FALSE), "counts")
  se <- filter_se()
  assays(se)$counts[1, 1] <- -1
  expect_error(filterSamples(se, shiny = FALSE), ">= 0")
  expect_error(filterSamples(filter_se()[, 1], shiny = FALSE), "ncol")
})

# filterSamplesPlot ----------------------------------------------------------------------

plot_se <- function() filter_se()

test_that("filterSamplesPlot returns a ggplot", {
  p <- suppressWarnings(filterSamplesPlot(plot_se(), min_in_peaks = 0.3,
                                          min_depth = 100, use_plotly = FALSE))
  expect_is(p, "ggplot")
})

test_that("filterSamplesPlot plots depth against the fraction of fragments in peaks", {
  p <- suppressWarnings(filterSamplesPlot(plot_se(), min_in_peaks = 0.3,
                                          min_depth = 100, use_plotly = FALSE))
  expect_equal(p$data$x, c(1000, 2000, 4000, 400, 8000))
  expect_equal(p$data$y, c(0.5, 0.5, 0.2, 0.5, 0.1))
  expect_equal(as.character(p$data$name), paste0("s", 1:5))
})

test_that("filterSamplesPlot flags the samples that pass the cutoffs", {
  se <- plot_se()
  p <- suppressWarnings(filterSamplesPlot(se, min_in_peaks = 0.3,
                                          min_depth = 100, use_plotly = FALSE))
  expect_equal(p$data$z, c(TRUE, TRUE, FALSE, TRUE, FALSE))
  kept <- filterSamples(se, min_in_peaks = 0.3, min_depth = 100, shiny = FALSE,
                        ix_return = TRUE)
  expect_equal(which(p$data$z), kept)
  p <- suppressWarnings(filterSamplesPlot(se, min_in_peaks = 0.5,
                                          min_depth = 1000, use_plotly = FALSE))
  expect_equal(p$data$z, c(TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("filterSamplesPlot draws the cutoffs as reference lines", {
  p <- suppressWarnings(filterSamplesPlot(plot_se(), min_in_peaks = 0.3,
                                          min_depth = 1234, use_plotly = FALSE))
  h <- unlist(lapply(p$layers, function(l) l$data$yintercept))
  v <- unlist(lapply(p$layers, function(l) l$data$xintercept))
  expect_equal(h, 0.3)
  expect_equal(v, 1234)
})

test_that("filterSamplesPlot labels its axes", {
  p <- suppressWarnings(filterSamplesPlot(plot_se(), min_in_peaks = 0.3,
                                          min_depth = 100, use_plotly = FALSE))
  expect_equal(p$labels$x, "Number of fragments")
  expect_equal(p$labels$y, "Proportion of fragments in peaks")
})

test_that("filterSamplesPlot estimates missing cutoffs like filterSamples", {
  se <- plot_se()
  msgs <- capture_messages(
    p <- suppressWarnings(filterSamplesPlot(se, use_plotly = FALSE)))
  expect_match(msgs, "min_in_peaks set to 0\\.25", all = FALSE)
  expect_match(msgs, "min_depth set to 500", all = FALSE)
  expect_equal(which(p$data$z), c(1L, 2L))
  h <- unlist(lapply(p$layers, function(l) l$data$yintercept))
  v <- unlist(lapply(p$layers, function(l) l$data$xintercept))
  expect_equal(h, 0.25)
  expect_equal(v, 500)
})

test_that("filterSamplesPlot ignores use_plotly outside interactive sessions", {
  skip_if(interactive())
  p <- suppressMessages(suppressWarnings(
    filterSamplesPlot(plot_se(), use_plotly = TRUE)))
  expect_is(p, "ggplot")
})

test_that("filterSamplesPlot requires a depth column and a valid counts object", {
  se <- plot_se()
  colData(se)$depth <- NULL
  expect_error(filterSamplesPlot(se),
               "colData for object must have column named depth")
  expect_error(filterSamplesPlot(1:3), "SummarizedExperiment")
  expect_error(filterSamplesPlot(plot_se()[, 1]), "ncol")
})

# filterPeaks ----------------------------------------------------------------------------

test_that("filterPeaks drops peaks below min_fragments_per_peak (inclusive)", {
  se <- peak_se("chr1", c(1, 100, 200, 300), c(50, 150, 250, 350),
                row_counts = c(0, 1, 2, 3))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2:4)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 2, ix_return = TRUE),
               3:4)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 3, ix_return = TRUE),
               4L)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 4, ix_return = TRUE),
               integer(0))
  expect_equal(filterPeaks(se, min_fragments_per_peak = 0, ix_return = TRUE),
               1:4)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 100, ix_return = TRUE),
               integer(0))
})

test_that("filterPeaks sums fragments over all samples", {
  m <- cbind(c(1, 0, 0), c(0, 1, 0), c(0, 1, 0))
  se <- SummarizedExperiment(
    assays = list(counts = m),
    rowRanges = GRanges("chr1", IRanges(c(1, 100, 200), width = 50)))
  expect_equal(filterPeaks(se, min_fragments_per_peak = 2, ix_return = TRUE), 2L)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 1, ix_return = TRUE),
               1:2)
})

test_that("filterPeaks returns the subsetted object unless ix_return is set", {
  se <- peak_se("chr1", c(1, 100, 200), c(50, 150, 250), row_counts = c(0, 3, 4))
  colData(se)$id <- c("a", "b")
  out <- filterPeaks(se)
  expect_is(out, "RangedSummarizedExperiment")
  expect_equal(dim(out), c(2, 2))
  expect_equal(start(rowRanges(out)), c(100, 200))
  expect_equal(unname(as.matrix(assays(out)$counts)), cbind(c(3, 4), 0))
  expect_equal(colData(out)$id, c("a", "b"))
})

test_that("filterPeaks keeps the higher-count peak of an overlapping pair", {
  se <- peak_se("chr1", c(1, 8), c(10, 20), row_counts = c(5, 9))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2L)
  se <- peak_se("chr1", c(1, 8), c(10, 20), row_counts = c(9, 5))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1L)
})

test_that("filterPeaks keeps overlapping peaks when non_overlapping = FALSE", {
  se <- peak_se("chr1", c(1, 8), c(10, 20), row_counts = c(5, 9))
  expect_equal(filterPeaks(se, non_overlapping = FALSE, ix_return = TRUE), 1:2)
})

test_that("filterPeaks resolves chains of overlaps", {
  se <- peak_se(c("chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2"),
                c(1, 8, 18, 50, 1, 5, 12), c(10, 20, 30, 60, 10, 15, 20),
                row_counts = c(5, 9, 3, 1, 2, 7, 4))
  expect_equal(filterPeaks(se, ix_return = TRUE), c(2L, 4L, 6L))
})

test_that("filterPeaks lets a long peak lose to or beat several short peaks", {
  se <- peak_se("chr1", c(1, 10, 30), c(100, 20, 40), row_counts = c(9, 5, 6))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1L)
  se <- peak_se("chr1", c(1, 10, 30), c(100, 20, 40), row_counts = c(2, 5, 6))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2:3)
})

test_that("filterPeaks treats peaks sharing one base as overlapping but not adjacent ones", {
  se <- peak_se("chr1", c(1, 10), c(10, 20), row_counts = c(5, 9))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2L)
  se <- peak_se("chr1", c(1, 11), c(10, 20), row_counts = c(5, 9))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1:2)
})

test_that("filterPeaks does not compare peaks on different chromosomes", {
  se <- peak_se(c("chr1", "chr2"), c(1, 1), c(10, 10), row_counts = c(5, 9))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1:2)
})

test_that("filterPeaks ignores strand when looking for overlaps", {
  se <- peak_se("chr1", c(1, 8), c(10, 20), row_counts = c(5, 9),
                strand = c("+", "-"))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2L)
})

test_that("filterPeaks removes zero-count peaks before resolving overlaps", {
  se <- peak_se("chr1", c(1, 8, 18), c(10, 20, 30), row_counts = c(0, 9, 3))
  expect_equal(filterPeaks(se, min_fragments_per_peak = 1, ix_return = TRUE),
               2L)
  se <- peak_se("chr1", c(1, 8, 100), c(10, 20, 110), row_counts = c(0, 9, 3))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2:3)
})

test_that("filterPeaks of a single peak or no passing peak is well defined", {
  se <- peak_se("chr1", c(1, 100), c(10, 110), row_counts = c(0, 3))
  expect_equal(filterPeaks(se, ix_return = TRUE), 2L)
  none <- filterPeaks(se, min_fragments_per_peak = 10)
  expect_equal(nrow(none), 0)
  expect_equal(filterPeaks(se, min_fragments_per_peak = 10,
                           non_overlapping = FALSE, ix_return = TRUE),
               integer(0))
})

test_that("filterPeaks leaves disjoint but unsorted peaks alone", {
  se <- peak_se(c("chr2", "chr1"), c(1, 1), c(10, 10), row_counts = c(4, 5))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1:2)
})

test_that("filterPeaks refuses to reduce overlapping peaks that are not sorted", {
  se <- peak_se("chr1", c(25, 20), c(35, 30), row_counts = c(4, 5))
  expect_error(filterPeaks(se), "Peaks must be sorted")
  expect_equal(filterPeaks(se, non_overlapping = FALSE, ix_return = TRUE), 1:2)
})

test_that("filterPeaks validates the counts object", {
  expect_error(filterPeaks(1:3), "SummarizedExperiment")
  se <- peak_se("chr1", 1, 10, row_counts = 3)
  expect_error(filterPeaks(se[, 1]), "ncol")
  assays(se)$counts[1, 1] <- -1
  expect_error(filterPeaks(se), ">= 0")
})

test_that("filterPeaks keeps the best peak of every overlap cluster and returns disjoint peaks", {
  local_seed(5)
  n <- 150
  chr <- rep(c("chr1", "chr2", "chr10"), each = n / 3)
  start <- unlist(lapply(1:3, function(i) sort(sample(1:3000, n / 3))))
  peaks <- GRanges(chr, IRanges(start, width = 60))
  peaks <- sortSeqlevels(peaks)
  peaks <- sort(peaks)
  m <- cbind(sample(1:1000, n), 0)
  se <- SummarizedExperiment(assays = list(counts = m), rowRanges = peaks)
  keep <- filterPeaks(se, ix_return = TRUE)
  expect_true(isDisjoint(peaks[keep]))
  expect_false(is.unsorted(keep))
  cluster <- S4Vectors::subjectHits(findOverlaps(peaks,
                                                  reduce(peaks, min.gapwidth = 0L)))
  best <- vapply(split(seq_len(n), cluster),
                 function(i) i[which.max(m[i, 1])], 0L)
  expect_true(all(best %in% keep))
  expect_true(length(keep) < n)
})

test_that("filterPeaks without overlaps in the data keeps every peak", {
  se <- peak_se("chr1", c(1, 100, 200), c(50, 150, 250), row_counts = c(1, 2, 3))
  expect_equal(filterPeaks(se, ix_return = TRUE), 1:3)
  expect_equal(dim(filterPeaks(se)), c(3, 2))
})

# bias_skew and upper_bias_limit_helper -------------------------------------------------------

test_that("bias_skew is zero when counts follow the expectation exactly", {
  a <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  b <- c(1, 2, 3, 4)
  peaks <- GRanges("chr1", IRanges(seq(1, 1000, by = 100), width = 5))
  peaks$bias <- a / 10
  se <- SummarizedExperiment(assays = list(counts = outer(a, b)),
                             rowRanges = peaks,
                             colData = DataFrame(depth = 100 * b))
  out <- chromVAR:::bias_skew(se, nbins = 2)
  expect_equal(dim(out), c(2, 4))
  expect_equal(unname(out), matrix(0, 2, 4), tolerance = 1e-8)
  by_counts <- chromVAR:::bias_skew(se, nbins = 5, what = "counts")
  expect_equal(dim(by_counts), c(5, 4))
  expect_equal(unname(by_counts), matrix(0, 5, 4), tolerance = 1e-8)
})

test_that("bias_skew requires every peak to have a fragment", {
  m <- outer(1:6, 1:3)
  m[3, ] <- 0
  se <- SummarizedExperiment(
    assays = list(counts = m),
    rowRanges = GRanges("chr1", IRanges(seq(1, 600, by = 100), width = 5)),
    colData = DataFrame(depth = c(10, 20, 30)))
  expect_error(chromVAR:::bias_skew(se, nbins = 2, what = "counts"),
               "All peaks must have at least one fragment")
})

test_that("bias_skew requires one expectation value per peak", {
  se <- SummarizedExperiment(
    assays = list(counts = outer(1:6, 1:3)),
    rowRanges = GRanges("chr1", IRanges(seq(1, 600, by = 100), width = 5)),
    colData = DataFrame(depth = c(10, 20, 30)))
  expect_error(chromVAR:::bias_skew(se, nbins = 2, expectation = 1:3,
                                    what = "counts"),
               "length\\(expectation\\) == nrow\\(object\\)")
})

test_that("upper_bias_limit_helper is the third quartile plus k interquartile ranges", {
  x <- 1:9
  expect_equal(unname(chromVAR:::upper_bias_limit_helper(x, 1.5)), 13)
  expect_equal(unname(chromVAR:::upper_bias_limit_helper(x, 0)), 7)
  expect_equal(unname(chromVAR:::upper_bias_limit_helper(c(x, NA), 1.5)), 13)
  expect_equal(unname(chromVAR:::upper_bias_limit_helper(rep(2, 10), 3)), 2)
})
