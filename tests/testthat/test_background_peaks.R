context("background peaks, permuted data and GC bias")

library(SummarizedExperiment)
library(Matrix)
library(GenomicRanges)
library(IRanges)
library(Biostrings)

local_serial_param()

make_gradient <- function(npeaks = 300, nsamples = 6) {
  i <- seq_len(npeaks)
  scale <- exp(seq(log(4), log(400), length.out = npeaks))
  depth <- seq(0.8, 1.2, length.out = nsamples)
  m <- round(outer(scale, depth)) + (i %% 3)
  bias <- 0.5 + 0.2 * sin(i * 12.9898)
  list(counts = m, bias = bias)
}

make_clusters <- function() {
  i <- seq_len(60)
  low <- 1 + (i %% 3)
  high <- 60 + 2 * (i %% 5)
  m <- cbind(c(low, high), c(low + 1, high + 2), c(low, high + 1),
             c(low + 2, high), c(low, high))
  bias <- c(0.30 + 0.004 * (i %% 7), 0.65 + 0.004 * (i %% 7))
  list(counts = m, bias = bias, low = 1:60, high = 61:120)
}

g <- make_gradient()
cl <- make_clusters()
intensity <- log10(rowSums(g$counts))

test_that("getBackgroundPeaks returns one row per peak and one column per iteration",
{
  for (ni in c(1, 7, 50)) {
    bg <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = ni)
    expect_equal(dim(bg), c(nrow(g$counts), ni))
  }
  bg <- getBackgroundPeaks(g$counts, bias = g$bias)
  expect_equal(dim(bg), c(nrow(g$counts), 50))
  bg_bs <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 5, bs = 10)
  expect_equal(dim(bg_bs), c(nrow(g$counts), 5))
})

test_that("getBackgroundPeaks returns valid 1-based peak indices",
{
  bg <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 20)
  expect_true(all(bg >= 1 & bg <= nrow(g$counts)))
  expect_equal(as.numeric(bg), round(as.numeric(bg)))
  expect_false(anyNA(bg))
})

test_that("getBackgroundPeaks gives the same result for every input class under the same seed",
{
  rd <- DataFrame(bias = g$bias)
  se <- SummarizedExperiment(assays = list(counts = g$counts), rowData = rd)
  rse <- SummarizedExperiment(
    assays = list(counts = g$counts),
    rowRanges = GRanges("chr1", IRanges(seq_len(nrow(g$counts)) * 1000,
                                        width = 100), bias = g$bias))
  local_seed(4)
  ref <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 10)
  local_seed(4)
  expect_identical(getBackgroundPeaks(Matrix::Matrix(g$counts), bias = g$bias,
                                      niterations = 10), ref)
  local_seed(4)
  expect_identical(getBackgroundPeaks(se, niterations = 10), ref)
  local_seed(4)
  expect_identical(getBackgroundPeaks(rse, niterations = 10), ref)
  local_seed(4)
  expect_identical(getBackgroundPeaks(se, bias = g$bias, niterations = 10), ref)
  se_dense <- SummarizedExperiment(
    assays = list(counts = Matrix::Matrix(g$counts, sparse = TRUE)), rowData = rd)
  local_seed(4)
  expect_identical(getBackgroundPeaks(se_dense, niterations = 10), ref)
})

test_that("getBackgroundPeaks is reproducible under set.seed and varies across seeds",
{
  local_seed(21)
  a <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 10)
  local_seed(21)
  b <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 10)
  local_seed(22)
  d <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 10)
  expect_identical(a, b)
  expect_false(identical(a, d))
})

test_that("background peaks of well separated peak groups come from the same group",
{
  local_seed(1)
  bg <- getBackgroundPeaks(cl$counts, bias = cl$bias, niterations = 30)
  expect_true(all(bg[cl$low, ] %in% cl$low))
  expect_true(all(bg[cl$high, ] %in% cl$high))
})

test_that("background peaks have similar fragment counts and bias to the peak",
{
  local_seed(2)
  bg <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 40)
  n <- nrow(bg)
  d_int_bg <- mean(abs(intensity[bg] - rep(intensity, ncol(bg))))
  d_bias_bg <- mean(abs(g$bias[bg] - rep(g$bias, ncol(bg))))
  rand <- matrix(sample.int(n, n * 40, replace = TRUE), nrow = n)
  d_int_rand <- mean(abs(intensity[rand] - rep(intensity, 40)))
  d_bias_rand <- mean(abs(g$bias[rand] - rep(g$bias, 40)))
  expect_lt(d_int_bg, 0.3 * d_int_rand)
  expect_lt(d_bias_bg, 0.25 * d_bias_rand)
})

test_that("a smaller w gives more similar background peaks than a larger w",
{
  local_seed(3)
  tight <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 20, w = 0.05)
  loose <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 20, w = 5)
  spread <- function(bg) mean(abs(intensity[bg] - rep(intensity, ncol(bg))))
  expect_lt(spread(tight), spread(loose))
})

test_that("background peaks are invariant to permuting peaks, up to relabeling",
{
  perm <- c(seq(300, 1, by = -1))
  local_seed(5)
  bg <- getBackgroundPeaks(g$counts[perm, ], bias = g$bias[perm],
                           niterations = 40)
  int_p <- intensity[perm]
  d_perm <- mean(abs(int_p[bg] - rep(int_p, ncol(bg))))
  local_seed(5)
  bg0 <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 40)
  d_orig <- mean(abs(intensity[bg0] - rep(intensity, ncol(bg0))))
  expect_equal(d_perm, d_orig, tolerance = 0.2)
  expect_true(all(bg >= 1 & bg <= 300))

  cperm <- c(61:120, 1:60)
  local_seed(6)
  bgc <- getBackgroundPeaks(cl$counts[cperm, ], bias = cl$bias[cperm],
                            niterations = 20)
  expect_true(all(bgc[1:60, ] %in% 1:60))
  expect_true(all(bgc[61:120, ] %in% 61:120))
})

test_that("background peaks do not depend on sample order",
{
  local_seed(7)
  a <- getBackgroundPeaks(g$counts, bias = g$bias, niterations = 10)
  local_seed(7)
  b <- getBackgroundPeaks(g$counts[, 6:1], bias = g$bias, niterations = 10)
  expect_identical(a, b)
})

test_that("getBackgroundPeaks rejects invalid input",
{
  expect_error(getBackgroundPeaks(g$counts, bias = g$bias[-1]),
               "length\\(bias\\) == length\\(fragments_per_peak\\)")
  zero <- g$counts
  zero[3, ] <- 0
  expect_error(getBackgroundPeaks(zero, bias = g$bias),
               "All peaks must have at least one fragment in one sample")
  se_nobias <- SummarizedExperiment(assays = list(counts = g$counts))
  expect_error(getBackgroundPeaks(se_nobias),
               "length\\(bias\\) == length\\(fragments_per_peak\\)")
  se_zero <- SummarizedExperiment(assays = list(counts = zero),
                                  rowData = DataFrame(bias = g$bias))
  expect_error(getBackgroundPeaks(se_zero),
               "All peaks must have at least one fragment in one sample")
  se_neg <- SummarizedExperiment(assays = list(counts = g$counts - 100),
                                 rowData = DataFrame(bias = g$bias))
  expect_error(getBackgroundPeaks(se_neg),
               "min\\(assays\\(object\\)\\$counts\\) >= 0")
  se_one <- SummarizedExperiment(assays = list(counts = g$counts[, 1, drop = FALSE]),
                                 rowData = DataFrame(bias = g$bias))
  expect_error(getBackgroundPeaks(se_one), "ncol\\(object\\) > 1")
  se_other <- SummarizedExperiment(assays = list(reads = g$counts),
                                   rowData = DataFrame(bias = g$bias))
  expect_error(getBackgroundPeaks(se_other),
               "\"counts\" %in% assayNames\\(object\\)")
})

test_that("getBackgroundPeaks works on the packaged mini_counts",
{
  data(mini_counts, package = "chromVAR")
  local_seed(1)
  bg <- getBackgroundPeaks(mini_counts, niterations = 25)
  expect_equal(dim(bg), c(1000, 25))
  expect_true(all(bg >= 1 & bg <= 1000))
  fr <- log10(getFragmentsPerPeak(mini_counts))
  b <- rowRanges(mini_counts)$bias
  expect_lt(mean(abs(fr[bg] - rep(fr, 25))),
            0.12 * mean(abs(fr - mean(fr))))
  expect_lt(mean(abs(b[bg] - rep(b, 25))), 0.45 * mean(abs(b - mean(b))))
})

test_that("getPermutedData adds one permuted counts assay per iteration",
{
  se <- SummarizedExperiment(assays = list(counts = cl$counts))
  rowData(se)$bias <- cl$bias
  local_seed(1)
  out <- getPermutedData(se, niterations = 3)
  expect_is(out, "SummarizedExperiment")
  expect_equal(assayNames(out), c("counts", paste0("perm_", 1:3)))
  expect_equal(assays(out)$counts, assays(se)$counts)
  for (k in 1:3) {
    p <- assays(out)[[paste0("perm_", k)]]
    expect_equal(dim(p), dim(cl$counts))
  }
  local_seed(1)
  one <- getPermutedData(se, niterations = 1)
  expect_equal(assayNames(one), c("counts", "perm_1"))
})

test_that("permuted counts replace each count by a count of a similar peak in the same sample",
{
  se <- SummarizedExperiment(assays = list(counts = cl$counts))
  rowData(se)$bias <- cl$bias
  local_seed(2)
  out <- getPermutedData(se, niterations = 2)
  for (k in 1:2) {
    p <- as.matrix(assays(out)[[paste0("perm_", k)]])
    for (j in seq_len(ncol(p))) {
      expect_true(all(p[cl$low, j] %in% cl$counts[cl$low, j]))
      expect_true(all(p[cl$high, j] %in% cl$counts[cl$high, j]))
    }
    expect_true(all(p[cl$low, ] <= 5))
    expect_true(all(p[cl$high, ] >= 60))
  }
})

test_that("getPermutedData is reproducible under a seeded BiocParallel backend and differs across iterations",
{
  se <- SummarizedExperiment(assays = list(counts = g$counts))
  rowData(se)$bias <- g$bias
  local_serial_param(BiocParallel::SerialParam(RNGseed = 9))
  a <- getPermutedData(se, niterations = 2)
  b <- getPermutedData(se, niterations = 2)
  expect_identical(assays(a)$perm_1, assays(b)$perm_1)
  expect_identical(assays(a)$perm_2, assays(b)$perm_2)
  expect_false(identical(assays(a)$perm_1, assays(a)$perm_2))
})

test_that("reorder_columns picks, in each column, the rows named by the matching index column",
{
  mat <- matrix(1:12, nrow = 4)
  ix <- cbind(c(4, 3, 2, 1), c(1, 1, 2, 2), c(3, 3, 3, 3))
  out <- chromVAR:::reorder_columns(mat, ix)
  expect_equal(out, cbind(mat[c(4, 3, 2, 1), 1], mat[c(1, 1, 2, 2), 2],
                          mat[c(3, 3, 3, 3), 3]), check.attributes = FALSE)
  expect_equal(dim(out), c(4, 3))
})

test_that("ProbSampleReplace samples 1-based indices with the requested probabilities",
{
  local_seed(1)
  n <- 40000
  p <- c(0.5, 0.25, 0.25, 0, 0)
  s <- chromVAR:::ProbSampleReplace(n, p)
  expect_equal(length(s), n)
  expect_true(all(s >= 1 & s <= 5))
  freq <- tabulate(s, nbins = 5) / n
  expect_equal(freq, p, tolerance = 0.025, scale = 1)
  expect_equal(freq[4:5], c(0, 0))
  s1 <- chromVAR:::ProbSampleReplace(10, c(0, 1, 0))
  expect_true(all(s1 == 2))
  local_seed(3)
  a <- chromVAR:::ProbSampleReplace(30, rep(0.2, 5))
  local_seed(3)
  b <- chromVAR:::ProbSampleReplace(30, rep(0.2, 5))
  expect_identical(a, b)
})

test_that("ProbSampleReplace with many skewed probabilities matches them",
{
  local_seed(2)
  p <- c(rep(0.001, 100), rep(0.9 / 10, 10))
  p <- p / sum(p)
  n <- 100000
  s <- chromVAR:::ProbSampleReplace(n, p)
  freq <- tabulate(s, nbins = length(p)) / n
  expect_equal(sum(freq[101:110]), sum(p[101:110]), tolerance = 0.01)
  expect_equal(freq[101:110], p[101:110], tolerance = 0.01, scale = 1)
})

test_that("bg_sample_helper draws, for each peak, from peaks weighted by bin probabilities",
{
  membership <- c(0, 0, 1, 1, 2)
  bin_p <- diag(3)
  density <- c(2, 2, 1)
  local_seed(1)
  out <- chromVAR:::bg_sample_helper(membership, bin_p, density, 20)
  expect_equal(dim(out), c(5, 20))
  expect_true(all(out[1:2, ] %in% 1:2))
  expect_true(all(out[3:4, ] %in% 3:4))
  expect_true(all(out[5, ] == 5))
  bin_p2 <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), 3)
  local_seed(1)
  out2 <- chromVAR:::bg_sample_helper(membership, bin_p2, density, 200)
  expect_true(all(out2[1:4, ] %in% 1:4))
  expect_true(all(1:4 %in% out2[1, ]))
  expect_true(all(out2[5, ] == 5))
})

test_that("bg_sample_helper weights each bin by its probability, not by its peak count",
{
  membership <- c(0, 0, 0, 1)
  bin_p <- matrix(1, 2, 2)
  density <- c(3, 1)
  local_seed(1)
  out <- chromVAR:::bg_sample_helper(membership, bin_p, density, 40000)
  freq <- tabulate(out[1, ], nbins = 4) / ncol(out)
  expect_equal(freq, c(1/6, 1/6, 1/6, 1/2), tolerance = 0.025, scale = 1)
})

test_that("euc_dist is the matrix of Euclidean distances between rows",
{
  x <- matrix(c(0, 3, 1, 4, 0, 4, 2, 2, 7), ncol = 3, byrow = TRUE)
  expect_equal(chromVAR:::euc_dist(x), as.matrix(dist(x)),
               check.attributes = FALSE)
  expect_equal(chromVAR:::euc_dist(x[1, , drop = FALSE]), matrix(0, 1, 1))
  y <- cbind(seq(0, 1, length.out = 17), sin(1:17))
  expect_equal(chromVAR:::euc_dist(y), as.matrix(dist(y)),
               check.attributes = FALSE)
  d <- chromVAR:::euc_dist(y)
  expect_equal(d, t(d))
  expect_equal(diag(d), rep(0, 17))
})

dss <- DNAStringSet(c(chr1 = "ACGTACGTNNGGGGCCCCAAAATTTT",
                      chr2 = "GGGGGGGGGGAAAAAAAAAA"))
peaks <- GRanges(c("chr1", "chr1", "chr2"), IRanges(c(1, 9, 5), c(8, 18, 14)))
expected_gc <- c(0.5, 1, 0.6)

test_that("addGCBias adds the GC fraction (A, C, G, T only) for a RangedSummarizedExperiment",
{
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)),
                             rowRanges = peaks)
  out <- addGCBias(se, genome = dss)
  expect_is(out, "RangedSummarizedExperiment")
  expect_equal(rowRanges(out)$bias, expected_gc)
  expect_equal(assays(out)$counts, assays(se)$counts)
  expect_equal(granges(rowRanges(out)), granges(peaks))
})

test_that("addGCBias adds a bias row-data column for a SummarizedExperiment with explicit peaks",
{
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)))
  out <- addGCBias(se, peaks = peaks, genome = dss)
  expect_equal(rowData(out)$bias, expected_gc)
  expect_equal(assays(out)$counts, assays(se)$counts)
})

test_that("addGCBias reads sequences from a FaFile",
{
  fa <- tempfile(fileext = ".fa")
  writeXStringSet(dss, fa)
  Rsamtools::indexFa(fa)
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)),
                             rowRanges = peaks)
  out <- addGCBias(se, genome = Rsamtools::FaFile(fa))
  expect_equal(rowRanges(out)$bias, expected_gc)
  se2 <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)))
  expect_equal(rowData(addGCBias(se2, peaks = peaks,
                                 genome = Rsamtools::FaFile(fa)))$bias,
               expected_gc)
})

test_that("addGCBias output is usable as the bias for getBackgroundPeaks",
{
  local_seed(1)
  n <- 120
  starts <- seq(1, by = 40, length.out = n)
  seqs <- vapply(seq_len(n), function(i) {
    paste(sample(c("A", "C", "G", "T"), 40, replace = TRUE,
                 prob = c(0.35, 0.15, 0.15, 0.35) + (i %% 4) * c(-0.07, 0.07, 0.07, -0.07)),
          collapse = "")
  }, "")
  chr <- DNAStringSet(c(chrA = paste(seqs, collapse = "")))
  gr <- GRanges("chrA", IRanges(starts, width = 40))
  cm <- matrix(rpois(n * 4, 20) + 1, nrow = n)
  se <- SummarizedExperiment(assays = list(counts = cm), rowRanges = gr)
  out <- addGCBias(se, genome = chr)
  gc_ref <- vapply(seqs, function(s) {
    ch <- strsplit(s, "")[[1]]
    mean(ch %in% c("C", "G"))
  }, 0, USE.NAMES = FALSE)
  expect_equal(rowRanges(out)$bias, gc_ref)
  bg <- getBackgroundPeaks(out, niterations = 3)
  expect_equal(dim(bg), c(n, 3))
})

test_that("addGCBias rejects missing or invalid genomes",
{
  se <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)),
                             rowRanges = peaks)
  expect_error(addGCBias(se), "No genome provided")
  se2 <- SummarizedExperiment(assays = list(counts = matrix(1:6, 3)))
  expect_error(addGCBias(se2, peaks = peaks), "No genome provided")
  expect_error(addGCBias(se, genome = 5),
               "genome input must be a BSgenome, DNAStringSet, or FaFile")
})
