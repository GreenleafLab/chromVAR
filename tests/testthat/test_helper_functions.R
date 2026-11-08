context("fragment count helpers")

library(SummarizedExperiment)
library(Matrix)

m <- matrix(c(1, 0, 3,
              2, 5, 0,
              0, 0, 4,
              7, 1, 1), nrow = 4, byrow = TRUE,
            dimnames = list(paste0("p", 1:4), paste0("s", 1:3)))

se_dense <- SummarizedExperiment(assays = list(counts = m))
se_sparse <- SummarizedExperiment(assays = list(counts = Matrix::Matrix(m, sparse = TRUE)))
se_rse <- SummarizedExperiment(
  assays = list(counts = m),
  rowRanges = GenomicRanges::GRanges("chr1", IRanges::IRanges(1:4 * 100, width = 50)))

test_that("getFragmentsPerPeak sums counts across samples for every input class",
{
  expected <- c(p1 = 4, p2 = 7, p3 = 4, p4 = 9)
  expect_equal(getFragmentsPerPeak(m), expected)
  expect_equal(getFragmentsPerPeak(Matrix::Matrix(m)), expected)
  expect_equal(getFragmentsPerPeak(Matrix::Matrix(m, sparse = TRUE)), expected)
  expect_equal(getFragmentsPerPeak(se_dense), expected)
  expect_equal(getFragmentsPerPeak(se_sparse), expected)
  expect_equal(getFragmentsPerPeak(se_rse), expected)
})

test_that("getFragmentsPerSample sums counts across peaks for every input class",
{
  expected <- c(s1 = 10, s2 = 6, s3 = 8)
  expect_equal(getFragmentsPerSample(m), expected)
  expect_equal(getFragmentsPerSample(Matrix::Matrix(m)), expected)
  expect_equal(getFragmentsPerSample(Matrix::Matrix(m, sparse = TRUE)), expected)
  expect_equal(getFragmentsPerSample(se_dense), expected)
  expect_equal(getFragmentsPerSample(se_sparse), expected)
  expect_equal(getFragmentsPerSample(se_rse), expected)
})

test_that("getTotalFragments sums all counts for every input class",
{
  expect_equal(getTotalFragments(m), 24)
  expect_equal(getTotalFragments(Matrix::Matrix(m)), 24)
  expect_equal(getTotalFragments(Matrix::Matrix(m, sparse = TRUE)), 24)
  expect_equal(getTotalFragments(se_dense), 24)
  expect_equal(getTotalFragments(se_sparse), 24)
  expect_equal(getTotalFragments(se_rse), 24)
})

test_that("the three summaries agree with each other and keep the right lengths",
{
  expect_equal(sum(getFragmentsPerPeak(m)), getTotalFragments(m))
  expect_equal(sum(getFragmentsPerSample(m)), getTotalFragments(m))
  big <- matrix(rpois(30 * 7, 4), 30)
  expect_length(getFragmentsPerPeak(big), 30)
  expect_length(getFragmentsPerSample(big), 7)
  expect_equal(unname(getFragmentsPerPeak(big)), rowSums(big))
  expect_equal(unname(getFragmentsPerSample(big)), colSums(big))
})

test_that("fragment helpers handle single rows, single columns and all-zero data",
{
  expect_equal(unname(getFragmentsPerPeak(m[1, , drop = FALSE])), 4)
  expect_equal(unname(getFragmentsPerSample(m[1, , drop = FALSE])), c(1, 0, 3))
  expect_equal(unname(getFragmentsPerPeak(m[, 1, drop = FALSE])), c(1, 2, 0, 7))
  expect_equal(unname(getFragmentsPerSample(m[, 1, drop = FALSE])), 10)
  z <- matrix(0, 3, 2)
  expect_equal(unname(getFragmentsPerPeak(z)), rep(0, 3))
  expect_equal(unname(getFragmentsPerSample(z)), rep(0, 2))
  expect_equal(getTotalFragments(z), 0)
  expect_equal(unname(getFragmentsPerPeak(Matrix::Matrix(z, sparse = TRUE))), rep(0, 3))
})

test_that("fragment helpers on packaged data match base R sums",
{
  data(mini_counts, package = "chromVAR")
  cm <- as.matrix(assays(mini_counts)$counts)
  expect_equal(unname(getFragmentsPerPeak(mini_counts)), rowSums(cm))
  expect_equal(getFragmentsPerSample(mini_counts), colSums(cm))
  expect_equal(getTotalFragments(mini_counts), sum(cm))
  expect_equal(names(getFragmentsPerSample(mini_counts)), colnames(mini_counts))
  data(example_counts, package = "chromVAR")
  ec <- as.matrix(assays(example_counts)$counts)
  expect_equal(getFragmentsPerSample(example_counts), colSums(ec))
  expect_equal(getTotalFragments(example_counts), sum(ec))
})

test_that("fragment helpers on a SummarizedExperiment require a counts assay",
{
  se <- SummarizedExperiment(assays = list(reads = m))
  expect_error(getFragmentsPerPeak(se), "\"counts\" %in% assayNames")
  expect_error(getFragmentsPerSample(se), "\"counts\" %in% assayNames")
  expect_error(getTotalFragments(se), "\"counts\" %in% assayNames")
})

test_that("quantile_helper returns quantiles of the values when no NA is present",
{
  v <- c(5, 1, 4, 2, 3)
  q <- chromVAR:::quantile_helper(v, c(0.25, 0.5, 0.75), na.rm = FALSE)
  expect_equal(unname(q), c(2, 3, 4))
  expect_equal(unname(chromVAR:::quantile_helper(v, c(0.25, 0.5, 0.75), TRUE)),
               c(2, 3, 4))
  expect_equal(unname(chromVAR:::quantile_helper(v, 0.5, FALSE)), 3)
})

test_that("quantile_helper drops NA with na.rm and returns NA otherwise",
{
  v <- c(5, 1, NA, 4, 2, 3)
  expect_equal(unname(chromVAR:::quantile_helper(v, c(0.25, 0.75), TRUE)),
               c(2, 4))
  out <- chromVAR:::quantile_helper(v, c(0.025, 0.5, 0.975), FALSE)
  expect_length(out, 3)
  expect_true(all(is.na(out)))
  out_nan <- chromVAR:::quantile_helper(c(1, NaN, 2), c(0.1, 0.9), FALSE)
  expect_true(all(is.na(out_nan)))
})
