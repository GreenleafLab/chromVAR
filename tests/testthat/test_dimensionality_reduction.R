context("dimensionality reduction")

skip_if_not_installed("Rtsne")

make_dev <- function(dev, z) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z))
  new("chromVARDeviations", se)
}

# z rows are standardized so that the row standard deviations equal `scales`
make_object <- function(n_rows, n_cols, scales, seed) {
  local_seed(seed)
  dev <- matrix(rnorm(n_rows * n_cols), n_rows,
                dimnames = list(paste0("r", seq_len(n_rows)),
                                paste0("s", seq_len(n_cols))))
  z <- t(scale(t(matrix(rnorm(n_rows * n_cols), n_rows)))) * scales
  dimnames(z) <- dimnames(dev)
  make_dev(dev, z)
}

scales16 <- rep(c(0.5, 2), each = 8)
obj <- make_object(16, 10, scales16, seed = 11)
dev <- SummarizedExperiment::assays(obj)$deviations
keep_hi <- which(scales16 > 1)

test_that("test data rows are not correlated above the removal cutoff", {
  rc <- cor(t(dev[keep_hi, ]))
  diag(rc) <- 0
  expect_true(all(rc < 0.9))
})

test_that("deviationsTsne for samples runs t-SNE on the variable rows and labels by sample", {
  local_seed(5)
  res <- deviationsTsne(obj, threshold = 1, perplexity = 2, max_iter = 100,
                        what = "samples")
  local_seed(5)
  ref <- Rtsne::Rtsne(t(dev[keep_hi, ]), perplexity = 2, max_iter = 100,
                      theta = 0.5)$Y
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(10, 2))
  expect_equal(rownames(res), colnames(obj))
  expect_equal(unname(res), ref)
  expect_true(all(is.finite(res)))
})

test_that("deviationsTsne for samples drops the lower-variability member of a correlated pair", {
  dev2 <- rbind(dev, r17 = 0.9 * dev["r9", ])
  z2 <- rbind(SummarizedExperiment::assays(obj)$z, r17 = rep(c(1.5, -1.5), 5) *
                1.1)
  obj2 <- make_dev(dev2, z2)
  expect_true(sd(z2["r17", ]) > 1 && sd(z2["r17", ]) < sd(z2["r9", ]))
  local_seed(5)
  res <- deviationsTsne(obj2, threshold = 1, perplexity = 2, max_iter = 100)
  local_seed(5)
  ref <- Rtsne::Rtsne(t(dev[keep_hi, ]), perplexity = 2, max_iter = 100)$Y
  expect_equal(unname(res), ref)
})

test_that("deviationsTsne passes theta and max_iter on to Rtsne", {
  local_seed(9)
  res <- deviationsTsne(obj, threshold = 1, perplexity = 2, max_iter = 50,
                        theta = 0)
  local_seed(9)
  ref <- Rtsne::Rtsne(t(dev[keep_hi, ]), perplexity = 2, max_iter = 50,
                      theta = 0)$Y
  expect_equal(unname(res), ref)
})

test_that("deviationsTsne for samples is reproducible under a fixed seed", {
  local_seed(1)
  a <- deviationsTsne(obj, threshold = 1, perplexity = 2, max_iter = 100)
  local_seed(1)
  b <- deviationsTsne(obj, threshold = 1, perplexity = 2, max_iter = 100)
  expect_equal(a, b)
})

test_that("deviationsTsne lowers perplexity when it is too high for the number of samples", {
  data(mini_dev, package = "chromVAR")
  local_seed(2)
  msgs <- testthat::capture_messages(
    res <- deviationsTsne(mini_dev, threshold = 0.8, max_iter = 100))
  expect_true(any(grepl("Perplexity given too high", msgs)))
  expect_true(any(grepl("Setting perplexity to 11", msgs)))
  devs <- SummarizedExperiment::assays(mini_dev)$deviations
  sds <- apply(SummarizedExperiment::assays(mini_dev)$z, 1, sd)
  local_seed(2)
  ref <- Rtsne::Rtsne(t(devs[sds > 0.8, ]), perplexity = 11, max_iter = 100)$Y
  expect_equal(unname(res), ref)
  expect_equal(rownames(res), colnames(mini_dev))
})

test_that("deviationsTsne for annotations returns a Dim1/Dim2 data.frame keyed by annotation", {
  local_seed(5)
  res <- deviationsTsne(obj, threshold = 1, perplexity = 2, max_iter = 100,
                        what = "annotations")
  local_seed(5)
  ref <- Rtsne::Rtsne(dev[keep_hi, ], perplexity = 2, max_iter = 100,
                      check_duplicates = FALSE)$Y
  expect_is(res, "data.frame")
  expect_equal(colnames(res), c("Dim1", "Dim2"))
  expect_equal(rownames(res), rownames(obj)[keep_hi])
  expect_equal(unname(as.matrix(res)), ref)
})

test_that("deviationsTsne for annotations keeps correlated rows and uses the documented default perplexity", {
  big <- make_object(40, 6, rep(2, 40), seed = 4)
  bigdev <- SummarizedExperiment::assays(big)$deviations
  bigdev[40, ] <- bigdev[39, ] * 2
  SummarizedExperiment::assays(big)$deviations <- bigdev
  local_seed(8)
  res <- deviationsTsne(big, threshold = 1, max_iter = 100, what = "annotations")
  local_seed(8)
  ref <- Rtsne::Rtsne(bigdev, perplexity = 8, max_iter = 100,
                      check_duplicates = FALSE)$Y
  expect_equal(nrow(res), 40)
  expect_equal(unname(as.matrix(res)), ref)
})

test_that("remove_correlated_helper removes rows correlated at or above the cutoff", {
	e1 <- c(1, -1, 1, -1, 1, -1) / sqrt(6)
	e2 <- c(1, 1, -1, -1, 0, 0) / 2
	with_cor <- function(r) rbind(a = e1, b = r * e1 + sqrt(1 - r^2) * e2)
	rem <- chromVAR:::remove_correlated_helper
	expect_equal(rem(with_cor(0.93), c(1, 2)), 2L)
	expect_equal(rem(with_cor(0.93), c(2, 1)), 1L)
	expect_equal(rem(with_cor(0.85), c(1, 2)), 1:2)
	expect_equal(rem(with_cor(0.85), c(1, 2), cutoff = 0.8), 2L)
})

test_that("deviationsTsne keeps rows strictly above the threshold only", {
	vars <- chromVAR:::row_sds(SummarizedExperiment::assays(obj)$z, FALSE)
	thr <- max(vars[scales16 == 0.5])
	local_seed(5)
	res <- deviationsTsne(obj, threshold = thr, perplexity = 2, max_iter = 50,
	                      what = "annotations")
	expect_equal(rownames(res), rownames(obj)[scales16 > 1])
	local_seed(5)
	res_s <- deviationsTsne(obj, threshold = thr, perplexity = 2, max_iter = 50)
	local_seed(5)
	ref_s <- Rtsne::Rtsne(t(dev[keep_hi, ]), perplexity = 2, max_iter = 50)$Y
	expect_equal(unname(res_s), ref_s)
	graded <- make_object(16, 10, seq(0.5, 2, length.out = 16), seed = 3)
	sds <- chromVAR:::row_sds(SummarizedExperiment::assays(graded)$z, FALSE)
	second <- sort(sds, decreasing = TRUE)[2]
	expect_error(deviationsTsne(graded, threshold = second, perplexity = 2),
	             "too few non/NA")
	expect_error(deviationsTsne(graded, threshold = second, perplexity = 2,
	                            what = "annotations"), "too few non/NA")
})

test_that("deviationsTsne rejects thresholds that leave too few rows", {
  expect_error(deviationsTsne(obj, threshold = 5, perplexity = 2),
               "Threshold too high")
  expect_error(deviationsTsne(obj, threshold = 5, perplexity = 2,
                              what = "annotations"), "threshold too high")
  graded <- make_object(16, 10, seq(0.5, 2, length.out = 16), seed = 3)
  one_left <- 1.95
  sds <- apply(SummarizedExperiment::assays(graded)$z, 1, sd)
  expect_equal(sum(sds > one_left), 1)
  expect_error(deviationsTsne(graded, threshold = one_left, perplexity = 2),
               "too few non/NA")
  expect_error(deviationsTsne(graded, threshold = one_left, perplexity = 2,
                              what = "annotations"), "too few non/NA")
})

test_that("deviationsTsne rejects invalid arguments", {
  expect_error(deviationsTsne(obj, what = "genes"), "should be one of")
  expect_error(deviationsTsne(dev, threshold = 1), "chromVARDeviations")
})
