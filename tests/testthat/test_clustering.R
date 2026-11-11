context("clustering")

make_dev <- function(dev, z = dev * 10) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z))
  new("chromVARDeviations", se)
}

# Rows r1 and r3 are orthogonal and centered; r2 = r1 / 2 is perfectly
# correlated with r1 (and has the smaller z sd); r4 = r1 + r3 has correlation
# sqrt(0.6) with r1 and sqrt(0.5) with r3, below the 0.9 removal cutoff;
# r5 = -r1 is perfectly anti-correlated.
r1 <- c(1, -1, 1, -1, 1, -1)
r3 <- c(1, 1, -1, -1, 0, 0)
samples <- paste0("s", 1:6)
dev_mat <- rbind(r1 = r1, r2 = r1 / 2, r3 = r3, r4 = r1 + r3, r5 = -r1)
colnames(dev_mat) <- samples
obj <- make_dev(dev_mat)

test_that("getSampleCorrelation removes the lower-variability duplicate row and correlates samples", {
  res <- getSampleCorrelation(obj, threshold = 0)
  expect_equal(res, cor(dev_mat[c("r1", "r3", "r4", "r5"), ]))
  expect_equal(dim(res), c(6, 6))
  expect_equal(dimnames(res), list(samples, samples))
})

test_that("getSampleCorrelation keeps the higher-variability member of a correlated pair", {
  z <- dev_mat * 10
  z["r2", ] <- z["r2", ] * 4
  res <- getSampleCorrelation(make_dev(dev_mat, z), threshold = 0)
  expect_equal(res, cor(dev_mat[c("r2", "r3", "r4", "r5"), ]))
})

test_that("getSampleCorrelation is a symmetric correlation matrix with unit diagonal", {
  res <- getSampleCorrelation(obj, threshold = 0)
  expect_equal(res, t(res))
  expect_equal(unname(diag(res)), rep(1, 6))
  expect_true(all(abs(res) <= 1 + 1e-12))
})

test_that("getSampleCorrelation threshold filters on the standard deviation of z", {
  z_sd <- apply(dev_mat * 10, 1, sd)
  cutoff <- (z_sd["r3"] + z_sd["r1"]) / 2
  expect_true(z_sd["r3"] < cutoff && z_sd["r1"] > cutoff && z_sd["r2"] < cutoff)
  res <- getSampleCorrelation(obj, threshold = cutoff)
  expect_equal(res, cor(dev_mat[c("r1", "r4", "r5"), ]))
})

test_that("getSampleCorrelation and getSampleDistance keep rows whose z sd equals the threshold", {
	z_sd <- chromVAR:::row_sds(deviationScores(obj), na_rm = TRUE)
	cutoff <- z_sd[3]
	expect_equal(getSampleCorrelation(obj, threshold = cutoff),
	             cor(dev_mat[c("r1", "r3", "r4", "r5"), ]))
	expect_equal(as.matrix(getSampleDistance(obj, threshold = cutoff)),
	             as.matrix(dist(t(dev_mat[c("r1", "r3", "r4", "r5"), ]))))
})

test_that("getSampleCorrelation uses z for filtering and deviations for the correlation", {
  z <- dev_mat * 10
  z["r4", ] <- z["r4", ] / 1e4
  res <- getSampleCorrelation(make_dev(dev_mat, z), threshold = 5)
  expect_equal(res, cor(dev_mat[c("r1", "r3", "r5"), ]))
})

test_that("getSampleCorrelation keeps anti-correlated rows", {
  res <- getSampleCorrelation(make_dev(dev_mat[c("r1", "r5"), ]), threshold = 0)
  expect_equal(res, cor(dev_mat[c("r1", "r5"), ]))
})

test_that("getSampleCorrelation uses pairwise complete observations", {
  local_seed(3)
  dev_na <- matrix(rnorm(8 * 6), nrow = 8,
                   dimnames = list(paste0("r", 1:8), samples))
  dev_na[3, 2] <- NA
  row_cor <- cor(t(dev_na), use = "pairwise.complete.obs")
  diag(row_cor) <- 0
  expect_true(all(row_cor < 0.9))
  res <- getSampleCorrelation(make_dev(dev_na, matrix(rnorm(48), 8)),
                              threshold = 0)
  expect_equal(res, cor(dev_na, use = "pairwise.complete.obs"))
  expect_false(anyNA(res))
})

test_that("getSampleCorrelation works on the mini_dev example", {
  data(mini_dev, package = "chromVAR")
  res <- getSampleCorrelation(mini_dev, threshold = 0.8)
  devs <- SummarizedExperiment::assays(mini_dev)$deviations
  sds <- apply(SummarizedExperiment::assays(mini_dev)$z, 1, sd)
  keep <- which(sds >= 0.8)
  expect_lt(cor(devs[keep[1], ], devs[keep[2], ]), 0.9)
  expect_equal(dim(res), c(35, 35))
  expect_equal(dimnames(res), list(colnames(mini_dev), colnames(mini_dev)))
  expect_equal(res, cor(devs[keep, ], use = "pairwise.complete.obs"))
})

test_that("getSampleCorrelation rejects objects that are not chromVARDeviations", {
  expect_error(getSampleCorrelation(dev_mat), "chromVARDeviations")
  expect_error(getSampleCorrelation(list(1, 2)), "chromVARDeviations")
})

# getSampleDistance -------------------------------------------------------------

kept <- t(dev_mat[c("r1", "r3", "r4", "r5"), ])

expect_dist_equal <- function(object, expected) {
  expect_is(object, "dist")
  expect_equal(as.matrix(object), as.matrix(expected))
}

test_that("getSampleDistance defaults to Euclidean distances between samples", {
  res <- getSampleDistance(obj, threshold = 0)
  expect_is(res, "dist")
  expect_dist_equal(res,
                   dist(kept))
  expect_equal(attr(res, "Labels"), samples)
  expect_equal(attr(res, "Size"), 6L)
})

test_that("getSampleDistance with PCA keeping every non-zero component preserves distances", {
  expect_dist_equal(getSampleDistance(obj, threshold = 0, initial_dims = 2),
                   dist(kept))
  expect_dist_equal(getSampleDistance(obj, threshold = 0, initial_dims = 3),
                   dist(kept))
})

test_that("getSampleDistance with one PCA dimension matches projection on the first singular vector", {
  res <- getSampleDistance(obj, threshold = 0, initial_dims = 1)
  sv <- svd(scale(kept, center = TRUE, scale = FALSE))
  score <- sv$u[, 1] * sv$d[1]
  expected <- dist(matrix(score, ncol = 1, dimnames = list(samples, NULL)))
  expect_equal(as.numeric(res), as.numeric(expected))
  expect_true(all(as.numeric(res) <= as.numeric(dist(kept)) + 1e-12))
  expect_equal(attr(res, "Labels"), samples)
})

test_that("getSampleDistance applies a custom distance_function", {
  manhattan <- function(m) dist(m, method = "manhattan")
  expect_dist_equal(getSampleDistance(obj, threshold = 0, distance_function = manhattan),
                   dist(kept, method = "manhattan"))
  cosine <- function(m) as.dist(1 - cor(t(m)))
  expect_dist_equal(getSampleDistance(obj, threshold = 0, distance_function = cosine),
                   as.dist(1 - cor(t(kept))))
})

test_that("getSampleDistance threshold selects rows by z standard deviation", {
  z_sd <- apply(dev_mat * 10, 1, sd)
  cutoff <- (z_sd["r3"] + z_sd["r1"]) / 2
  expect_dist_equal(getSampleDistance(obj, threshold = cutoff),
                   dist(t(dev_mat[c("r1", "r4", "r5"), ])))
})

test_that("getSampleDistance works on the mini_dev example", {
  data(mini_dev, package = "chromVAR")
  res <- getSampleDistance(mini_dev, threshold = 0.8)
  devs <- SummarizedExperiment::assays(mini_dev)$deviations
  sds <- apply(SummarizedExperiment::assays(mini_dev)$z, 1, sd)
  expect_equal(attr(res, "Size"), 35L)
  expect_equal(attr(res, "Labels"), colnames(mini_dev))
  expect_dist_equal(res,
                   dist(t(devs[sds >= 0.8, ])))
})

test_that("getSampleDistance rejects invalid input", {
  expect_error(getSampleDistance(dev_mat), "chromVARDeviations")
  expect_error(getSampleDistance(obj, initial_dims = 0), "initial_dims >= 1")
})
