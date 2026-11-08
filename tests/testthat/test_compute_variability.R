context("computeVariability")

library(SummarizedExperiment)
library(Matrix)

local_serial_param()

make_dev <- function(z, rownames = paste0("set", seq_len(nrow(z))),
                     name = NULL) {
  rownames(z) <- rownames
  se <- SummarizedExperiment(assays = list(deviations = z, z = z))
  if (!is.null(name)) rowData(se)$name <- name
  new("chromVARDeviations", se)
}

local_seed(31)
zmat <- matrix(rnorm(6 * 12, sd = rep(c(0.5, 1, 2), 2)), nrow = 6)
dev <- make_dev(zmat)

ref_p <- function(z) {
  n <- ncol(z)
  s <- apply(z, 1, sd)
  pchisq((n - 1) * s^2, df = n - 1, lower.tail = FALSE)
}

test_that("computeVariability without bootstrap equals sd, chi-square p-value and BH adjustment",
{
  out <- computeVariability(dev, bootstrap_error = FALSE)
  s <- apply(zmat, 1, sd)
  expect_is(out, "data.frame")
  expect_equal(colnames(out), c("name", "variability", "p_value", "p_value_adj"))
  expect_equal(out$variability, s)
  expect_equal(out$p_value, ref_p(zmat))
  expect_equal(out$p_value_adj, p.adjust(ref_p(zmat), method = "BH"))
  expect_equal(rownames(out), rownames(dev))
  expect_equal(out$name, rownames(dev))
  expect_equal(nrow(out), 6)
})

test_that("p-values are small for large variability and near one for small variability",
{
  z <- rbind(rep(c(-3, 3), 10), rep(c(-0.01, 0.01), 10))
  out <- computeVariability(make_dev(z), bootstrap_error = FALSE)
  expect_lt(out$p_value[1], 1e-10)
  expect_gt(out$p_value[2], 1 - 1e-10)
  expect_equal(out$variability, c(sd(z[1, ]), sd(z[2, ])))
})

test_that("computeVariability uses the name row-data column when present",
{
  d <- make_dev(zmat, name = paste0("TF", 1:6))
  out <- computeVariability(d, bootstrap_error = FALSE)
  expect_equal(out$name, paste0("TF", 1:6))
  expect_equal(rownames(out), rownames(d))
})

test_that("computeVariability bootstrap output has the documented columns and ordered bounds",
{
  out <- computeVariability(dev, bootstrap_samples = 300)
  expect_equal(colnames(out),
               c("name", "variability", "bootstrap_lower_bound",
                 "bootstrap_upper_bound", "p_value", "p_value_adj"))
  expect_equal(out$variability, apply(zmat, 1, sd))
  expect_equal(out$p_value, ref_p(zmat))
  expect_true(all(out$bootstrap_lower_bound >= 0))
  expect_true(all(out$bootstrap_lower_bound <= out$bootstrap_upper_bound))
  expect_equal(rownames(out), rownames(dev))
})

test_that("bootstrap bounds for two samples are the two possible resampled standard deviations",
{
  z <- rbind(c(0, 2), c(1, 5), c(-4, -1))
  d <- make_dev(z)
  s <- c(2, 4, 3) / sqrt(2)
  out <- computeVariability(d, bootstrap_samples = 1000)
  expect_equal(out$variability, s)
  expect_equal(out$bootstrap_lower_bound, rep(0, 3))
  expect_equal(out$bootstrap_upper_bound, s)
  out2 <- computeVariability(d, bootstrap_samples = 600,
                             bootstrap_quantiles = c(0.1, 0.9))
  expect_equal(out2$bootstrap_lower_bound, rep(0, 3))
  expect_equal(out2$bootstrap_upper_bound, s)
})

test_that("bootstrap interval brackets the variability for mini_dev",
{
  data(mini_dev, package = "chromVAR")
  out <- computeVariability(mini_dev)
  z <- assays(mini_dev)$z
  expect_equal(unname(out$variability), unname(apply(z, 1, sd)))
  expect_equal(unname(out$p_value), unname(ref_p(z)))
  expect_equal(rownames(out), rownames(mini_dev))
  expect_equal(unname(out$name), unname(rowData(mini_dev)$name))
  expect_true(all(out$bootstrap_lower_bound <= out$variability))
  expect_true(all(out$bootstrap_upper_bound >= out$variability))
  expect_true(all(out$p_value_adj >= out$p_value))
})

test_that("a constant z row has zero variability and a degenerate bootstrap interval",
{
  z <- rbind(rep(1.5, 10), c(-2, 2, 1, -1, 0, 3, -3, 1, 2, -2))
  out <- computeVariability(make_dev(z), bootstrap_samples = 100)
  expect_equal(out$variability[1], 0)
  expect_equal(out$bootstrap_lower_bound[1], 0)
  expect_equal(out$bootstrap_upper_bound[1], 0)
  expect_equal(out$p_value[1], 1)
})

test_that("variability is unchanged by sample order and by row order",
{
  base <- computeVariability(dev, bootstrap_error = FALSE)
  perm <- c(5, 3, 12, 1, 8, 2, 10, 4, 11, 6, 9, 7)
  out <- computeVariability(make_dev(zmat[, perm]), bootstrap_error = FALSE)
  expect_equal(out$variability, base$variability)
  rp <- 6:1
  out_r <- computeVariability(make_dev(zmat[rp, ], rownames = rownames(dev)[rp]),
                              bootstrap_error = FALSE)
  expect_equal(out_r$variability, base$variability[rp])
  expect_equal(rownames(out_r), rownames(dev)[rp])
})

test_that("computeVariability handles NA z-scores according to na.rm",
{
  z <- zmat
  z[2, 3] <- NA
  d <- make_dev(z)
  keep <- computeVariability(d, bootstrap_error = FALSE, na.rm = TRUE)
  expect_equal(keep$variability[2], sd(z[2, ], na.rm = TRUE))
  expect_equal(keep$variability[-2], apply(zmat, 1, sd)[-2])
  drop <- computeVariability(d, bootstrap_error = FALSE, na.rm = FALSE)
  expect_true(is.na(drop$variability[2]))
  expect_equal(drop$variability[-2], apply(zmat, 1, sd)[-2])
  boot_keep <- computeVariability(d, bootstrap_samples = 50, na.rm = TRUE)
  expect_false(is.na(boot_keep$bootstrap_lower_bound[2]))
  boot_drop <- computeVariability(d, bootstrap_samples = 50, na.rm = FALSE)
  expect_true(is.na(boot_drop$bootstrap_lower_bound[2]))
  expect_true(is.na(boot_drop$bootstrap_upper_bound[2]))
  expect_false(anyNA(boot_drop$bootstrap_lower_bound[-2]))
})

test_that("a row with at most one finite z-score has undefined variability under na.rm",
{
  z <- rbind(c(1, NA, NA, NA), c(NA, NA, NA, NA), c(1, 2, 3, 4))
  out <- computeVariability(make_dev(z), bootstrap_error = FALSE)
  expect_true(is.na(out$variability[1]))
  expect_true(is.na(out$variability[2]))
  expect_equal(out$variability[3], sd(1:4))
})

test_that("computeVariability accepts a SummarizedExperiment with deviations and z assays",
{
  se <- SummarizedExperiment(assays = list(deviations = zmat, z = zmat))
  rownames(se) <- rownames(dev)
  out <- computeVariability(se, bootstrap_error = FALSE)
  expect_equal(out$variability, apply(zmat, 1, sd))
})

test_that("computeVariability on computeDeviations output agrees with the sd of the z-scores",
{
  local_seed(2)
  cm <- matrix(rpois(50 * 10, 6) + 1, 50)
  bg <- matrix(sample.int(50, 50 * 15, replace = TRUE), 50)
  d <- computeDeviations(cm, list(a = 1:5, b = 6:20, c = 21:22), bg)
  out <- computeVariability(d, bootstrap_error = FALSE)
  expect_equal(unname(out$variability), unname(apply(assays(d)$z, 1, sd)))
  expect_equal(rownames(out), c("a", "b", "c"))
})

test_that("computeVariability rejects invalid input",
{
  expect_error(computeVariability(matrix(1, 2, 2)),
               "inherits\\(object, \"chromVARDeviations\"\\)")
  expect_error(computeVariability(dev, bootstrap_samples = 0),
               "bootstrap_samples > 0")
  expect_error(computeVariability(dev, bootstrap_samples = -5),
               "bootstrap_samples > 0")
  expect_error(computeVariability(dev, bootstrap_samples = 2.5),
               "all_whole\\(bootstrap_samples\\)")
  expect_error(computeVariability(dev, bootstrap_samples = 5,
                                  bootstrap_quantiles = c(0, 0.9)),
               "bootstrap_quantiles > 0")
  expect_error(computeVariability(dev, bootstrap_samples = 5,
                                  bootstrap_quantiles = c(0.1, 1)),
               "bootstrap_quantiles < 1")
  expect_error(computeVariability(dev, bootstrap_samples = 5,
                                  bootstrap_quantiles = c(0.9, 0.1)),
               "bootstrap_quantiles\\[2\\] > bootstrap_quantiles\\[1\\]")
  expect_silent(computeVariability(dev, bootstrap_error = FALSE,
                                   bootstrap_samples = -1))
})

test_that("row_sds computes the sample standard deviation of each row",
{
  x <- matrix(c(1, 2, 4, 10, 10, 10, -1, 0, 1), nrow = 3, byrow = TRUE)
  expect_equal(chromVAR:::row_sds(x), apply(x, 1, sd))
  expect_equal(chromVAR:::row_sds(x, TRUE), apply(x, 1, sd))
  expect_equal(chromVAR:::row_sds(matrix(1:6, nrow = 3)), rep(sqrt(2) * 1.5, 3))
  y <- matrix(c(1, NA, 3, 5, NA, 6, 7, 8, 9, 12), nrow = 2, byrow = TRUE)
  expect_equal(chromVAR:::row_sds(y, TRUE)[1], sd(c(1, 3, 5)))
  expect_equal(chromVAR:::row_sds(y, TRUE)[2], sd(c(6, 7, 8, 9, 12)))
  expect_true(is.na(chromVAR:::row_sds(y, FALSE)[1]))
  expect_equal(chromVAR:::row_sds(y, FALSE)[2], sd(c(6, 7, 8, 9, 12)))
  expect_true(is.na(chromVAR:::row_sds(matrix(c(1, NA), nrow = 1), TRUE)))
  expect_equal(length(chromVAR:::row_sds(matrix(1:4, nrow = 2), FALSE)), 2)
  big <- matrix(c(1e8 + 1, 1e8 + 2, 1e8 + 3), nrow = 1)
  expect_equal(chromVAR:::row_sds(big), 1)
})

test_that("row_sds_perm bootstraps columns jointly for all rows",
{
  x <- rbind(1:9, (1:9) * 2, rep(5, 9))
  local_seed(1)
  s <- chromVAR:::row_sds_perm(x)
  expect_equal(length(s), 3)
  expect_equal(s[2], 2 * s[1])
  expect_equal(s[3], 0)
  two <- rbind(c(0, 3), c(1, 1))
  sds <- vapply(seq_len(2000), function(i) chromVAR:::row_sds_perm(two),
                c(0, 0))
  expect_true(all(abs(sds[1, ]) < 1e-12 | abs(sds[1, ] - 3 / sqrt(2)) < 1e-12))
  expect_true(all(sds[2, ] == 0))
  expect_equal(mean(sds[1, ] > 1e-12), 0.5, tolerance = 0.15, scale = 1)
})

test_that("row_sds_perm gives the average of resampled sds near the sd for a large sample",
{
  local_seed(3)
  x <- matrix(rnorm(2 * 400, sd = c(1, 3)), nrow = 2)
  boot <- t(vapply(seq_len(300), function(i) chromVAR:::row_sds_perm(x),
                   c(0, 0)))
  expect_equal(colMeans(boot), apply(x, 1, sd), tolerance = 0.03)
})

test_that("row_sds_perm honours na_rm",
{
  x <- rbind(c(1, NA, 3, 4, 8, 9), c(1, 2, 3, 4, 8, 9))
  local_seed(5)
  with_na <- chromVAR:::row_sds_perm(x, na_rm = FALSE)
  expect_equal(length(with_na), 2)
  local_seed(5)
  s_rm <- chromVAR:::row_sds_perm(x, na_rm = TRUE)
  expect_equal(length(s_rm), 2)
  expect_false(is.na(s_rm[2]))
})

test_that("compute_variability_single equals the sd of the z-scores of the deviations",
{
  local_seed(4)
  cm <- matrix(rpois(40 * 9, 7) + 1, 40)
  bg <- matrix(sample.int(40, 40 * 10, replace = TRUE), 40)
  e <- rowSums(cm) / sum(cm)
  v <- chromVAR:::compute_variability_single(c(2, 3, 4, 5), cm, bg, e)
  d <- computeDeviations(cm, list(x = c(2, 3, 4, 5)), bg)
  expect_equal(v, sd(assays(d)$z["x", ], na.rm = TRUE))
})
