context("variability synergy")

library(SummarizedExperiment)
library(Matrix)

local_serial_param()

vs_np <- 60
vs_ns <- 12
local_seed(2)
vs_lambda <- matrix(5, vs_np, vs_ns)
vs_lambda[1:10, 1:6] <- 12
vs_lambda[11:20, c(2, 4, 6, 8)] <- 14
vs_lambda[21:30, 7:12] <- 12
vs_counts <- matrix(rpois(vs_np * vs_ns, vs_lambda) + 1, vs_np, vs_ns)
vs_bg <- vapply(1:15, function(k) ((seq_len(vs_np) + k * 2 - 1) %% vs_np) + 1,
                rep(0, vs_np))
vs_exp <- rowSums(vs_counts) / sum(vs_counts)
vs_se <- SummarizedExperiment(assays = list(counts = vs_counts),
                              rowData = DataFrame(bias = runif(vs_np, 0.3, 0.7)))

vs_list <- list(A = 1:10, B = 11:20, C = 5:14)
vs_sparse <- sparseMatrix(i = unlist(vs_list),
                          j = rep(1:3, lengths(vs_list)),
                          x = TRUE, dims = c(vs_np, 3),
                          dimnames = list(NULL, names(vs_list)))
vs_se_ann <- SummarizedExperiment(assays = list(matches = vs_sparse))

# Deviation of a peak set: (observed - expected) / expected per sample, minus the
# mean of the same quantity over the background sets; z divides by their sd.
ref_dev <- function(peak_set, counts, bg, expectation) {
  depth <- colSums(counts)
  dev_of <- function(p) {
    obs <- colSums(counts[p, , drop = FALSE])
    ex <- sum(expectation[p]) * depth
    (obs - ex) / ex
  }
  obs_dev <- dev_of(peak_set)
  bg_dev <- vapply(seq_len(ncol(bg)), function(k) dev_of(bg[peak_set, k]),
                   rep(0, ncol(counts)))
  list(dev = obs_dev - rowMeans(bg_dev),
       z = (obs_dev - rowMeans(bg_dev)) / apply(bg_dev, 1, sd))
}

ref_variability <- function(peak_set) {
  sd(ref_dev(peak_set, vs_counts, vs_bg, vs_exp)$z)
}

ref_correlation <- function(sets) {
  vars <- vapply(sets, ref_variability, 0)
  sets <- sets[order(vars, decreasing = TRUE)]
  n <- length(sets)
  out <- matrix(1, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      s1 <- setdiff(sets[[i]], sets[[j]])
      s2 <- setdiff(sets[[j]], sets[[i]])
      out[i, j] <- out[j, i] <-
        cor(ref_dev(s1, vs_counts, vs_bg, vs_exp)$dev,
            ref_dev(s2, vs_counts, vs_bg, vs_exp)$dev)
    }
  }
  out
}

test_that("getAnnotationCorrelation matches the definition on overlapping sets",
          {
            res <- getAnnotationCorrelation(vs_counts, vs_list,
                                            background_peaks = vs_bg,
                                            expectation = vs_exp)
            expect_equal(res, ref_correlation(vs_list))
          })

test_that("getAnnotationCorrelation is symmetric with unit diagonal and in [-1, 1]",
          {
            res <- getAnnotationCorrelation(vs_counts, vs_list,
                                            background_peaks = vs_bg,
                                            expectation = vs_exp)
            expect_equal(dim(res), c(3, 3))
            expect_equal(res, t(res))
            expect_equal(unname(diag(res)), rep(1, 3))
            expect_true(all(res >= -1 - 1e-12 & res <= 1 + 1e-12))
          })

test_that("getAnnotationCorrelation orders annotations by decreasing variability",
          {
            res <- getAnnotationCorrelation(vs_counts, vs_list,
                                            background_peaks = vs_bg,
                                            expectation = vs_exp)
            vars <- vapply(vs_list, ref_variability, 0)
            expect_equal(rownames(res), names(sort(vars, decreasing = TRUE)))
            expect_equal(colnames(res), rownames(res))
          })

test_that("getAnnotationCorrelation uses supplied variabilities for the order",
          {
            res <- getAnnotationCorrelation(
              vs_counts, vs_list, background_peaks = vs_bg, expectation = vs_exp,
              variabilities = data.frame(variability = c(1, 3, 2)))
            expect_equal(rownames(res), c("B", "C", "A"))
            ref <- ref_correlation(vs_list)
            expect_equal(res, ref[c("B", "C", "A"), c("B", "C", "A")])
          })

test_that("getAnnotationCorrelation gives NA when one set is contained in the other",
          {
            sets <- list(A = 1:10, A2 = 1:10, B = 11:20, S = 3:5)
            res <- getAnnotationCorrelation(vs_counts, sets,
                                            background_peaks = vs_bg,
                                            expectation = vs_exp)
            expect_true(is.na(res["A", "A2"]))
            expect_true(is.na(res["A2", "A"]))
            expect_true(is.na(res["A", "S"]))
            expect_false(is.na(res["A", "B"]))
            expect_false(is.na(res["B", "S"]))
          })

test_that("getAnnotationCorrelation names unnamed annotations by position",
          {
            res <- getAnnotationCorrelation(vs_counts, unname(vs_list),
                                            background_peaks = vs_bg,
                                            expectation = vs_exp)
            ref <- ref_correlation(setNames(unname(vs_list), c("1", "2", "3")))
            expect_equal(res, ref)
          })

test_that("getAnnotationCorrelation default expectation is the fraction of reads per peak",
          {
            explicit <- getAnnotationCorrelation(vs_counts, vs_list,
                                                 background_peaks = vs_bg,
                                                 expectation = vs_exp)
            default <- getAnnotationCorrelation(vs_counts, vs_list,
                                                background_peaks = vs_bg)
            expect_equal(default, explicit)
          })

test_that("getAnnotationCorrelation gives the same answer for every input type",
          {
            ref <- ref_correlation(vs_list)
            objs <- list(matrix = vs_counts, Matrix = Matrix::Matrix(vs_counts),
                         se = vs_se)
            annos <- list(list = vs_list, sparse = vs_sparse,
                          dense = as.matrix(vs_sparse), se = vs_se_ann)
            for (o in names(objs)) {
              for (a in names(annos)) {
                res <- getAnnotationCorrelation(objs[[o]], annos[[a]],
                                                background_peaks = vs_bg,
                                                expectation = vs_exp)
                expect_equal(res, ref, info = paste(o, a))
              }
            }
          })

test_that("getAnnotationCorrelation with default background peaks returns a valid correlation matrix",
          {
            res <- getAnnotationCorrelation(vs_se, vs_list)
            expect_equal(dim(res), c(3, 3))
            expect_equal(sort(rownames(res)), c("A", "B", "C"))
            expect_equal(res, t(res))
            expect_equal(unname(diag(res)), rep(1, 3))
            expect_true(all(abs(res) <= 1 + 1e-12))
          })

test_that("getAnnotationCorrelation rejects invalid inputs",
          {
            f <- function(...) getAnnotationCorrelation(vs_counts, ...,
                                                        background_peaks = vs_bg,
                                                        expectation = vs_exp)
            expect_error(f(list(A = 1:10, B = 0:3)), "Annotations are not valid")
            expect_error(f(list(A = 1:10, B = 55:61)), "Annotations are not valid")
            expect_error(f(list(A = 1:10)), "l >= 2")
            expect_error(f(vs_list, variabilities = 1:3),
                         "Incorrect variabilities input")
            expect_error(f(vs_list, variabilities = data.frame(variability = 1:2)),
                         "length\\(variabilities\\)")
            expect_error(getAnnotationCorrelation(vs_counts, vs_list,
                                                  background_peaks = vs_bg[-1, ],
                                                  expectation = vs_exp),
                         "nrow\\(background_peaks\\)")
            expect_error(getAnnotationCorrelation(vs_counts, vs_list,
                                                  background_peaks = vs_bg,
                                                  expectation = vs_exp[-1]),
                         "length\\(expectation\\)")
            zero <- vs_counts
            zero[3, ] <- 0
            expect_error(getAnnotationCorrelation(zero, vs_list,
                                                  background_peaks = vs_bg,
                                                  expectation = vs_exp),
                         "All peaks must have at least one fragment")
            expect_error(getAnnotationCorrelation(vs_counts,
                                                  matrix("a", vs_np, 2),
                                                  background_peaks = vs_bg,
                                                  expectation = vs_exp))
            no_matches <- SummarizedExperiment(assays = list(foo = as.matrix(vs_sparse)))
            expect_error(getAnnotationCorrelation(vs_counts, no_matches,
                                                  background_peaks = vs_bg,
                                                  expectation = vs_exp),
                         "matches")
            expect_error(getAnnotationCorrelation(vs_se[, 1], vs_list,
                                                  background_peaks = vs_bg),
                         "ncol")
          })

test_that("getAnnotationSynergy returns a symmetric matrix with zero diagonal",
          {
            res <- getAnnotationSynergy(vs_counts, vs_list,
                                        background_peaks = vs_bg,
                                        expectation = vs_exp, nbg = 5)
            expect_equal(dim(res), c(3, 3))
            expect_equal(res, t(res))
            expect_equal(unname(diag(res)), rep(0, 3))
            vars <- vapply(vs_list, ref_variability, 0)
            expect_equal(rownames(res), names(sort(vars, decreasing = TRUE)))
            expect_equal(colnames(res), rownames(res))
          })

test_that("getAnnotationSynergy is NA for sets that share no peaks",
          {
            res <- getAnnotationSynergy(vs_counts, vs_list,
                                        background_peaks = vs_bg,
                                        expectation = vs_exp, nbg = 5)
            expect_true(is.na(res["A", "B"]))
            expect_true(is.na(res["B", "A"]))
            expect_false(is.na(res["A", "C"]))
            expect_false(is.na(res["B", "C"]))
          })

test_that("getAnnotationSynergy orders by supplied variabilities",
          {
            res <- getAnnotationSynergy(
              vs_counts, vs_list, background_peaks = vs_bg, expectation = vs_exp,
              variabilities = data.frame(variability = c(1, 3, 2)), nbg = 5)
            expect_equal(rownames(res), c("B", "C", "A"))
            expect_equal(colnames(res), c("B", "C", "A"))
          })

test_that("getAnnotationSynergy names unnamed annotations by position",
          {
            res <- getAnnotationSynergy(
              vs_counts, unname(vs_list), background_peaks = vs_bg,
              expectation = vs_exp,
              variabilities = data.frame(variability = c(1, 3, 2)), nbg = 5)
            expect_equal(rownames(res), c("2", "3", "1"))
          })

test_that("getAnnotationSynergy is the z-score of the overlap variability against random subsets",
          {
            p <- c(1, 2, 11, 12)
            q <- c(1, 2, 11, 40, 41)
            vals <- vapply(combn(p, 3, simplify = FALSE), ref_variability, 0)
            target <- ref_variability(intersect(q, p))
            mu <- mean(vals)
            sigma <- sqrt(mean((vals - mu)^2))
            expected <- (target - mu) / sigma
            local_seed(1)
            res <- getAnnotationSynergy(
              vs_counts, list(P = p, Q = q), background_peaks = vs_bg,
              expectation = vs_exp, nbg = 300,
              variabilities = data.frame(variability = c(2, 1)))
            expect_equal(rownames(res), c("P", "Q"))
            expect_lt(abs(res["P", "Q"] - expected), 0.4)
            expect_equal(res["P", "Q"], res["Q", "P"])
          })

test_that("getAnnotationSynergy works for every input type",
          {
            objs <- list(matrix = vs_counts, Matrix = Matrix::Matrix(vs_counts),
                         se = vs_se)
            annos <- list(list = vs_list, sparse = vs_sparse,
                          dense = as.matrix(vs_sparse), se = vs_se_ann)
            for (o in names(objs)) {
              for (a in names(annos)) {
                res <- getAnnotationSynergy(objs[[o]], annos[[a]],
                                            background_peaks = vs_bg,
                                            expectation = vs_exp, nbg = 3)
                info <- paste(o, a)
                expect_equal(dim(res), c(3, 3), info = info)
                expect_equal(sort(rownames(res)), c("A", "B", "C"), info = info)
                expect_equal(res, t(res), info = info)
                expect_equal(unname(diag(res)), rep(0, 3), info = info)
              }
            }
          })

test_that("getAnnotationSynergy with default background peaks and expectation runs",
          {
            res <- getAnnotationSynergy(vs_se, vs_list, nbg = 3)
            expect_equal(dim(res), c(3, 3))
            expect_equal(res, t(res))
          })

test_that("getAnnotationSynergy rejects invalid inputs",
          {
            f <- function(...) getAnnotationSynergy(vs_counts, ...,
                                                    background_peaks = vs_bg,
                                                    expectation = vs_exp,
                                                    nbg = 3)
            expect_error(f(list(A = 1:10, B = 0:3)), "Annotations are not valid")
            expect_error(f(list(A = 1:10, B = 55:61)), "Annotations are not valid")
            expect_error(f(list(A = 1:10)), "l >= 2")
            expect_error(f(vs_list, variabilities = 1:3),
                         "Incorrect variabilities input")
            expect_error(f(vs_list, variabilities = data.frame(variability = 1:2)),
                         "length\\(variabilities\\)")
            expect_error(getAnnotationSynergy(vs_counts, vs_list,
                                              background_peaks = vs_bg[-1, ],
                                              expectation = vs_exp),
                         "nrow\\(background_peaks\\)")
            expect_error(getAnnotationSynergy(vs_counts, vs_list,
                                              background_peaks = vs_bg,
                                              expectation = vs_exp[-1]),
                         "length\\(expectation\\)")
            zero <- vs_counts
            zero[3, ] <- 0
            expect_error(getAnnotationSynergy(zero, vs_list,
                                              background_peaks = vs_bg,
                                              expectation = vs_exp),
                         "All peaks must have at least one fragment")
            expect_error(getAnnotationSynergy(vs_counts, matrix("a", vs_np, 2),
                                              background_peaks = vs_bg,
                                              expectation = vs_exp))
          })

test_that("deviationsCovariability divides covariances by the variance of the row",
          {
            z <- rbind(c(1, 2, 3, 4), c(2, 4, 6, 8), c(4, 1, 3, 2))
            se <- SummarizedExperiment(assays = list(z = z))
            res <- deviationsCovariability(se)
            expect_equal(res[1:2, 1:2], matrix(c(1, 0.5, 2, 1), 2))
            expect_equal(dim(res), c(3, 3))
            expect_equal(unname(diag(res)), rep(1, 3))
            for (i in 1:3) {
              for (j in 1:3) {
                expect_equal(res[i, j], cov(z[i, ], z[j, ]) / var(z[i, ]))
              }
            }
          })

test_that("deviationsCovariability keeps the deviation names and scales covariance by row variance",
          {
            data(mini_dev, package = "chromVAR")
            res <- deviationsCovariability(mini_dev)
            z <- assays(mini_dev)$z
            expect_equal(dim(res), c(3, 3))
            expect_equal(rownames(res), rownames(mini_dev))
            expect_equal(colnames(res), rownames(mini_dev))
            expect_equal(res * apply(z, 1, var), t(res * apply(z, 1, var)))
            expect_equal(res[2, 3], cov(z[2, ], z[3, ]) / var(z[2, ]))
          })
