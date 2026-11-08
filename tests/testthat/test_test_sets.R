context("test sets")

library(SummarizedExperiment)
library(GenomicRanges)
library(IRanges)
library(Matrix)

local_serial_param()

np <- 60
ns <- 8
local_seed(11)
ts_counts <- matrix(rpois(np * ns, 5) + 1, nrow = np)
ts_bias <- runif(np, 0.3, 0.7)

ts_se <- SummarizedExperiment(assays = list(counts = ts_counts),
                              rowData = DataFrame(bias = ts_bias))
ts_rse <- SummarizedExperiment(assays = list(counts = ts_counts),
                               rowRanges = GRanges("chr1",
                                                   IRanges(start = seq_len(np) * 10,
                                                           width = 5),
                                                   bias = ts_bias))

singletons <- as.list(seq_len(np))
names(singletons) <- paste0("s", seq_len(np))

ann_list <- list(a = 1:30, b = 20:50, c = c(5, 9, 33))
ann_sparse <- sparseMatrix(i = unlist(ann_list),
                           j = rep(1:3, lengths(ann_list)),
                           x = TRUE, dims = c(np, 3),
                           dimnames = list(NULL, names(ann_list)))
ann_dense <- as.matrix(ann_sparse)
ann_se <- SummarizedExperiment(assays = list(matches = ann_sparse),
                               colData = DataFrame(name = c("A", "B", "C"),
                                                   row.names = names(ann_list)))

whiten <- function(x) t(forwardsolve(t(chol(cov(x))), t(x)))

neighbor_rank <- function(counts, bias) {
  w <- whiten(cbind(rowSums(counts), bias))
  d <- as.matrix(dist(w))
  t(apply(d, 1, rank))
}

chosen_peak <- function(res) {
  vapply(seq_len(ncol(res)), function(j) which(as.logical(assay(res)[, j])), 0L)
}

test_that("makePermutedSets output has one logical column per annotation",
          {
            res <- makePermutedSets(ts_se, ann_list, window = 4)
            expect_is(res, "SummarizedExperiment")
            expect_equal(assayNames(res), "matches")
            expect_equal(dim(res), c(np, 3))
            expect_equal(colnames(res), paste0("permuted.", names(ann_list)))
            expect_is(assay(res), "lMatrix")
            expect_equal(unname(res$name), paste0("permuted.", names(ann_list)))
            expect_equal(rownames(colData(res)), colnames(res))
          })

test_that("permuted singleton sets pick one of the window nearest neighbors",
          {
            ranks <- neighbor_rank(ts_counts, ts_bias)
            for (window in c(2, 5, 10)) {
              res <- makePermutedSets(ts_se, singletons, window = window)
              expect_equal(unname(colSums(assay(res))), rep(1, np))
              q <- chosen_peak(res)
              expect_false(any(q == seq_len(np)))
              picked_rank <- ranks[cbind(seq_len(np), q)] - 1
              expect_true(all(picked_rank >= 1 & picked_rank <= window))
              expect_gt(max(picked_rank), window / 2)
            }
          })

test_that("default window is 10 nearest neighbors",
          {
            res <- makePermutedSets(ts_se, singletons)
            ranks <- neighbor_rank(ts_counts, ts_bias)
            picked_rank <- ranks[cbind(seq_len(np), chosen_peak(res))] - 1
            expect_true(all(picked_rank <= 10))
            expect_gt(max(picked_rank), 5)
          })

test_that("explicit bias argument replaces the bias stored in the object",
          {
            other_bias <- rev(ts_bias)
            res <- makePermutedSets(ts_se, singletons, bias = other_bias,
                                    window = 3)
            ranks <- neighbor_rank(ts_counts, other_bias)
            picked_rank <- ranks[cbind(seq_len(np), chosen_peak(res))] - 1
            expect_true(all(picked_rank <= 3))
            res2 <- makePermutedSets(ts_counts, singletons, bias = other_bias,
                                     window = 3)
            picked_rank <- ranks[cbind(seq_len(np), chosen_peak(res2))] - 1
            expect_true(all(picked_rank <= 3))
          })

test_that("RangedSummarizedExperiment takes bias from rowRanges",
          {
            res <- makePermutedSets(ts_rse, singletons, window = 3)
            ranks <- neighbor_rank(ts_counts, ts_bias)
            picked_rank <- ranks[cbind(seq_len(np), chosen_peak(res))] - 1
            expect_true(all(picked_rank <= 3))
          })

test_that("permuted sets are drawn from the peaks and never grow",
          {
            objects <- list(se = ts_se, rse = ts_rse, matrix = ts_counts,
                            Matrix = Matrix::Matrix(ts_counts))
            annos <- list(list = ann_list, sparse = ann_sparse,
                          dense = ann_dense, se = ann_se)
            for (o in names(objects)) {
              for (a in names(annos)) {
                if (o %in% c("matrix", "Matrix") && a == "se") next
                if (o %in% c("matrix", "Matrix")) {
                  res <- makePermutedSets(objects[[o]], annos[[a]],
                                          bias = ts_bias, window = 5)
                } else {
                  res <- makePermutedSets(objects[[o]], annos[[a]],
                                          window = 5)
                }
                info <- paste(o, a)
                expect_equal(dim(res), c(np, 3), info = info)
                expect_equal(colnames(res), paste0("permuted.", names(ann_list)),
                             info = info)
                sizes <- unname(colSums(assay(res)))
                expect_true(all(sizes >= 1 & sizes <= lengths(ann_list)),
                            info = info)
              }
            }
          })

test_that("annotation SummarizedExperiment colData is carried with permuted names",
          {
            res <- makePermutedSets(ts_se, ann_se, window = 4)
            expect_equal(unname(res$name), paste0("permuted.", c("A", "B", "C")))
            expect_equal(rownames(colData(res)), paste0("permuted.", names(ann_list)))
            noname <- SummarizedExperiment(assays = list(matches = ann_sparse))
            res2 <- makePermutedSets(ts_se, noname, window = 4)
            expect_equal(unname(res2$name), paste0("permuted.", names(ann_list)))
          })

test_that("numeric annotation matches in a SummarizedExperiment are converted with a warning",
          {
            numeric_se <- SummarizedExperiment(assays = list(matches =
                                                  as.matrix(ann_sparse) * 1))
            expect_warning(res <- makePermutedSets(ts_se, numeric_se, window = 4),
                           "Annotation object matches converted to logical")
            expect_equal(dim(res), c(np, 3))
          })

test_that("makePermutedSets rejects peaks without fragments",
          {
            zero_counts <- ts_counts
            zero_counts[7, ] <- 0
            expect_error(makePermutedSets(zero_counts, ann_list, bias = ts_bias),
                         "All peaks must have at least one fragment")
            zero_se <- SummarizedExperiment(assays = list(counts = zero_counts),
                                            rowData = DataFrame(bias = ts_bias))
            expect_error(makePermutedSets(zero_se, ann_list),
                         "All peaks must have at least one fragment")
          })

test_that("makePermutedSets validates its SummarizedExperiment inputs",
          {
            no_counts <- SummarizedExperiment(assays = list(foo = ts_counts),
                                              rowData = DataFrame(bias = ts_bias))
            expect_error(makePermutedSets(no_counts, ann_list), "counts")
            neg <- SummarizedExperiment(assays = list(counts = -ts_counts),
                                        rowData = DataFrame(bias = ts_bias))
            expect_error(makePermutedSets(neg, ann_list), "counts")
            expect_error(makePermutedSets(ts_se[, 1], ann_list), "ncol")
            no_matches <- SummarizedExperiment(assays = list(foo = ann_dense))
            expect_error(makePermutedSets(ts_se, no_matches), "matches")
          })

test_that("getPermutedData adds named permuted assays and keeps the counts",
          {
            res <- getPermutedData(ts_se, niterations = 3)
            expect_equal(assayNames(res), c("counts", paste0("perm_", 1:3)))
            expect_equal(assay(res, "counts"), ts_counts)
            for (k in 1:3) {
              expect_equal(dim(assay(res, paste0("perm_", k))), dim(ts_counts))
            }
            expect_equal(rowData(res), rowData(ts_se))
            expect_equal(colnames(res), colnames(ts_se))
          })

test_that("getPermutedData default makes 10 permutations",
          {
            res <- getPermutedData(ts_se)
            expect_equal(assayNames(res), c("counts", paste0("perm_", 1:10)))
          })

test_that("getPermutedData replaces each count by a count of a similar peak in the same sample",
          {
            npk <- 300
            nsm <- 6
            local_seed(5)
            a <- sample(2:40, npk, replace = TRUE)
            cm <- matrix(1000 * a, npk, nsm) + seq_len(npk)
            bias <- runif(npk, 0.3, 0.7)
            rse <- SummarizedExperiment(assays = list(counts = cm),
                                        rowRanges = GRanges("chr1",
                                                            IRanges(seq_len(npk) * 10,
                                                                    width = 5),
                                                            bias = bias))
            res <- getPermutedData(rse, niterations = 2)
            expect_equal(rowRanges(res), rowRanges(rse))
            w <- whiten(cbind(log10(rowSums(cm)), bias))
            for (k in 1:2) {
              perm <- as.matrix(assay(res, paste0("perm_", k)))
              source_peak <- perm %% 1000
              expect_true(all(source_peak >= 1 & source_peak <= npk))
              for (j in seq_len(nsm)) {
                expect_equal(perm[, j], cm[source_peak[, j], j])
              }
              dists <- sqrt(rowSums((w - w[source_peak[, 1], ])^2))
              expect_true(all(dists < 1))
            }
          })

test_that("getPermutedData accepts w and bs",
          {
            res <- getPermutedData(ts_se, niterations = 2, w = 0.5, bs = 20)
            expect_equal(assayNames(res), c("counts", "perm_1", "perm_2"))
            for (k in 1:2) {
              perm <- as.matrix(assay(res, paste0("perm_", k)))
              for (j in seq_len(ns)) {
                expect_true(all(perm[, j] %in% ts_counts[, j]))
              }
            }
          })

test_that("getPermutedData rejects peaks without fragments",
          {
            zero_counts <- ts_counts
            zero_counts[3, ] <- 0
            zero_se <- SummarizedExperiment(assays = list(counts = zero_counts),
                                            rowData = DataFrame(bias = ts_bias))
            expect_error(getPermutedData(zero_se, niterations = 1),
                         "All peaks must have at least one fragment")
          })
