context("computeDeviations and computeExpectations, definitions")

library(SummarizedExperiment)
library(Matrix)
library(GenomicRanges)
library(IRanges)

local_serial_param()

make_counts <- function(npeaks = 60, nsamples = 8, seed = 11) {
  local_seed(seed)
  m <- matrix(rpois(npeaks * nsamples, rep(c(3, 8, 15), length.out = npeaks)),
              nrow = npeaks) + 1L
  storage.mode(m) <- "double"
  colnames(m) <- paste0("s", seq_len(nsamples))
  m
}

make_bg <- function(npeaks, niter, seed = 5) {
  local_seed(seed)
  matrix(sample.int(npeaks, npeaks * niter, replace = TRUE), nrow = npeaks)
}

ref_deviations <- function(counts, peaks, bg, expectation, threshold = 1) {
  counts <- as.matrix(counts)
  if (length(peaks) == 0) {
    return(list(z = rep(NA_real_, ncol(counts)),
                dev = rep(NA_real_, ncol(counts)), matches = 0,
                overlap = NA_real_))
  }
  depth <- colSums(counts)
  observed <- colSums(counts[peaks, , drop = FALSE])
  expected <- sum(expectation[peaks]) * depth
  obs_dev <- (observed - expected) / expected
  niter <- ncol(bg)
  samp_dev <- matrix(NA_real_, niter, ncol(counts))
  for (k in seq_len(niter)) {
    b <- bg[peaks, k]
    s <- colSums(counts[b, , drop = FALSE])
    se <- sum(expectation[b]) * depth
    samp_dev[k, ] <- (s - se) / se
  }
  dev <- obs_dev - colMeans(samp_dev)
  z <- dev / apply(samp_dev, 2, sd)
  low <- expected < threshold
  z[low] <- NA
  dev[low] <- NA
  overlap <- mean(vapply(seq_len(niter),
                         function(k) sum(bg[peaks, k] %in% peaks),
                         0)) / length(peaks)
  list(z = z, dev = dev, matches = length(peaks) / nrow(counts),
       overlap = overlap)
}

ref_matrix <- function(counts, anno, bg, expectation) {
  res <- lapply(anno, function(p) ref_deviations(counts, p, bg, expectation))
  list(z = t(vapply(res, function(x) x$z, rep(0, ncol(counts)))),
       dev = t(vapply(res, function(x) x$dev, rep(0, ncol(counts)))),
       matches = vapply(res, function(x) x$matches, 0),
       overlap = vapply(res, function(x) x$overlap, 0))
}

cnt <- make_counts()
bgm <- make_bg(nrow(cnt), 12)
expc <- rowSums(cnt) / sum(cnt)
anno <- list(a = c(1, 5, 9, 20, 33), b = c(2, 3, 4), c = seq(10, 60, by = 5),
             d = 7)
anno_lmat <- sparseMatrix(i = unlist(anno), j = rep(seq_along(anno), lengths(anno)),
                          x = TRUE, dims = c(nrow(cnt), length(anno)),
                          dimnames = list(NULL, names(anno)))
ref <- ref_matrix(cnt, anno, bgm, expc)

check_against_ref <- function(dev, ref, rn = names(anno)) {
  expect_is(dev, "chromVARDeviations")
  expect_equal(dim(dev), c(length(rn), ncol(cnt)))
  expect_equal(rownames(dev), rn)
  expect_equal(unname(assays(dev)$z), unname(ref$z))
  expect_equal(unname(assays(dev)$deviations), unname(ref$dev))
  expect_equal(unname(rowData(dev)$fractionMatches), unname(ref$matches))
  expect_equal(unname(rowData(dev)$fractionBackgroundOverlap),
               unname(ref$overlap))
}

test_that("computeDeviations matches a base-R reference for a hand-worked example",
{
  m <- rbind(c(10, 30), c(20, 20), c(30, 10), c(40, 40))
  bg <- cbind(c(2, 2, 2, 2), c(3, 3, 3, 3))
  dev <- computeDeviations(m, list(first = 1), bg)
  s <- 0.5 / sqrt(2)
  expect_equal(unname(deviations(dev)["first", ]), c(-0.75, 0.75))
  expect_equal(unname(deviationScores(dev)["first", ]),
               c(-0.75 / s, 0.75 / s))
  expect_equal(unname(rowData(dev)$fractionMatches), 0.25)
  expect_equal(unname(rowData(dev)$fractionBackgroundOverlap), 0)
})

test_that("computeDeviations equals the reference definition for every object/annotation signature",
{
  se <- SummarizedExperiment(assays = list(counts = cnt),
                             colData = DataFrame(grp = rep(c("x", "y"), 4),
                                                 row.names = colnames(cnt)))
  anno_se <- SummarizedExperiment(assays = list(matches = anno_lmat),
                                  colData = DataFrame(name = names(anno),
                                                      row.rn = names(anno)))
  check_against_ref(computeDeviations(cnt, anno, bgm), ref)
  check_against_ref(computeDeviations(cnt, anno_lmat, bgm), ref)
  check_against_ref(computeDeviations(cnt, anno_se, bgm), ref)
  check_against_ref(computeDeviations(Matrix::Matrix(cnt), anno, bgm), ref)
  check_against_ref(computeDeviations(se, anno, bgm), ref)
  check_against_ref(computeDeviations(se, anno_lmat, bgm), ref)
  check_against_ref(computeDeviations(se, anno_se, bgm), ref)
  expect_equal(rowData(computeDeviations(se, anno_se, bgm))$name,
               names(anno))
  expect_equal(colData(computeDeviations(se, anno, bgm)), colData(se))
})

test_that("computeDeviations accepts dense logical and numeric 0/1 annotation matrices",
{
  dense <- as.matrix(anno_lmat)
  check_against_ref(computeDeviations(cnt, dense, bgm), ref, rn = names(anno))
  check_against_ref(computeDeviations(cnt, dense * 1, bgm), ref,
                    rn = names(anno))
})

test_that("computeDeviations with a single sample-annotation combination, single peak, and all peaks",
{
  one <- computeDeviations(cnt, list(only = c(4, 8, 15)), bgm)
  r1 <- ref_matrix(cnt, list(only = c(4, 8, 15)), bgm, expc)
  check_against_ref(one, r1, rn = "only")
  single_peak <- computeDeviations(cnt, list(p = 42), bgm)
  check_against_ref(single_peak, ref_matrix(cnt, list(p = 42), bgm, expc),
                    rn = "p")
  all_peaks <- computeDeviations(cnt, list(all = seq_len(nrow(cnt))), bgm)
  check_against_ref(all_peaks,
                    ref_matrix(cnt, list(all = seq_len(nrow(cnt))), bgm, expc),
                    rn = "all")
  expect_equal(unname(rowData(all_peaks)$fractionMatches), 1)
  expect_equal(unname(rowData(all_peaks)$fractionBackgroundOverlap), 1)
})

test_that("computeDeviations with a user-supplied expectation uses it in place of the default",
{
  e2 <- seq_len(nrow(cnt))
  e2 <- e2 / sum(e2)
  dev <- computeDeviations(cnt, anno, bgm, expectation = e2)
  check_against_ref(dev, ref_matrix(cnt, anno, bgm, e2))
  expect_false(isTRUE(all.equal(assays(dev)$z,
                                assays(computeDeviations(cnt, anno, bgm))$z)))
})

test_that("computeDeviations is invariant to the order of the annotations and of the samples",
{
  dev <- computeDeviations(cnt, anno, bgm)
  rev_dev <- computeDeviations(cnt, rev(anno), bgm)
  expect_equal(assays(rev_dev)$z[rev(seq_along(anno)), ], assays(dev)$z)
  perm <- c(3, 1, 8, 2, 7, 4, 6, 5)
  pdev <- computeDeviations(cnt[, perm], anno, bgm)
  expect_equal(assays(pdev)$z, assays(dev)$z[, perm])
  expect_equal(assays(pdev)$deviations, assays(dev)$deviations[, perm])
  expect_equal(colnames(pdev), colnames(cnt)[perm])
})

test_that("computeDeviations equals the reference when peaks are permuted together with background and expectation",
{
  perm <- c(rev(seq_len(30)), 31:60)
  inv <- order(perm)
  cnt_p <- cnt[perm, ]
  bg_p <- matrix(inv[bgm[perm, ]], nrow = nrow(bgm))
  anno_p <- lapply(anno, function(x) inv[x])
  dev_p <- computeDeviations(cnt_p, anno_p, bg_p)
  dev <- computeDeviations(cnt, anno, bgm)
  expect_equal(assays(dev_p)$z, assays(dev)$z)
  expect_equal(assays(dev_p)$deviations, assays(dev)$deviations)
})

test_that("computeDeviations scales observed counts consistently: deviations are (obs - exp)/exp relative to background",
{
  dev <- computeDeviations(cnt, anno, bgm)
  dev2 <- computeDeviations(cnt * 3, anno, bgm)
  expect_equal(assays(dev2)$z, assays(dev)$z)
  expect_equal(assays(dev2)$deviations, assays(dev)$deviations)
})

test_that("computeDeviations returns NA where the expected count is below one",
{
  m <- rbind(c(1, 1, 100), c(0, 1, 100), c(1, 0, 100), c(0, 1, 100),
             c(1, 0, 100), c(0, 0, 100))
  m <- cbind(m, m[, 3] * 0 + c(3, 4, 5, 6, 7, 8))
  bg <- make_bg(nrow(m), 10)
  e <- rowSums(m) / sum(m)
  an <- list(a = c(1, 2), b = c(3, 4, 5))
  dev <- computeDeviations(m, an, bg)
  r <- ref_matrix(m, an, bg, e)
  expect_true(any(is.na(r$z)))
  expect_equal(is.na(assays(dev)$z), is.na(r$z), check.attributes = FALSE)
  expect_equal(is.na(assays(dev)$deviations), is.na(r$dev),
               check.attributes = FALSE)
  expect_equal(unname(assays(dev)$z), unname(r$z))
})

test_that("an expected count of exactly one is kept and below one is NA",
{
  m <- rbind(c(1, 1), c(49, 49), c(50, 50))
  bg <- cbind(c(2, 2, 2), c(3, 3, 3))
  dev <- computeDeviations(m, list(at_one = 1), bg)
  expect_equal(unname(deviations(dev)["at_one", ]), c(0, 0))
  m2 <- rbind(c(1, 0), c(50, 50), c(49, 49))
  dev2 <- computeDeviations(m2, list(below = 1), bg,
                            expectation = c(0.005, 0.5, 0.495))
  expect_true(all(is.na(deviations(dev2))))
})

test_that("single-peak background overlap is the fraction of background draws equal to the peak",
{
  local_seed(1)
  m <- matrix(rpois(10 * 4, 6) + 1, 10)
  bg <- matrix(rep(c(2, 3), length.out = 10 * 8), 10)
  bg[4, ] <- c(4, 4, 4, 5, 6, 7, 8, 9)
  dev <- computeDeviations(m, list(p = 4), bg)
  expect_equal(unname(rowData(dev)$fractionBackgroundOverlap), 3 / 8)
})

test_that("computeDeviations handles an annotation without matching peaks",
{
  an <- list(a = c(1, 2, 3), empty = integer(0), b = 10)
  dev <- computeDeviations(cnt, an, bgm)
  expect_true(all(is.na(assays(dev)$z["empty", ])))
  expect_true(all(is.na(assays(dev)$deviations["empty", ])))
  expect_equal(unname(rowData(dev)$fractionMatches[2]), 0)
  expect_true(is.na(rowData(dev)$fractionBackgroundOverlap[2]))
  expect_equal(unname(assays(dev)$z[c("a", "b"), ]),
               unname(ref_matrix(cnt, an[c("a", "b")], bgm, expc)$z))

  empty_col <- anno_lmat
  empty_col[, 4] <- FALSE
  dev2 <- computeDeviations(cnt, empty_col, bgm)
  expect_true(all(is.na(assays(dev2)$z["d", ])))
  expect_equal(unname(rowData(dev2)$fractionMatches[4]), 0)
})

test_that("computeDeviations without annotations runs on each peak and messages",
{
  small <- cnt[1:10, ]
  bg <- make_bg(10, 6)
  expect_message(dev <- computeDeviations(small, background_peaks = bg),
                 "Annotations not provided")
  expect_equal(dim(dev), c(10, ncol(small)))
  expect_equal(rownames(dev), as.character(1:10))
  ref_single <- ref_matrix(small, as.list(1:10), bg, rowSums(small) / sum(small))
  expect_equal(unname(assays(dev)$z), unname(ref_single$z))
  expect_message(computeDeviations(small, NULL, bg),
                 "Annotations not provided")
  se <- SummarizedExperiment(assays = list(counts = small))
  expect_message(dev_se <- computeDeviations(se, background_peaks = bg),
                 "Annotations not provided")
  expect_equal(assays(dev_se)$z, assays(dev)$z)
})

test_that("computeDeviations names unnamed annotation lists and matrices 1..k and keeps sample names",
{
  dev <- computeDeviations(cnt, unname(anno), bgm)
  expect_equal(rownames(dev), as.character(seq_along(anno)))
  expect_equal(colnames(dev), colnames(cnt))
  nonames <- anno_lmat
  colnames(nonames) <- NULL
  expect_equal(rownames(computeDeviations(cnt, nonames, bgm)),
               as.character(seq_along(anno)))
})

test_that("computeDeviations computes background peaks and expectations by default",
{
  se <- SummarizedExperiment(assays = list(counts = cnt),
                             rowData = DataFrame(bias = seq(0.2, 0.7, length.out = nrow(cnt)) +
                                                   0.1 * sin(seq_len(nrow(cnt)))))
  local_seed(3)
  dev <- computeDeviations(se, anno)
  local_seed(3)
  bg <- getBackgroundPeaks(se)
  check_against_ref(dev, ref_matrix(cnt, anno, bg, expc))
  local_seed(3)
  dev_again <- computeDeviations(se, anno)
  expect_equal(assays(dev_again)$z, assays(dev)$z)
})

test_that("computeDeviations on mini_counts and mini_ix agrees with the reference",
{
  data(mini_counts, package = "chromVAR")
  data(mini_ix, package = "chromVAR")
  local_seed(8)
  bg <- getBackgroundPeaks(mini_counts, niterations = 20)
  dev <- computeDeviations(mini_counts, mini_ix, bg)
  cm <- as.matrix(counts(mini_counts))
  an <- lapply(seq_len(ncol(mini_ix)),
                               function(j) which(as.matrix(assays(mini_ix)[[1]][, j])))
  names(an) <- colnames(mini_ix)
  r <- ref_matrix(cm, an, bg, rowSums(cm) / sum(cm))
  expect_equal(dim(dev), c(3, 35))
  expect_equal(rownames(dev), colnames(mini_ix))
  expect_equal(unname(assays(dev)$z), unname(r$z))
  expect_equal(unname(assays(dev)$deviations), unname(r$dev))
  expect_equal(colData(dev), colData(mini_counts))
})

test_that("computeDeviations converts non-logical annotation matches with a warning",
{
  num_se <- SummarizedExperiment(
    assays = list(matches = Matrix::Matrix(as.matrix(anno_lmat) * 1, sparse = TRUE)))
  expect_warning(dev <- computeDeviations(cnt, num_se, bgm),
                 "converted to logical")
  check_against_ref(dev, ref)
  lgl_se <- SummarizedExperiment(assays = list(matches = anno_lmat))
  expect_warning(computeDeviations(cnt, lgl_se, bgm), NA)
})

test_that("computeDeviations gives the same result under SerialParam and a two-worker MulticoreParam",
{
  skip_on_os("windows")
  skip_if(nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_")) &&
            !identical(Sys.getenv("NOT_CRAN"), "true"),
          "no forked workers inside R CMD check")
  skip_if(parallel::detectCores() < 2, "fewer than two cores")
  serial <- computeDeviations(cnt, anno, bgm)
  local_serial_param(BiocParallel::MulticoreParam(workers = 2,
                                                  progressbar = FALSE))
  multi <- computeDeviations(cnt, anno, bgm)
  expect_equal(assays(multi)$z, assays(serial)$z)
  expect_equal(assays(multi)$deviations, assays(serial)$deviations)
  expect_equal(rowData(multi), rowData(serial))
})

test_that("computeDeviations rejects invalid input",
{
  zero_row <- cnt
  zero_row[5, ] <- 0
  expect_error(computeDeviations(zero_row, anno, bgm,
                                 expectation = rep(1 / 60, 60)),
               "All peaks must have at least one fragment in one sample")
  expect_error(computeDeviations(cnt, anno, bgm[-1, ]),
               "nrow\\(counts_mat\\) == nrow\\(background_peaks\\)")
  expect_error(computeDeviations(cnt, anno, bgm, expectation = expc[-1]),
               "length\\(expectation\\) == nrow\\(counts_mat\\)")
  expect_error(computeDeviations(cnt, list(c(1, 61)), bgm),
               "Annotations are not valid")
  expect_error(computeDeviations(cnt, list(c(0, 2)), bgm),
               "Annotations are not valid")
  expect_error(computeDeviations(cnt, list(), bgm), "Annotations are not valid")
  expect_error(computeDeviations(cnt, list(c(1.5, 2)), bgm))
  se_neg <- SummarizedExperiment(assays = list(counts = cnt - 3))
  expect_error(computeDeviations(se_neg, anno, bgm),
               "min\\(assays\\(object\\)\\$counts\\) >= 0")
  se_one <- SummarizedExperiment(assays = list(counts = cnt[, 1, drop = FALSE]))
  expect_error(computeDeviations(se_one, anno, bgm), "ncol\\(object\\) > 1")
  se_nocount <- SummarizedExperiment(assays = list(other = cnt))
  expect_error(computeDeviations(se_nocount, anno, bgm),
               "\"counts\" %in% assayNames\\(object\\)")
  no_matches <- SummarizedExperiment(assays = list(other = anno_lmat))
  expect_error(computeDeviations(cnt, no_matches, bgm), "matches\" %in% assayNames")
})

test_that("computeExpectations without options is the fraction of all fragments per peak",
{
  e <- computeExpectations(cnt)
  expect_equal(e, rowSums(cnt) / sum(cnt))
  expect_equal(sum(e), 1)
  expect_equal(computeExpectations(Matrix::Matrix(cnt)), e)
  expect_equal(computeExpectations(
    Matrix::Matrix(cnt, sparse = TRUE)), e)
  se <- SummarizedExperiment(assays = list(counts = cnt))
  expect_equal(computeExpectations(se), e)
  expect_length(e, nrow(cnt))
})

test_that("computeExpectations with norm = TRUE weights samples equally (shape of the expectation)",
{
  depth <- colSums(cnt)
  ref_e <- rowSums(sweep(cnt, 2, depth, "/"))
  e <- computeExpectations(cnt, norm = TRUE)
  expect_equal(e / sum(e), ref_e / sum(ref_e))
  expect_equal(computeExpectations(Matrix::Matrix(cnt), norm = TRUE), e)
  se <- SummarizedExperiment(assays = list(counts = cnt))
  expect_equal(computeExpectations(se, norm = TRUE), e)
  skewed <- cnt
  skewed[, 1] <- skewed[, 1] * 50
  expect_equal(computeExpectations(skewed, norm = TRUE) /
                 sum(computeExpectations(skewed, norm = TRUE)),
               rowSums(sweep(skewed, 2, colSums(skewed), "/")) /
                 sum(sweep(skewed, 2, colSums(skewed), "/")))
})

test_that("computeExpectations by group averages the per-group fractions",
{
  grp <- rep(c("g1", "g2", "g3", "g1"), 2)
  ref_e <- rowMeans(vapply(c("g1", "g2", "g3"), function(g) {
    sub <- cnt[, grp == g, drop = FALSE]
    rowSums(sub) / sum(sub)
  }, rep(0, nrow(cnt))))
  expect_equal(computeExpectations(cnt, group = grp), ref_e)
  expect_equal(computeExpectations(cnt, group = factor(grp)), ref_e)
  expect_equal(computeExpectations(Matrix::Matrix(cnt), group = grp), ref_e)
  se <- SummarizedExperiment(assays = list(counts = cnt),
                             colData = DataFrame(celltype = grp))
  expect_equal(computeExpectations(se, group = grp), ref_e)
  expect_equal(computeExpectations(se, group = "celltype"), ref_e)
  expect_equal(sum(computeExpectations(cnt, group = grp)), 1)
  one_group <- computeExpectations(cnt, group = rep("a", ncol(cnt)))
  expect_equal(one_group, computeExpectations(cnt))
})

test_that("computeExpectations by group with norm = TRUE and equal-sized groups has the per-sample-fraction shape",
{
  grp <- rep(c("a", "b"), each = 4)
  fr <- sweep(cnt, 2, colSums(cnt), "/")
  ref_e <- rowMeans(cbind(rowSums(fr[, grp == "a"]), rowSums(fr[, grp == "b"])))
  e <- computeExpectations(cnt, norm = TRUE, group = grp)
  expect_equal(e / sum(e), ref_e / sum(ref_e))
})

test_that("computeExpectations rejects a group of the wrong length",
{
  expect_error(computeExpectations(cnt, group = c(1, 2)),
               "Group must be vector of length of columns of object")
  se <- SummarizedExperiment(assays = list(counts = cnt))
  expect_error(computeExpectations(se, group = c("a", "b")),
               "Group must be vector of length of columns of object")
  expect_error(computeExpectations(se, group = "nocolumn"),
               "Group must be vector of length of columns of object")
  se_neg <- SummarizedExperiment(assays = list(counts = cnt - 3))
  expect_error(computeExpectations(se_neg), "min\\(assays\\(object\\)\\$counts\\) >= 0")
})

test_that("accessors return the assays of a chromVARDeviations object",
{
  dev <- computeDeviations(cnt, anno, bgm)
  expect_identical(deviations(dev), assays(dev)$deviations)
  expect_identical(deviationScores(dev), assays(dev)$z)
  expect_equal(assayNames(dev), c("deviations", "z"))
  data(mini_dev, package = "chromVAR")
  expect_equal(dim(deviations(mini_dev)), c(3, 35))
  expect_equal(dim(deviationScores(mini_dev)), c(3, 35))
})

test_that("rbind stacks chromVARDeviations results and keeps shared row data",
{
  data(mini_dev, package = "chromVAR")
  dd <- rbind(mini_dev, mini_dev)
  expect_is(dd, "chromVARDeviations")
  expect_equal(dim(dd), c(6, 35))
  expect_equal(unname(deviations(dd)), rbind(deviations(mini_dev),
                                              deviations(mini_dev)),
               check.attributes = FALSE)
  expect_equal(unname(deviationScores(dd)[4:6, ]),
               unname(deviationScores(mini_dev)))
  expect_equal(colData(dd), colData(mini_dev))
  expect_equal(colnames(rowData(dd)), colnames(rowData(mini_dev)))

  a <- computeDeviations(cnt, anno[1:2], bgm)
  b <- computeDeviations(cnt, anno[3:4], bgm)
  rowData(b)$extra <- 1
  ab <- rbind(a, b)
  expect_equal(rownames(ab), names(anno))
  expect_false("extra" %in% colnames(rowData(ab)))
  expect_equal(assays(ab)$z, assays(computeDeviations(cnt, anno, bgm))$z)
})

test_that("cbind on chromVARDeviations is an error",
{
  data(mini_dev, package = "chromVAR")
  expect_error(cbind(mini_dev, mini_dev),
               "Can't concatenate chromVARDeviations horizontally")
})

test_that("a SummarizedExperiment with deviations and z assays coerces to chromVARDeviations",
{
  se <- SummarizedExperiment(assays = list(deviations = matrix(0, 2, 3),
                                           z = matrix(1, 2, 3)))
  out <- as(se, "chromVARDeviations")
  expect_is(out, "chromVARDeviations")
  expect_equal(dim(out), c(2, 3))
  bad <- SummarizedExperiment(assays = list(z = matrix(1, 2, 3)))
  expect_error(as(bad, "chromVARDeviations"),
               "assays slot must contain 'deviations' and 'z'")
})
