context("getAnnotations and annotationMatches")

library(GenomicRanges)
library(SummarizedExperiment)
library(Matrix)

peaks3 <- GRanges("chr1", IRanges(start = c(101, 201, 301), width = 50))

lmat <- function(x) as.matrix(annotationMatches(x))

anno_grl <- GRangesList(
  a = GRanges("chr1", IRanges(120, 130)),
  b = GRanges("chr1", IRanges(c(250, 400), c(260, 410))))

expected_grl <- matrix(c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE), nrow = 3,
                       dimnames = list(NULL, c("a", "b")))

test_that("getAnnotations from a GRangesList marks overlapping peaks", {
  res <- getAnnotations(anno_grl, rowRanges = peaks3)
  expect_true(is(res, "RangedSummarizedExperiment"))
  expect_equal(assayNames(res), "annotationMatches")
  expect_equal(lmat(res), expected_grl)
  expect_equal(colnames(res), c("a", "b"))
  expect_equal(dim(res), c(3, 2))
  expect_equal(start(rowRanges(res)), start(peaks3))
  expect_equal(end(rowRanges(res)), end(peaks3))
  expect_true(is(annotationMatches(res), "lMatrix"))
})

test_that("getAnnotations from a GRangesList treats range ends as inclusive", {
  anno <- GRangesList(
    ends_before = GRanges("chr1", IRanges(50, 100)),
    starts_at = GRanges("chr1", IRanges(50, 101)),
    ends_at = GRanges("chr1", IRanges(150, 160)),
    starts_after = GRanges("chr1", IRanges(151, 160)),
    other_chr = GRanges("chr2", IRanges(101, 150)))
  res <- getAnnotations(anno, rowRanges = peaks3)
  expect_equal(unname(lmat(res)[1, ]), c(FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_false(any(lmat(res)[2:3, ]))
})

test_that("getAnnotations from a GRangesList accepts a RangedSummarizedExperiment", {
  rse <- SummarizedExperiment(assays = list(counts = matrix(0, 3, 2)),
                              rowRanges = peaks3)
  res <- getAnnotations(anno_grl, rowRanges = rse)
  expect_equal(lmat(res), expected_grl)
  expect_equal(start(rowRanges(res)), start(peaks3))
})

test_that("getAnnotations from a GRangesList handles one annotation or none", {
  res <- getAnnotations(anno_grl[2], rowRanges = peaks3)
  expect_equal(dim(res), c(3, 1))
  expect_equal(unname(lmat(res)[, 1]), c(FALSE, TRUE, FALSE))
  expect_equal(colnames(res), "b")
  res <- getAnnotations(GRangesList(), rowRanges = peaks3)
  expect_equal(dim(res), c(3, 0))
  res <- getAnnotations(anno_grl[1], rowRanges = peaks3[1])
  expect_equal(dim(res), c(1, 1))
  expect_true(lmat(res)[1, 1])
})

test_that("getAnnotations passes additional arguments to SummarizedExperiment", {
  res <- getAnnotations(anno_grl, rowRanges = peaks3,
                        colData = DataFrame(group = c("x", "y"),
                                            row.names = c("a", "b")))
  expect_equal(as.character(colData(res)$group), c("x", "y"))
  res <- getAnnotations(expected_grl, rowRanges = peaks3,
                        metadata = list(note = "n"))
  expect_equal(metadata(res), list(note = "n"))
})

test_that("getAnnotations from a logical matrix keeps its values and names", {
  res <- getAnnotations(expected_grl, rowRanges = peaks3)
  expect_true(is(res, "RangedSummarizedExperiment"))
  expect_equal(lmat(res), expected_grl)
  expect_equal(colnames(res), c("a", "b"))
  expect_true(is(annotationMatches(res), "lMatrix"))
})

test_that("getAnnotations coerces numeric matrices to logical", {
  num <- matrix(c(1, 0, 0, 0, 2, 0), nrow = 3,
                dimnames = list(NULL, c("a", "b")))
  res <- getAnnotations(num, rowRanges = peaks3)
  expect_true(is(annotationMatches(res), "lMatrix"))
  expect_equal(lmat(res), expected_grl)
})

test_that("getAnnotations accepts sparse Matrix objects", {
  sp <- Matrix::Matrix(expected_grl, sparse = TRUE)
  expect_true(is(sp, "lgCMatrix"))
  res <- getAnnotations(sp, rowRanges = peaks3)
  expect_equal(lmat(res), expected_grl)
  expect_true(is(annotationMatches(res), "lMatrix"))
  dense <- Matrix::Matrix(expected_grl, sparse = FALSE)
  res <- getAnnotations(dense, rowRanges = peaks3)
  expect_equal(lmat(res), expected_grl)
})

test_that("getAnnotations accepts a data.frame of logicals", {
  df <- data.frame(a = c(TRUE, FALSE, FALSE), b = c(FALSE, TRUE, FALSE))
  res <- getAnnotations(df, rowRanges = peaks3)
  expect_equal(lmat(res), expected_grl)
  expect_equal(colnames(res), c("a", "b"))
  expect_equal(dim(res), c(3, 2))
})

test_that("getAnnotations from a list of indices builds a sparse matrix", {
  res <- getAnnotations(list(a = c(1L, 3L), b = 2L), npeaks = 3)
  expect_equal(lmat(res), matrix(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE), 3,
                                 dimnames = list(NULL, c("a", "b"))))
  expect_true(is(annotationMatches(res), "lMatrix"))
  expect_false(is(res, "RangedSummarizedExperiment"))
  res <- getAnnotations(list(c(2, 3), integer(0)), npeaks = 4)
  expect_equal(unname(lmat(res)),
               matrix(c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
                      4))
  expect_null(colnames(res))
})

test_that("getAnnotations from a list takes npeaks from rowRanges or rowData", {
  res <- getAnnotations(list(1, 2), rowRanges = peaks3)
  expect_equal(dim(res), c(3, 2))
  expect_true(is(res, "RangedSummarizedExperiment"))
  res <- getAnnotations(list(1, 2), rowData = DataFrame(x = 1:4))
  expect_equal(dim(res), c(4, 2))
  expect_equal(rowData(res)$x, 1:4)
  res <- getAnnotations(list(1, 2), npeaks = 3, rowRanges = peaks3)
  expect_equal(dim(res), c(3, 2))
})

test_that("getAnnotations from a list needs the number of peaks", {
  expect_error(getAnnotations(list(1, 2)),
               "Must provide npeaks, rowData, or rowRanges")
})

bed_peaks <- GRanges("chr1", IRanges(start = c(100, 101, 110, 111), width = 1))

write_bed <- function(lines) {
  f <- tempfile(fileext = ".bed")
  writeLines(lines, f)
  f
}

test_that("BED files are read as zero-based half-open intervals", {
  f <- write_bed("chr1\t100\t110")
  res <- getAnnotations(f, rowRanges = bed_peaks)
  expect_equal(unname(lmat(res)[, 1]), c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(dim(res), c(4, 1))
})

test_that("each BED file in a vector is one annotation", {
  f1 <- write_bed("chr1\t100\t110")
  f2 <- write_bed(c("chr1\t110\t111", "chr2\t0\t500"))
  res <- getAnnotations(c(f1, f2), rowRanges = bed_peaks)
  expect_equal(unname(lmat(res)),
               cbind(c(FALSE, TRUE, TRUE, FALSE), c(FALSE, FALSE, FALSE, TRUE)))
  expect_equal(dim(res), c(4, 2))
  expect_true(is(res, "RangedSummarizedExperiment"))
})

test_that("a BED file with a name column is split into one annotation per name", {
  f <- write_bed(c("chr1\t100\t101\tB", "chr1\t109\t110\tA",
                   "chr1\t110\t111\tB"))
  res <- getAnnotations(f, rowRanges = bed_peaks, column = 4)
  expect_equal(colnames(res), c("A", "B"))
  expect_equal(lmat(res), cbind(A = c(FALSE, FALSE, TRUE, FALSE),
                                B = c(FALSE, TRUE, FALSE, TRUE)))
})

test_that("a single BED file without column is one annotation", {
  f <- write_bed(c("chr1\t100\t101", "chr1\t110\t111"))
  res <- getAnnotations(f, rowRanges = bed_peaks)
  expect_equal(dim(res), c(4, 1))
  expect_equal(unname(lmat(res)[, 1]), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("BED files are read the same way without readr", {
  f <- write_bed(c("chr1\t100\t101\tB", "chr1\t109\t110\tA",
                   "chr1\t110\t111\tB"))
  f2 <- write_bed("chr1\t100\t110")
  with_readr <- getAnnotations(f, rowRanges = bed_peaks, column = 4)
  skip_if(utils::packageVersion("testthat") < "3.2.0",
          "local_mocked_bindings needs testthat 3.2.0")
  testthat::local_mocked_bindings(is.installed = function(pkg) FALSE,
                                  .package = "chromVAR")
  without_readr <- getAnnotations(f, rowRanges = bed_peaks, column = 4)
  expect_equal(lmat(without_readr), lmat(with_readr))
  multi <- getAnnotations(c(f2, f2), rowRanges = bed_peaks)
  expect_equal(unname(lmat(multi)[, 2]), c(FALSE, TRUE, TRUE, FALSE))
})

annotation_se <- function(name, value = matrix(c(TRUE, FALSE, TRUE, TRUE), 2)) {
  SummarizedExperiment(assays = setNames(list(Matrix::Matrix(value)), name))
}

test_that("annotationMatches reads each supported assay name", {
  value <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2)
  for (nm in c("annotationMatches", "annotation_matches", "motifMatches",
               "matches")) {
    expect_equal(unname(lmat(annotation_se(nm, value))), value)
  }
})

test_that("annotationMatches warns for the deprecated motif_matches assay", {
  value <- matrix(c(TRUE, FALSE, TRUE, TRUE), 2)
  expect_warning(res <- annotationMatches(annotation_se("motif_matches", value)),
                 "motif_matches assay deprecated")
  expect_equal(unname(as.matrix(res)), value)
})

test_that("annotationMatches picks one of several recognised assays by the current precedence (characterisation, not in the docs)", {
  a <- matrix(c(TRUE, FALSE, FALSE, FALSE), 2)
  b <- matrix(c(FALSE, TRUE, TRUE, TRUE), 2)
  se <- SummarizedExperiment(assays = list(matches = Matrix::Matrix(b),
                                           motifMatches = Matrix::Matrix(b),
                                           annotationMatches = Matrix::Matrix(a)))
  expect_equal(unname(lmat(se)), a)
  se <- SummarizedExperiment(assays = list(matches = Matrix::Matrix(b),
                                           motifMatches = Matrix::Matrix(a)))
  expect_equal(unname(lmat(se)), a)
  se <- SummarizedExperiment(assays = list(matches = Matrix::Matrix(b),
                                           annotation_matches = Matrix::Matrix(a)))
  expect_equal(unname(lmat(se)), a)
})

test_that("annotationMatches fails without a matches assay", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1, 2, 2)))
  expect_error(annotationMatches(se), "No appropriately named assay")
})

test_that("annotationMatches<- adds an annotationMatches assay", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1, 2, 2)))
  value <- Matrix::Matrix(matrix(c(TRUE, FALSE, TRUE, TRUE), 2))
  annotationMatches(se) <- value
  expect_equal(assayNames(se), c("counts", "annotationMatches"))
  expect_equal(unname(lmat(se)), as.matrix(value))
})

test_that("annotationMatches<- replaces the assay that is already there", {
  value <- matrix(c(FALSE, TRUE, FALSE, FALSE), 2)
  for (nm in c("annotationMatches", "annotation_matches", "motifMatches",
               "matches")) {
    se <- annotation_se(nm)
    annotationMatches(se) <- Matrix::Matrix(value)
    expect_equal(assayNames(se), nm)
    expect_equal(unname(lmat(se)), value)
  }
  se <- annotation_se("motif_matches")
  expect_warning(annotationMatches(se) <- Matrix::Matrix(value),
                 "motif_matches assay deprecated")
  expect_equal(assayNames(se), "motif_matches")
})

test_that("annotationMatches<- converts non-logical values with a warning", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1, 2, 2)))
  expect_warning(annotationMatches(se) <- matrix(c(1, 0, 1, 1), 2),
                 "converted to logical")
  expect_true(is.logical(as.matrix(annotationMatches(se))))
  expect_equal(unname(lmat(se)), matrix(c(TRUE, FALSE, TRUE, TRUE), 2))
})

test_that("annotationMatches<- accepts logical Matrix values without a warning", {
  se <- SummarizedExperiment(assays = list(counts = matrix(1, 2, 2)))
  expect_silent(annotationMatches(se) <- Matrix::Matrix(matrix(c(TRUE, FALSE, TRUE,
                                                         TRUE), 2)))
})

test_that("annotationMatches returns what getAnnotations stored", {
  res <- getAnnotations(list(a = c(1, 3), b = 2), npeaks = 3)
  copy <- res
  annotationMatches(copy) <- annotationMatches(res)
  expect_equal(lmat(copy), lmat(res))
})
