context("utils")

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(SummarizedExperiment)
})

counts_check <- chromVAR:::counts_check
matches_check <- chromVAR:::matches_check
`%ni%` <- chromVAR:::`%ni%`
mean_smooth <- chromVAR:::mean_smooth
tabulate2 <- chromVAR:::tabulate2
mxsort <- chromVAR:::mxsort

# counts_check --------------------------------------------------------------------

counts_se <- function(m = matrix(c(1, 0, 2, 3, 4, 5), 3, 2)) {
  SummarizedExperiment(assays = list(counts = m))
}

test_that("counts_check converts a base matrix of counts to a Matrix", {
  out <- counts_check(counts_se())
  expect_true(is(assays(out)$counts, "Matrix"))
  expect_equal(as.matrix(assays(out)$counts), matrix(c(1, 0, 2, 3, 4, 5), 3, 2))
})

test_that("counts_check leaves sparse counts as they are", {
  se <- counts_se(Matrix::Matrix(matrix(c(1, 0, 2, 0, 4, 0), 3, 2),
                                 sparse = TRUE))
  out <- counts_check(se)
  expect_true(is(assays(out)$counts, "sparseMatrix"))
  expect_equal(as.matrix(assays(out)$counts), matrix(c(1, 0, 2, 0, 4, 0), 3, 2))
})

test_that("counts_check accepts counts of exactly zero", {
  expect_is(counts_check(counts_se(matrix(0, 3, 2))), "SummarizedExperiment")
})

test_that("counts_check rejects objects that are not SummarizedExperiments", {
  expect_error(counts_check(matrix(1:6, 3, 2)), "SummarizedExperiment")
  expect_error(counts_check(NULL), "SummarizedExperiment")
})

test_that("counts_check requires an assay named counts", {
  se <- SummarizedExperiment(assays = list(other = matrix(1:6, 3, 2)))
  expect_error(counts_check(se), "counts")
})

test_that("counts_check rejects negative counts", {
  expect_error(counts_check(counts_se(matrix(c(1, -1, 2, 3, 4, 5), 3, 2))),
               ">= 0")
})

test_that("counts_check requires more than one sample", {
  expect_error(counts_check(counts_se(matrix(1:3, 3, 1))), "ncol")
})

test_that("counts_check rejects counts that cannot be a Matrix", {
  expect_error(counts_check(counts_se(matrix(c("a", "b", "c", "d"), 2, 2))))
})

# matches_check -------------------------------------------------------------------------

test_that("matches_check accepts a logical sparse matrix without a warning", {
  m <- Matrix::Matrix(matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE), 3, 2),
                      sparse = TRUE)
  se <- SummarizedExperiment(assays = list(matches = m))
  expect_warning(out <- matches_check(se), NA)
  expect_equal(as.matrix(annotationMatches(out)), as.matrix(m))
})

test_that("matches_check converts non-logical matches to logical with a warning", {
  m <- matrix(c(1, 0, 1, 1, 0, 0), 3, 2)
  se <- SummarizedExperiment(assays = list(annotationMatches = m))
  expect_warning(out <- matches_check(se),
                 "Annotation object matches converted to logical")
  converted <- annotationMatches(out)
  expect_true(is(converted, "lMatrix"))
  expect_equal(as.matrix(converted) * 1, m, check.attributes = FALSE)
})

test_that("matches_check requires an assay with a recognized name", {
  se <- SummarizedExperiment(assays = list(foo = matrix(TRUE, 2, 2)))
  expect_error(matches_check(se), "assayNames")
})

test_that("matches_check rejects objects that are not SummarizedExperiments", {
  expect_error(matches_check(matrix(TRUE, 2, 2)), "SummarizedExperiment")
})

# Small vector helpers -------------------------------------------------------------------

test_that("%ni% is the negation of %in%", {
  expect_equal(c(1, 2, 3, NA) %ni% c(2, NA), c(TRUE, FALSE, TRUE, FALSE))
  expect_equal(character(0) %ni% "a", logical(0))
})

test_that("last_element returns the final element", {
  expect_equal(chromVAR:::last_element(c(4, 5, 6)), 6)
  expect_equal(chromVAR:::last_element("a"), "a")
  expect_equal(chromVAR:::last_element(list(1, "z")), "z", check.attributes = FALSE)
})

test_that("is.installed tells installed from missing packages", {
  expect_true(chromVAR:::is.installed("stats"))
  expect_true(chromVAR:::is.installed("Matrix"))
  expect_false(chromVAR:::is.installed("noSuchPackageForChromVAR"))
})

test_that("zero_or_null is TRUE only for NULL and zero", {
  expect_true(chromVAR:::zero_or_null(NULL))
  expect_true(chromVAR:::zero_or_null(0))
  expect_false(chromVAR:::zero_or_null(3))
  expect_false(chromVAR:::zero_or_null(-1))
})

test_that("all_true and all_false test logical vectors", {
  expect_true(chromVAR:::all_true(c(TRUE, TRUE)))
  expect_false(chromVAR:::all_true(c(TRUE, FALSE)))
  expect_true(chromVAR:::all_true(logical(0)))
  expect_true(chromVAR:::all_false(c(FALSE, FALSE)))
  expect_false(chromVAR:::all_false(c(FALSE, TRUE)))
})

test_that("all_true and all_false reject non-logical input", {
  expect_error(chromVAR:::all_true(1), "logical")
  expect_error(chromVAR:::all_false("a"), "logical")
})

test_that("all_equal tests whether all elements equal the first", {
  expect_true(chromVAR:::all_equal(c(2, 2, 2)))
  expect_false(chromVAR:::all_equal(c(2, 2, 3)))
  expect_true(chromVAR:::all_equal(c("a", "a")))
  expect_error(chromVAR:::all_equal(matrix(1:4, 2)), "is.vector")
})

test_that("all_whole tests for integer-valued numbers within a tolerance", {
  expect_true(chromVAR:::all_whole(c(1, 2, 3)))
  expect_true(chromVAR:::all_whole(c(1, 2 + 1e-10)))
  expect_false(chromVAR:::all_whole(c(1, 1.5)))
  expect_true(chromVAR:::all_whole(c(1, 1.1), tol = 0.2))
  expect_true(chromVAR:::all_whole(c(-3, 0, 4)))
})

test_that("getstrand returns strands as characters", {
  gr <- GRanges("chr1", IRanges(1:3, width = 5), strand = c("+", "-", "*"))
  expect_equal(chromVAR:::getstrand(gr), c("+", "-", "*"))
})

# Mathematical helpers ----------------------------------------------------------------------

test_that("rms is the root mean square", {
  expect_equal(chromVAR:::rms(c(3, 4)), sqrt(12.5))
  expect_equal(chromVAR:::rms(c(-2, 2)), 2)
  expect_equal(chromVAR:::rms(rep(0, 4)), 0)
})

test_that("means.along averages over the chosen dimension", {
  local_seed(3)
  a <- array(rnorm(60), dim = c(3, 4, 5))
  expect_equal(chromVAR:::means.along(a, 3), apply(a, c(1, 2), mean))
  expect_equal(chromVAR:::means.along(a, 1), apply(a, c(2, 3), mean))
  expect_equal(chromVAR:::means.along(a, 2), apply(a, c(1, 3), mean))
})

test_that("sd.along is the sample standard deviation over the chosen dimension", {
  local_seed(4)
  a <- array(rnorm(60), dim = c(3, 4, 5))
  expect_equal(chromVAR:::sd.along(a, 3), apply(a, c(1, 2), sd))
  expect_equal(chromVAR:::sd.along(a, 1), apply(a, c(2, 3), sd))
  expect_equal(chromVAR:::sd.along(a, 2), apply(a, c(1, 3), sd))
})

test_that("mean_smooth averages a window with the edges reflected", {
  reference <- function(x, w) {
    n <- length(x)
    offs <- seq(-(w %/% 2), (w - 1) %/% 2)
    vapply(seq_len(n), function(i) {
      j <- i + offs
      j <- ifelse(j < 1, 1 - j, ifelse(j > n, 2 * n + 1 - j, j))
      mean(x[j])
    }, 0)
  }
  local_seed(8)
  x <- rnorm(12)
  for (w in 2:7) {
    expect_equal(mean_smooth(x, w), reference(x, w))
  }
  expect_equal(mean_smooth(1:10, 3)[c(1, 2, 10)], c(4 / 3, 2, 29 / 3))
  expect_equal(mean_smooth(1:10, 4)[c(1, 2, 10)], c(1.5, 1.75, 9.25))
})

test_that("mean_smooth keeps the length and constants", {
  expect_equal(length(mean_smooth(runif(20), 5)), 20)
  expect_equal(mean_smooth(rep(3, 15), 4), rep(3, 15))
})

test_that("mean_smooth validates the window", {
  x <- 1:10
  expect_error(mean_smooth(x, 1), "window must be an integer between 2")
  expect_error(mean_smooth(x, 10), "window must be an integer between 2")
  expect_error(mean_smooth(x, 25), "window must be an integer between 2")
  expect_error(mean_smooth(x, NA), "window must be an integer between 2")
})

test_that("get_jaccard gives pairwise Jaccard indices between columns", {
  local_seed(12)
  x <- matrix(rbinom(100, 1, 0.4), 20, 5)
  expected <- outer(seq_len(5), seq_len(5), Vectorize(function(i, j) {
    sum(x[, i] & x[, j]) / sum(x[, i] | x[, j])
  }))
  expect_equal(unname(chromVAR:::get_jaccard(x)), expected)
  expect_equal(unname(chromVAR:::get_jaccard(Matrix::Matrix(x, sparse = TRUE))),
               expected)
  expect_equal(unname(diag(chromVAR:::get_jaccard(x))), rep(1, 5))
})

test_that("get_jaccard is 0 for disjoint columns and 1 for identical ones", {
  x <- cbind(c(1, 1, 0, 0), c(0, 0, 1, 1), c(1, 1, 0, 0))
  j <- chromVAR:::get_jaccard(x)
  expect_equal(unname(j[1, 2]), 0)
  expect_equal(unname(j[1, 3]), 1)
})

# tabulate2 ---------------------------------------------------------------------------------------

test_that("tabulate2 counts values over a positive range", {
  expect_equal(tabulate2(c(1, 2, 2, 5), 1, 5),
               c("1" = 1, "2" = 2, "3" = 0, "4" = 0, "5" = 1))
  expect_equal(tabulate2(c(3, 4), 2, 6),
               c("2" = 0, "3" = 1, "4" = 1, "5" = 0, "6" = 0))
  expect_equal(tabulate2(c(1, 10), 2, 4), c("2" = 0, "3" = 0, "4" = 0))
})

test_that("tabulate2 counts values over a range starting at zero", {
  expect_equal(tabulate2(c(0, 1, 1), 0, 3),
               c("0" = 1, "1" = 2, "2" = 0, "3" = 0))
  expect_equal(tabulate2(integer(0), 0, 2), c("0" = 0, "1" = 0, "2" = 0))
})

test_that("tabulate2 counts values over a range spanning zero", {
  expect_equal(tabulate2(c(-1, 0, 0, 2), -2, 2),
               c("-2" = 0, "-1" = 1, "0" = 2, "1" = 0, "2" = 1))
  expect_equal(tabulate2(c(-5, 1), -2, 2),
               c("-2" = 0, "-1" = 0, "0" = 0, "1" = 1, "2" = 0))
})

test_that("tabulate2 counts values over a range ending at zero", {
  expect_equal(tabulate2(c(-1, -1, 0), -3, 0),
               c("-3" = 0, "-2" = 0, "-1" = 2, "0" = 1))
})

test_that("tabulate2 counts values over a negative range", {
  expect_equal(tabulate2(c(-1, -3), -3, -1),
               c("-3" = 1, "-2" = 0, "-1" = 1))
})

test_that("tabulate2 requires max_val to exceed min_val", {
  expect_error(tabulate2(1, 3, 3), "max_val must be greater than min_val")
  expect_error(tabulate2(1, 4, 2), "max_val must be greater than min_val")
})

# merge_lists -----------------------------------------------------------------------------------

test_that("merge_lists by order concatenates elements at the same position", {
  out <- chromVAR:::merge_lists(list(1, 2), list(3, 4))
  expect_equal(out, list(c(1, 3), c(2, 4)))
  out <- chromVAR:::merge_lists(list(c(1, 2), "a"), list(3, "b"), list(5, "c"))
  expect_equal(out, list(c(1, 2, 3, 5), c("a", "b", "c")))
})

test_that("merge_lists by name concatenates elements with the same name", {
  out <- chromVAR:::merge_lists(list(a = 1, b = 2), list(b = 3, c = 4),
                                by = "name")
  expect_equal(out, list(b = c(2, 3), a = 1, c = 4))
})

test_that("merge_lists by name merges repeated names within one list", {
  out <- chromVAR:::merge_lists(list(a = 1, a = 2, b = 3), by = "name")
  expect_equal(out, list(a = c(1, 2), b = 3))
  unique_names <- list(a = 1, b = 2)
  expect_equal(chromVAR:::merge_lists(unique_names, by = "name"), unique_names)
})

test_that("merge_lists by name handles three lists", {
  out <- chromVAR:::merge_lists(list(a = 1), list(a = 2, b = 3),
                                list(b = 4, c = 5), by = "name")
  expect_equal(sort(names(out)), c("a", "b", "c"))
  expect_equal(out$a, c(1, 2))
  expect_equal(out$b, c(3, 4))
  expect_equal(out$c, 5)
})

test_that("merge_lists rejects an unknown by argument", {
  expect_error(chromVAR:::merge_lists(list(1), by = "position"),
               "should be one of")
})

# remove_overlap and remove_nonoverlap -----------------------------------------------------------

test_that("remove_overlap drops the given items from every list element", {
  out <- chromVAR:::remove_overlap(list(a = 1:5, b = 3:7), c(3, 4))
  expect_equal(out, list(a = c(1L, 2L, 5L), b = 5:7))
})

test_that("remove_nonoverlap keeps only the given items in every list element", {
  out <- chromVAR:::remove_nonoverlap(list(a = 1:5, b = 3:7), c(3, 4))
  expect_equal(out, list(a = 3:4, b = 3:4))
  expect_equal(chromVAR:::remove_nonoverlap(list(a = 1:2), 9),
               list(a = integer(0)))
})

test_that("remove_overlap and remove_nonoverlap reject input that is not a list", {
  expect_error(chromVAR:::remove_overlap(sum, 1), "x must be list or vector")
  expect_error(chromVAR:::remove_nonoverlap(sum, 1), "x must be list or vector")
})

# split_alpha_numeric and mxsort ----------------------------------------------------------------------

test_that("split_alpha_numeric splits runs of letters and digits", {
  split <- chromVAR:::split_alpha_numeric
  expect_equal(split("abc12de"), c("abc", "12", "de"))
  expect_equal(split("a1b22c"), c("a", "1", "b", "22", "c"))
  expect_equal(split("ab12"), c("ab", "12"))
  expect_equal(split("abc"), "abc")
})

test_that("mxsort orders strings so that embedded numbers sort numerically", {
  expect_equal(mxsort(c("g10", "g2", "g1")), c("g1", "g2", "g10"))
  expect_equal(mxsort(c("a2", "a10", "b1", "a1")), c("a1", "a2", "a10", "b1"))
  expect_equal(mxsort(c("x-10-a", "x-2-b", "x-2-a")),
               c("x-2-a", "x-2-b", "x-10-a"))
  expect_equal(mxsort(c("chr10", "chr2", "chr1")), c("chr1", "chr2", "chr10"))
})

test_that("mxsort orders strings without digits alphabetically", {
  expect_equal(mxsort(c("b", "c", "a")), c("a", "b", "c"))
})

test_that("mxsort sorts non-character input numerically", {
  expect_equal(mxsort(c(10, 2, 1)), c(1, 2, 10))
})

test_that("mxsort keeps all elements and is a permutation of the input", {
  x <- c("RG-12", "RG-3", "RG-3b", "RG-100", "RG-7")
  out <- mxsort(x)
  expect_equal(sort(out), sort(x))
  expect_equal(out[1:2], c("RG-3", "RG-3b"))
  expect_equal(out[length(out)], "RG-100")
})

# Index lists and matrices -----------------------------------------------------------------------------

test_that("convert_to_ix_list gives the row indices of each column", {
  m <- matrix(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE), 3,
              dimnames = list(NULL, c("a", "b")))
  expected <- list(a = c(1L, 3L), b = 3L)
  expect_equal(chromVAR:::convert_to_ix_list(m), expected)
  expect_equal(chromVAR:::convert_to_ix_list(Matrix::Matrix(m, sparse = TRUE)),
               expected)
  expect_equal(chromVAR:::convert_to_ix_list(Matrix::Matrix(m, sparse = FALSE)),
               expected)
})

test_that("convert_to_ix_list gives an empty entry for an all-FALSE column", {
  m <- Matrix::sparseMatrix(i = c(1, 3), j = c(1, 1), x = TRUE, dims = c(3, 2))
  out <- chromVAR:::convert_to_ix_list(m)
  expect_equal(length(out), 2)
  expect_equal(as.integer(out[[1]]), c(1L, 3L))
  expect_equal(length(out[[2]]), 0)
  out <- chromVAR:::convert_to_ix_list(matrix(FALSE, 3, 2))
  expect_equal(lengths(out), c(0L, 0L), check.attributes = FALSE)
})

test_that("convert_to_ix_list rejects input that is not a matrix", {
  expect_error(chromVAR:::convert_to_ix_list(1:3), "Matrix")
  expect_error(chromVAR:::convert_to_ix_list(list(1)), "Matrix")
})

test_that("convert_from_ix_list builds a logical sparse matrix", {
  m <- chromVAR:::convert_from_ix_list(list(a = c(1, 3), b = 3), 4)
  expect_equal(dim(m), c(4, 2))
  expect_equal(colnames(m), c("a", "b"))
  expect_equal(unname(as.matrix(m)),
               matrix(c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
                      4, 2))
  expect_true(is(m, "lMatrix"))
})

test_that("convert_from_ix_list handles empty index vectors", {
  m <- chromVAR:::convert_from_ix_list(list(a = integer(0), b = c(2, 3)), 3)
  expect_equal(unname(as.matrix(m)),
               matrix(c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE), 3, 2))
})

test_that("convert_from_ix_list rejects input that is not a list", {
  expect_error(chromVAR:::convert_from_ix_list(1:3, 4), "list")
})

test_that("convert_to_ix_list undoes convert_from_ix_list", {
  ix <- list(a = c(1L, 4L, 6L), b = c(2L, 3L), c = 6L)
  expect_equal(chromVAR:::convert_to_ix_list(
    chromVAR:::convert_from_ix_list(ix, 6)), ix, check.attributes = FALSE)
})

# validate_genome_input --------------------------------------------------------------------------------------

test_that("validate_genome_input returns DNAStringSet genomes unchanged", {
  g <- Biostrings::DNAStringSet(c(chr1 = "ACGT", chr2 = "TTGA"))
  expect_identical(chromVAR:::validate_genome_input(g), g)
})

test_that("validate_genome_input returns FaFile genomes unchanged", {
  fa <- Rsamtools::FaFile(tempfile(fileext = ".fa"))
  expect_identical(chromVAR:::validate_genome_input(fa), fa)
})

test_that("validate_genome_input looks up a genome name with BSgenome", {
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg19")
  g <- chromVAR:::validate_genome_input("hg19")
  expect_is(g, "BSgenome")
  expect_is(chromVAR:::validate_genome_input(c("hg19", "hg19")), "BSgenome")
})

test_that("validate_genome_input rejects a missing genome name", {
  expect_error(chromVAR:::validate_genome_input(NA_character_),
               "No genome provided")
})

test_that("validate_genome_input rejects several different genome names", {
  expect_error(chromVAR:::validate_genome_input(c("hg19", "mm10")),
               "genome == genome")
})

test_that("validate_genome_input rejects other types", {
  expect_error(chromVAR:::validate_genome_input(1),
               "genome input must be a BSgenome, DNAStringSet, or FaFile")
  expect_error(chromVAR:::validate_genome_input(NULL),
               "genome input must be a BSgenome")
  expect_error(chromVAR:::validate_genome_input(list("hg19")),
               "genome input must be a BSgenome")
})
