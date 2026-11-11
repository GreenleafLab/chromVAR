context("classes and accessors")

library(SummarizedExperiment)

dev_a <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2,
                dimnames = list(c("a1", "a2"), c("s1", "s2", "s3")))
z_a <- dev_a * 10
dev_b <- matrix(c(7, 8, 9, 10, 11, 12, 13, 14, 15), nrow = 3,
                dimnames = list(c("b1", "b2", "b3"), c("s1", "s2", "s3")))
z_b <- dev_b * 10
cd <- S4Vectors::DataFrame(group = c("x", "y", "x"), depth = c(10, 20, 30),
                           row.names = c("s1", "s2", "s3"))

make_obj <- function(dev, z, rd) {
  se <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                             colData = cd, rowData = rd)
  new("chromVARDeviations", se)
}

obj_a <- make_obj(dev_a, z_a, S4Vectors::DataFrame(
  name = c("A1", "A2"), only_a = c(TRUE, FALSE), shared = c(1, 2),
  row.names = rownames(dev_a)))
obj_b <- make_obj(dev_b, z_b, S4Vectors::DataFrame(
  name = c("B1", "B2", "B3"), only_b = c("p", "q", "r"), shared = c(3, 4, 5),
  row.names = rownames(dev_b)))

# chromVARDeviations class -------------------------------------------------------

test_that("chromVARDeviations extends SummarizedExperiment", {
  expect_true(extends("chromVARDeviations", "SummarizedExperiment"))
  expect_true(is(obj_a, "SummarizedExperiment"))
  expect_silent(validObject(obj_a))
})

test_that("chromVARDeviations validity requires both deviations and z assays", {
  only_dev <- SummarizedExperiment(assays = list(deviations = dev_a))
  only_z <- SummarizedExperiment(assays = list(z = z_a))
  neither <- SummarizedExperiment(assays = list(counts = dev_a))
  msg <- "The assays slot must contain 'deviations' and 'z'"
  expect_error(new("chromVARDeviations", only_dev), msg)
  expect_error(new("chromVARDeviations", only_z), msg)
  expect_error(new("chromVARDeviations", neither), msg)
})

test_that("chromVARDeviations validity is re-checked when assays are removed", {
  obj <- obj_a
  assays(obj) <- assays(obj)["z"]
  expect_error(validObject(obj), "The assays slot must contain")
})

test_that("chromVARDeviations accepts extra assays", {
  se <- SummarizedExperiment(assays = list(deviations = dev_a, z = z_a,
                                           extra = dev_a + 1))
  obj <- new("chromVARDeviations", se)
  expect_equal(assayNames(obj), c("deviations", "z", "extra"))
})

test_that("class unions accept their members", {
  expect_true(is(matrix(1), "MatrixOrmatrix"))
  expect_true(is(Matrix::Matrix(1), "MatrixOrmatrix"))
  expect_false(is(data.frame(a = 1), "MatrixOrmatrix"))
  expect_true(is(NULL, "missingOrNULL"))
  expect_true(extends("missing", "missingOrNULL"))
  expect_false(is(1, "missingOrNULL"))
})

# accessors --------------------------------------------------------------------

test_that("deviations and deviationScores return the matching assays", {
  expect_equal(deviations(obj_a), dev_a)
  expect_equal(deviationScores(obj_a), z_a)
  expect_true(is.matrix(deviations(obj_a)))
})

test_that("deviations and deviationScores keep the dimension names", {
  expect_equal(dimnames(deviations(obj_b)), dimnames(dev_b))
  expect_equal(dimnames(deviationScores(obj_b)), dimnames(z_b))
})

test_that("accessors reject objects that are not chromVARDeviations", {
  se <- SummarizedExperiment(assays = list(deviations = dev_a, z = z_a))
  expect_error(deviations(se))
  expect_error(deviationScores(se))
  expect_error(deviations(dev_a))
})

test_that("subsetting a chromVARDeviations keeps the class, assays and annotations", {
  sub <- obj_b[c(3, 1), c("s3", "s1")]
  expect_is(sub, "chromVARDeviations")
  expect_equal(deviations(sub), dev_b[c(3, 1), c(3, 1)])
  expect_equal(deviationScores(sub), z_b[c(3, 1), c(3, 1)])
  expect_equal(rowData(sub)$name, c("B3", "B1"))
  expect_equal(sub$group, c("x", "x"))
  expect_equal(sub$depth, c(30, 10))
})

test_that("subsetting down to one row or column keeps matrix assays", {
  one_row <- obj_b["b2", ]
  expect_equal(dim(deviations(one_row)), c(1, 3))
  one_col <- obj_b[, "s2"]
  expect_equal(dim(deviationScores(one_col)), c(3, 1))
  expect_is(one_col, "chromVARDeviations")
})

test_that("the mini_dev example is a valid chromVARDeviations object", {
  data(mini_dev, package = "chromVAR")
  expect_is(mini_dev, "chromVARDeviations")
  expect_silent(validObject(mini_dev))
  expect_equal(dim(mini_dev), c(3, 35))
  expect_equal(assayNames(mini_dev), c("deviations", "z"))
  expect_equal(dim(deviations(mini_dev)), dim(deviationScores(mini_dev)))
  expect_true(all(c("name", "fractionMatches", "fractionBackgroundOverlap") %in%
                    colnames(rowData(mini_dev))))
})

# rbind / cbind -------------------------------------------------------------------

test_that("rbind stacks deviations and z in input order", {
  res <- rbind(obj_a, obj_b)
  expect_is(res, "chromVARDeviations")
  expect_equal(deviations(res), rbind(dev_a, dev_b))
  expect_equal(deviationScores(res), rbind(z_a, z_b))
  expect_equal(rownames(res), c("a1", "a2", "b1", "b2", "b3"))
  expect_equal(colnames(res), colnames(obj_a))
})

test_that("rbind keeps only rowData columns shared by all inputs", {
  res <- rbind(obj_a, obj_b)
  expect_equal(sort(colnames(rowData(res))), c("name", "shared"))
  expect_equal(rowData(res)$name, c("A1", "A2", "B1", "B2", "B3"))
  expect_equal(rowData(res)$shared, c(1, 2, 3, 4, 5))
})

test_that("rbind takes colData from the first object", {
  other <- obj_b
  other$group <- c("q", "q", "q")
  res <- rbind(obj_a, other)
  expect_equal(res$group, c("x", "y", "x"))
  expect_equal(res$depth, c(10, 20, 30))
})

test_that("rbind of three objects and of an object with itself", {
  res3 <- rbind(obj_a, obj_b, obj_a)
  expect_equal(dim(res3), c(7, 3))
  expect_equal(deviations(res3), rbind(dev_a, dev_b, dev_a))
  data(mini_dev, package = "chromVAR")
  doubled <- rbind(mini_dev, mini_dev)
  expect_equal(dim(doubled), c(6, 35))
  expect_equal(deviations(doubled), rbind(deviations(mini_dev), deviations(mini_dev)))
  expect_equal(deviationScores(doubled)[4:6, ], deviationScores(mini_dev))
})

test_that("cbind on chromVARDeviations is refused", {
  expect_error(cbind(obj_a, obj_a), "Can't concatenate chromVARDeviations horizontally")
})

# annotationMatches -------------------------------------------------------------

match_mat <- matrix(c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE), nrow = 3,
                    dimnames = list(c("p1", "p2", "p3"), c("m1", "m2")))

make_matches <- function(assay_name, mat = match_mat) {
  assays <- list(mat)
  names(assays) <- assay_name
  SummarizedExperiment(assays = assays)
}

test_that("annotationMatches reads each accepted assay name", {
  for (nm in c("annotationMatches", "annotation_matches", "motifMatches",
               "matches")) {
    expect_equal(annotationMatches(make_matches(nm)), match_mat)
  }
})

test_that("annotationMatches warns when reading the deprecated motif_matches assay", {
  expect_warning(res <- annotationMatches(make_matches("motif_matches")),
                 "motif_matches assay deprecated")
  expect_equal(res, match_mat)
})

test_that("annotationMatches picks one of several recognised assays by the current precedence (characterisation, not in the docs)", {
  other <- !match_mat
  se <- SummarizedExperiment(assays = list(matches = other,
                                           motifMatches = other,
                                           annotation_matches = match_mat))
  expect_equal(annotationMatches(se), match_mat)
  se2 <- SummarizedExperiment(assays = list(matches = other,
                                            motifMatches = match_mat))
  expect_equal(annotationMatches(se2), match_mat)
  se3 <- SummarizedExperiment(assays = list(matches = other,
                                            annotationMatches = match_mat,
                                            annotation_matches = other))
  expect_equal(annotationMatches(se3), match_mat)
})

test_that("annotationMatches errors without a recognised assay", {
  expect_error(annotationMatches(make_matches("counts")),
               "No appropriately named assay")
})

test_that("annotationMatches works on the mini_ix example", {
  data(mini_ix, package = "chromVAR")
  res <- annotationMatches(mini_ix)
  expect_equal(dim(res), dim(mini_ix))
  expect_true(is(res, "lMatrix") || is.logical(res))
})

test_that("annotationMatches<- replaces the assay that is already present", {
  for (nm in c("annotationMatches", "annotation_matches", "motifMatches",
               "matches")) {
    se <- make_matches(nm)
    annotationMatches(se) <- Matrix::Matrix(!match_mat)
    expect_equal(assayNames(se), nm)
    expect_equal(as.matrix(annotationMatches(se)), !match_mat)
  }
})

test_that("annotationMatches<- creates an annotationMatches assay when none exists", {
  se <- make_matches("counts", matrix(1:6, 3, dimnames = dimnames(match_mat)))
  annotationMatches(se) <- Matrix::Matrix(match_mat)
  expect_equal(assayNames(se), c("counts", "annotationMatches"))
  expect_equal(as.matrix(annotationMatches(se)), match_mat)
})

test_that("annotationMatches<- warns when writing the deprecated motif_matches assay", {
  se <- make_matches("motif_matches")
  expect_warning(annotationMatches(se) <- Matrix::Matrix(!match_mat),
                 "motif_matches assay deprecated")
  expect_equal(as.matrix(assays(se)$motif_matches), !match_mat)
})

test_that("annotationMatches<- converts non-logical input with a warning", {
  se <- make_matches("matches")
  numeric_mat <- matrix(c(1, 0, 1, 0, 0, 1), 3, dimnames = dimnames(match_mat))
  expect_warning(annotationMatches(se) <- numeric_mat, "converted to logical")
  expect_equal(as.matrix(annotationMatches(se)),
               matrix(c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE), 3,
                      dimnames = dimnames(match_mat)))
})

test_that("annotationMatches<- accepts a sparse logical Matrix without a warning", {
  se <- make_matches("matches")
  sparse <- Matrix::Matrix(match_mat, sparse = TRUE)
  expect_warning(annotationMatches(se) <- sparse, NA)
  expect_equal(as.matrix(annotationMatches(se)), match_mat)
})

# counts ---------------------------------------------------------------------------

test_that("counts returns the counts assay and errors without one", {
  cnt <- matrix(1:6, 3, dimnames = list(c("p1", "p2", "p3"), c("s1", "s2")))
  se <- SummarizedExperiment(assays = list(counts = cnt))
  expect_equal(counts(se), cnt)
  expect_error(counts(SummarizedExperiment(assays = list(other = cnt))),
               "counts")
})

test_that("counts<- sets, replaces and validates the counts assay", {
  se <- SummarizedExperiment(assays = list(other = matrix(0, 3, 2)))
  new_counts <- matrix(1:6, 3)
  counts(se) <- new_counts
  expect_equal(assayNames(se), c("other", "counts"))
  expect_equal(unname(counts(se)), new_counts)
  counts(se) <- new_counts * 2L
  expect_equal(unname(counts(se)), new_counts * 2L)
  expect_equal(assayNames(se), c("other", "counts"))
})

test_that("counts<- accepts a sparse Matrix", {
  se <- SummarizedExperiment(assays = list(counts = matrix(0, 3, 2)))
  counts(se) <- Matrix::Matrix(matrix(c(0, 1, 0, 2, 0, 0), 3), sparse = TRUE)
  expect_true(is(counts(se), "Matrix"))
  expect_equal(as.matrix(counts(se))[2, 1], 1)
})

test_that("counts<- rejects values of the wrong dimension or type", {
  se <- SummarizedExperiment(assays = list(counts = matrix(0, 3, 2)))
  expect_error(counts(se) <- matrix(1, 4, 2))
  expect_error(counts(se) <- data.frame(a = 1:3, b = 1:3))
})

test_that("counts works on the mini_counts example", {
  data(mini_counts, package = "chromVAR")
  cnt <- counts(mini_counts)
  expect_equal(dim(cnt), dim(mini_counts))
  expect_true(all(as.matrix(cnt) >= 0))
})

# package-level -----------------------------------------------------------------

test_that("the package exports the documented functions", {
  exports <- getNamespaceExports("chromVAR")
  expected <- c("differentialDeviations", "differentialVariability",
                "deviationsTsne", "plotDeviationsTsne", "getSampleCorrelation",
                "getSampleDistance", "plotVariability", "chromVAR_theme",
                "annotationMatches", "deviations", "deviationScores", "counts")
  expect_true(all(expected %in% exports))
  for (fn in expected) {
    expect_true(is.function(get(fn, envir = asNamespace("chromVAR"))))
  }
})

test_that("the example datasets have mutually consistent dimensions", {
  data(mini_counts, package = "chromVAR")
  data(mini_ix, package = "chromVAR")
  data(mini_dev, package = "chromVAR")
  data(example_counts, package = "chromVAR")
  expect_is(mini_counts, "RangedSummarizedExperiment")
  expect_is(example_counts, "RangedSummarizedExperiment")
  expect_equal(nrow(mini_ix), nrow(mini_counts))
  expect_equal(colnames(mini_dev), colnames(mini_counts))
  expect_equal(nrow(mini_dev), ncol(annotationMatches(mini_ix)))
  expect_true("counts" %in% assayNames(example_counts))
})
