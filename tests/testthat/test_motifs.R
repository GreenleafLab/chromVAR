context("getJasparMotifs")

library(TFBSTools)

nucs <- c("A", "C", "G", "T")

make_pwm <- function(id, nm) {
  PWMatrix(ID = id, name = nm,
           profileMatrix = matrix(log2(c(0.7, 0.1, 0.1, 0.1) / 0.25), 4,
                                  dimnames = list(nucs, NULL)))
}

fake_motif_set <- function() {
  out <- PWMatrixList(make_pwm("MA0001.1", "AGL3"), make_pwm("MA0002.2", "RUNX1"))
  names(out) <- c("MA0001.1", "MA0002.2")
  out
}

skip_if_no_mocking <- function() {
  skip_if(utils::packageVersion("testthat") < "3.2.0",
          "local_mocked_bindings needs testthat 3.2.0")
}

test_that("getJasparMotifs queries human CORE motifs by default", {
  skip_if_no_mocking()
  seen <- list()
  testthat::local_mocked_bindings(
    getMatrixSet = function(x, opts) {
      seen <<- list(db = x, opts = opts)
      fake_motif_set()
    },
    .package = "TFBSTools")
  getJasparMotifs(jaspar_db = "fake_db")
  expect_equal(seen$db, "fake_db")
  expect_equal(seen$opts, list(species = "Homo sapiens", collection = "CORE"))
})

test_that("getJasparMotifs passes species and collection on", {
  skip_if_no_mocking()
  seen <- list()
  testthat::local_mocked_bindings(
    getMatrixSet = function(x, opts) {
      seen <<- opts
      fake_motif_set()
    },
    .package = "TFBSTools")
  getJasparMotifs(species = "Mus musculus", collection = "CNE",
                  jaspar_db = "fake_db")
  expect_equal(seen, list(species = "Mus musculus", collection = "CNE"))
  getJasparMotifs(species = 9606, collection = "PBM", jaspar_db = "fake_db")
  expect_equal(seen, list(species = 9606, collection = "PBM"))
})

test_that("getJasparMotifs appends additional options", {
  skip_if_no_mocking()
  seen <- list()
  testthat::local_mocked_bindings(
    getMatrixSet = function(x, opts) {
      seen <<- opts
      fake_motif_set()
    },
    .package = "TFBSTools")
  getJasparMotifs(jaspar_db = "fake_db", all_versions = TRUE, min_ic = 10)
  expect_equal(seen, list(species = "Homo sapiens", collection = "CORE",
                          all_versions = TRUE, min_ic = 10))
})

test_that("getJasparMotifs names motifs by ID and name", {
  skip_if_no_mocking()
  testthat::local_mocked_bindings(
    getMatrixSet = function(x, opts) fake_motif_set(),
    .package = "TFBSTools")
  res <- getJasparMotifs(jaspar_db = "fake_db")
  expect_equal(names(res), c("MA0001.1_AGL3", "MA0002.2_RUNX1"))
  expect_equal(length(res), 2)
  expect_true(is(res, "PWMatrixList"))
  expect_equal(unname(TFBSTools::ID(res)), c("MA0001.1", "MA0002.2"))
  expect_equal(unname(TFBSTools::name(res)), c("AGL3", "RUNX1"))
})

test_that("getJasparMotifs returns the matrices unchanged", {
  skip_if_no_mocking()
  fake <- fake_motif_set()
  testthat::local_mocked_bindings(
    getMatrixSet = function(x, opts) fake,
    .package = "TFBSTools")
  res <- getJasparMotifs(jaspar_db = "fake_db")
  expect_equal(as.matrix(res[[1]]), as.matrix(fake[[1]]))
  expect_equal(as.matrix(res[[2]]), as.matrix(fake[[2]]))
})

test_that("getJasparMotifs reads motifs from the JASPAR2016 database", {
  skip_if_not_installed("JASPAR2016")
  res <- getJasparMotifs()
  expect_true(length(res) > 100)
  expect_true(all(grepl("_", names(res), fixed = TRUE)))
  expect_equal(names(res), paste(TFBSTools::ID(res), TFBSTools::name(res),
                                 sep = "_"), ignore_attr = TRUE)
})
