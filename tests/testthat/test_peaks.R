context("peaks")

# Test reading in peaks---------------------------------------------------------

peaks_file <- system.file("extdata", "test_bed.txt", package = "chromVAR")

test_that("can read peaks file", {
  peaks <- get_peaks(peaks_file)
  expect_is(peaks, "GenomicRanges")
  expect_equal_to_reference(peaks, "test_peaks.rds")
})

test_that("can read peaks file with one extra column", {
  peaks <- get_peaks(peaks_file, extra_cols = c("name" = 4))
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), "name") 
})

test_that("can read peaks file with several extra columns", {
  peaks <- get_peaks(peaks_file, extra_cols = c("name" = 4, "score" = 5))
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), c("name","score")) 
})