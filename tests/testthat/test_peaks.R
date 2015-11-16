context("peaks")

# Test reading in peaks---------------------------------------------------------

peaks_file <- system.file("extdata", "test_bed.txt", package = "tfdev")

test_that("can read peaks file", {
  peaks <- read_peaks(peaks_file)
  expect_is(peaks, "GenomicRanges")
  expect_equal_to_reference(peaks, "test_peaks.rds")
})

test_that("can read peaks file with one extra column", {
  peaks <- read_peaks(peaks_file, extra_cols = 4, extra_names = "name")
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), "name") 
})

test_that("can read peaks file with several extra columns", {
  peaks <- read_peaks(peaks_file, extra_cols = c(4,5), 
                      extra_names = c("name","score"))
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), c("name","score")) 
})