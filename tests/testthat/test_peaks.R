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

test_that("reduce peaks works as expected"){
  tmp = GenomicRanges::GRanges(c("chr1","chr1","chr1","chr1","chr2","chr2"), 
                               IRanges::IRanges(start = c(1,10,20,30,10,11), end = c(9,23,33,42,18,14)))
  tmp2 = matrix(c(1,3,2,4,7,3), ncol = 1)
  expect_equal(filter_peaks(tmp2, tmp), c(1,2,4,5))
}