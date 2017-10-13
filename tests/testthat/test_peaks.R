context("getPeaks")

# Test reading in peaks---------------------------------------------------------

peaks_file <- system.file("extdata", "test_bed.txt", package = "chromVAR")

test_that("can read peaks file", {
  peaks <- getPeaks(peaks_file, sort_peaks = TRUE)
  expect_is(peaks, "GenomicRanges")
  expect_equal(start(peaks), 
               c(143913189,64977984,176032858,130506305,100183786,25315858,
                 76585873,6255789,45973111,42772928,45996001))
  expect_equal(width(peaks) , rep(501, length(peaks)))
  expect_equal(as.character(seqnames(peaks)), 
               c("chr1", "chr2", "chr2", "chr5", "chr7", "chr8", "chr10",
                 "chr11", "chr17", "chr19", "chr19"))
})

test_that("can read peaks file with one extra column", {
  peaks <- getPeaks(peaks_file, extra_cols = c("name" = 4), sort_peaks = TRUE)
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), "name")
})

test_that("can read peaks file with several extra columns", {
  peaks <- getPeaks(peaks_file, extra_cols = c("name" = 4, "score" = 5), 
                    sort_peaks = TRUE)
  expect_is(peaks, "GenomicRanges")
  expect_equal(colnames(mcols(peaks)), c("name","score"))
})

test_that("reduce peaks works as expected", {
  tmp <- GenomicRanges::GRanges(c("chr1","chr1","chr1","chr1","chr2","chr2"),
                               IRanges::IRanges(start = 
                                                  c(1,10,20,30,10,11), 
                                                end = c(9,23,33,42,18,14)))
  tmp2 <- 
    SummarizedExperiment(assays = list(counts  = matrix(rep(c(1,3,2,4,7,3),2),
                                                        ncol = 2)), 
                         rowRanges = tmp)
  expect_equal(filterPeaks(tmp2, ix_return = TRUE), c(1,2,4,5))
})
