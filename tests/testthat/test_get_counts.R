context('fragmentCounts')

# Read in data needed for tests ------------------------------------------------

test_rg <- system.file("extdata", "test_RG.bam", package = "chromVAR")
test_bam1 <- system.file("extdata", "test_single1.bam", package = "chromVAR")
test_bam2 <- system.file("extdata", "test_single2.bam", package = "chromVAR")
test_bam3 <- system.file("extdata", "test_single3.bam", package = "chromVAR")
peaks_file <- system.file("extdata", "test_bed.txt", package = "chromVAR")
test_peaks <- get_peaks(peaks_file)


# Test fragment counts with RG ___________--------------------------------------

test_that("can count fragments using RG tags", {
  counts <- get_counts(test_rg, test_peaks, by_rg = TRUE, paired = TRUE)
  expect_is(counts, "SummarizedExperiment")
  expect_equal_to_reference(counts, "counts_rg.rds")
})


# Test fragment counts with multiple bam --------------------------------------

test_that("can count fragments with multiple bams", {
  counts <- get_counts(c(test_bam1, test_bam2, test_bam3), test_peaks, by_rg = FALSE, paired = TRUE)
  expect_is(counts, "SummarizedExperiment")
  expect_equal_to_reference(counts, "counts_multiple_bam.rds")
})









