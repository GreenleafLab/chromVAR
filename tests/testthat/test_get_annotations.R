context("getAnnotations")

anno_file1 <- system.file("extdata", "test_anno1.bed", package = "chromVAR")
anno_file2 <- system.file("extdata", "test_anno2.bed", package = "chromVAR")
anno_file3 <- system.file("extdata", "test_anno3.bed", package = "chromVAR")

peaks <- 
  GenomicRanges::GRanges(
    seqnames = c("chr1","chr2","chr2"),
    ranges = IRanges::IRanges(start = c(76585873,42772928,100183786),
                              width = 500))

anno_ranges <- 
  GenomicRanges::GRangesList(
    GenomicRanges::GRanges(seqnames = c("chr1"),
                           ranges = IRanges::IRanges(start = 76585933,
                                                     end = 76585942)),
    GenomicRanges::GRanges(seqnames = c("chr2"),
                           ranges = 
                             IRanges::IRanges(start = c(42773087,100183852,
                                                        100183852),
                                              end = c(42773097,100183862,
                                                      100183862))),
    GenomicRanges::GRanges())


anno_mat <- matrix(c(TRUE,FALSE,FALSE,FALSE,TRUE,TRUE), nrow = 3, byrow = FALSE)

anno_list <- list(c(1),c(2,3),c())

test_that("Can read annotations from multiple bed file",{
  res <- getAnnotations(c(anno_file1, anno_file2), rowRanges = peaks)
  expect_equal(as.matrix(annotationMatches(res))[,1], c(TRUE, FALSE, FALSE))
  expect_equal(as.matrix(annotationMatches(res))[,2], c(FALSE, TRUE, TRUE))
  expect_equal(dim(res), c(3,2))
  expect_is(res, "RangedSummarizedExperiment")
})

test_that("Can read annotations from single bed file",{
  res <- getAnnotations(anno_file3, rowRanges = peaks, column = 4)
  expect_equal(as.matrix(annotationMatches(res))[,1], c(TRUE, FALSE, FALSE))
  expect_equal(as.matrix(annotationMatches(res))[,2], c(FALSE, TRUE, TRUE))
  expect_equal(dim(res), c(3,2))
  expect_is(res, "RangedSummarizedExperiment")
})

test_that("Can get annotations from genomic ranges",{
  res <- getAnnotations(anno_ranges, rowRanges = peaks)
  expect_equal(as.matrix(annotationMatches(res))[,1], c(TRUE, FALSE, FALSE))
  expect_equal(as.matrix(annotationMatches(res))[,2], c(FALSE, TRUE, TRUE))
  expect_equal(as.matrix(annotationMatches(res))[,3], c(FALSE, FALSE, FALSE))
  expect_equal(dim(res), c(3,3))
  expect_is(res, "RangedSummarizedExperiment")
})

test_that("Can get annotations from matrix",{
  res <- getAnnotations(anno_mat, peaks)
  expect_equal(as.matrix(annotationMatches(res))[,1], c(TRUE, FALSE, FALSE))
  expect_equal(as.matrix(annotationMatches(res))[,2], c(FALSE, TRUE, TRUE))
  expect_equal(dim(res), c(3,2))
  expect_is(res, "RangedSummarizedExperiment")
})

test_that("Can get annotations from list",{
  res <- getAnnotations(anno_list, npeaks = 3)
  expect_equal(as.matrix(annotationMatches(res))[,1], c(TRUE, FALSE, FALSE))
  expect_equal(as.matrix(annotationMatches(res))[,2], c(FALSE, TRUE, TRUE))
  expect_equal(as.matrix(annotationMatches(res))[,3], c(FALSE, FALSE, FALSE))
  expect_equal(dim(res), c(3,3))
  expect_is(res, "SummarizedExperiment")
})

