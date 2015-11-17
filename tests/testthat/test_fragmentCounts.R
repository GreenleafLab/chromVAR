context('fragmentCounts')

# Read in data needed for tests ------------------------------------------------

test_rg <- system.file("extdata", "test_RG.bam", package = "tfdev")
test_bam1 <- system.file("extdata", "test_single1.bam", package = "tfdev")
test_bam2 <- system.file("extdata", "test_single2.bam", package = "tfdev")
test_bam3 <- system.file("extdata", "test_single3.bam", package = "tfdev")
test_peaks <- readRDS("test_peaks.rds") 

# Test reading with RG ---------------------------------------------------------


test_that("can read bam with RG", {
  fragments <- bamToFragmentsByRG(bamfile =  test_rg)
  expect_is(fragments, "list")
  expect_equal_to_reference(fragments, "rg_fragments.rds")
})

# Test fragment counts with RG ___________--------------------------------------

test_that("can count fragments using RG tags", {
  counts <- getFragmentCountsByRG(test_rg, test_peaks)
  expect_is(counts, "fragmentCounts")
  expect_equal_to_reference(counts, "counts_rg.rds")
})

# Test reading  bam without RG -------------------------------------------------

test_that("can read bam without RG", {
  fragments <- bamToFragments(bamfile = test_bam1)
  expect_is(fragments, "GenomicRangesList")
  expect_equal_to_reference(fragments, "fragments1.rds")
})

# Test fragment counts with multiple bam --------------------------------------

test_that("can count fragments with multiple bams", {
  counts <- getFragmentCounts(c(test_bam1, test_bam2, test_bam3), test_peaks)
  expect_is(counts, "fragmentCounts")
  expect_equal_to_reference(counts, "counts_multiple_bam.rds")
})

# Test subsetting of fragment counts -------------------------------------------

test_that("can subset fragmentCounts by numeric j", {
  counts <- readRDS("counts_rg.rds") 
  counts <- counts[,1:3]
  expect_is(counts, "fragmentCounts")
  expect_equal(counts@nsample, 3)
  expect_equal(ncol(counts@counts),3)
})

test_that("can subset fragmentCounts by numeric i", {
  counts <- readRDS("counts_rg.rds") 
  counts <- counts[1:3,]
  expect_is(counts, "fragmentCounts")
  expect_equal(counts@npeak, 3)
  expect_equal(nrow(counts@counts),3)
})

test_that("can subset fragmentCounts by numeric i,j", {
  counts <- readRDS("counts_rg.rds") 
  counts <- counts[1:3,1:3]
  expect_is(counts, "fragmentCounts")
  expect_equal(counts@npeak, 3)
  expect_equal(nrow(counts@counts),3)
  expect_equal(counts@nsample, 3)
  expect_equal(ncol(counts@counts),3)
})








