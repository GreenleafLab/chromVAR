context("computeDeviations")

library(SummarizedExperiment)
library(IRanges)
library(Matrix)

mat <- matrix(c(0,1,2,1,1,0,2,0,1,0,1,1), nrow = 3)

counts <- SummarizedExperiment( assays= list(counts = Matrix(mat)),
  rowRanges = GRanges(seqnames = c("chr1","chr2","chr2"),
                                ranges = IRanges(start = c(76585873,
                                                           42772928,
                                                           100183786),
                                                          width = 500)))


anno_list <- list(c(1,2,3),c(1,2),c(2,3))

anno_mat <- sparseMatrix(j = c(rep(1,3),rep(2,2),rep(3,2)),
                         i = c(1,2,3,1,2,2,3),
                         x = TRUE)

anno_se <- SummarizedExperiment(assays = list(matches = anno_mat))

e <- rowSums(mat) / sum(mat) 

bg <- matrix(c(2,2,2,1,1,3,2,2,1), nrow = 3, byrow = TRUE)

deviations_test <- readRDS(system.file("extdata", "deviations_test", 
                                       package = "chromVAR"))



test_that("computeDeviations works with SummarizedExperiment for both args",
          {
            dev <- computeDeviations(counts, anno_se, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with SummarizedExperiment & list",
          {
            dev <- computeDeviations(counts, anno_list, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with SummarizedExperiment & Matrix",
          {
            dev <- computeDeviations(counts, anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with matrix & SummarizedExperient",
          {
            dev <- computeDeviations(mat, anno_se, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with matrix & list",
          {
            dev <- computeDeviations(mat, anno_list, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with matrix & Matrix",
          {
            dev <- computeDeviations(mat, anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("computeDeviations works with Matrix & Matrix",
          {
            dev <- computeDeviations(Matrix(mat), anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })