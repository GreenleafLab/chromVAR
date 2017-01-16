context("compute_deviations")

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

e = rowSums(mat) / sum(mat) 

bg = matrix(c(2,2,2,1,1,3,2,2,1), nrow = 3, byrow = TRUE)

deviations_test <- readRDS(system.file("extdata", "deviations_test", 
                                       package = "chromVAR"))



test_that("compute_deviations works with SummarizedExperiment & SummarizedExperient",
          {
            dev <- compute_deviations(counts, anno_se, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with SummarizedExperiment & list",
          {
            dev <- compute_deviations(counts, anno_list, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with SummarizedExperiment & Matrix",
          {
            dev <- compute_deviations(counts, anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with matrix & SummarizedExperient",
          {
            dev <- compute_deviations(mat, anno_se, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with matrix & list",
          {
            dev <- compute_deviations(mat, anno_list, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with matrix & Matrix",
          {
            dev <- compute_deviations(mat, anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })

test_that("compute_deviations works with Matrix & Matrix",
          {
            dev <- compute_deviations(Matrix(mat), anno_mat, bg)
            expect_equal(assays(dev)$deviations, deviations_test$deviations)
            expect_equal(assays(dev)$z, deviations_test$z)
            expect_is(dev, "chromVARDeviations")
          })