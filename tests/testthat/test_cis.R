context("cis groups")

library(SummarizedExperiment)
library(GenomicRanges)
library(IRanges)
library(Matrix)

local_serial_param()

group_members <- function(res) {
  lapply(seq_len(ncol(res)), function(j) which(as.logical(assay(res)[, j])))
}

make_peaks <- function(counts_per_chr) {
  chrs <- names(counts_per_chr)
  sort(GRanges(rep(chrs, counts_per_chr),
               IRanges(start = unlist(lapply(counts_per_chr, seq_len)) * 100,
                       width = 50)))
}

test_that("getCisGroups makes overlapping windows of grpsize peaks every stepsize peaks",
          {
            peaks <- make_peaks(c(chr1 = 12))
            res <- getCisGroups(peaks, grpsize = 4, stepsize = 3)
            expect_equal(group_members(res),
                         list(1:4, 4:7, 7:10, 10:12))
            expect_equal(colnames(res), paste0("chr1_", 1:4))
          })

test_that("getCisGroups uses grpsize 25 and stepsize 10 by default",
          {
            peaks <- make_peaks(c(chr1 = 60))
            res <- getCisGroups(peaks)
            expect_equal(group_members(res),
                         list(1:25, 11:35, 21:45, 31:55, 41:60, 51:60))
            expect_equal(colnames(res), paste0("chr1_", 1:6))
          })

test_that("getCisGroups never mixes chromosomes and numbers groups per chromosome",
          {
            peaks <- make_peaks(c(chr1 = 9, chr2 = 7, chr3 = 8))
            res <- getCisGroups(peaks, grpsize = 3, stepsize = 2)
            expect_equal(group_members(res),
                         list(1:3, 3:5, 5:7, 7:9,
                              10:12, 12:14, 14:16,
                              17:19, 19:21, 21:23, 23:24))
            expect_equal(colnames(res),
                         c(paste0("chr1_", 1:4), paste0("chr2_", 1:3),
                           paste0("chr3_", 1:4)))
            chrs <- as.character(seqnames(peaks))
            for (ix in group_members(res)) {
              expect_equal(length(unique(chrs[ix])), 1)
            }
          })

test_that("getCisGroups returns the annotationMatches assay and keeps the peaks",
          {
            peaks <- make_peaks(c(chr1 = 20, chr2 = 15))
            res <- getCisGroups(peaks, grpsize = 5, stepsize = 5)
            expect_is(res, "RangedSummarizedExperiment")
            expect_equal(assayNames(res), "annotationMatches")
            expect_is(assay(res), "lMatrix")
            expect_equal(nrow(res), 35)
            expect_equal(rowRanges(res), peaks)
            expect_equal(ncol(res), 4 + 3)
          })

test_that("getCisGroups orders groups by the chromosome order of the peaks",
          {
            peaks <- GRanges(rep(c("chr2", "chr1"), c(6, 6)),
                             IRanges(start = rep(seq_len(6) * 100, 2), width = 50))
            seqlevels(peaks) <- c("chr2", "chr1")
            peaks <- sort(peaks)
            res <- getCisGroups(peaks, grpsize = 3, stepsize = 3)
            expect_equal(colnames(res), c("chr2_1", "chr2_2", "chr1_1", "chr1_2"))
            expect_equal(group_members(res), list(1:3, 4:6, 7:9, 10:12))
          })

test_that("getCisGroups on a RangedSummarizedExperiment equals GenomicRanges input",
          {
            peaks <- make_peaks(c(chr1 = 30, chr2 = 18))
            rse <- SummarizedExperiment(assays = list(counts = matrix(1, 48, 3)),
                                        rowRanges = peaks)
            res_gr <- getCisGroups(peaks, grpsize = 6, stepsize = 4)
            res_rse <- getCisGroups(rse, grpsize = 6, stepsize = 4)
            expect_equal(assay(res_rse), assay(res_gr))
            expect_equal(colnames(res_rse), colnames(res_gr))
            expect_equal(rowRanges(res_rse), rowRanges(res_gr))
            expect_equal(dim(getCisGroups(rse)), dim(getCisGroups(peaks)))
          })

test_that("getCisGroups ignores seqlevels without peaks",
          {
            peaks <- GRanges("chr1", IRanges(start = seq_len(8) * 100, width = 50))
            seqlevels(peaks) <- c("chr1", "chrUnused")
            res <- getCisGroups(peaks, grpsize = 4, stepsize = 2)
            expect_equal(colnames(res), paste0("chr1_", 1:4))
            expect_equal(group_members(res), list(1:4, 3:6, 5:8, 7:8))
          })

test_that("getCisGroups on sorted mini_counts gives contiguous single-chromosome groups",
          {
            data(mini_counts, package = "chromVAR")
            sorted <- sort(mini_counts)
            res <- getCisGroups(sorted)
            expect_equal(nrow(res), nrow(sorted))
            chrs <- as.character(seqnames(rowRanges(sorted)))
            members <- group_members(res)
            expect_true(all(lengths(members) >= 1 & lengths(members) <= 25))
            for (ix in members) {
              expect_equal(ix, seq(min(ix), max(ix)))
              expect_equal(length(unique(chrs[ix])), 1)
            }
            per_chr <- table(factor(chrs, levels = seqlevels(sorted)))
            expected_groups <- sum(per_chr[per_chr > 10] %/% 10)
            expect_equal(ncol(res), expected_groups)
          })

test_that("getCisGroups skips chromosomes with at most stepsize peaks and groups the others",
          {
            peaks <- make_peaks(c(chr1 = 5, chr2 = 12))
            res <- getCisGroups(peaks, grpsize = 4, stepsize = 5)
            expect_equal(colnames(res), c("chr2_1", "chr2_2"))
            expect_equal(group_members(res), list(6:9, 11:14))
            none <- getCisGroups(peaks, grpsize = 4, stepsize = 50)
            expect_equal(ncol(none), 0)
            expect_equal(nrow(none), 17)
          })

test_that("getCisGroups refuses unsorted peaks",
          {
            peaks <- make_peaks(c(chr1 = 10))
            expect_error(getCisGroups(peaks[c(5, 1, 3, 2, 4, 6:10)]),
                         "peaks must be sorted")
            rse <- SummarizedExperiment(assays = list(counts = matrix(1, 10, 2)),
                                        rowRanges = peaks[10:1])
            expect_error(getCisGroups(rse), "peaks must be sorted")
          })
