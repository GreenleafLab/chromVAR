context("computeDeviations rejects invalid inputs")

library(SummarizedExperiment)
library(Matrix)

local_serial_param()

data(mini_counts, package = "chromVAR")

validation_input <- function() {
	local_seed(1)
	n <- nrow(mini_counts)
	bg <- getBackgroundPeaks(mini_counts, niterations = 5)
	anno <- matrix(FALSE, n, 3)
	for (a in 1:3) anno[sample(n, 20 + 10 * a), a] <- TRUE
	list(n = n, bg = bg, anno = anno, expectation = computeExpectations(mini_counts))
}

test_that("valid inputs run", {
	x <- validation_input()
	out <- computeDeviations(mini_counts, x$anno, x$bg)
	expect_equal(dim(deviations(out)), c(3L, ncol(mini_counts)))
})

test_that("background_peaks values outside 1..nrow(counts) are an R error", {
	x <- validation_input()
	for (bad in list(0L, -1L, x$n + 1L, 1000000L, NA_integer_)) {
		bg <- x$bg
		bg[, 1] <- bad
		expect_error(computeDeviations(mini_counts, x$anno, bg))
	}
	for (bad in list(NaN, Inf, -Inf)) {
		bg <- x$bg
		storage.mode(bg) <- "double"
		bg[, 1] <- bad
		expect_error(computeDeviations(mini_counts, x$anno, bg))
	}
})

test_that("background_peaks with a different number of rows than counts is an R error", {
	x <- validation_input()
	expect_error(computeDeviations(mini_counts, x$anno, x$bg[-x$n, ]))
	expect_error(computeDeviations(mini_counts, x$anno, rbind(x$bg, x$bg[1, ])))
})

test_that("background_peaks that is not a matrix is an R error", {
	x <- validation_input()
	expect_error(computeDeviations(mini_counts, x$anno, x$bg[, 1]))
	expect_error(computeDeviations(mini_counts, x$anno, NULL))
})

test_that("annotation indices outside 1..nrow(counts) are an R error", {
	x <- validation_input()
	counts <- counts(mini_counts)
	for (bad in list(0L, -1L, x$n + 1L, NA_integer_, 7.5)) {
		anno <- list(a = c(5L, 7L, 9L), b = c(bad, 11L, 13L))
		expect_error(computeDeviations(counts, anno, x$bg, x$expectation))
	}
	expect_error(computeDeviations(counts, list(), x$bg, x$expectation))
})

test_that("an annotation matrix with more rows than counts is an R error", {
	x <- validation_input()
	expect_error(computeDeviations(mini_counts, rbind(x$anno, TRUE), x$bg))
})

test_that("an expectation vector of the wrong length is an R error", {
	x <- validation_input()
	expect_error(computeDeviations(mini_counts, x$anno, x$bg, x$expectation[-1]))
	expect_error(computeDeviations(mini_counts, x$anno, x$bg, c(x$expectation, 1)))
	expect_error(computeDeviations(mini_counts, x$anno, x$bg, NULL))
})

test_that("counts with a different number of peaks than background_peaks is an R error", {
	x <- validation_input()
	expect_error(computeDeviations(mini_counts[-x$n, ], x$anno, x$bg))
	expect_error(computeDeviations(mini_counts[0, ], x$anno, x$bg))
})

test_that("counts that are negative, NA or NaN or have no cells are an R error", {
	x <- validation_input()
	for (bad in list(-2, NA_real_, NaN)) {
		m <- as.matrix(counts(mini_counts))
		m[which(x$anno[, 1])[1], 1] <- bad
		se <- SummarizedExperiment(assays = list(counts = m))
		expect_error(computeDeviations(se, x$anno, x$bg))
	}
	expect_error(computeDeviations(mini_counts[, 0], x$anno, x$bg))
})
