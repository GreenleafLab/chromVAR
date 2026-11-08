# Shared helpers, sourced by testthat before the test files.

# Seed the RNG for the calling scope (a test_that block, a function or a test
# file) and put the previous RNG state back when that scope ends.
local_seed <- function(seed, .local_envir = parent.frame()) {
	withr::local_seed(seed, .local_envir = .local_envir)
}

# Register a SerialParam as the default BiocParallel back-end for the calling
# scope and restore the previously registered back-end afterwards.
local_serial_param <- function(param = BiocParallel::SerialParam(),
                               .local_envir = parent.frame()) {
	old <- BiocParallel::bpparam()
	BiocParallel::register(param)
	withr::defer(BiocParallel::register(old), envir = .local_envir)
	invisible(old)
}
