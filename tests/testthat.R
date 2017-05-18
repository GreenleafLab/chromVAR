library(testthat)
library(chromVAR)
BiocParallel::register(BiocParallel::SerialParam())
test_check("chromVAR")
