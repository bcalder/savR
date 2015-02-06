require(savR)
context("Data Accessor Test")

# use skip_on_cran() for long running tests

# test for avaliablitily of example data in extdata folder
miseq_available <- function() {
  if (!file.exists(system.file("extdata", "MiSeq", package="savR"))) {
    skip("Example data not available")
  }
}

# data folder not included in built package, skip if missing
data_available <- function() {
  if (!file.exists("../data")) {
    skip("Test data folder not available")
  }
}

test_that("data accessors", {
  miseq_available()
  fc <- savR(system.file("extdata", "MiSeq", package="savR"))
  expect_equal(dim(correctedIntensities(fc)), c(3116, 18), label = "correctedIntensities")
  expect_equal(dim(qualityMetrics(fc)), c(3116, 53), label = "qualityMetrics")
  expect_equal(dim(tileMetrics(fc)), c(648, 4), label = "tileMetrics")
  expect_equal(dim(extractionMetrics(fc)), c(3116, 11), label = "extractionMetrics")
})

test_that("errorMetrics", {
  skip_on_cran()
  data_available()
  fc <- savR("../data/AAF39")
  expect_equal(dim(errorMetrics(fc)), c(700, 9), label = "errorMetrics")
})