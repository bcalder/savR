require(savR)
context("Load InterOp File Test")

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

test_that("Load sample MiSeq data and check values", {
  miseq_available()
  fc <- savR(system.file("extdata", "MiSeq", package="savR"))
  expect_equivalent(run(fc)["Id"], "131030_M01243_0072_000000000-A58WM", label = "FCID check")
  expect_equal(cycles(fc), 82, label = "cycle count")
  expect_equal(directions(fc), 1, label = "directions")
  expect_equal(clusters(fc,1L),35470731, label="clusters")
})

test_that("Load sample HiSeq RapidRun (QMetricsOut.bin v5) data and check values", {
  skip_on_cran()
  data_available()
  expect_warning(fc <- savR("../data/BHF5GNADXX"), regexp = "incomplete row", label = "load RapidRun with truncated QMetricsOut.bin v5")
  expect_equivalent(run(fc)["Id"], "150115_SN7001401_0253_BHF5GNADXX", label = "FCID check")
  expect_equal(cycles(fc), 209, label = "cycle count")
  expect_equal(directions(fc), 2, label = "directions")
  expect_equal(clusters(fc,1L), 184790371, label="clusters")
})

test_that("Load sample MiSeq with ErrorMetricsOut.bin", {
  skip_on_cran()
  data_available()
  fc <- savR("../data/AAF39")
  expect_equal(dim(errorMetrics(fc)), c(700, 9))
})

test_that("Load problematic flowcell", {
  skip_on_cran()
  data_available()
  expect_warning(fc <- savR("../data/AC5J04ACXX"), regexp = "Unable to parse binary data: number of lanes in data file",
                 label="load truncated ErrorMetrics")
})


