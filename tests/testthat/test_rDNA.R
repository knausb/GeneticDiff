
library(GeneticDiff)
#
library(testthat)
context("rDNA")

##### ##### ##### ##### #####

test_that("rDNA returns a matrix", {
  myDNA <- rDNA()
  expect_is(myDNA, "matrix")
  expect_equal(nrow(myDNA), 1)
  expect_equal(ncol(myDNA), 100)
})


