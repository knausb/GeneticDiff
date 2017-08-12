
library(GeneticDiff)
#
library(testthat)
context("genotype_split")

##### ##### ##### ##### #####

test_that("genotype_split returns a vector", {
  myAlleles <- genotype_split("01/1/1|12")
  expect_true(is.vector(myAlleles))
  expect_equal(length(myAlleles), 4)
})


