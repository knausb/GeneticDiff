
library(GeneticDiff)
#library(testthat)
context("rgt")

##### ##### ##### ##### #####

test_that("rgt returns a matrix", {
  myGT <- rgt()
  expect_is(myGT, "matrix")
})

test_that("rgt, all unphased", {
  myGT <- rgt(pphased = 0)
  expect_equal( length(grep("/", myGT)), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all phased", {
  myGT <- rgt(pphased = 1)
  expect_equal( length(grep("\\|", myGT)), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all haploid", {
  myGT <- rgt(pploid = c(1))
  expect_equal( sum(nchar(myGT) == 1), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all diploid", {
  myGT <- rgt(pploid = c(0,1))
  expect_equal( sum(nchar(myGT) == 3), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all triploid", {
  myGT <- rgt(pploid = c(0,0,1))
  expect_equal( sum(nchar(myGT) == 5), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all tetraploid", {
  myGT <- rgt(pploid = c(0,0,0,1))
  expect_equal( sum(nchar(myGT) == 7), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all 0", {
  myGT <- rgt(pploid = c(1), pallele=c(1))
  expect_equal(sum(myGT == 0), nrow(myGT) * ncol(myGT) )
})

test_that("rgt, all 1", {
  myGT <- rgt(pploid = c(1), pallele=c(0,1))
#  expect_equal(sum(myGT == 1), nrow(myGT) * ncol(myGT) )
})

##### ##### ##### ##### #####
# EOF.