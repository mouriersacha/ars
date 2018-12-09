## TESTS

library(testthat)
library(ars)
context("tests")
library(stats)

test_that("h_x_d_d is two derivation after take log", {
  myfun <- function(x) exp(-x*x/2)
  x <- 3
  res <- h_x_d_d(myfun)(x)
  expect_lte(abs(res+1), 1e-6)
})

test_that("initialize the breakpoints and z", {
  mylist <- initialize(dnorm, -Inf, 5)
  expect_equal(mylist, list(c(-1,2.5), c(-Inf, 0.75, 5)))
  mylist <- initialize(dnorm, -5, Inf)
  expect_equal(mylist, list(c(-2.5,1), c(-5, -0.75, Inf)))
})

test_that("add point into breakpoints", {
  breakpoints=c(-2.5,2.5)
  z=c(-5,0,5)
  mylist <- add(dnorm, breakpoints, z, 1, -5, 5)
  expect_equal(mylist, list(c(-2.5, 1, 2.5),c(-5, -0.75, 1.75, 5)))
})

test_that("generate sample in the right section", {
  breakpoints <- c(-1.0, 2.5)
  z <- c(-Inf, 0.75, 5)
  mylist <- sample_one(dnorm, breakpoints, z)
  group <- mylist[[1]]
  candidate <- mylist[[2]]
  expect_lte(candidate, z[group+1])
  expect_gte(candidate, z[group])
})

test_that("unit test for normal distribution", {
  x <- dnorm
  samples <- ars(x, 200,-Inf, Inf)
  p <- ks.test(samples, "pnorm")
  expect_gte(p$p.value, 0.05)
})

test_that("unit test for beta(2,2) distribution", {
  x <- function(x) dbeta(x, 2, 2)
  samples <- ars(x, 200, 0.01, 0.99)
  p <- ks.test(samples, 'pbeta', 2, 2)
  expect_gte(p$p.value, 0.05)
})

test_that("unit test for gamma(2,2) distribution", {
  x <- function(x) dgamma(x, 2, 2)
  samples <- ars(x, 200, 0.01, Inf)
  p <- ks.test(samples, 'pgamma', 2, 2)
  expect_gte(p$p.value, 0.05)
})
