#' ---
#' title: "Test: Coverage"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Coverage}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ include = FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#'
#+
library(testthat)
library(ramR)
#'
#+ E
A <- Omega <- matrix(
  0,
  ncol = 3,
  nrow = 2
)
test_that("E: square.", {
  expect_error(
    E(
      A = A
    )
  )
})
test_that("Sigmatheta: Omega is a square matrix.", {
  expect_error(
    Sigmatheta(
      A = A,
      Omega = Omega
    )
  )
})
A <- matrix(
  0,
  ncol = 3,
  nrow = 3
)
Omega <- matrix(
  0,
  ncol = 2,
  nrow = 2
)
test_that("Sigmatheta: A and Omega have identical dimensions.", {
  expect_error(
    Sigmatheta(
      A = A,
      Omega = Omega
    )
  )
})
A <- matrix(
  data = c(
    0,
    -1,
    -1,
    0
  ),
  nrow = 2,
  ncol = 2
)
test_that("E: inverse.", {
  expect_error(
    E(
      A = A
    )
  )
})
test_that("Omega_linreg: symmetric SigmaX.", {
  expect_error(
    Omega_linreg(
      sigma2 = 0,
      SigmaX = matrix(
        data = 0,
        ncol = 3,
        nrow = 2
      )
    )
  )
})
beta <- c(
  1
)
sigma2 <- 1
SigmaX <- 1
muX <- 1
test_that("ram_linreg: wrong dimensions beta.", {
  expect_error(
    ramR::ram_linreg(
      beta = beta,
      sigma2 = sigma2,
      SigmaX = SigmaX,
      muX = muX
    )
  )
})
beta <- c(
  1,
  1,
  1
)
sigma2 <- 1
SigmaX <- matrix(
  data = 1,
  nrow = 2,
  ncol = 2
)
muX <- 1
test_that("ram_linreg: wrong dimensions mux.", {
  expect_error(
    ramR::ram_linreg(
      beta = beta,
      sigma2 = sigma2,
      SigmaX = SigmaX,
      muX = muX
    )
  )
})
